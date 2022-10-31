from geovoronoi import voronoi_regions_from_coords
import networkx as nx
import pickle
import numpy as np
import pandas as pd
import math
from shapely.geometry import shape as Shape
from shapely.geometry import Point, Polygon, MultiPolygon, mapping

import json
from sortedcontainers import SortedSet
import itertools
from scipy.spatial import ConvexHull
from shapely.geometry import MultiPoint
from shapely import affinity
import mapply
import lzma

mapply.init(
    n_workers=40,
    chunk_size=20,
    progressbar=False
)

def earth_distance(origin, destination):
    lon1, lat1 = origin
    lon2, lat2 = destination
    radius = 6371000  # meters!!!!

    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) * math.sin(dlat / 2) +
         math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) *
         math.sin(dlon / 2) * math.sin(dlon / 2))
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
    d = radius * c

    return round(d, 2)

class VoronoiBoost:

    def __init__(self, sites, border, model_path):
        
        self.taus = [.25, .5, .75, .85, .95]
        self.taus_scales = {}
        
        self.border = border
        self.sites = np.array(sites)
        self.n = len(sites)
        self.ids = list(range(len(sites)))
        self.lats = [site[0] for site in sites]
        self.lons = [site[1] for site in sites]
        self.model_path = model_path

        self.df_bs = pd.DataFrame(data={
                        'id': self.ids, 
                        'lon': self.lons,
                        'lat': self.lats,
                        })

        self.g_delaunay = nx.Graph()
        self.model = None

        self.selected_features = [
            'mean_d_neighbors_6_1_div',
            'mean_d_neighbors_8_1_div',
            'mean_d_neighbors_10_1_div',
            'mean_distance_6_1_div',
            'mean_distance_8_1_div',
            'mean_distance_10_1_div',
            'min_d_neighbors_8',
            'convex_hull_width_8_1_div',
            'mean_area_neighbors_4_1_div',
            'mean_area_neighbors_6_1_div',
            'mean_area_neighbors_8_1_div',
            'mean_d_t_barc_v_10_1_div',
            'd_v_max',
            'v_diameter_v_width_div',
            ]
        
    
    def load_model(self):
        print('Loading model 📦...')
        fd = lzma.open(self.model_path, 'rb')
        self.model = pickle.load(fd)
        fd.close()
        return self.model

    def compute_voronoi_tessellation(self):
        print('Computing Voronoi Tessellation 💻...')
        region_polys, region_pts = voronoi_regions_from_coords(self.sites, self.border)

        # Check that the Voronoi polygons are valid
        if len(region_pts) != self.n:
            # show a sample of the points that were assigned to more than one polygon
            print(list(filter(lambda k_v: len(k_v[1]) > 1, region_pts.items()))[0:10])
            raise Exception('Number of sites and assignments do not match 😨')
        else:
            print(f'Voronoi Tessellation successful 🤩.')
        voronois = [0 for i in range(self.n)]

        for voronoi_index, pts_index,  in region_pts.items():
            pts_index = pts_index[0] # only one point per polygon
            voronois[pts_index] = region_polys[voronoi_index]

        # if multiple polygons are returned, keep only the one that contains the BS
        for index, (lat, lon, voronoi) in enumerate(zip(self.lats, self.lons, voronois)):
            if voronoi.type == 'Polygon':
                continue
            # else is a MultiPolygon
            for polygon in voronoi:
                if Point(lat, lon).within(polygon):
                    voronois[index] = polygon

        self.df_bs['voronoi'] = voronois

        return list(self.df_bs['voronoi'])


    def compute_delaunay(self):
        print('Computing Delaunay Graph 💻...')
        map_node_coords = {}
        vertexs = {}

        # set of voronois that share a vertex
        for bs in self.df_bs.to_dict(orient='records'):
            bs_id = bs['id']
            lat, lon = bs['lat'], bs['lon']
            voronoi = bs['voronoi']
            
            voronoi_lats = voronoi.exterior.coords.xy[0]
            voronoi_lons = voronoi.exterior.coords.xy[1]
                
            for lat, lon in zip(voronoi_lats, voronoi_lons):
                vertex = (lat, lon)
                if vertex not in vertexs:
                    vertexs[vertex] = SortedSet()
                vertexs[vertex].add(bs_id)

        # add nodes to the graph
        for bs in self.df_bs.to_dict(orient='records'):
            bs_id = bs['id']
            lat, lon = bs['lat'], bs['lon']
            self.g_delaunay.add_node(bs_id, lat=lat, lon=lon)

        # add edges to the graph
        for vertex in vertexs:
            # all pair of BSs that share a vertex
            # given that self.g_delaunay is undirected, we only need to add one edge
            for pair in itertools.combinations(vertexs[vertex], r=2):
                node_1 = pair[0]
                node_2 = pair[1]

                node_1_coord = self.sites[node_1]
                node_2_coord = self.sites[node_2]
                distance = earth_distance(node_1_coord, node_2_coord)
                self.g_delaunay.add_edge(node_1, node_2, weight=distance)
        
        print(f'Delaunay Graph computed 🤩.')
        return self.g_delaunay


    def get_average_d_neighbors(self, node_id, level):
        neighbors_level_down = list(nx.single_source_shortest_path_length(self.g_delaunay, node_id, cutoff=level-1).keys())
        neighbors_level = list(nx.single_source_shortest_path_length(self.g_delaunay, node_id, cutoff=level).keys())
        neighbors_level_exclusive = set(neighbors_level) - set(neighbors_level_down)

        if len(neighbors_level_exclusive) == 0:

            if len(neighbors_level) == self.n:
                min_distance, mean_distance = self.get_average_d_neighbors(node_id, level-1)
                return min_distance, mean_distance
            else:
                min_distance, mean_distance = self.get_average_d_neighbors(node_id, level+1)
                return min_distance, mean_distance

        site = self.sites[node_id]
        distances = []
        for node_id_neighbor in neighbors_level_exclusive:
            site_neighbor = self.sites[node_id_neighbor]
            distance = earth_distance(site, site_neighbor)
            distances.append(distance)

        min_distance = min(distances)
        mean_distance = np.mean(distances)

        return min_distance, mean_distance


    def get_distance_between_neighbors_level(self, node_id, level):
        neighbors_level_down = list(nx.single_source_shortest_path_length(self.g_delaunay, node_id, cutoff=level-1).keys())
        neighbors_level = list(nx.single_source_shortest_path_length(self.g_delaunay, node_id, cutoff=level).keys())
        # set of exclusive neighbors in level i
        neighbors_level_exclusive = set(neighbors_level) - set(neighbors_level_down)

        # if in this level the number of neighboords dont allow to compute the feature, go to the next level
        if len(neighbors_level_exclusive) < 2:
            if len(neighbors_level) == self.n:
                mean_distance = self.get_distance_between_neighbors_level(node_id, level-1)
                return mean_distance
            else:
                mean_distance = self.get_distance_between_neighbors_level(node_id, level+1)
                return mean_distance

        # set of all possible pairs of neighbors
        neighbors_level_pairs = list(itertools.combinations(neighbors_level_exclusive, 2))

        distances = []
        for node_id_1, node_id_2 in neighbors_level_pairs:
            site_1 = self.sites[node_id_1]
            # if node_id_2 in set(self.g_delaunay.neighbors(node_id_1)): # if node_id_2 is neighbor of node_id_1 and viceversa
            site_2 = self.sites[node_id_2]
            distance = earth_distance(site_1, site_2)
            distances.append(distance)

        if len(distances) == 0:
            raise Exception('No pair of neighbors in level i 😨')

        mean_distance = np.mean(distances)

        return mean_distance


    def get_d_vk(self, node_id): # distance to vertexs k (k=3)
        site = self.sites[node_id]

        bs = self.df_bs[self.df_bs['id'] == node_id].iloc[0]
        voronoi = bs['voronoi']        
        voronoi_lats, voronoi_lons = voronoi.exterior.coords.xy
        vertexs = list(zip(voronoi_lats, voronoi_lons))

        distances = [earth_distance(site, vertex) for vertex in vertexs]
    
        v_max = max(distances)
        v_mean = np.mean(distances)

        return v_max


    def get_voronoi_diameter_width(self, node_id):
        bs = self.df_bs[self.df_bs['id'] == node_id].iloc[0]
        voronoi = bs['voronoi']        
        
        minimum_rotated_rectangle = voronoi.minimum_rotated_rectangle
        lats, lons = minimum_rotated_rectangle.exterior.coords.xy

        side_1 = earth_distance((lats[0], lons[0]), (lats[1], lons[1]))
        side_2 = earth_distance((lats[1], lons[1]), (lats[2], lons[2]))
        
        voronoi_diameter = max([side_1, side_2]) # diameter
        voronoi_width = min([side_1, side_2]) # width
    
        return voronoi_diameter, voronoi_width


    def get_average_area_neighbors(self, node_id, level):
        # TODO: train model with latlon area
        neighbors_level = list(nx.single_source_shortest_path_length(self.g_delaunay, node_id, cutoff=level).keys())

        voronoi_areas = []
        for node_id_neighbor in neighbors_level:
            bs = self.df_bs[self.df_bs['id'] == node_id_neighbor].iloc[0]
            voronoi = bs['voronoi']
            voronoi_areas.append(voronoi.area)

        mean_area = np.mean(voronoi_areas)
        return mean_area


    def get_convex_hull_area_perimeter_diameter_width(self, node_id, level):

        neighbors_level = list(nx.single_source_shortest_path_length(self.g_delaunay, node_id, cutoff=level).keys())
        points = []
        for node_id_neighbor in neighbors_level:
            points.append(self.sites[node_id_neighbor])

        # if in this level the number of neighboords dont allow to compute the feature, go to the next level
        if len(points) < 3:
            hull_width = self.get_convex_hull_area_perimeter_diameter_width(node_id, level+1)
            return hull_width

        hull = ConvexHull(points)
        hull = MultiPoint(hull.points).convex_hull
   
        minimum_rotated_rectangle = hull.minimum_rotated_rectangle
        lons, lats = minimum_rotated_rectangle.exterior.coords.xy

        side_1 = earth_distance((lats[0], lons[0]), (lats[1], lons[1]))
        side_2 = earth_distance((lats[1], lons[1]), (lats[2], lons[2]))
        
        hull_diameter = max([side_1, side_2]) # diameter
        hull_width = min([side_1, side_2]) # width

        return hull_width


    def get_average_d_t_barc_v(self, node_id, level):
        neighbors_level = list(nx.single_source_shortest_path_length(self.g_delaunay, node_id, cutoff=level).keys())
        
        distances_to_centroid = []
        for node_id_neighbor in neighbors_level:
            site = self.sites[node_id_neighbor]

            bs = self.df_bs[self.df_bs['id'] == node_id_neighbor].iloc[0]
            voronoi = bs['voronoi']     
            voronoi_centroid = np.array(voronoi.centroid.coords[0])

            distance_to_centroid = earth_distance(site, voronoi_centroid)
            distances_to_centroid.append(distance_to_centroid)

        mean_distance_to_centroid = np.mean(distances_to_centroid)
        return mean_distance_to_centroid

    def get_features(self):
        print('Computing features 🤖...')
        self.df_bs[['min_d_neighbors_1', 'mean_d_neighbors_1']] = self.df_bs.mapply(lambda row: self.get_average_d_neighbors(row['id'], 1), axis=1, result_type='expand')
        self.df_bs[['min_d_neighbors_6', 'mean_d_neighbors_6']] = self.df_bs.mapply(lambda row: self.get_average_d_neighbors(row['id'], 6), axis=1, result_type='expand')
        self.df_bs[['min_d_neighbors_8', 'mean_d_neighbors_8']] = self.df_bs.mapply(lambda row: self.get_average_d_neighbors(row['id'], 8), axis=1, result_type='expand')
        self.df_bs[['min_d_neighbors_10', 'mean_d_neighbors_10']] = self.df_bs.mapply(lambda row: self.get_average_d_neighbors(row['id'], 10), axis=1, result_type='expand')
        
        self.df_bs['mean_d_neighbors_6_1_div'] = self.df_bs['mean_d_neighbors_6']/self.df_bs['mean_d_neighbors_1']
        self.df_bs['mean_d_neighbors_8_1_div'] = self.df_bs['mean_d_neighbors_8']/self.df_bs['mean_d_neighbors_1']
        self.df_bs['mean_d_neighbors_10_1_div'] = self.df_bs['mean_d_neighbors_10']/self.df_bs['mean_d_neighbors_1']

        self.df_bs['mean_distance_1'] = self.df_bs.mapply(lambda row: self.get_distance_between_neighbors_level(row['id'], 1), axis=1)
        self.df_bs['mean_distance_6'] = self.df_bs.mapply(lambda row: self.get_distance_between_neighbors_level(row['id'], 6), axis=1)
        self.df_bs['mean_distance_8'] = self.df_bs.mapply(lambda row: self.get_distance_between_neighbors_level(row['id'], 8), axis=1)
        self.df_bs['mean_distance_10'] = self.df_bs.mapply(lambda row: self.get_distance_between_neighbors_level(row['id'], 10), axis=1)
        
        self.df_bs['mean_distance_6_1_div'] = self.df_bs['mean_distance_6']/self.df_bs['mean_distance_1']
        self.df_bs['mean_distance_8_1_div'] = self.df_bs['mean_distance_8']/self.df_bs['mean_distance_1']
        self.df_bs['mean_distance_10_1_div'] = self.df_bs['mean_distance_10']/self.df_bs['mean_distance_1']

        self.df_bs['d_v_max'] = self.df_bs.mapply(lambda row: self.get_d_vk(row['id']), axis=1)

        self.df_bs[['v_diameter', 'v_width']] = self.df_bs.mapply(lambda row: self.get_voronoi_diameter_width(row['id']), axis=1, result_type='expand')
        self.df_bs['v_diameter_v_width_div'] = self.df_bs.mapply(lambda row: row['v_diameter']/row['v_width'], axis=1)


        self.df_bs['mean_area_neighbors_1'] = self.df_bs.mapply(lambda row: self.get_average_area_neighbors(row['id'], 1), axis=1)
        self.df_bs['mean_area_neighbors_4'] = self.df_bs.mapply(lambda row: self.get_average_area_neighbors(row['id'], 4), axis=1)
        self.df_bs['mean_area_neighbors_6'] = self.df_bs.mapply(lambda row: self.get_average_area_neighbors(row['id'], 6), axis=1)
        self.df_bs['mean_area_neighbors_8'] = self.df_bs.mapply(lambda row: self.get_average_area_neighbors(row['id'], 8), axis=1)


        self.df_bs['mean_area_neighbors_4_1_div'] = self.df_bs['mean_area_neighbors_4'] / self.df_bs['mean_area_neighbors_1']
        self.df_bs['mean_area_neighbors_6_1_div'] = self.df_bs['mean_area_neighbors_6'] / self.df_bs['mean_area_neighbors_1']
        self.df_bs['mean_area_neighbors_8_1_div'] = self.df_bs['mean_area_neighbors_8'] / self.df_bs['mean_area_neighbors_1']


        self.df_bs['convex_hull_width_1'] = self.df_bs.mapply(lambda row: self.get_convex_hull_area_perimeter_diameter_width(row['id'], 1), axis=1)
        self.df_bs['convex_hull_width_8'] = self.df_bs.mapply(lambda row: self.get_convex_hull_area_perimeter_diameter_width(row['id'], 8), axis=1)
        self.df_bs['convex_hull_width_8_1_div'] = self.df_bs['convex_hull_width_8']/self.df_bs[f'convex_hull_width_1']


        self.df_bs['mean_d_t_barc_v_1'] = self.df_bs.mapply(lambda row: self.get_average_d_t_barc_v(row['id'], 1), axis=1)
        self.df_bs['mean_d_t_barc_v_10'] = self.df_bs.mapply(lambda row: self.get_average_d_t_barc_v(row['id'], 10), axis=1)
        self.df_bs['mean_d_t_barc_v_10_1_div'] = self.df_bs['mean_d_t_barc_v_10']/self.df_bs['mean_d_t_barc_v_1']


        print('Features computed 🤓.')
        self.df_bs = self.df_bs[ ['id', 'lon', 'lat', 'voronoi'] + self.selected_features ].copy()
        return self.df_bs

    def get_prediction(self):
        print('Predicting 🤖...')
        for tau in self.taus:
            tau_column = f'tau_{tau}'
            self.df_bs[tau_column] = tau

            input = self.df_bs[ [f'tau_{tau}'] + 
                            self.selected_features ].to_numpy()

            self.df_bs[f'pred_opt_scale_{tau}'] = self.model.predict(input)
        
        return self.df_bs

    def _get_correct_scales(self, row):
        
        taus_scales = []
        for tau in self.taus:
            pred_opt_scale_tau = row[f'pred_opt_scale_{tau}']
            pred_opt_scale_tau = round(pred_opt_scale_tau, 3)
            taus_scales.append((tau, pred_opt_scale_tau))
        
        #reverse 
        taus_scales = taus_scales[::-1]

        taus_scales_decrease = [taus_scales[0]]

        for index in range(len(taus_scales)-1):
            
            if taus_scales[index+1][1] < taus_scales[index][1]:
                taus_scales_decrease.append(taus_scales[index+1])
        
        # reverse again
        taus_scales = taus_scales_decrease[::-1]
        all_taus = len(self.taus) == len(taus_scales)

        return taus_scales, all_taus
    
    def get_corrected_scales(self):
        print('Correcting scales ➕...')

        self.df_bs[['taus_scales', 'all_taus']] = self.df_bs.mapply(lambda row: self._get_correct_scales(row), 
                                                        axis=1, 
                                                        result_type='expand')

        print('Scales corrected ✅.')

        return self.df_bs

    def _get_voronois_overlap(self, row):
        voronoi = row['voronoi']
        taus_scales = row['taus_scales']
        taus = [tau_scale[0] for tau_scale in taus_scales]

        voronois_scaled = []
        for tau, scale in taus_scales:
            voronoi_scaled = affinity.scale(voronoi, scale, scale)
            voronois_scaled.append(voronoi_scaled)

        return list(zip(voronois_scaled, taus))


    def get_voronois_overlap(self):
        print('Computing Voronois overlap 🤓...')
        
        self.df_bs['voronoi_boost'] = self.df_bs.mapply(lambda row: self._get_voronois_overlap(row), axis=1)
        
        print('Voronois overlap computed ✅.')
        return self.df_bs

    def compute_voronoiBoost(self):
        print('Computing VoronoiBoost 🤓...')
        self.compute_voronoi_tessellation()
        self.compute_delaunay()

        self.get_features()
        self.load_model()
        self.get_prediction()

        self.get_corrected_scales()
        self.get_voronois_overlap()

        self.df_bs = self.df_bs[['id', 'lat', 'lon', 'voronoi', 'voronoi_boost']].copy()

        print('VoronoiBoost computed ✅.')
        return self.df_bs