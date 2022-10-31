# VoronoiBoost

*VoronoiBoost*, is a data-driven model that scales Voronoi cells to match the probabilistic distribution of
users associated to each base station. 
*VoronoiBoost* relies on the same input as traditional Voronoi decompositions, but provides a richer and more accurate rendering of where users are located.

For more details please refer to our paper 'VoronoiBoost: Data-driven Probabilistic Spatial Mapping of Mobile Network Metadata' 
waiting for publication on Secon 2022 Conference proceeding (https://secon2022.ieee-secon.org/program/)


| Ground Truth                             | Voronoi                                    | VoronoiBoost                             |
| ---------------------------------------- | ------------------------------------------ | ---------------------------------------- |
| ![alt text](images/PARIS_810_p_l_t.png)  | ![alt text](images/PARIS_810_voronoi.png)  | ![alt text](images/PARIS_810_model.png)  |
| ![alt text](images/PARIS_1648_p_l_t.png) | ![alt text](images/PARIS_1648_voronoi.png) | ![alt text](images/PARIS_1648_model.png) |

## Installation

Clone this repository and install the requirements:

```bash
git clone https://github.com/nds-group/voronoiBoost.git

pip install -r requirements.txt
```

## Usage
First is need to import the VoronoiBoost class from the voronoiBoost.py file:

```python
import pandas as pd

from shapely.geometry import shape as Shape
from shapely.geometry import Polygon

from voronoiBoost import VoronoiBoost
```

**VoronoiBoost** use the same input as a standard voronoi tesselation, 
* set of points (base stations) 
* a border area that defines the area of interest

then, and instance can be lunch by provinding also the path to the trained model file.

```python
voronoiBoost = VoronoiBoost(sites, city_shape, model_path)

df_bs = voronoiBoost.compute_voronoiBoost()
df_bs.head(2)
```

The output is a pandas dataframe with the following columns:
``lon``, ``lat``, ``voronoi`` & ``voronois_scaled_overlap``

where ``voronoi`` is the voronoi cell associated to the base station and ``voronois_scaled_overlap`` is the scaled voronoi cell.

The following images show the legacy voronoi decomposition and the voronoiBoost decomposition for the same base stations. Other base stations are hidden for clarity.

| Legacy tessallation | VoronoiBoost |
| ------------------- | ------------ |
| ![alt text](images/Voronoi_Orange_Paris_same_color.png)  | ![alt text](images/VoronoiBoost_Orange_Paris_same_color.png) 