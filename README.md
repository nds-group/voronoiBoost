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
To load the model and predict the VoronoiBoost you can use the following code:

```python
import pickle

filename = f'./model/voronoiBoost.sav'
model = pickle.load(open(filename, 'rb'))
model
```