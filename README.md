# Seareport Models
Concept: Build reproducible models in 3 steps:
 1. `mesh`: reproduce the meshing
 2. `model`: reproduce the model
 3. `results`: reproduce the results data extraction (2D and 1D Time Series)

Example: Using `scripts/model_v2.2.py`, you can produce this 3km global mesh and model metocean data
![world_2.2](./assets/v2.2.png)

## current meshes/models available

| Mesh version | mesher used| resolution | hole in the north pole | bathy gradient | other              | model application | solver supported | Number of nodes | Number of elements |
|--------------|------------|------------|------------------------|----------------|--------------------|-------------------|------------------|-----------------|--------------------|
| `v0.0`       | oceanmesh  | 50km       |         :x:            |      yes       |         -          |        2D         | SCHISM/TELEMAC2D |      84,689     |     164,735        |
| `v0.1`       | oceanmesh  | 20km       |         :x:            |      yes       |         -          |        2D         | SCHISM/TELEMAC2D |     365,494     |     716,158        |
| `v0.2`       | oceanmesh  | 20km       |         :x:            |      yes       |  max depth at -20m |        2D         | SCHISM/TELEMAC2D |     365,494     |     716,158        |
| `v0.3`       | oceanmesh  | 30km       |         yes            |      :x:       |         -          |      waves        | TELEMAC-TOMAWAC  |      77,669     |     146,413        |
| `v0.4`       | oceanmesh  | 20km       |         yes            |      yes       |         -          |        3D         |  TELEMAC3D       |     369,029     |     723,382        |
| `v1.2`       | oceanmesh  | 6km        |         :x:            |      yes       |         -          |        2D         |  SCHISM/TELEMAC  |   1,613,172     |   3,154,713        |
| `v1.3`       | oceanmesh  | 6km        |         yes            |      :x:       |         -          |      waves        | TELEMAC-TOMAWAC  |     591,858     |   1,114,831        |
| `v1.4`       | oceanmesh  | 6km        |         yes            |      yes       |         -          |        3D         |  TELEMAC3D       |   1,538,748     |   3,006,764        |
| `v1.5`       | JIGSAW     | 6km        |         yes            |      yes       |         -          |        2D         | SCHISM/TELEMAC2D |   1,657,207     |   3,247,598        |
| `v2.2`       | oceanmesh  | 3km        |         :x:            |      yes       |         -          |        2D         | SCHISM/TELEMAC2D |   3,028,101     |   5,887,129        |
| `v2.3`       | JIGSAW     | 3km        |         :x:            |      yes       |         -          |        2D         |  TELEMAC2D       |   2,827,515     |   5,464,377        |
| `v3.0`       | JIGSAW     | 2km        |         :x:            |      yes       |         -          |        2D         |  TELEMAC2D       |   4,208,101     |   8,096,114        |
| `v3.1`       | JIGSAW     | 1km        |         :x:            |      yes       |         -          |        2D         |  TELEMAC2D       |   7,304,154     |  13,873,525        |
| `v3.2`       | oceanmesh  | 1km        |         :x:            |      yes       |         -          |        2D         | SCHISM/TELEMAC2D |   7,990,779     |  15,312,900        |

## Install
first build the binaries: 
```
cd v0.0
mamba env create -n v0 -f binary-om-telemac-openmpi-p3.11.yml
```

then build the correct python libraries 
```
python -mvenv .venv
source .venv/bin/activate
poetry install
```
then run the meshing script: 
```
cd ..
python scripts/mesh_v0.0.py
```

## Mesh, model and/or export results
By default, all 3 steps are activated, to deactivate meshing:

change:
```python
if __name__ == "__main__":
    main()
```
into:
```python
if __name__ == "__main__":
    main(mesh = False)
```

in `scripts/mesh_v0.0.py`