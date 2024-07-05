# Seareport Models
Concept: Build reproducible models in 3 steps:
 1. `mesh`: reproduce the meshing
 2. `model`: reproduce the model
 3. `results`: reproduce the results data extraction (2D and 1D Time Series)

## 1 - Installl
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

## 2 - Reproduce
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