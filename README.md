## example for v0 build

to recreate the environment to build the mesh: 

first build the binaries: 
```
cd v0
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
python scripts/mesh_v0.py
```
