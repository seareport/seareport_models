import pyposeidon.model as pm
import pyposeidon.mesh as pmesh
import pandas as pd
import xarray as xr
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main():
    # general meshing settings (oceanmesh settings below)
    # Folders and files for mesh generation
    seaset_path = 'v0.1/seaset_full.csv'
    seaset_full = pd.read_csv(
        'https://raw.githubusercontent.com/tomsail/seaset/2176a61197d136121878a5412f39e3351d646af5/Notebooks/catalog_full.csv', 
        index_col=0
    ).to_csv(seaset_path)
    # 
    model = {
        # "type": "tri2d",
        # model
        "solver_name": 'telemac',
        "start_date": "2023-7-1 0:0:0",
        "end_date": "2023-7-1 23:0:0",
        # "meteo_source": wind,  # path to meteo files
        "meteo_input360": True,  # if meteo files longitudes go from from 0 to 360
        "monitor": True,  # get time series for observation points
        # "obs": seaset_path,
        "update": ["dem"],
        "parameters": {
            "dt": 400,
            "chezy": 30,
            "rnday": 30,
            "hotout": 0,
            "ihot": 0,
            "nspool": 9,
            "ihfskip": 36,
            "hotout_write": 108,
        },  # set param.nml components
    }

    mesh_file_0p1 = "meshes/global-v0.1.gr3"
    mesh_file_0p2 = "meshes/global-v0.2.gr3"
    mesh = pmesh.set(type='tri2d', mesh_file=mesh_file_0p1)
    z = mesh.Dataset.depth
    z[z < -20] = -20
    mesh.Dataset["depth"] = xr.Variable(('nSCHISM_hgrid_node'), z)
    mesh.to_file(mesh_file_0p2)

    # second step --- run model
    if model:
        for solver in ["schism", "telemac"]:
            rpath = 'v0.2/' + solver
            model["mesh_file"] = mesh_file_0p2
            model["rpath"] = rpath
            model["solver_name"] = solver
            model["coastlines"] = None # skip this step, we don't need coastlines anymore
            model['obs'] = seaset_path # apply obs only to export "station.in" file 
            model["update"]= ["model"]
            b = pm.set(**model)
            b.create()
            b.output()
            b.save()
            b.set_obs()

if __name__ == "__main__":
    main()