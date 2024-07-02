import pyposeidon.model as pm
import geopandas as gp
import pandas as pd
import logging
import xarray as xr
from pyposeidon.telemac import flip
import numpy as np

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def dist(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)


def center_pole(x, y):
    idx = np.where(dist(x, y, -180, 90) < 10e-6)
    x[idx] = 0
    return x


def is_overlapping(tris, meshx):
    PIR = 180
    x1, x2, x3 = meshx[tris].T
    condition = (
        # fmt: off
          ((x1 * x2 < 0) & (np.abs(x1 - x2) >= PIR))
        | ((x1 * x3 < 0) & (np.abs(x1 - x3) >= PIR))
        | ((x2 * x3 < 0) & (np.abs(x2 - x3) >= PIR))
        # fmt: on
    )
    return condition


def flip(tris):
    return np.column_stack((tris[:, 2], tris[:, 1], tris[:, 0]))


def fix_mesh(b, corrections = None):
    tri = b.mesh.Dataset.SCHISM_hgrid_face_nodes
    x = b.mesh.Dataset.SCHISM_hgrid_node_x
    y = b.mesh.Dataset.SCHISM_hgrid_node_y
    x = center_pole(x, y)
    m = is_overlapping(tri, x)
    tri[m] = flip(tri[m])
    rem_mask = np.zeros(len(tri), dtype=bool)
    # manual additional corrections
    if corrections is not None:
        for rev in corrections["reverse"]:
            tri[rev : rev + 1] = flip(tri[rev : rev + 1])
        for rem in corrections["remove"]:
            rem_mask[rem] = True
    b.mesh.Dataset["SCHISM_hgrid_face_nodes"] = xr.Variable(
        ("nSCHISM_hgrid_face", "nMaxSCHISM_hgrid_face_nodes"), tri[~rem_mask, :]
    )


def main(mesh: bool = True, model: bool = True):
    # general meshing settings (oceanmesh settings below)
    resolution = "f"
    res_min = 0.06
    res_max = 2
    cbuffer = res_min / 10
    # Folders and files for mesh generation
    seaset_path = "v1.2/seaset_full.csv"
    seaset_full = pd.read_csv(
        "https://raw.githubusercontent.com/tomsail/seaset/2176a61197d136121878a5412f39e3351d646af5/Notebooks/catalog_full.csv",
        index_col=0,
    ).to_csv(seaset_path)
    gshhs_folder = "../coastlines/out/"
    fdem = "00_bathy/etopo/ETOPO_0.03.nc"
    coasts = gp.read_parquet(gshhs_folder + f"gshhg_{resolution}.parquet")
    #
    rpath = "v1.2"
    model = {
        # "type": "tri2d",
        "coastlines": coasts,
        "mesh_generator": "oceanmesh",
        "geometry": "global",
        "cbuffer": cbuffer,
        "dem_source": fdem,
        "bgmesh": "om",
        "resolution_min": res_min,
        "resolution_max": res_max,
        # oceanmesh settings
        "grad": 0.12,
        "iterations": 200,
        "bathy_gradient": True,
        "alpha_wavelength": 100,  # number of element to resolve WL
        "alpha_slope": 10,  # number of element to resolve bathy gradient
        "plot": False,
        "rpath": rpath,
        # model
        "solver_name": "telemac",
        "start_date": "2023-7-1 0:0:0",
        "end_date": "2023-7-1 23:0:0",
        "monitor": True,  # get time series for observation points
        # "obs": seaset_path,
        "update": ["dem"],
        "parameters": {
            "dt": 300,
            "chezy": 30,
            "rnday": 30,
            "hotout": 0,
            "ihot": 0,
            "nspool": 9,
            "ihfskip": 36,
            "hotout_write": 108,
        },  # set param.nml components
    }

    # first step --- create mesh
    mesh_file = rpath + f"/GSHHS_{resolution}_{res_min}.gr3"
    if mesh:
        b = pm.set(**model)
        b.create()
        b.mesh.to_file(mesh_file)

    # second step --- run model
    if model:
        for solver in ["schism", "telemac"]:
            rpath = "v1.2/" + solver
            model["mesh_file"] = mesh_file
            model["rpath"] = rpath
            model["solver_name"] = solver
            model["coastlines"] = None  # skip this step, we don't need coastlines
            model["obs"] = seaset_path  # apply obs only to export "station.in" file
            model["update"] = ["model"]
            b = pm.set(**model)
            b.create()
            b.output()
            corr = {"reverse": [4], "remove": []}
            fix_mesh(b, corrections=corr)
            b.mesh.to_file("v1.2/schism/hgrid.gr3")
            b.save()
            b.set_obs()


if __name__ == "__main__":
    main(mesh=True)
