import os
import glob
import geopandas as gp
import pandas as pd
import numpy as np
import searvey
import logging
import xarray as xr
from pyposeidon.telemac import flip
from pyposeidon.utils import data
import pyposeidon.model as pm
import pyposeidon.meteo as pmeteo

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


def fix_mesh(b, corrections=None):
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


def get_meta() -> gp.GeoDataFrame:
    meta_web = searvey.get_ioc_stations().drop(columns=["lon", "lat"])
    meta_api = (
        pd.read_json(
            "http://www.ioc-sealevelmonitoring.org/service.php?query=stationlist&showall=all"
        )
        .drop_duplicates()
        .drop(columns=["lon", "lat"])
        .rename(columns={"Code": "ioc_code", "Lon": "lon", "Lat": "lat"})
    )
    merged = pd.merge(
        meta_web,
        meta_api[["ioc_code", "lon", "lat"]].drop_duplicates(),
        on=["ioc_code"],
    )
    return merged.drop(columns=["geometry"])


def ioc_subset_from_files_in_folder(
    df: pd.DataFrame, folder: str, ext: str = ".parquet"
):
    """this function return a subset of the ioc database from all the files (json or parquet)
    present in a folder
    """
    list_files = []
    for file in os.listdir(folder):
        name = file.split(ext)[0]
        if file.endswith(ext):
            list_files.append(name)
    return df[df.ioc_code.isin(list_files)]


YEAR = 2023
PROJECT = ""
WIND = glob.glob(PROJECT + f"02_meteo/era5/lon_lat/netcdf/{YEAR}*.nc")
WIND += glob.glob(PROJECT + f"02_meteo/era5/lon_lat/netcdf/{YEAR+1}*.nc")


def main(mesh: bool = True, model: bool = True, results=True):
    # general meshing settings (oceanmesh settings below)
    resolution = "f"
    res_min = 0.06
    res_max = 2
    cbuffer = res_min / 10
    # obs
    obs_path = "v1.2/ioc_.csv"
    ioc_ = get_meta()
    ioc_cleanup = ioc_subset_from_files_in_folder(ioc_, "../ioc_cleanup/clean/")
    ioc_cleanup.to_csv(obs_path)
    # Folders and files for mesh generation
    gshhs_folder = "../coastlines/out/"
    fdem = "00_bathy/etopo/ETOPO_0.03.nc"
    coasts = gp.read_parquet(gshhs_folder + f"gshhg_{resolution}.parquet")
    #
    rpath = "v1.2"
    MODEL = {
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
        "end_date": "2023-7-02 23:0:0",
        "meteo_input360": True,  # if meteo files longitudes go from from 0 to 360
        "meteo_source": None,  # path to meteo files
        "monitor": True,  # get time series for observation points
        # "obs": seaset_path,
        "update": ["dem"],
        "fortran": "./scripts/temp_fortran/out_history.F90", # can be a file or a folder
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

    # first step --- create mesh
    mesh_file = rpath + f"/GSHHS_{resolution}_{res_min}.gr3"
    if mesh:
        b = pm.set(**MODEL)
        b.create()
        b.mesh.to_file(mesh_file)

    # second step --- run model
    if model:
        for solver in ["schism", "telemac"]:
            rpath = f"v1.2/{solver}"
            MODEL["mesh_file"] = mesh_file
            MODEL["rpath"] = rpath
            MODEL["solver_name"] = solver
            MODEL["coastlines"] = None  # skip this step, we don't need coastlines
            MODEL["obs"] = obs_path  # apply obs only to export "station.in" file
            MODEL["update"] = ["model"]
            meteo = pmeteo.Meteo(WIND)
            meteo.Dataset = meteo.Dataset.sel(
                time=slice(MODEL["start_date"], MODEL["end_date"])
            )
            MODEL["meteo_source"] = meteo
            b = pm.set(**MODEL)
            b.create()
            b.output()
            if solver == "schism":
                corr = {"reverse": [2], "remove": []}  # this fix is for schism
            else:
                corr = None
            fix_mesh(b, corrections=corr)
            b.mesh.to_file(f"v1.2/{solver}/hgrid.gr3")
            b.save()
            b.set_obs()
            b.run()

    # third step --- extract results
    if results:
        for solver in ["telemac"]:  # schism not tested
            MODEL["rpath"] = f"v1.2/{solver}"
            MODEL["solver_name"] = solver
            MODEL["result_type"] = "2D"
            MODEL["convert_results"] = True
            MODEL["extract_TS"] = True
            res = data.get_output(**MODEL)


if __name__ == "__main__":
    main()
