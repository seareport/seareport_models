import searvey
import glob
import pandas as pd
import geopandas as gp
from pyposeidon.utils.obs import serialize_stations
from pyposeidon.utils.cpoint import find_nearest_nodes
from pyposeidon.utils.cfl import parse_hgrid

def get_stofs2d_meta():
    stofs2d = pd.read_csv(
        "https://polar.ncep.noaa.gov/stofs/data/stofs_2d_glo_elev_stat_v2_1_0",
        names=["coords", "name"],
        sep="!",
        header=None,
        skiprows=1
    )
    stofs2d = stofs2d.assign(
        lon=stofs2d.coords.str.split("\t", n=1).str[0].astype(float),
        lat=stofs2d.coords.str.strip().str.rsplit("\t", n=1).str[1].astype(float),
        stofs2d_name=stofs2d.name.str.strip(),
    ).drop(columns=["coords", "name"])
    return stofs2d


def get_ioc_meta() -> gp.GeoDataFrame:
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
    return merged


def merge_ioc_and_stofs(ioc: pd.DataFrame, stofs2d: pd.DataFrame) -> pd.DataFrame:
    stations = pd.concat((ioc, stofs2d), ignore_index=True)
    stations = stations.assign(unique_id=stations.ioc_code.combine_first(stations.stofs2d_name))
    return stations

def merge_stofs_ioc(obs_path:str):
    ioc = get_ioc_meta()
    stofs2d = get_stofs2d_meta()
    m = merge_ioc_and_stofs(ioc=ioc, stofs2d=stofs2d)
    m.to_csv(obs_path)
    return m

def export_station_in(stations: pd.DataFrame):
    for vmodel in sorted(glob.glob("v*.*/")):
        print(vmodel)
        for model in ["schism", "telemac"]:
            if len(glob.glob(f"{vmodel}{model}/*hgrid.gr3"))>0:
                hgrid=glob.glob(f"{vmodel}{model}/*hgrid.gr3")[0]
            elif len(glob.glob(f"{vmodel}{model}/*/hgrid.gr3")):
                hgrid=glob.glob(f"{vmodel}{model}/*/hgrid.gr3")[0]
            else:
                hgrid = None
            if hgrid:
                mesh = parse_hgrid(hgrid)
                mesh_nodes = find_nearest_nodes(mesh["nodes"],stations)
                if model == "telemac":
                    mesh_nodes.index += 1
                    mesh_nodes["ind"] = mesh_nodes.index
                serialize_stations(mesh_nodes.iloc[:-4], f"{vmodel}/{model}/station.in", format=model, duration = 3600*24*365, timestep=400)
                print(f"  > exported {vmodel}/{model}/station.in")

if __name__ == "__main__":
    # s = merge_stofs_ioc("scripts/ioc_stofs.csv")
    s = pd.read_csv("scripts/ioc_cleanup_2024.csv", index_col=0)
    s['unique_id'] = s['ioc_code']
    export_station_in(s)
    