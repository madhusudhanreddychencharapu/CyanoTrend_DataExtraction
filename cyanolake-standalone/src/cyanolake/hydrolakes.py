"""Standalone hydrolakes implementation."""

from __future__ import annotations

import math
import os
import zipfile
from pathlib import Path

import geopandas as gpd
import netCDF4
import numpy as np
import pandas as pd
import pyogrio
import requests
from shapely import box, contains_xy, union_all
from tqdm.auto import tqdm

from . import science as _m_science
from . import settings as _m_settings


# Reference cell 29, lines 81-106.
def download_stream(url: str, destination: Path) -> Path:
    """Resumable streaming download used only when the user requests HydroLAKES."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    part = destination.with_suffix(destination.suffix + ".part")
    start = part.stat().st_size if part.exists() else 0
    headers = {"Range": f"bytes={start}-"} if start else {}
    response = requests.get(url, headers=headers, stream=True, timeout=(60, 900))
    if response.status_code == 416 and part.exists() and zipfile.is_zipfile(part):
        os.replace(part, destination)
        return destination
    if start and response.status_code != 206:
        start = 0
    response.raise_for_status()
    mode = "ab" if start and response.status_code == 206 else "wb"
    total = int(response.headers.get("Content-Length", 0)) + (start if mode == "ab" else 0)
    with (
        part.open(mode) as handle,
        tqdm(
            total=total or None,
            initial=start if mode == "ab" else 0,
            unit="B",
            unit_scale=True,
            desc=destination.name,
        ) as bar,
    ):
        for chunk in response.iter_content(8 * 1024 * 1024):
            if chunk:
                handle.write(chunk)
                bar.update(len(chunk))
    os.replace(part, destination)
    return destination


# Reference cell 29, lines 109-114.
def _zip_vector_path(zip_path: Path, suffix: str = ".shp") -> str:
    with zipfile.ZipFile(zip_path) as archive:
        candidates = [name for name in archive.namelist() if name.lower().endswith(suffix)]
    if len(candidates) != 1:
        raise RuntimeError(f"Expected one {suffix} in {zip_path.name}; found {len(candidates)}")
    return f"/vsizip/{zip_path.resolve()}/{candidates[0]}"


# Reference cell 29, lines 117-139.
def _normalize_hydrolakes(frame: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    aliases = {column.lower(): column for column in frame.columns}
    required = {"hylak_id", "lake_area", "geometry"}
    missing = sorted(required - set(aliases))
    if missing:
        raise ValueError(f"HydroLAKES file lacks columns: {missing}")
    rename = {aliases["hylak_id"]: "Hylak_id", aliases["lake_area"]: "Lake_area"}
    for optional, canonical in (
        ("lake_name", "Lake_name"),
        ("country", "Country"),
        ("lake_type", "Lake_type"),
        ("continent", "Continent"),
    ):
        if optional in aliases:
            rename[aliases[optional]] = canonical
    frame = frame.rename(columns=rename).copy()
    frame["Hylak_id"] = pd.to_numeric(frame["Hylak_id"], errors="raise").astype("int64")
    frame["Lake_area"] = pd.to_numeric(frame["Lake_area"], errors="coerce")
    frame = frame.set_crs(4326) if frame.crs is None else frame.to_crs(4326)
    frame = frame[frame.geometry.notna() & ~frame.geometry.is_empty].copy()
    try:
        frame.geometry = frame.geometry.make_valid()
    except AttributeError:
        frame.geometry = frame.geometry.buffer(0)
    return frame


# Reference cell 29, lines 142-183.
def prepare_hydrolakes(
    min_area_km2: float = 10.0, region: str = "USA", delete_global_archive: bool = True
) -> tuple[gpd.GeoDataFrame, Path]:
    """Download/cache an official HydroLAKES subset only on explicit request."""
    region = str(region).upper()
    region_slug = "USA" if region == "USA" else "WORLD"
    cache_path = (
        _m_settings.HYDRO_ROOT / f"HydroLAKES_{region_slug}_FULLPOLY_v14_ge{min_area_km2:g}km2.gpkg"
    )
    if cache_path.exists():
        return (_normalize_hydrolakes(gpd.read_file(cache_path)), cache_path)
    archive = _m_settings.HYDRO_ROOT / "HydroLAKES_polys_v10_shp.zip"
    if not archive.exists():
        print("Downloading the official global HydroLAKES archive (large; cached temporarily).")
        download_stream(_m_settings.HYDROLAKES_URL, archive)
    columns = ["Hylak_id", "Lake_name", "Country", "Continent", "Lake_type", "Lake_area"]
    lakes = pyogrio.read_dataframe(
        _zip_vector_path(archive), columns=columns, where=f"Lake_area >= {float(min_area_km2)}"
    )
    lakes = _normalize_hydrolakes(gpd.GeoDataFrame(lakes, geometry="geometry", crs=4326))
    if region == "USA":
        states_zip = _m_settings.HYDRO_ROOT / "cb_2025_us_state_5m.zip"
        if not states_zip.exists():
            download_stream(_m_settings.US_STATES_URL, states_zip)
        states = gpd.read_file(f"zip://{states_zip}").to_crs(4326)
        usa = union_all(states.geometry.values)
        lakes = lakes[lakes.intersects(usa)].copy()
        lakes = lakes[~lakes.geometry.is_empty].copy()
        lakes["boundary_scope"] = "full polygon; intersects USA"
        lakes["analysis_area_km2"] = lakes.to_crs(6933).geometry.area / 1000000.0
    else:
        lakes["boundary_scope"] = "full polygon; worldwide catalogue"
        lakes["analysis_area_km2"] = lakes.to_crs(6933).geometry.area / 1000000.0
    lakes.to_file(cache_path, driver="GPKG", layer="hydrolakes", index=False)
    if delete_global_archive:
        archive.unlink(missing_ok=True)
    print(f"Cached {len(lakes):,} lakes: {cache_path}")
    return (lakes.reset_index(drop=True), cache_path)


# Reference cell 29, lines 186-193.
def load_hydrolakes_file(
    path: str | Path, min_area_km2: float = 0.0
) -> tuple[gpd.GeoDataFrame, Path]:
    path = Path(path).expanduser().resolve()
    frame = gpd.read_file(f"zip://{path}" if path.suffix.lower() == ".zip" else path)
    frame = _normalize_hydrolakes(frame)
    frame = frame[frame["Lake_area"] >= float(min_area_km2)].copy()
    cache = _m_settings.HYDRO_ROOT / f"uploaded_HydroLAKES_ge{min_area_km2:g}km2.gpkg"
    frame.to_file(cache, driver="GPKG", layer="hydrolakes", index=False)
    return (frame.reset_index(drop=True), cache)


# Reference cell 29, lines 196-211.
def save_uploaded_vector(upload_widget) -> Path:
    value = upload_widget.value
    if not value:
        raise ValueError("Choose a GeoPackage, GeoJSON, or zipped shapefile first")
    if isinstance(value, dict):
        item = next(iter(value.values()))
        name = item.get("metadata", {}).get("name") or item.get("name") or "hydrolakes_upload"
    else:
        item = value[0]
        name = item.get("name", "hydrolakes_upload")
    content = item["content"]
    if hasattr(content, "tobytes"):
        content = content.tobytes()
    destination = _m_settings.HYDRO_ROOT / Path(name).name
    destination.write_bytes(bytes(content))
    return destination


# Reference cell 29, lines 214-219.
def scene_geolocation(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with netCDF4.Dataset(path) as root:
        if "navigation_data" in root.groups:
            nav = root.groups["navigation_data"]
            return (
                _m_science.to_float(nav.variables["latitude"][:]),
                _m_science.to_float(nav.variables["longitude"][:]),
            )
        return (
            _m_science.to_float(root.variables["latitude"][:]),
            _m_science.to_float(root.variables["longitude"][:]),
        )


# Reference cell 29, lines 222-232.
def subset_lakes_to_scene(path: Path, lakes: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    latitude, longitude = scene_geolocation(path)
    stride = max(1, int(math.sqrt(max(1, latitude.size / 350000))))
    lat = latitude[::stride, ::stride]
    lon = longitude[::stride, ::stride]
    valid = np.isfinite(lon) & np.isfinite(lat) & (np.abs(lon) <= 180) & (np.abs(lat) <= 90)
    if not valid.any():
        raise RuntimeError("Scene has no valid geolocation sample")
    footprint = box(
        float(lon[valid].min()),
        float(lat[valid].min()),
        float(lon[valid].max()),
        float(lat[valid].max()),
    )
    return lakes[lakes.intersects(footprint)].copy().reset_index(drop=True)


# Reference cell 29, lines 235-243.
def erode_lake_geometries(lakes: gpd.GeoDataFrame, distance_m: float) -> gpd.GeoDataFrame:
    result = lakes.copy()
    if distance_m <= 0 or result.empty:
        return result
    local_crs = result.estimate_utm_crs() or "EPSG:6933"
    projected = result.to_crs(local_crs)
    projected.geometry = projected.geometry.buffer(-float(distance_m))
    projected = projected[~projected.geometry.is_empty].copy()
    return projected.to_crs(4326)


# Reference cell 29, lines 246-259.
def hydrolakes_point_mask(
    latitude: np.ndarray,
    longitude: np.ndarray,
    lakes: gpd.GeoDataFrame | None,
    shore_buffer_m: float = 0.0,
) -> np.ndarray | None:
    if lakes is None or lakes.empty:
        return None
    buffered = erode_lake_geometries(lakes, float(shore_buffer_m))
    if buffered.empty:
        return np.zeros(latitude.shape, dtype=bool)
    geometry = union_all(buffered.geometry.values)
    geolocated = np.isfinite(latitude) & np.isfinite(longitude)
    result = np.zeros(latitude.shape, dtype=bool)
    result[geolocated] = contains_xy(geometry, longitude[geolocated], latitude[geolocated])
    return result


# Reference cell 29, lines 262-268.
def selected_lake_bbox(lake: gpd.GeoDataFrame | None, pad_fraction: float = 0.12):
    if lake is None or lake.empty:
        return None
    west, south, east, north = map(float, lake.total_bounds)
    dx, dy = (max(east - west, 0.01), max(north - south, 0.01))
    return (
        west - dx * pad_fraction,
        south - dy * pad_fraction,
        east + dx * pad_fraction,
        north + dy * pad_fraction,
    )
