"""Standalone admin implementation."""

from __future__ import annotations

import contextlib
import json
import re

import geopandas as gpd
import gradio as gr
import numpy as np
import pandas as pd
import requests

from . import settings as _m_settings
from . import state as _m_state
from . import ui_helpers as _m_ui_helpers


# Reference cell 50, lines 65-67.
def _safe_admin_slug(text):
    text = re.sub("[^A-Za-z0-9._-]+", "_", str(text or "").strip())
    return text.strip("_") or "region"


# Reference cell 50, lines 70-78.
def _adm1_name_column(gdf):
    for name in ("shapeName", "NAME_1", "name", "Name", "NAME", "admin1Name", "ADM1_EN"):
        if name in gdf.columns:
            return name
    for name in gdf.columns:
        if name != gdf.geometry.name and gdf[name].dtype == object:
            return name
    raise ValueError("Could not identify an ADM1 name column in the boundary file")


# Reference cell 50, lines 81-127.
def _load_geoboundaries_adm1(iso3):
    iso3 = str(iso3 or "").upper().strip()
    if len(iso3) != 3:
        raise ValueError("Choose a country first")
    _m_settings.ADMIN_ROOT.mkdir(parents=True, exist_ok=True)
    meta_path = _m_settings.ADMIN_ROOT / f"geoBoundaries_{iso3}_ADM1_metadata.json"
    geo_path = _m_settings.ADMIN_ROOT / f"geoBoundaries_{iso3}_ADM1_simplified.geojson"
    metadata = None
    if meta_path.exists() and geo_path.exists():
        with contextlib.suppress(Exception):
            metadata = json.loads(meta_path.read_text())
    if metadata is None:
        url = _m_settings.GEOB_API.format(iso3=iso3)
        r = requests.get(url, timeout=90)
        r.raise_for_status()
        metadata = r.json()
        gj = metadata.get("simplifiedGeometryGeoJSON") or metadata.get("gjDownloadURL")
        if not gj:
            raise RuntimeError(f"geoBoundaries returned no ADM1 GeoJSON URL for {iso3}")
        rr = requests.get(gj, timeout=180)
        rr.raise_for_status()
        geo_path.write_bytes(rr.content)
        meta_path.write_text(json.dumps(metadata, indent=2))
    gdf = gpd.read_file(geo_path)
    if gdf.crs is None:
        gdf = gdf.set_crs(4326)
    else:
        gdf = gdf.to_crs(4326)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy().reset_index(drop=True)
    if gdf.empty:
        raise RuntimeError(f"ADM1 boundary file for {iso3} contains no usable geometries")
    name_col = _adm1_name_column(gdf)
    lookup = {}
    choices = []
    for idx, row in gdf.iterrows():
        name = str(row.get(name_col, f"ADM1 {idx + 1}") or f"ADM1 {idx + 1}")
        code = ""
        for candidate in ("shapeISO", "HASC_1", "iso_3166_2", "ISO_3166_2"):
            if candidate in gdf.columns and pd.notna(row.get(candidate)):
                code = str(row.get(candidate))
                break
        shape_id = str(row.get("shapeID", idx)) if "shapeID" in gdf.columns else str(idx)
        key = f"{idx}:{shape_id}"
        label = f"{name} ({code})" if code and code.lower() != "nan" else name
        lookup[key] = int(idx)
        choices.append((label, key))
    _m_state.ADM1_STATE.update(country_iso3=iso3, gdf=gdf, lookup=lookup, metadata=metadata)
    return (choices, metadata)


# Reference cell 50, lines 130-143.
def load_adm1_callback(country_iso3):
    try:
        choices, meta = _load_geoboundaries_adm1(country_iso3)
        value = choices[0][1] if choices else None
        detail = f"{len(choices)} ADM1 units loaded for {meta.get('boundaryName', country_iso3)}; year={meta.get('boundaryYearRepresented', 'unknown')}; source={meta.get('boundarySource', 'geoBoundaries')}; license={meta.get('boundaryLicense', 'see geoBoundaries metadata')}."
        return (
            gr.Dropdown(choices=choices, value=value),
            _m_ui_helpers._app_status("State/province list ready", detail),
        )
    except Exception as exc:
        _m_state.ADM1_STATE.update(country_iso3=None, gdf=None, lookup={}, metadata={})
        return (
            gr.Dropdown(choices=[], value=None),
            _m_ui_helpers._app_status("ADM1 loading failed", str(exc), False),
        )


# Reference cell 50, lines 146-162.
def _selected_adm1_geometry(country_iso3, admin_key):
    iso3 = str(country_iso3 or "").upper().strip()
    if _m_state.ADM1_STATE.get("gdf") is None or _m_state.ADM1_STATE.get("country_iso3") != iso3:
        _load_geoboundaries_adm1(iso3)
    key = str(admin_key or "")
    if key not in _m_state.ADM1_STATE.get("lookup", {}):
        raise ValueError("Choose a state/province after loading the country ADM1 list")
    idx = _m_state.ADM1_STATE["lookup"][key]
    row = _m_state.ADM1_STATE["gdf"].iloc[idx]
    name_col = _adm1_name_column(_m_state.ADM1_STATE["gdf"])
    name = str(row.get(name_col, key))
    geom = row.geometry
    if geom is None or geom.is_empty:
        raise ValueError("Selected ADM1 geometry is empty")
    if not geom.is_valid:
        geom = geom.buffer(0)
    return (geom, name)


# Reference cell 50, lines 165-175.
def _lakes_intersecting_geometry(lakes, geom):
    if lakes is None or lakes.empty:
        return lakes.iloc[0:0].copy()
    try:
        idx = lakes.sindex.query(geom, predicate="intersects")
        out = lakes.iloc[np.asarray(idx, dtype=int)].copy()
    except Exception:
        out = lakes[lakes.intersects(geom)].copy()
    if out.empty:
        return out
    return out[out.intersects(geom)].copy().reset_index(drop=True)
