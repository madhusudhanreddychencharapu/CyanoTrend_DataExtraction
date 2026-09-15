"""Standalone archive implementation."""

from __future__ import annotations

import contextlib
import json
import math
import zipfile
from pathlib import Path

import geopandas as gpd
import gradio as gr
import netCDF4
import numpy as np
import pandas as pd
from shapely.geometry import box as shapely_box

from . import catalogue as _m_catalogue
from . import registry as _m_registry
from . import science as _m_science
from . import settings as _m_settings
from . import state as _m_state
from . import statistics as _m_statistics
from . import ui_helpers as _m_ui_helpers
from . import workspace as _m_workspace


# Reference cell 50, lines 346-350.
def _scene_choices():
    t = _m_registry.registry_table(_m_state.GLOBAL_STATE.get("config_hash"))
    if t.empty:
        return []
    t = t[t.status == "done"].sort_values("acquisition_start")
    return [
        (f"{r.get('acquisition_start', '')} · {r.get('scene_name', '')}", str(r.scene_id))
        for _, r in t.iterrows()
    ]


# Reference cell 58, lines 39-50.
def refresh_scene_choices_callback():
    c = _scene_choices()
    v = c[0][1] if c else None
    n_lakes = 0
    if v:
        with contextlib.suppress(Exception):
            n_lakes = int(_load_scene_pixels(v, "scene", None).Hylak_id.nunique())
    return (
        gr.Dropdown(choices=c, value=v),
        _scene_lake_update(v),
        _m_ui_helpers._app_status(
            "Scene list refreshed",
            f"{len(c):,} completed scenes; selected scene contains {n_lakes:,} stored lakes.",
        ),
    )


# Reference cell 58, lines 16-32.
def _scene_lake_update(scene_id):
    """Dropdown derived from the actual compact scene pixels, never the global HydroLAKES table."""
    if not scene_id:
        return gr.Dropdown(choices=[], value=None)
    try:
        f = _load_scene_pixels(scene_id, "scene", None)
    except Exception:
        return gr.Dropdown(choices=[], value=None)
    ids = np.sort(pd.to_numeric(f["Hylak_id"], errors="coerce").dropna().astype(np.int64).unique())
    name_map = {}
    lakes = _m_state.GLOBAL_STATE.get("lakes")
    if lakes is not None and (not lakes.empty) and ("Hylak_id" in lakes.columns):
        subset = lakes[lakes.Hylak_id.astype(np.int64).isin(ids)]
        for _, r in subset.iterrows():
            name_map[int(r.Hylak_id)] = str(r.get("Lake_name", "") or "Unnamed lake")
    choices = [
        (f"{name_map.get(int(lid), 'Unnamed lake')} · Hylak_id {int(lid)}", str(int(lid)))
        for lid in ids
    ]
    return gr.Dropdown(choices=choices, value=choices[0][1] if choices else None)


# Reference cell 58, lines 35-36.
def scene_changed_callback(scene_id):
    return _scene_lake_update(scene_id)


# Reference cell 50, lines 371-380.
def _scene_compact_files(scene_id):
    s = _m_statistics.all_stats(_m_state.GLOBAL_STATE.get("config_hash"))
    t = s[s.scene_id.astype(str) == str(scene_id)] if not s.empty else pd.DataFrame()
    files = []
    for p in t.get("compact_netcdf", pd.Series(dtype=str)).dropna().astype(str).unique():
        q = Path(p)
        if q.exists():
            files.append(q)
    if not files:
        cfg = _m_state.GLOBAL_STATE.get("config_hash")
        if not cfg:
            raise ValueError("Select a processing configuration before loading scene pixels")
        files = sorted(_m_settings.COMPACT_NC_DIR.rglob(f"{scene_id}_*_{cfg}_lakepixels.nc"))
    return files


# Reference cell 57, lines 326-352.
def _load_scene_pixels(scene_id, scope="scene", lake_id=None):
    frames = []
    wanted = None if scope == "scene" else int(lake_id)
    names = [
        "latitude",
        "longitude",
        "Hylak_id",
        "l2_flags",
        "valid_water_mask",
        "cyan_strict_valid_mask",
        "bloom_rescue_mask",
        "cldice_mask",
        "hisatzen_mask",
        "navfail_mask",
        "cloud_excluded_mask",
        "ci_valid_mask",
        "ci_candidate_mask",
        "ci_detection_mask",
        "ndci_valid_mask",
        "mph_valid_mask",
        "fai_valid_mask",
        "MPH_peak_nm",
        *_m_settings.SCALAR_METRICS,
    ]
    for p in _scene_compact_files(scene_id):
        with netCDF4.Dataset(p) as ds:
            ids = np.asarray(ds.variables["Hylak_id"][:], dtype=np.int64)
            sel = np.ones(len(ids), bool) if wanted is None else ids == wanted
            if not sel.any():
                continue
            data = {"Hylak_id": ids[sel]}
            for n in names:
                if n == "Hylak_id":
                    continue
                if n == "l2_flags" and n in ds.variables:
                    data[n] = np.asarray(ds.variables[n][:], dtype=np.uint32)[sel]
                elif n in ds.variables:
                    data[n] = _m_science.to_float(ds.variables[n][:])[sel]
                else:
                    data[n] = np.full(sel.sum(), np.nan)
            frames.append(pd.DataFrame(data))
    if not frames:
        raise ValueError("No compact lake pixels found for selected scope")
    f = pd.concat(frames, ignore_index=True)
    f["_latr"] = f.latitude.round(6)
    f["_lonr"] = f.longitude.round(6)
    f = f.drop_duplicates(["Hylak_id", "_latr", "_lonr"]).drop(columns=["_latr", "_lonr"])
    return f


# Reference cell 50, lines 531-534.
def _lake_boundaries_for_ids(ids):
    try:
        lakes = _m_workspace.load_lakes()[0]
        ids = set(map(int, ids))
        return lakes[lakes.Hylak_id.astype(int).isin(ids)].copy()
    except Exception:
        return gpd.GeoDataFrame()


# Reference cell 50, lines 842-849.
def _compact_download(scene_id):
    files = _scene_compact_files(scene_id)
    if not files:
        return None
    if len(files) == 1:
        return str(files[0])
    d = _m_settings.EXPORT_DIR / "scene_bundles"
    d.mkdir(parents=True, exist_ok=True)
    z = d / f"{scene_id}_compact_windows.zip"
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as arc:
        for p in files:
            arc.write(p, arcname=p.name)
    return str(z)


# Reference cell 50, lines 864-865.
def _grid_label(west, south, size=5):
    return f"G{int(size):02d}_{('N' if south >= 0 else 'S')}{abs(int(south)):02d}_{('E' if west >= 0 else 'W')}{abs(int(west)):03d}"


# Reference cell 50, lines 867-871.
def grid_id_to_bbox(grid_id):
    m = _m_settings.GRID_ID_RE.match(str(grid_id or "").strip())
    if not m:
        raise ValueError("Grid ID example: G05_N40_W085")
    size = int(m.group("size"))
    south = int(m.group("lat")) * (1 if m.group("ns").upper() == "N" else -1)
    west = int(m.group("lon")) * (1 if m.group("ew").upper() == "E" else -1)
    return (float(west), float(south), float(west + size), float(south + size))


# Reference cell 50, lines 873-885.
def scene_grid_ids(scene, size=5):
    geom = _m_catalogue.scene_geometry(scene)
    if geom is None or geom.is_empty:
        return []
    minx, miny, maxx, maxy = geom.bounds
    size = int(size)
    ids = []
    w0 = math.floor(minx / size) * size
    s0 = math.floor(miny / size) * size
    w = w0
    while w < maxx:
        s = s0
        while s < maxy:
            if geom.intersects(shapely_box(w, s, w + size, s + size)):
                ids.append(_grid_label(w, s, size))
            s += size
        w += size
    return sorted(set(ids))


# Reference cell 50, lines 887-908.
def build_archive_index(config_hash=None):
    t = _m_registry.registry_table(config_hash)
    rows = []
    if t.empty:
        pd.DataFrame().to_csv(_m_settings.ARCHIVE_INDEX_CSV, index=False)
        return pd.DataFrame()
    for _, r in t[t.status == "done"].iterrows():
        try:
            scene = json.loads(r.scene_json) if isinstance(r.scene_json, str) else {}
        except Exception:
            scene = {}
        try:
            meta = _m_registry.parse_s3_frame_metadata(r.scene_name)
        except Exception:
            meta = {
                "platform": str(r.scene_name)[:3],
                "relative_orbit": np.nan,
                "frame": np.nan,
                "native_frame_id": "",
            }
        cfg = str(r.config_hash)
        stats_p = _m_settings.STATS_SCENE_DIR / f"{r.scene_id}_{cfg}.parquet"
        stats_c = stats_p.with_suffix(".csv")
        manifest = _m_settings.EXPORT_DIR / str(r.scene_id) / f"manifest_{cfg}.json"
        comp = _scene_compact_files(str(r.scene_id))
        rows.append(
            {
                "scene_id": str(r.scene_id),
                "scene_name": r.scene_name,
                "acquisition_start": r.acquisition_start,
                "platform": meta.get("platform"),
                "relative_orbit": meta.get("relative_orbit"),
                "frame": meta.get("frame"),
                "native_frame_id": meta.get("native_frame_id"),
                "grid_ids": ";".join(scene_grid_ids(scene, 5)),
                "config_hash": cfg,
                "candidate_lakes": r.candidate_lakes,
                "stats_rows": r.stats_rows,
                "compact_MB": sum((p.stat().st_size for p in comp)) / 1024**2 if comp else 0.0,
                "compact_files": "|".join(map(str, comp)),
                "stats_csv": str(stats_c) if stats_c.exists() else "",
                "stats_parquet": str(stats_p) if stats_p.exists() else "",
                "manifest": str(manifest) if manifest.exists() else "",
            }
        )
    out = (
        pd.DataFrame(rows).sort_values(["acquisition_start", "scene_name"])
        if rows
        else pd.DataFrame()
    )
    _m_settings.ARCHIVE_INDEX_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(_m_settings.ARCHIVE_INDEX_CSV, index=False)
    return out


# Reference cell 50, lines 910-917.
def archive_refresh_callback(grid_filter="", orbit_filter=""):
    try:
        t = build_archive_index(_m_state.GLOBAL_STATE.get("config_hash"))
        if str(grid_filter or "").strip():
            t = t[
                t.grid_ids.astype(str).str.contains(str(grid_filter).strip(), case=False, na=False)
            ]
        if str(orbit_filter or "").strip():
            orbit = int(float(orbit_filter))
            t = t[pd.to_numeric(t.relative_orbit, errors="coerce") == orbit]
        return (
            t.head(1000),
            str(_m_settings.ARCHIVE_INDEX_CSV),
            _m_ui_helpers._app_status(
                "Archive index ready",
                f"{len(t):,} completed scene records. Grid IDs are search cross-references; scene ID remains the deduplication key.",
            ),
        )
    except Exception as exc:
        return (
            pd.DataFrame(),
            None,
            _m_ui_helpers._app_status("Archive index failed", str(exc), False),
        )


# Reference cell 57, lines 417-433.
def _compact_qa_metadata(scene_id):
    files = _scene_compact_files(scene_id)
    meta = {
        "primary_qa_profile": "unknown",
        "primary_hard_exclude_flags": [],
        "cyan_strict_exclude_flags": [],
        "conditional_cldice_bloom_recovery": False,
        "require_l2_nonland": False,
    }
    if not files:
        return meta
    with netCDF4.Dataset(files[0]) as ds:
        meta.update(
            {
                "primary_qa_profile": str(
                    getattr(ds, "primary_qa_profile", getattr(ds, "schema", "unknown"))
                ),
                "primary_hard_exclude_flags": [
                    x for x in str(getattr(ds, "primary_hard_exclude_flags", "")).split(",") if x
                ],
                "cyan_strict_exclude_flags": [
                    x for x in str(getattr(ds, "cyan_strict_exclude_flags", "")).split(",") if x
                ],
                "conditional_cldice_bloom_recovery": bool(
                    int(getattr(ds, "conditional_cldice_bloom_recovery", 0))
                ),
                "bloom_rescue_definition": str(getattr(ds, "bloom_rescue_definition", "")),
                "bloom_rescue_method_note": str(getattr(ds, "bloom_rescue_method_note", "")),
                "require_l2_nonland": bool(int(getattr(ds, "require_l2_nonland", 1))),
                "shore_buffer_m": float(getattr(ds, "shore_buffer_m", 0.0)),
                "ci_detection_limit": float(
                    getattr(ds, "ci_detection_limit", _m_settings.CYAN_CI_DETECTION_LIMIT)
                ),
            }
        )
    return meta
