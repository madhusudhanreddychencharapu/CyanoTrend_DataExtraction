"""Standalone callbacks implementation."""

from __future__ import annotations

import geopandas as gpd
import gradio as gr
import pandas as pd
import plotly.express as px
from shapely.geometry import box as shapely_box
from tqdm.auto import tqdm

from . import admin as _m_admin
from . import archive as _m_archive
from . import catalogue as _m_catalogue
from . import hydrolakes as _m_hydrolakes
from . import pipeline as _m_pipeline
from . import registry as _m_registry
from . import settings as _m_settings
from . import state as _m_state
from . import statistics as _m_statistics
from . import ui_helpers as _m_ui_helpers
from . import workspace as workspace


# Reference cell 50, lines 178-217.
def _plan_catalogue_shard(
    start_date,
    end_date,
    all_lakes,
    min_area_km2,
    config_region,
    bbox,
    shard_geometry,
    shard_label,
    max_products=None,
):
    """Plan one discovery shard while preserving scene-level global deduplication.

    Only scenes containing target lakes inside the selected shard are queued. Processing
    still receives the canonical loaded HydroLAKES universe, so the first time a frame is
    processed it can retain all target lakes in that frame. Adjacent states therefore reuse
    the same completed Sentinel-3 scene instead of repeating L2Gen.
    """
    if all_lakes is None or all_lakes.empty:
        raise ValueError("Load HydroLAKES before planning scenes")
    candidate_pool = (
        _m_admin._lakes_intersecting_geometry(all_lakes, shard_geometry)
        if shard_geometry is not None
        else all_lakes
    )
    if candidate_pool.empty:
        raise ValueError(f"No target HydroLAKES polygons intersect planning region: {shard_label}")
    ch = _m_registry.processing_config_hash(float(min_area_km2), str(config_region))
    catalog = _m_catalogue.search_olci_l1_catalog(
        start_date, end_date, bbox=bbox, max_products=max_products
    )
    planned_ids = []
    for scene in tqdm(catalog, desc=f"Planning {shard_label}"):
        sg = _m_catalogue.scene_geometry(scene)
        if shard_geometry is not None and (sg is None or not sg.intersects(shard_geometry)):
            continue
        subset = _m_catalogue.candidate_lakes_for_scene(scene, candidate_pool)
        if subset.empty:
            continue
        rec = dict(scene)
        rec["candidate_lakes"] = int(len(subset))
        rec["planning_shard"] = str(shard_label)
        _m_registry.registry_upsert_scene(rec, ch, len(subset))
        planned_ids.append(str(rec["id"]))
    table = _m_registry.registry_table(ch)
    if planned_ids:
        table = table[table.scene_id.astype(str).isin(set(planned_ids))].copy()
    else:
        table = table.iloc[0:0].copy()
    plan_dir = _m_settings.CATALOG_DIR / "plans" / ch
    plan_dir.mkdir(parents=True, exist_ok=True)
    stamp = f"{pd.Timestamp(start_date):%Y%m%d}_{pd.Timestamp(end_date):%Y%m%d}"
    stamp += "_" + pd.Timestamp.now(tz="UTC").strftime("%Y%m%dT%H%M%S%f")
    plan_csv = plan_dir / f"plan_{_m_admin._safe_admin_slug(shard_label)}_{stamp}.csv"
    table.to_csv(plan_csv, index=False)
    return (table, ch, planned_ids, plan_csv, int(len(candidate_pool)))


# Reference cell 50, lines 228-232.
def _load_target_lakes_if_needed():
    if _m_state.GLOBAL_STATE["lakes"] is not None:
        return _m_state.GLOBAL_STATE["lakes"]
    if _m_settings.TARGET_LAKES_GPKG.exists():
        _m_state.GLOBAL_STATE["lakes"] = _m_hydrolakes._normalize_hydrolakes(
            gpd.read_file(_m_settings.TARGET_LAKES_GPKG)
        )
        return _m_state.GLOBAL_STATE["lakes"]
    raise RuntimeError("Load/prepare HydroLAKES first")


# Reference cell 50, lines 235-252.
def load_lakes_callback(region, min_area, existing_path):
    try:
        min_area = float(min_area)
        existing_path = str(existing_path or "").strip()
        if existing_path:
            lakes, path = _m_hydrolakes.load_hydrolakes_file(existing_path, min_area)
        else:
            base_region = "USA" if str(region).upper() == "USA" else "WORLD"
            lakes, path = _m_hydrolakes.prepare_hydrolakes(
                min_area, base_region, delete_global_archive=True
            )
            if str(region).upper() not in {"WORLD", "USA"}:
                if "Continent" not in lakes.columns:
                    raise ValueError("HydroLAKES subset lacks Continent column")
                target = str(region).replace("_", " ").lower()
                lakes = lakes[lakes["Continent"].astype(str).str.lower() == target].copy()
        _m_catalogue.save_target_lakes(lakes)
        _m_state.GLOBAL_STATE.update(
            lakes=lakes,
            min_area_km2=min_area,
            region=str(region).upper(),
            planned_scene_ids=[],
            planning_shard="",
            planning_geometry=None,
        )
        workspace.select_configuration(lakes)
        preview = lakes[
            [
                c
                for c in ["Hylak_id", "Lake_name", "Lake_area", "Country", "Continent"]
                if c in lakes.columns
            ]
        ].head(50)
        note = (
            " WORLD is recommended for country/state-by-state global production so one canonical lake universe is used for deduplication."
            if str(region).upper() != "WORLD"
            else ""
        )
        return (
            preview,
            _m_ui_helpers._app_status(
                "HydroLAKES ready", f"{len(lakes):,} target lakes cached.{note}"
            ),
        )
    except Exception as exc:
        return (pd.DataFrame(), _m_ui_helpers._app_status("HydroLAKES failed", str(exc), False))


# Reference cell 50, lines 255-261.
def _queue_view(table):
    if table is None or table.empty:
        return pd.DataFrame()
    v = table.copy()
    if "size_bytes" in v:
        v["L1_size_GB"] = pd.to_numeric(v["size_bytes"], errors="coerce") / 1024**3
    if "compact_bytes" in v:
        v["compact_MB"] = pd.to_numeric(v["compact_bytes"], errors="coerce") / 1024**2
    cols = [
        "scene_name",
        "acquisition_start",
        "candidate_lakes",
        "L1_size_GB",
        "status",
        "attempts",
        "processing_strategy",
        "window_count",
        "download_method",
        "compact_MB",
        "stats_rows",
        "last_error",
    ]
    return v[[c for c in cols if c in v.columns]]


# Reference cell 50, lines 264-294.
def plan_callback(
    start_date, end_date, planning_mode, country_iso3, admin_key, bbox_text, grid_id, max_products
):
    try:
        lakes = _load_target_lakes_if_needed()
        mode = str(planning_mode or "adm1")
        max_products = None if not max_products or int(max_products) <= 0 else int(max_products)
        geom = None
        bbox = None
        shard = "WORLD"
        if mode == "adm1":
            geom, name = _m_admin._selected_adm1_geometry(country_iso3, admin_key)
            bbox = tuple(map(float, geom.bounds))
            country_name = (_m_state.ADM1_STATE.get("metadata") or {}).get(
                "boundaryName", str(country_iso3)
            )
            shard = f"ADM1 · {country_name} · {name}"
        elif mode == "bbox":
            bbox = _m_ui_helpers._parse_bbox(bbox_text)
            if bbox is None:
                raise ValueError("Enter bbox west,south,east,north")
            geom = shapely_box(*bbox)
            shard = f"BBox · {','.join((f'{x:.4f}' for x in bbox))}"
        elif mode == "grid":
            bbox = _m_archive.grid_id_to_bbox(grid_id)
            geom = shapely_box(*bbox)
            shard = f"Grid · {str(grid_id).strip().upper()}"
        else:
            raise ValueError(f"Unknown planning mode: {mode}")
        table, ch, ids, plan_csv, n_lakes = _plan_catalogue_shard(
            start_date,
            end_date,
            lakes,
            _m_state.GLOBAL_STATE["min_area_km2"],
            _m_state.GLOBAL_STATE["region"],
            bbox,
            geom,
            shard,
            max_products=max_products,
        )
        _m_state.GLOBAL_STATE.update(
            config_hash=ch, planned_scene_ids=ids, planning_shard=shard, planning_geometry=geom
        )
        total = (
            float(pd.to_numeric(table.get("size_bytes", 0), errors="coerce").fillna(0).sum())
            / 1024**3
            if len(table)
            else 0
        )
        done = int((table.status == "done").sum()) if len(table) else 0
        detail = f"{len(table):,} scene/config records for {shard}; {n_lakes:,} target lakes in shard; {done:,} already complete; L1 catalogue volume ≈ {total:,.1f} GiB; config={ch}. Plan snapshot: {plan_csv.name}"
        return (
            _queue_view(table),
            str(plan_csv),
            _m_ui_helpers._app_status("Scene queue planned", detail),
        )
    except Exception as exc:
        return (pd.DataFrame(), None, _m_ui_helpers._app_status("Planning failed", str(exc), False))


# Reference cell 50, lines 297-303.
def _current_plan_table():
    ch = _m_state.GLOBAL_STATE.get("config_hash")
    if not ch:
        return pd.DataFrame()
    t = _m_registry.registry_table(ch)
    ids = {str(x) for x in _m_state.GLOBAL_STATE.get("planned_scene_ids") or []}
    if _m_state.GLOBAL_STATE.get("planning_shard"):
        t = t[t.scene_id.astype(str).isin(ids)].copy()
    return t


# Reference cell 50, lines 306-309.
def refresh_queue_callback():
    try:
        t = _current_plan_table()
        return (_queue_view(t), str(_m_registry.registry_export_csv()))
    except Exception:
        return (pd.DataFrame(), None)


# Reference cell 50, lines 312-321.
def run_batch_callback(username, password, n_scenes, method, retry_failed):
    try:
        if not str(username or "").strip() or not str(password or ""):
            raise ValueError("Enter Copernicus Data Space username and password")
        lakes = _load_target_lakes_if_needed()
        ch = _m_state.GLOBAL_STATE.get("config_hash")
        if not ch:
            raise RuntimeError("Plan the scene queue first")
        from .runtime import preflight

        report = preflight()
        if not report["ok"]:
            raise RuntimeError("; ".join(report["errors"]))
        before = _current_plan_table()
        todo = (
            before.status.isin(["queued", "running"] + (["failed"] if retry_failed else [])).sum()
            if not before.empty
            else 0
        )
        if todo == 0:
            return (
                _queue_view(before),
                _m_ui_helpers._app_status(
                    "Nothing queued",
                    f"Current shard complete: {_m_state.GLOBAL_STATE.get('planning_shard', '')}",
                ),
                _m_pipeline.storage_summary(),
            )
        _m_pipeline.process_next_queued(
            ch,
            lakes,
            username,
            password,
            n_scenes=int(n_scenes),
            download_method=str(method).lower(),
            retry_failed=bool(retry_failed),
            scene_ids=_m_state.GLOBAL_STATE.get("planned_scene_ids"),
        )
        after = _current_plan_table()
        return (
            _queue_view(after),
            _m_ui_helpers._app_status(
                "Batch finished",
                f"{_m_state.GLOBAL_STATE.get('planning_shard', '')} · queue status: {after.status.value_counts().to_dict()}",
            ),
            _m_pipeline.storage_summary(),
        )
    except Exception as exc:
        return (
            pd.DataFrame(),
            _m_ui_helpers._app_status("Batch failed", str(exc), False),
            _m_pipeline.storage_summary(),
        )


# Reference cell 50, lines 324-324.
def _stats():
    return _m_statistics.all_stats(_m_state.GLOBAL_STATE.get("config_hash"))


# Reference cell 50, lines 326-330.
def _lake_choices():
    s = _stats()
    if s.empty:
        return []
    cols = ["Hylak_id"] + (["Lake_name"] if "Lake_name" in s.columns else [])
    u = s[cols].drop_duplicates().sort_values("Hylak_id")
    return [
        (
            f"{str(r.get('Lake_name', '') or 'Unnamed lake')} · Hylak_id {int(r.Hylak_id)}",
            str(int(r.Hylak_id)),
        )
        for _, r in u.iterrows()
    ]


# Reference cell 50, lines 333-334.
def refresh_lakes_callback():
    c = _lake_choices()
    return (
        gr.Dropdown(choices=c, value=c[0][1] if c else None),
        _m_ui_helpers._app_status(
            "Lake list refreshed", f"{len(c):,} lakes have processed observations"
        ),
    )


# Reference cell 50, lines 337-343.
def lake_timeseries_callback(lake_id, metric):
    try:
        s = _stats()
        lake_id = int(lake_id)
        t = s[s.Hylak_id == lake_id].copy()
        t["time"] = pd.to_datetime(t["acquisition_start"], errors="coerce")
        t = t.sort_values("time")
        if metric not in t:
            raise KeyError(metric)
        fig = px.scatter(
            t,
            x="time",
            y=metric,
            hover_data=[
                c
                for c in ["scene_name", "window_id", "retained_pixels", "ci_input_valid_pixels"]
                if c in t
            ],
            title=f"Hylak_id {lake_id}: {metric}",
        )
        fig.update_traces(mode="lines+markers")
        return (fig, t)
    except Exception:
        return (None, pd.DataFrame())


# Reference cell 50, lines 1163-1163.
def storage_callback():
    return (_m_pipeline.storage_summary(), str(_m_registry.registry_export_csv()))


# Reference cell 50, lines 1164-1166.
def clear_scratch_callback():
    try:
        return (
            _m_pipeline.storage_summary(),
            _m_ui_helpers._app_status("Scratch cleared", _m_pipeline.clear_scratch()),
        )
    except Exception as exc:
        return (
            _m_pipeline.storage_summary(),
            _m_ui_helpers._app_status("Scratch cleanup failed", str(exc), False),
        )
