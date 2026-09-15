"""Standalone comparison implementation."""

from __future__ import annotations

import datetime as dt
import netrc as netrc_module
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import rasterio
import requests
import xarray as xr
from pyproj import Transformer
from scipy import stats as scipy_stats
from scipy.spatial import cKDTree
from shapely import contains_xy

from . import archive as _m_archive
from . import callbacks as _m_callbacks
from . import quality as _m_quality
from . import registry as _m_registry
from . import science as _m_science
from . import settings as _m_settings
from . import state as _m_state
from . import ui_helpers as _m_ui_helpers
from . import utils as _m_utils


# Reference cell 50, lines 1044-1049.
def _scene_utc_date_global(scene_id, explicit=""):
    if str(explicit or "").strip():
        return dt.date.fromisoformat(str(explicit).strip())
    t = _m_registry.registry_table(_m_state.GLOBAL_STATE.get("config_hash"))
    r = t[t.scene_id.astype(str) == str(scene_id)]
    if r.empty:
        raise KeyError(scene_id)
    value = str(r.iloc[0].acquisition_start)
    return pd.to_datetime(value, utc=True).date()


# Reference cell 50, lines 1052-1053.
def _obpg_filename_global(acquisition_date):
    return f"L{acquisition_date.year}{acquisition_date.timetuple().tm_yday:03d}.L3m_DAY_CYAN_CI_cyano_CYAN_CONUS_300m.tif"


# Reference cell 50, lines 1056-1082.
def _download_obpg_cyano_global(acquisition_date):
    filename = _obpg_filename_global(acquisition_date)
    destination = _m_settings.PERSISTENT_ROOT / "obpg_cyan" / filename
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and destination.stat().st_size > 1024:
        return destination
    netrc_path = Path.home() / ".netrc"
    if not netrc_path.exists():
        raise RuntimeError(
            "Earthdata credentials are not configured. Run cyanolake configure-earthdata."
        )
    credentials = netrc_module.netrc(str(netrc_path)).authenticators("urs.earthdata.nasa.gov")
    if not credentials:
        raise RuntimeError("No urs.earthdata.nasa.gov entry exists in ~/.netrc")

    class EarthdataSession(requests.Session):
        AUTH_HOST = "urs.earthdata.nasa.gov"

        def rebuild_auth(self, prepared_request, response):
            headers = prepared_request.headers
            if "Authorization" not in headers:
                return
            original = urlparse(response.request.url).hostname
            redirect = urlparse(prepared_request.url).hostname
            if original != redirect and original != self.AUTH_HOST and (redirect != self.AUTH_HOST):
                del headers["Authorization"]

    url = "https://oceandata.sci.gsfc.nasa.gov/getfile/" + filename
    temporary = destination.with_suffix(".tif.part")
    session = EarthdataSession()
    session.auth = (credentials[0], credentials[2])
    with session.get(url, stream=True, timeout=(30, 180), allow_redirects=True) as response:
        response.raise_for_status()
        content_type = response.headers.get("content-type", "").lower()
        first = b""
        with temporary.open("wb") as fh:
            for chunk in response.iter_content(1024 * 1024):
                if not chunk:
                    continue
                if not first:
                    first = chunk[:16]
                fh.write(chunk)
        if "html" in content_type or first.lstrip().startswith(b"<"):
            temporary.unlink(missing_ok=True)
            raise RuntimeError("NASA returned an Earthdata sign-in page instead of GeoTIFF")
        if first[:4] not in (b"II*\x00", b"MM\x00*"):
            temporary.unlink(missing_ok=True)
            raise RuntimeError("OBPG response is not a TIFF file")
    temporary.replace(destination)
    return destination


# Reference cell 50, lines 1085-1086.
def _obpg_ci_from_dn_global(raw_dn, valid_mask):
    raw = np.asarray(raw_dn)
    science = np.asarray(valid_mask, bool) & np.isfinite(raw) & (raw >= 0) & (raw <= 250)
    ci = np.full(raw.shape, np.nan, np.float32)
    ci[science] = np.power(10.0, 3.0 / 250.0 * raw[science].astype(float) - 4.2)
    return (ci, science)


# Reference cell 50, lines 1089-1092.
def _selected_target_lake(lake_id):
    lakes = _m_callbacks._load_target_lakes_if_needed()
    selected = lakes[pd.to_numeric(lakes.Hylak_id, errors="coerce") == int(lake_id)].copy()
    if selected.empty:
        raise KeyError(f"Hylak_id {lake_id} is not in the current target lake catalogue")
    return selected.to_crs(4326)


# Reference cell 50, lines 1095-1104.
def _obpg_lake_window_global(obpg_path, lake_frame):
    with rasterio.open(obpg_path) as source:
        if source.crs is None:
            raise ValueError("OBPG GeoTIFF has no CRS")
        projected = lake_frame.to_crs(source.crs)
        west, south, east, north = projected.total_bounds
        window = (
            rasterio.windows.from_bounds(west, south, east, north, source.transform)
            .round_offsets()
            .round_lengths()
        )
        full = rasterio.windows.Window(0, 0, source.width, source.height)
        try:
            window = window.intersection(full)
        except rasterio.errors.WindowError as exc:
            raise ValueError("Selected lake does not overlap this OBPG raster") from exc
        if window.width < 1 or window.height < 1:
            raise ValueError("Selected lake does not overlap this OBPG raster")
        raw = source.read(1, window=window, masked=False)
        valid = source.read_masks(1, window=window) > 0
        transform = rasterio.windows.transform(window, source.transform)
        rows, cols = np.indices(raw.shape)
        xs, ys = rasterio.transform.xy(transform, rows, cols, offset="center")
        xs = np.asarray(xs)
        ys = np.asarray(ys)
        tr = Transformer.from_crs(source.crs, 4326, always_xy=True)
        lon, lat = tr.transform(xs, ys)
        lon = np.asarray(lon).reshape(raw.shape)
        lat = np.asarray(lat).reshape(raw.shape)
        profile = {"crs": source.crs.to_string(), "transform": tuple(transform)}
    geom = lake_frame.geometry.union_all()
    inside = contains_xy(geom, lon, lat)
    return (raw, valid, lat, lon, inside, profile)


# Reference cell 50, lines 1107-1114.
def _cyan_metrics(obpg, local):
    paired = (
        np.isfinite(obpg)
        & np.isfinite(local)
        & (obpg > _m_settings.CYAN_CI_DETECTION_LIMIT)
        & (local > _m_settings.CYAN_CI_DETECTION_LIMIT)
    )
    x = np.asarray(obpg, float)[paired]
    y = np.asarray(local, float)[paired]
    out = {
        "paired_above_detection_limit_pixels": int(len(x)),
        "pearson_r": np.nan,
        "spearman_rho": np.nan,
        "mean_bias_local_minus_obpg": np.nan,
        "mae": np.nan,
        "rmse": np.nan,
    }
    if len(x):
        d = y - x
        out.update(
            mean_bias_local_minus_obpg=float(np.mean(d)),
            mae=float(np.mean(np.abs(d))),
            rmse=float(np.sqrt(np.mean(d * d))),
        )
    if len(x) >= 2 and np.nanstd(x) > 0 and (np.nanstd(y) > 0):
        out["pearson_r"] = float(scipy_stats.pearsonr(x, y).statistic)
        out["spearman_rho"] = float(scipy_stats.spearmanr(x, y).statistic)
    return (out, paired)


# Reference cell 50, lines 1117-1128.
def prepare_cyan_source_callback(scene_id, explicit_date, source_mode, uploaded_tif):
    try:
        date = _scene_utc_date_global(scene_id, explicit_date)
        if source_mode == "upload":
            if not uploaded_tif:
                raise ValueError("Upload an OBPG daily CIcyano GeoTIFF")
            path = Path(str(uploaded_tif))
        else:
            path = _download_obpg_cyano_global(date)
        with rasterio.open(path) as ds:
            if ds.count < 1:
                raise ValueError("GeoTIFF has no raster band")
        _m_state.CYAN_STATE.update(source_path=Path(path), source_date=date.isoformat())
        return (
            str(path),
            _m_ui_helpers._app_status(
                "NASA CyAN source ready", f"{date.isoformat()} · {path.name}"
            ),
        )
    except Exception as exc:
        return (None, _m_ui_helpers._app_status("NASA CyAN source failed", str(exc), False))


# Reference cell 50, lines 1131-1160.
def compare_cyan_callback(
    scene_id, lake_id, explicit_date, source_mode, uploaded_tif, max_distance_m
):
    try:
        if not scene_id or not lake_id:
            raise ValueError("Choose a completed scene and lake")
        date = _scene_utc_date_global(scene_id, explicit_date)
        if source_mode == "upload":
            if not uploaded_tif:
                raise ValueError("Upload the NASA OBPG GeoTIFF")
            obpg_path = Path(str(uploaded_tif))
        elif (
            _m_state.CYAN_STATE.get("source_path")
            and _m_state.CYAN_STATE.get("source_date") == date.isoformat()
        ):
            obpg_path = Path(_m_state.CYAN_STATE["source_path"])
        else:
            obpg_path = _download_obpg_cyano_global(date)
        lake = _selected_target_lake(lake_id)
        raw, raster_valid, lat, lon, inside, profile = _obpg_lake_window_global(obpg_path, lake)
        obpg_ci, science = _obpg_ci_from_dn_global(raw, raster_valid)
        local = _m_archive._load_scene_pixels(scene_id, "lake", lake_id)
        source = local[
            _m_quality._cyan_comparison_selector(local)
            & np.isfinite(local.latitude)
            & np.isfinite(local.longitude)
        ].reset_index(drop=True)
        if source.empty:
            raise ValueError("No compact native pixels exist for selected lake")
        target_lat = lat[inside]
        target_lon = lon[inside]
        tree = cKDTree(_m_science._unit_sphere_xyz(source.latitude, source.longitude))
        chord, idx = tree.query(_m_science._unit_sphere_xyz(target_lat, target_lon), k=1)
        dist = 6371000.0 * (2 * np.arcsin(np.minimum(1.0, np.asarray(chord, float) / 2)))
        matched = dist <= float(max_distance_m)
        local_val = np.full(len(target_lat), np.nan, np.float32)
        local_ci_valid = np.zeros(len(target_lat), bool)
        local_candidate = np.zeros(len(target_lat), bool)
        local_detected = np.zeros(len(target_lat), bool)
        if matched.any():
            q = source.iloc[np.asarray(idx[matched], int)]
            local_val[matched] = pd.to_numeric(q.CI_cyano, errors="coerce").to_numpy(np.float32)
            local_ci_valid[matched] = (
                pd.to_numeric(q.get("ci_valid_mask", 0), errors="coerce").fillna(0).to_numpy(float)
                > 0
            )
            local_candidate[matched] = (
                pd.to_numeric(q.get("ci_candidate_mask", 0), errors="coerce")
                .fillna(0)
                .to_numpy(float)
                > 0
            )
            local_detected[matched] = (
                pd.to_numeric(q.get("ci_detection_mask", 0), errors="coerce")
                .fillna(0)
                .to_numpy(float)
                > 0
            )
        raw_inside = raw[inside]
        science_inside = science[inside]
        obpg_inside = obpg_ci[inside]
        obpg_det = (
            science_inside
            & np.isfinite(obpg_inside)
            & (obpg_inside > _m_settings.CYAN_CI_DETECTION_LIMIT)
        )
        local_det = (
            matched
            & local_detected
            & np.isfinite(local_val)
            & (local_val > _m_settings.CYAN_CI_DETECTION_LIMIT)
        )
        metrics, paired = _cyan_metrics(obpg_inside, np.where(local_det, local_val, np.nan))
        table = pd.DataFrame(
            {
                "latitude": target_lat,
                "longitude": target_lon,
                "obpg_raw_dn": raw_inside.astype(np.int16),
                "obpg_science_valid": science_inside,
                "obpg_CI_cyano": obpg_inside,
                "obpg_detected_gt_0_0001": obpg_det,
                "local_match_found": matched,
                "local_match_distance_m": np.where(matched, dist, np.nan),
                "local_ci_input_valid": local_ci_valid,
                "local_ci_candidate": local_candidate,
                "local_ci_cyano_detected": local_det,
                "local_CI_cyano": local_val,
                "paired_above_0_0001_for_correlation": paired,
            }
        )
        coverage = pd.DataFrame(
            [
                {
                    "date_utc": date.isoformat(),
                    "Hylak_id": int(lake_id),
                    "lake_name": str(lake.iloc[0].get("Lake_name", "") or ""),
                    "obpg_grid_cells_inside_lake": int(inside.sum()),
                    "obpg_science_valid_cells": int(science_inside.sum()),
                    "obpg_detection_cells_gt_0_0001": int(obpg_det.sum()),
                    "local_nearest_matches": int(matched.sum()),
                    "local_detection_matches_gt_0_0001": int(local_det.sum()),
                    "paired_detection_pixels": metrics["paired_above_detection_limit_pixels"],
                    "max_nearest_distance_m": float(max_distance_m),
                    "comparison_basis": "same-day OBPG merged 300 m cells versus nearest native OLCI pixel passing the stored CyAN-strict QA mask",
                }
            ]
        )
        metric_table = pd.DataFrame([metrics])
        d = _m_settings.EXPORT_DIR / "cyan_comparison" / str(scene_id)
        d.mkdir(parents=True, exist_ok=True)
        stem = f"{scene_id}_Hylak{int(lake_id)}_{date.isoformat()}_CyAN"
        csv = d / f"{stem}_pairs.csv"
        met = d / f"{stem}_metrics.csv"
        nc = d / f"{stem}_comparison.nc"
        table.to_csv(csv, index=False)
        pd.concat([coverage, metric_table], axis=1).to_csv(met, index=False)
        shape = raw.shape

        def grid(values, fill=np.nan, dtype=np.float32):
            out = np.full(shape, fill, dtype=dtype)
            out[inside] = np.asarray(values, dtype=dtype)
            return out

        ds = xr.Dataset(
            {
                "latitude": (("y", "x"), lat.astype(np.float32)),
                "longitude": (("y", "x"), lon.astype(np.float32)),
                "inside_lake": (("y", "x"), inside.astype(np.uint8)),
                "obpg_raw_dn": (("y", "x"), raw.astype(np.int16)),
                "obpg_CI_cyano": (("y", "x"), np.where(inside, obpg_ci, np.nan).astype(np.float32)),
                "local_CI_cyano_nearest": (("y", "x"), grid(local_val)),
                "local_match_distance_m": (("y", "x"), grid(np.where(matched, dist, np.nan))),
                "paired_detection_mask": (("y", "x"), grid(paired, 0, np.uint8)),
            },
            attrs={
                "title": "Local L2Gen CyAN-strict CIcyano versus NASA OBPG Merged-S3-CYAN comparison",
                "created_utc": _m_utils.utc_now(),
                "scene_id": str(scene_id),
                "acquisition_date_utc": date.isoformat(),
                "selected_hylak_id": int(lake_id),
                "source_obpg_geotiff": str(obpg_path),
                "pairing": "nearest compact native OLCI lake pixel to each in-lake OBPG grid cell centre",
                "max_nearest_distance_m": float(max_distance_m),
                "ci_cyano_detection_limit": float(_m_settings.CYAN_CI_DETECTION_LIMIT),
                "note": "OBPG is same-day merged S3A/S3B Level-3; local source is one OLCI scene.",
            },
        )
        ds.to_netcdf(
            nc, engine="netcdf4", encoding={n: {"zlib": True, "complevel": 4} for n in ds.data_vars}
        )
        ds.close()
        finite = table[table.paired_above_0_0001_for_correlation]
        fig = go.Figure()
        fig.add_trace(
            go.Scattergl(
                x=finite.obpg_CI_cyano,
                y=finite.local_CI_cyano,
                mode="markers",
                name="Paired detections",
                marker={"size": 5, "opacity": 0.55},
            )
        )
        if not finite.empty:
            lo = float(min(finite.obpg_CI_cyano.min(), finite.local_CI_cyano.min()))
            hi = float(max(finite.obpg_CI_cyano.max(), finite.local_CI_cyano.max()))
            fig.add_trace(
                go.Scatter(
                    x=[lo, hi],
                    y=[lo, hi],
                    mode="lines",
                    name="1:1",
                    line={"dash": "dash", "color": "black"},
                )
            )
        fig.update_layout(
            title=f"{coverage.iloc[0].lake_name} · local scene versus same-day NASA CyAN",
            xaxis={"title": "NASA OBPG CIcyano", "type": "log"},
            yaxis={"title": "Local L2Gen CIcyano", "type": "log"},
            template="plotly_white",
            height=600,
        )
        return (
            fig,
            coverage.round(6),
            metric_table.round(8),
            str(csv),
            str(met),
            str(nc),
            str(obpg_path),
            _m_ui_helpers._app_status(
                "NASA CyAN comparison complete",
                f"Paired detections >0.0001: {metrics['paired_above_detection_limit_pixels']:,}. This is same-day merged Level-3 versus one scene, not exact-overpass independent validation.",
            ),
        )
    except Exception as exc:
        return (
            None,
            pd.DataFrame(),
            pd.DataFrame(),
            None,
            None,
            None,
            None,
            _m_ui_helpers._app_status("NASA CyAN comparison failed", str(exc), False),
        )


# Reference cell 58, lines 617-629.
def _cyan_metrics_v265(obpg, local):
    paired = (
        np.isfinite(obpg)
        & np.isfinite(local)
        & (obpg > _m_settings.CYAN_CI_DETECTION_LIMIT)
        & (local > _m_settings.CYAN_CI_DETECTION_LIMIT)
    )
    x = np.asarray(obpg, float)[paired]
    y = np.asarray(local, float)[paired]
    out = {
        "paired_above_detection_limit_pixels": int(len(x)),
        "pearson_r": np.nan,
        "spearman_rho": np.nan,
        "mean_bias_local_minus_obpg": np.nan,
        "mae": np.nan,
        "rmse": np.nan,
        "linear_slope": np.nan,
        "linear_intercept": np.nan,
        "median_local_over_obpg": np.nan,
        "median_ratio_obpg_ge_0_01": np.nan,
        "mean_bias_obpg_ge_0_01": np.nan,
    }
    if len(x):
        d = y - x
        out.update(
            mean_bias_local_minus_obpg=float(np.mean(d)),
            mae=float(np.mean(np.abs(d))),
            rmse=float(np.sqrt(np.mean(d * d))),
            median_local_over_obpg=float(np.nanmedian(y / x)),
        )
        hi = x >= 0.01
        if hi.any():
            out["median_ratio_obpg_ge_0_01"] = float(np.nanmedian(y[hi] / x[hi]))
            out["mean_bias_obpg_ge_0_01"] = float(np.nanmean(y[hi] - x[hi]))
    if len(x) >= 2 and np.nanstd(x) > 0 and (np.nanstd(y) > 0):
        out["pearson_r"] = float(scipy_stats.pearsonr(x, y).statistic)
        out["spearman_rho"] = float(scipy_stats.spearmanr(x, y).statistic)
        slope, intercept = np.polyfit(x, y, 1)
        out["linear_slope"] = float(slope)
        out["linear_intercept"] = float(intercept)
    return (out, paired)


# Reference cell 60, lines 19-453.
def compare_cyan_callback_v265(
    scene_id, lake_id, explicit_date, source_mode, uploaded_tif, max_distance_m, qa_basis="primary"
):
    """Compare NASA Merged-S3-CYAN to one local OLCI scene.

    Pairing semantics are intentionally native-first:
      NASA cell centre -> nearest same-lake native OLCI pixel -> inspect QA.

    The KD-tree is built from all finite stored native observations for the selected
    Hylak_id, not from a QA-prefiltered subset. Therefore a QA-rejected bright/cloudy
    pixel is never replaced by a farther "clear" pixel.
    """
    try:
        if not scene_id or not lake_id:
            raise ValueError("Choose a completed scene and lake")
        date = _scene_utc_date_global(scene_id, explicit_date)
        if source_mode == "upload":
            if not uploaded_tif:
                raise ValueError("Upload the NASA OBPG GeoTIFF")
            obpg_path = Path(str(uploaded_tif))
        elif (
            _m_state.CYAN_STATE.get("source_path")
            and _m_state.CYAN_STATE.get("source_date") == date.isoformat()
        ):
            obpg_path = Path(_m_state.CYAN_STATE["source_path"])
        else:
            obpg_path = _download_obpg_cyano_global(date)
        lake = _selected_target_lake(lake_id)
        raw, raster_valid, lat, lon, inside, profile = _obpg_lake_window_global(obpg_path, lake)
        obpg_ci, science = _obpg_ci_from_dn_global(raw, raster_valid)
        local = _m_archive._load_scene_pixels(scene_id, "lake", lake_id).reset_index(drop=True)
        geo = np.isfinite(
            pd.to_numeric(local.latitude, errors="coerce").to_numpy(float)
        ) & np.isfinite(pd.to_numeric(local.longitude, errors="coerce").to_numpy(float))
        source = local.loc[geo].reset_index(drop=True)
        if source.empty:
            raise ValueError("No finite compact native OLCI pixels exist for the selected lake")
        target_lat = lat[inside]
        target_lon = lon[inside]
        if not len(target_lat):
            raise ValueError("No NASA CyAN cell centres fall inside the selected lake polygon")
        tree = cKDTree(_m_science._unit_sphere_xyz(source.latitude, source.longitude))
        chord, idx = tree.query(
            _m_science._unit_sphere_xyz(target_lat, target_lon), k=1, workers=-1
        )
        dist = 6371000.0 * (2 * np.arcsin(np.minimum(1.0, np.asarray(chord, float) / 2)))
        matched = np.isfinite(dist) & (dist <= float(max_distance_m))
        n = len(target_lat)
        local_val = np.full(n, np.nan, np.float32)
        local_ci_valid = np.zeros(n, bool)
        local_candidate = np.zeros(n, bool)
        local_detected_raw = np.zeros(n, bool)
        nearest_primary_valid = np.zeros(n, bool)
        nearest_strict_valid = np.zeros(n, bool)
        nearest_selected_qa_valid = np.zeros(n, bool)
        nearest_bloom_rescued = np.zeros(n, bool)
        nearest_cldice = np.zeros(n, bool)
        nearest_hisatzen = np.zeros(n, bool)
        nearest_navfail = np.zeros(n, bool)
        nearest_lat = np.full(n, np.nan, float)
        nearest_lon = np.full(n, np.nan, float)
        if matched.any():
            pos = np.flatnonzero(matched)
            q = source.iloc[np.asarray(idx[matched], dtype=int)].reset_index(drop=True)
            nearest_lat[pos] = pd.to_numeric(q.latitude, errors="coerce").to_numpy(float)
            nearest_lon[pos] = pd.to_numeric(q.longitude, errors="coerce").to_numpy(float)
            primary_valid_q = _m_quality._primary_selector(q)
            strict_valid_q = _m_quality._cyan_comparison_selector(q)
            if str(qa_basis).lower() == "strict":
                selected_valid_q = strict_valid_q
                qa_label = "generic_CLDICE_excluded_sensitivity"
            else:
                selected_valid_q = primary_valid_q
                qa_label = "primary_bloom_aware"
            nearest_primary_valid[pos] = primary_valid_q
            nearest_strict_valid[pos] = strict_valid_q
            nearest_selected_qa_valid[pos] = selected_valid_q

            def _mask_col(frame, name):
                if name not in frame.columns:
                    return np.zeros(len(frame), bool)
                return pd.to_numeric(frame[name], errors="coerce").fillna(0).to_numpy(float) > 0

            nearest_bloom_rescued[pos] = _mask_col(q, "bloom_rescue_mask")
            nearest_cldice[pos] = _mask_col(q, "cldice_mask")
            nearest_hisatzen[pos] = _mask_col(q, "hisatzen_mask")
            nearest_navfail[pos] = _mask_col(q, "navfail_mask")
            q_val = pd.to_numeric(q.CI_cyano, errors="coerce").to_numpy(np.float32)
            q_ci_valid = _mask_col(q, "ci_valid_mask")
            q_candidate = _mask_col(q, "ci_candidate_mask")
            q_detected = _mask_col(q, "ci_detection_mask")
            local_val[pos] = q_val
            local_ci_valid[pos] = q_ci_valid
            local_candidate[pos] = q_candidate
            local_detected_raw[pos] = q_detected
        else:
            qa_label = (
                "generic_CLDICE_excluded_sensitivity"
                if str(qa_basis).lower() == "strict"
                else "primary_bloom_aware"
            )
        raw_inside = raw[inside]
        science_inside = science[inside]
        obpg_inside = obpg_ci[inside]
        obpg_det = (
            science_inside
            & np.isfinite(obpg_inside)
            & (obpg_inside > _m_settings.CYAN_CI_DETECTION_LIMIT)
        )
        local_det = (
            matched
            & nearest_selected_qa_valid
            & local_detected_raw
            & np.isfinite(local_val)
            & (local_val > _m_settings.CYAN_CI_DETECTION_LIMIT)
        )
        metrics, paired = _cyan_metrics_v265(obpg_inside, np.where(local_det, local_val, np.nan))
        table = pd.DataFrame(
            {
                "latitude": target_lat,
                "longitude": target_lon,
                "obpg_raw_dn": raw_inside.astype(np.int16),
                "obpg_science_valid": science_inside,
                "obpg_CI_cyano": obpg_inside,
                "obpg_detected_gt_0_0001": obpg_det,
                "local_match_found": matched,
                "local_match_distance_m": np.where(matched, dist, np.nan),
                "nearest_native_latitude": nearest_lat,
                "nearest_native_longitude": nearest_lon,
                "nearest_primary_bloom_aware_valid": nearest_primary_valid,
                "nearest_generic_strict_valid": nearest_strict_valid,
                "nearest_selected_QA_valid": nearest_selected_qa_valid,
                "nearest_CLDICE": nearest_cldice,
                "nearest_HISATZEN": nearest_hisatzen,
                "nearest_NAVFAIL": nearest_navfail,
                "nearest_bloom_rescued": nearest_bloom_rescued,
                "local_ci_input_valid": local_ci_valid,
                "local_ci_candidate": local_candidate,
                "local_ci_detection_before_selected_QA": local_detected_raw,
                "local_ci_cyano_detected": local_det,
                "local_CI_cyano": local_val,
                "paired_above_0_0001_for_correlation": paired,
            }
        )
        coverage = pd.DataFrame(
            [
                {
                    "date_utc": date.isoformat(),
                    "Hylak_id": int(lake_id),
                    "lake_name": str(lake.iloc[0].get("Lake_name", "") or ""),
                    "obpg_grid_cells_inside_lake": int(inside.sum()),
                    "obpg_science_valid_cells": int(science_inside.sum()),
                    "obpg_detection_cells_gt_0_0001": int(obpg_det.sum()),
                    "local_nearest_same_lake_matches": int(matched.sum()),
                    "nearest_primary_valid": int((matched & nearest_primary_valid).sum()),
                    "nearest_selected_QA_valid": int((matched & nearest_selected_qa_valid).sum()),
                    "nearest_CLDICE_pixels": int((matched & nearest_cldice).sum()),
                    "nearest_bloom_rescued_pixels": int((matched & nearest_bloom_rescued).sum()),
                    "local_detection_matches_gt_0_0001": int(local_det.sum()),
                    "paired_detection_pixels": metrics["paired_above_detection_limit_pixels"],
                    "max_nearest_distance_m": float(max_distance_m),
                    "local_QA_basis": qa_label,
                    "pairing_order": "nearest same-lake native OLCI observation FIRST; selected QA inspected SECOND; no farther clear-pixel substitution",
                    "comparison_basis": "same-day NASA Merged-S3-CYAN daily maximum versus one local OLCI scene",
                }
            ]
        )
        metric_table = pd.DataFrame([metrics])
        d = _m_settings.EXPORT_DIR / "cyan_comparison" / str(scene_id)
        d.mkdir(parents=True, exist_ok=True)
        stem = f"{scene_id}_Hylak{int(lake_id)}_{date.isoformat()}_CyAN_native_first_{qa_label}"
        csv = d / f"{stem}_pairs.csv"
        met = d / f"{stem}_metrics.csv"
        nc = d / f"{stem}_comparison.nc"
        table.to_csv(csv, index=False)
        pd.concat([coverage, metric_table], axis=1).to_csv(met, index=False)
        shape = raw.shape

        def grid(values, fill=np.nan, dtype=np.float32):
            out = np.full(shape, fill, dtype=dtype)
            out[inside] = np.asarray(values, dtype=dtype)
            return out

        xds = xr.Dataset(
            {
                "latitude": (("y", "x"), lat.astype(np.float32)),
                "longitude": (("y", "x"), lon.astype(np.float32)),
                "inside_lake": (("y", "x"), inside.astype(np.uint8)),
                "obpg_raw_dn": (("y", "x"), raw.astype(np.int16)),
                "obpg_CI_cyano": (("y", "x"), np.where(inside, obpg_ci, np.nan).astype(np.float32)),
                "local_CI_cyano_nearest": (("y", "x"), grid(local_val)),
                "local_match_distance_m": (("y", "x"), grid(np.where(matched, dist, np.nan))),
                "nearest_selected_QA_valid": (
                    ("y", "x"),
                    grid(nearest_selected_qa_valid, 0, np.uint8),
                ),
                "nearest_CLDICE": (("y", "x"), grid(nearest_cldice, 0, np.uint8)),
                "nearest_bloom_rescued": (("y", "x"), grid(nearest_bloom_rescued, 0, np.uint8)),
                "paired_detection_mask": (("y", "x"), grid(paired, 0, np.uint8)),
            },
            attrs={
                "title": "Native-first local L2Gen CIcyano versus NASA OBPG Merged-S3-CYAN comparison",
                "schema": _m_settings.V268_SCHEMA,
                "created_utc": _m_utils.utc_now(),
                "scene_id": str(scene_id),
                "acquisition_date_utc": date.isoformat(),
                "selected_hylak_id": int(lake_id),
                "source_obpg_geotiff": str(obpg_path),
                "local_QA_basis": qa_label,
                "pairing": "nearest same-lake native OLCI lake pixel selected before QA; QA-rejected nearest pixel is not replaced by a farther pixel",
                "max_nearest_distance_m": float(max_distance_m),
                "ci_cyano_detection_limit": float(_m_settings.CYAN_CI_DETECTION_LIMIT),
                "note": "NASA Merged-S3-CYAN is the daily per-pixel maximum across S3A/S3B; the local source is one scene.",
            },
        )
        xds.to_netcdf(
            nc,
            engine="netcdf4",
            encoding={n: {"zlib": True, "complevel": 4} for n in xds.data_vars},
        )
        xds.close()
        finite_pairs = table[table.paired_above_0_0001_for_correlation].copy()
        fig = go.Figure()
        fig.add_trace(
            go.Scattergl(
                x=finite_pairs.obpg_CI_cyano,
                y=finite_pairs.local_CI_cyano,
                mode="markers",
                name="Paired detections",
                marker={"size": 5, "opacity": 0.48},
            )
        )
        if not finite_pairs.empty:
            lo = float(min(finite_pairs.obpg_CI_cyano.min(), finite_pairs.local_CI_cyano.min()))
            hi = float(max(finite_pairs.obpg_CI_cyano.max(), finite_pairs.local_CI_cyano.max()))
            fig.add_trace(
                go.Scatter(
                    x=[lo, hi],
                    y=[lo, hi],
                    mode="lines",
                    name="1:1",
                    line={"dash": "dash", "color": "black"},
                )
            )
            xmin = max(float(finite_pairs.obpg_CI_cyano.min()), _m_settings.CYAN_CI_DETECTION_LIMIT)
            xmax = float(finite_pairs.obpg_CI_cyano.max())
            if xmax > xmin:
                edges = np.geomspace(xmin, xmax, 9)
                medx, medy = ([], [])
                for a, b in zip(edges[:-1], edges[1:]):
                    upper = b if b < edges[-1] else b * 1.000001
                    qbin = finite_pairs[
                        (finite_pairs.obpg_CI_cyano >= a) & (finite_pairs.obpg_CI_cyano < upper)
                    ]
                    if len(qbin) >= 5:
                        medx.append(float(np.median(qbin.obpg_CI_cyano)))
                        medy.append(float(np.median(qbin.local_CI_cyano)))
                if medx:
                    fig.add_trace(
                        go.Scatter(
                            x=medx,
                            y=medy,
                            mode="lines+markers",
                            name="Binned median local",
                            line={"width": 3},
                        )
                    )
        fig.update_layout(
            title=f"{coverage.iloc[0].lake_name} · native-first local scene versus same-day NASA CyAN",
            xaxis={"title": "NASA OBPG Merged-S3-CYAN CIcyano", "type": "log"},
            yaxis={"title": "Local L2Gen CIcyano", "type": "log"},
            template="plotly_white",
            height=600,
        )
        note = f"Paired detections >0.0001: {metrics['paired_above_detection_limit_pixels']:,}. Nearest same-lake native OLCI observation is selected before QA; a rejected nearest pixel is never replaced by a farther clear pixel. Default comparison QA is the bloom-aware science mask. NASA Merged-S3-CYAN is a same-day S3A/S3B daily maximum, while the local product is one scene."
        return (
            fig,
            coverage.round(6),
            metric_table.round(8),
            str(csv),
            str(met),
            str(nc),
            str(obpg_path),
            _m_ui_helpers._app_status("NASA CyAN comparison complete", note),
        )
    except Exception as exc:
        return (
            None,
            pd.DataFrame(),
            pd.DataFrame(),
            None,
            None,
            None,
            None,
            _m_ui_helpers._app_status("NASA CyAN comparison failed", str(exc), False),
        )
