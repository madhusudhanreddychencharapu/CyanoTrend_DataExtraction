"""Standalone statistics implementation."""

from __future__ import annotations

import contextlib
from pathlib import Path

import geopandas as gpd
import netCDF4
import numpy as np
import pandas as pd

from . import science as _m_science
from . import settings as _m_settings


# Reference cell 57, lines 288-323.
def compact_stats(nc_path: Path, lake_metadata: gpd.GeoDataFrame | None = None) -> pd.DataFrame:
    """Primary bloom-aware lake stats plus explicit CyAN-strict/rescue QA counts."""
    with netCDF4.Dataset(nc_path) as ds:
        ids = np.asarray(ds.variables["Hylak_id"][:], dtype=np.int64)
        n = len(ids)
        data = {
            "Hylak_id": ids,
            "valid": _mask_values_from_ds(ds, "valid_water_mask", n),
            "strict": _mask_values_from_ds(ds, "cyan_strict_valid_mask", n),
            "rescue": _mask_values_from_ds(ds, "bloom_rescue_mask", n),
            "cldice": _mask_values_from_ds(ds, "cldice_mask", n),
            "ci_valid": _mask_values_from_ds(ds, "ci_valid_mask", n),
            "candidate": _mask_values_from_ds(ds, "ci_candidate_mask", n),
            "detected": _mask_values_from_ds(ds, "ci_detection_mask", n),
        }
        for name in ("CI", "CI_cyano", "NDCI", "MPH", "FAI"):
            data[name] = (
                _m_science.to_float(ds.variables[name][:])
                if name in ds.variables
                else np.full(n, np.nan)
            )
        frame = pd.DataFrame(data)
        attrs = {
            k: getattr(ds, k, None)
            for k in ("scene_id", "scene_name", "acquisition_start", "window_id")
        }
    meta = {}
    if lake_metadata is not None and (not lake_metadata.empty):
        for _, row in lake_metadata.iterrows():
            meta[int(row["Hylak_id"])] = {
                "Lake_name": row.get("Lake_name", "") or "",
                "Lake_area_km2": float(row.get("Lake_area", np.nan)),
                "Country": row.get("Country", "") if "Country" in row else "",
                "Continent": row.get("Continent", "") if "Continent" in row else "",
            }

    def vals(series):
        a = pd.to_numeric(series, errors="coerce").to_numpy(float)
        return a[np.isfinite(a)]

    def summary(prefix, a):
        if not len(a):
            return {
                f"{prefix}_mean_retained": np.nan,
                f"{prefix}_median_retained": np.nan,
                f"{prefix}_p10_retained": np.nan,
                f"{prefix}_p90_retained": np.nan,
                f"{prefix}_max_retained": np.nan,
            }
        return {
            f"{prefix}_mean_retained": float(np.nanmean(a)),
            f"{prefix}_median_retained": float(np.nanmedian(a)),
            f"{prefix}_p10_retained": float(np.nanpercentile(a, 10)),
            f"{prefix}_p90_retained": float(np.nanpercentile(a, 90)),
            f"{prefix}_max_retained": float(np.nanmax(a)),
        }

    records = []
    for lake_id, g in frame.groupby("Hylak_id", sort=False):
        ci = vals(g.loc[g["candidate"], "CI"])
        cy = vals(g.loc[g["detected"], "CI_cyano"])
        nd = vals(g.loc[g["valid"], "NDCI"])
        valid_n = int(g["valid"].sum())
        strict_n = int(g["strict"].sum())
        ci_valid_n = int(g["ci_valid"].sum())
        detected_n = int(g["detected"].sum())
        rec = {
            **attrs,
            "Hylak_id": int(lake_id),
            **meta.get(int(lake_id), {}),
            "native_polygon_pixels": int(len(g)),
            "retained_pixels": valid_n,
            "retained_percent": 100.0 * valid_n / len(g) if len(g) else np.nan,
            "cyan_strict_pixels": strict_n,
            "bloom_rescued_pixels": int(g["rescue"].sum()),
            "CLDICE_pixels": int(g["cldice"].sum()),
            "ci_input_valid_pixels": ci_valid_n,
            "ci_candidate_pixels": int(g["candidate"].sum()),
            "ci_cyano_detection_pixels": detected_n,
            "ci_cyano_detection_fraction_of_ci_valid": detected_n / ci_valid_n
            if ci_valid_n
            else np.nan,
            "ci_cyano_detection_fraction_of_retained": detected_n / valid_n if valid_n else np.nan,
            "CI_mean_candidates": float(np.nanmean(ci)) if len(ci) else np.nan,
            "CI_median_candidates": float(np.nanmedian(ci)) if len(ci) else np.nan,
            "CI_p10_candidates": float(np.nanpercentile(ci, 10)) if len(ci) else np.nan,
            "CI_p90_candidates": float(np.nanpercentile(ci, 90)) if len(ci) else np.nan,
            "CI_max_candidates": float(np.nanmax(ci)) if len(ci) else np.nan,
            "CIcyano_mean_detections": float(np.nanmean(cy)) if len(cy) else np.nan,
            "CIcyano_median_detections": float(np.nanmedian(cy)) if len(cy) else np.nan,
            "CIcyano_p90_detections": float(np.nanpercentile(cy, 90)) if len(cy) else np.nan,
            "NDCI_mean_retained": float(np.nanmean(nd)) if len(nd) else np.nan,
            "NDCI_median_retained": float(np.nanmedian(nd)) if len(nd) else np.nan,
            "NDCI_p90_retained": float(np.nanpercentile(nd, 90)) if len(nd) else np.nan,
            "compact_netcdf": str(nc_path),
        }
        for k in _m_settings.ADDITIONAL_INDEX_NAMES:
            rec.update(summary(k, vals(g.loc[g["valid"], k])))
        records.append(rec)
    return pd.DataFrame(records)


# Reference cell 46, lines 170-178.
def all_stats(config_hash: str | None = None) -> pd.DataFrame:
    files = sorted(_m_settings.STATS_SCENE_DIR.glob("*.parquet"))
    if config_hash:
        files = [p for p in files if p.stem.endswith("_" + config_hash)]
    frames = []
    for p in files:
        with contextlib.suppress(Exception):
            frames.append(pd.read_parquet(p))
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


# Reference cell 57, lines 284-285.
def _mask_values_from_ds(ds, name, n):
    return np.asarray(ds.variables[name][:]) == 1 if name in ds.variables else np.zeros(n, bool)
