"""Standalone exports implementation."""

from __future__ import annotations

import contextlib
import json
import math
import zipfile
from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd

from . import archive as _m_archive
from . import mapping as _m_mapping
from . import quality as _m_quality
from . import registry as _m_registry
from . import settings as _m_settings
from . import state as _m_state
from . import ui_helpers as _m_ui_helpers
from . import utils as _m_utils


# Reference cell 59, lines 226-254.
def export_snap_gridded_netcdf(scene_id, scope, lake_id, metric, frame, resolution_m):
    """Fast displayed-scope EPSG:4326 export using same-lake nearest-native values."""
    if frame is None or len(frame) == 0:
        raise ValueError("No displayed pixels are available for quick export")
    names = ["rhos_665", "rhos_560", "rhos_490"] if metric == "true_color" else [metric]
    g = _m_mapping._science_geographic_setup_v265(frame, float(resolution_m))
    ncell = int(g["ncell"])
    if ncell > 35000000:
        raise RuntimeError(
            f"Exact {g['resolution_m']:.0f} m quick export would require {ncell / 1000000.0:.1f} million cells. Use the planning-region/selected-lake scope, or the sparse 3-file share exporter."
        )
    clip_geom = _m_state.GLOBAL_STATE.get("planning_geometry") if str(scope) == "planning" else None
    rows, cols, src, dist = _m_mapping._native_assignments_geographic_v267(
        frame, g, scene_id, clip_geom=clip_geom
    )
    matched = src >= 0
    source_values = {
        n: pd.to_numeric(frame[n], errors="coerce").to_numpy(float)
        for n in names
        if n in frame.columns
    }
    SNAP_DIR = _m_settings.EXPORT_DIR / "snap_gridded"
    SNAP_DIR.mkdir(parents=True, exist_ok=True)
    out = (
        SNAP_DIR
        / f"S3_OLCI_{scene_id}_{scope}_{lake_id or 'ALL'}_{metric}_{int(round(g['resolution_m']))}m_CF.nc"
    )
    from pyproj import CRS

    wkt = CRS.from_epsg(4326).to_wkt()
    chunk = 128
    chunks = (min(chunk, g["ny"]), min(chunk, g["nx"]))
    with netCDF4.Dataset(out, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("latitude", g["ny"])
        ds.createDimension("longitude", g["nx"])
        ds.setncatts(
            {
                "Conventions": "CF-1.8",
                "title": "Sentinel-3 OLCI displayed-scope same-lake native-nearest export",
                "scene_id": str(scene_id),
                "scope": str(scope),
                "nominal_grid_resolution_m": float(g["resolution_m"]),
                "science_grid_method": _m_settings.V267_GRID_METHOD,
                "nearest_native_radius_m": max(
                    _m_settings.V267_NEAREST_RADIUS_M, 2.0 * float(g["resolution_m"])
                ),
                "resampling": "same-lake nearest-native assignment; no averaging/value interpolation or gap filling",
                "qa_basis": "primary bloom-aware mask already encoded in displayed science values",
            }
        )
        y = ds.createVariable("latitude", "f8", ("latitude",))
        y[:] = g["latitude"]
        y.setncatts({"standard_name": "latitude", "units": "degrees_north", "axis": "Y"})
        x = ds.createVariable("longitude", "f8", ("longitude",))
        x[:] = g["longitude"]
        x.setncatts({"standard_name": "longitude", "units": "degrees_east", "axis": "X"})
        crs = ds.createVariable("crs", "i4")
        crs.assignValue(0)
        crs.setncatts(
            {
                "grid_mapping_name": "latitude_longitude",
                "epsg_code": "EPSG:4326",
                "spatial_ref": wkt,
                "crs_wkt": wkt,
                "semi_major_axis": 6378137.0,
                "inverse_flattening": 298.257223563,
            }
        )
        vars_ = {}
        for n in source_values:
            on = "rhos_884" if n == "rhos_885" else n
            v = ds.createVariable(
                on,
                "f4",
                ("latitude", "longitude"),
                zlib=True,
                complevel=5,
                shuffle=True,
                chunksizes=chunks,
                fill_value=np.float32(-9999),
            )
            v.setncatts(
                {
                    "coordinates": "latitude longitude",
                    "grid_mapping": "crs",
                    "units": "1",
                    "long_name": on,
                }
            )
            vars_[n] = v
        if len(rows):
            nxc = int(math.ceil(g["nx"] / chunk))
            cid = rows // chunk * nxc + cols // chunk
            order = np.argsort(cid, kind="mergesort")
            cids = cid[order]
            cuts = np.flatnonzero(np.r_[True, cids[1:] != cids[:-1], True])
            for a, b in zip(cuts[:-1], cuts[1:]):
                pos = order[a:b]
                r0 = int(rows[pos[0]] // chunk * chunk)
                c0 = int(cols[pos[0]] // chunk * chunk)
                r1 = min(r0 + chunk, g["ny"])
                c1 = min(c0 + chunk, g["nx"])
                rr = rows[pos] - r0
                cc = cols[pos] - c0
                for n, v in vars_.items():
                    block = np.full((r1 - r0, c1 - c0), np.float32(-9999), np.float32)
                    good = src[pos] >= 0
                    if good.any():
                        vals = source_values[n][src[pos][good]]
                        ok = np.isfinite(vals)
                        block[rr[good][ok], cc[good][ok]] = vals[ok].astype(np.float32)
                    v[r0:r1, c0:c1] = np.ma.masked_equal(block, np.float32(-9999))
        if metric == "true_color":
            ds.setncattr("true_color_mapping", "red=rhos_665; green=rhos_560; blue=rhos_490")
    return str(out)


# Reference cell 50, lines 668-676.
def _science_grid_setup(
    frame, requested_resolution_m=_m_settings.SCIENCE_EXPORT_DEFAULT_RESOLUTION_M
):
    """Prepare one regular EPSG:4326 grid for a complete SNAP science product."""
    lat = frame.latitude.to_numpy(float)
    lon = frame.longitude.to_numpy(float)
    g = _m_mapping._grid_geometry_geographic(
        lat, lon, float(requested_resolution_m), max_cells=_m_settings.SCIENCE_EXPORT_MAX_CELLS
    )
    rr = np.floor((g["north"] - lat) / g["dlat"]).astype(np.int64)
    cc = np.floor((lon - g["west"]) / g["dlon"]).astype(np.int64)
    inside = (
        np.isfinite(lat)
        & np.isfinite(lon)
        & (rr >= 0)
        & (rr < g["ny"])
        & (cc >= 0)
        & (cc < g["nx"])
    )
    flat = rr[inside] * g["nx"] + cc[inside]
    return (g, inside, flat)


# Reference cell 50, lines 679-688.
def _mean_grid_from_flat(values, inside, flat, ncell):
    values = np.asarray(values, float)[inside]
    good = np.isfinite(values)
    arr = np.full(int(ncell), np.nan, dtype=np.float32)
    if good.any():
        sums = np.bincount(flat[good], weights=values[good], minlength=int(ncell))
        counts = np.bincount(flat[good], minlength=int(ncell))
        nz = counts > 0
        arr[nz] = (sums[nz] / counts[nz]).astype(np.float32)
    return arr


# Reference cell 50, lines 691-697.
def _mask_grid_from_flat(values, inside, flat, ncell):
    values = np.asarray(values, float)[inside]
    arr = np.zeros(int(ncell), dtype=np.uint8)
    good = np.isfinite(values) & (values > 0)
    if good.any():
        np.maximum.at(arr, flat[good], np.uint8(1))
    return arr


# Reference cell 59, lines 257-296.
def export_complete_scene_science_netcdf(
    scene_id, requested_resolution_m=_m_settings.SCIENCE_EXPORT_DEFAULT_RESOLUTION_M, force=False
):
    """Production share raster: exact nominal grid, same-lake nearest-native remap, sparse writes."""
    if not scene_id:
        raise ValueError("Choose a completed scene")
    frame = _m_archive._load_scene_pixels(scene_id, "scene", None).copy()
    g = _m_mapping._science_geographic_setup_v265(frame, float(requested_resolution_m))
    export_dir = _m_settings.EXPORT_DIR / "snap_complete"
    export_dir.mkdir(parents=True, exist_ok=True)
    out = export_dir / f"S3_OLCI_{scene_id}_ALL_SCIENCE_{int(round(g['resolution_m']))}m_CF.nc"
    source_files = _m_archive._scene_compact_files(scene_id)
    newest = max((p.stat().st_mtime_ns for p in source_files), default=0)
    if out.exists() and (not force) and (out.stat().st_mtime_ns >= newest):
        with contextlib.suppress(Exception):
            with netCDF4.Dataset(out) as old:
                if getattr(
                    old, "science_grid_method", ""
                ) == _m_settings.V267_GRID_METHOD and float(
                    getattr(old, "nominal_grid_resolution_m", -1)
                ) == float(g["resolution_m"]):
                    return (str(out), g)
    rows, cols, src, dist = _m_mapping._native_assignments_geographic_v267(frame, g, scene_id)
    matched = src >= 0
    primary = _m_quality._primary_selector(frame)
    qa = _m_archive._compact_qa_metadata(scene_id)
    source_values = {
        n: pd.to_numeric(frame[n], errors="coerce").to_numpy(float)
        for n in _m_settings.SCIENCE_EXPORT_VARIABLES
        if n in frame.columns
    }
    source_hylak = pd.to_numeric(frame.Hylak_id, errors="coerce").fillna(-1).to_numpy(np.int32)
    mask_input_names = [
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
        "mph_valid_mask",
        "fai_valid_mask",
    ]
    source_masks = {
        n: pd.to_numeric(frame[n], errors="coerce").fillna(0).to_numpy(float) > 0
        for n in mask_input_names
        if n in frame.columns
    }
    mask_output_name = {"cyan_strict_valid_mask": "generic_flag_exclusion_sensitivity_mask"}
    tmp = out.with_suffix(".nc.tmp")
    tmp.unlink(missing_ok=True)
    chunk = int(_m_settings.V267_EXPORT_CHUNK)
    with netCDF4.Dataset(tmp, "w", format="NETCDF4_CLASSIC") as ds:
        ds.createDimension("latitude", g["ny"])
        ds.createDimension("longitude", g["nx"])
        ds.setncatts(
            {
                "Conventions": "CF-1.8",
                "title": "Sentinel-3 OLCI HydroLAKES bloom-aware all-science efficient 300 m product",
                "scene_id": str(scene_id),
                "science_grid_method": _m_settings.V267_GRID_METHOD,
                "nominal_grid_resolution_m": float(g["resolution_m"]),
                "actual_grid_resolution_m": float(g["resolution_m"]),
                "cell_height_m_nominal": float(g["resolution_m"]),
                "cell_width_m_at_reference_latitude": float(g["resolution_m"]),
                "reference_latitude_for_longitude_spacing": float(g["reference_latitude"]),
                "geospatial_lat_min": float(g["south"]),
                "geospatial_lat_max": float(g["north"]),
                "geospatial_lon_min": float(g["west"]),
                "geospatial_lon_max": float(g["east"]),
                "latitude_step_degrees": float(g["dlat"]),
                "longitude_step_degrees": float(g["dlon"]),
                "nearest_native_radius_m": max(
                    _m_settings.V267_NEAREST_RADIUS_M, 2.0 * float(g["resolution_m"])
                ),
                "mapping": "each target lake cell maps only to the nearest stored native OLCI observation from the same Hylak_id",
                "qa_rule": "nearest same-lake source is chosen before QA; if that native source is not primary-valid, science value remains fill/NaN",
                "resampling": "nearest-native assignment; no averaging/value interpolation or gap filling",
                "striping_policy": "target lake support is populated from same-lake nearest native observations instead of center-only binning",
                "primary_qa_mask": "valid_water_mask (bloom-aware)",
                "generic_CLDICE_sensitivity": "diagnostic only; not claimed equivalent to NASA CyAN upstream cloud processing",
                "NASA_area_weighting": "not used in production v2.6.7; retained as optional validation concept only",
                "statistics_source": "native compact OLCI lake observations, not mapped raster",
                "true_color_mapping": "red=rhos_665; green=rhos_560; blue=rhos_490",
                "CI_cyano_detection_limit": float(_m_settings.CYAN_CI_DETECTION_LIMIT),
                "qa_policy_json": json.dumps(qa, default=str),
            }
        )
        yv = ds.createVariable("latitude", "f8", ("latitude",))
        yv[:] = g["latitude"]
        yv.setncatts({"standard_name": "latitude", "units": "degrees_north", "axis": "Y"})
        xv = ds.createVariable("longitude", "f8", ("longitude",))
        xv[:] = g["longitude"]
        xv.setncatts({"standard_name": "longitude", "units": "degrees_east", "axis": "X"})
        from pyproj import CRS

        wkt4326 = CRS.from_epsg(4326).to_wkt()
        crsv = ds.createVariable("crs", "i4")
        crsv.assignValue(0)
        crsv.setncatts(
            {
                "grid_mapping_name": "latitude_longitude",
                "epsg_code": "EPSG:4326",
                "spatial_ref": wkt4326,
                "crs_wkt": wkt4326,
                "semi_major_axis": 6378137.0,
                "inverse_flattening": 298.257223563,
            }
        )
        chunks = (min(chunk, g["ny"]), min(chunk, g["nx"]))
        vars_float = {}
        long_names = {
            "CI": "Cyanobacteria Index candidate magnitude",
            "CI_cyano": "cyanobacteria index detections",
            "NDCI": "Normalized Difference Chlorophyll Index",
            "MPH": "Maximum Peak Height",
            "FAI": "OLCI-adapted Floating Algae Index (AFAI form)",
        }
        for sn, on in [
            (n, "rhos_884" if n == "rhos_885" else n)
            for n in _m_settings.SCIENCE_EXPORT_VARIABLES
            if n in source_values
        ]:
            v = ds.createVariable(
                on,
                "f4",
                ("latitude", "longitude"),
                zlib=True,
                complevel=6,
                shuffle=True,
                chunksizes=chunks,
                fill_value=np.float32(-9999),
            )
            v.setncatts(
                {
                    "coordinates": "latitude longitude",
                    "grid_mapping": "crs",
                    "units": "1",
                    "long_name": long_names.get(on, f"Rayleigh-corrected OLCI reflectance {on}"),
                    "qa_basis": "primary bloom-aware valid_water_mask",
                }
            )
            vars_float[sn] = v
        vars_mask = {}
        for sn in source_masks:
            on = mask_output_name.get(sn, sn)
            v = ds.createVariable(
                on,
                "i1",
                ("latitude", "longitude"),
                zlib=True,
                complevel=6,
                shuffle=True,
                chunksizes=chunks,
                fill_value=np.int8(-1),
            )
            v.setncatts(
                {
                    "coordinates": "latitude longitude",
                    "grid_mapping": "crs",
                    "long_name": on.replace("_", " "),
                    "flag_values": np.array([0, 1], dtype=np.int8),
                    "flag_meanings": "false true",
                    "valid_range": np.array([0, 1], dtype=np.int8),
                }
            )
            vars_mask[sn] = v
        support_var = ds.createVariable(
            "lake_support_mask",
            "i1",
            ("latitude", "longitude"),
            zlib=True,
            complevel=5,
            shuffle=True,
            chunksizes=chunks,
            fill_value=np.int8(-1),
        )
        support_var.setncatts(
            {
                "coordinates": "latitude longitude",
                "grid_mapping": "crs",
                "long_name": "target cell lies inside the analysis HydroLAKES polygon after shoreline buffer",
                "flag_values": np.array([0, 1], dtype=np.int8),
                "flag_meanings": "false true",
            }
        )
        match_var = ds.createVariable(
            "native_source_match_mask",
            "i1",
            ("latitude", "longitude"),
            zlib=True,
            complevel=5,
            shuffle=True,
            chunksizes=chunks,
            fill_value=np.int8(-1),
        )
        match_var.setncatts(
            {
                "coordinates": "latitude longitude",
                "grid_mapping": "crs",
                "long_name": "target lake cell has a same-lake native OLCI source within search radius",
                "flag_values": np.array([0, 1], dtype=np.int8),
                "flag_meanings": "false true",
            }
        )
        dist_var = ds.createVariable(
            "native_match_distance_m",
            "f4",
            ("latitude", "longitude"),
            zlib=True,
            complevel=5,
            shuffle=True,
            chunksizes=chunks,
            fill_value=np.float32(-9999),
        )
        dist_var.setncatts(
            {
                "coordinates": "latitude longitude",
                "grid_mapping": "crs",
                "units": "m",
                "long_name": "distance to nearest same-lake stored native OLCI observation",
            }
        )
        hylak_var = ds.createVariable(
            "nearest_Hylak_id",
            "i4",
            ("latitude", "longitude"),
            zlib=True,
            complevel=5,
            shuffle=True,
            chunksizes=chunks,
            fill_value=np.int32(-1),
        )
        hylak_var.setncatts(
            {
                "coordinates": "latitude longitude",
                "grid_mapping": "crs",
                "long_name": "HydroLAKES identifier of same-lake native source",
            }
        )
        if len(rows):
            nxc = int(math.ceil(g["nx"] / chunk))
            cid = rows // chunk * nxc + cols // chunk
            order = np.argsort(cid, kind="mergesort")
            cids = cid[order]
            cuts = np.flatnonzero(np.r_[True, cids[1:] != cids[:-1], True])
            for a, b in zip(cuts[:-1], cuts[1:]):
                pos = order[a:b]
                r0 = int(rows[pos[0]] // chunk * chunk)
                c0 = int(cols[pos[0]] // chunk * chunk)
                r1 = min(r0 + chunk, g["ny"])
                c1 = min(c0 + chunk, g["nx"])
                rr = rows[pos] - r0
                cc = cols[pos] - c0
                spos = src[pos]
                m = spos >= 0
                sb = np.full((r1 - r0, c1 - c0), np.int8(-1), np.int8)
                sb[rr, cc] = 1
                support_var[r0:r1, c0:c1] = np.ma.masked_equal(sb, np.int8(-1))
                mb = np.full_like(sb, np.int8(-1))
                mb[rr, cc] = m.astype(np.int8)
                match_var[r0:r1, c0:c1] = np.ma.masked_equal(mb, np.int8(-1))
                db = np.full((r1 - r0, c1 - c0), np.float32(-9999), np.float32)
                valid_d = m & np.isfinite(dist[pos])
                db[rr[valid_d], cc[valid_d]] = dist[pos][valid_d].astype(np.float32)
                dist_var[r0:r1, c0:c1] = np.ma.masked_equal(db, np.float32(-9999))
                hb = np.full((r1 - r0, c1 - c0), np.int32(-1), np.int32)
                hb[rr[m], cc[m]] = source_hylak[spos[m]]
                hylak_var[r0:r1, c0:c1] = np.ma.masked_equal(hb, np.int32(-1))
                for sn, v in vars_float.items():
                    block = np.full((r1 - r0, c1 - c0), np.float32(-9999), np.float32)
                    if m.any():
                        vv = source_values[sn][spos[m]]
                        ok = np.isfinite(vv) & primary[spos[m]]
                        block[rr[m][ok], cc[m][ok]] = vv[ok].astype(np.float32)
                    v[r0:r1, c0:c1] = np.ma.masked_equal(block, np.float32(-9999))
                for sn, v in vars_mask.items():
                    block = np.full((r1 - r0, c1 - c0), np.int8(-1), np.int8)
                    if m.any():
                        block[rr[m], cc[m]] = source_masks[sn][spos[m]].astype(np.int8)
                    v[r0:r1, c0:c1] = np.ma.masked_equal(block, np.int8(-1))
    tmp.replace(out)
    g.update(
        {
            "science_grid_method": _m_settings.V267_GRID_METHOD,
            "target_lake_cells": int(len(rows)),
            "matched_target_cells": int(matched.sum()),
            "nearest_radius_m": max(
                _m_settings.V267_NEAREST_RADIUS_M, 2.0 * float(g["resolution_m"])
            ),
        }
    )
    return (str(out), g)


# Reference cell 59, lines 325-335.
def create_scene_share_bundle(
    scene_id, science_resolution_m=_m_settings.SCIENCE_EXPORT_DEFAULT_RESOLUTION_M
):
    """Exactly three files from the efficient production branch."""
    if not scene_id:
        raise ValueError("Choose a completed scene")
    t = _m_registry.registry_table(_m_state.GLOBAL_STATE.get("config_hash"))
    r = t[t.scene_id.astype(str) == str(scene_id)]
    if r.empty:
        raise KeyError(scene_id)
    row = r.iloc[0]
    cfg = str(row.config_hash)
    d = _m_settings.EXPORT_DIR / "scene_bundles"
    d.mkdir(parents=True, exist_ok=True)
    complete_nc, g = export_complete_scene_science_netcdf(
        scene_id, float(science_resolution_m), force=False
    )
    stats_csv, stats_table = _share_all_index_stats(scene_id, cfg, row)
    qa = _m_archive._compact_qa_metadata(scene_id)
    frame_meta = _m_registry.parse_s3_frame_metadata(row.scene_name)
    meta = {
        "schema": "global-efficient-share-v2.6.7",
        "created_utc": _m_utils.utc_now(),
        "scene_id": str(scene_id),
        "scene_name": str(row.scene_name),
        "acquisition_start": str(row.acquisition_start),
        "config_hash": cfg,
        **frame_meta,
        "files": {
            "snap_all_science_netcdf": Path(complete_nc).name,
            "all_index_stats_csv": Path(stats_csv).name,
        },
        "snap_variables": [
            "rhos_490",
            "rhos_560",
            "rhos_620",
            "rhos_665",
            "rhos_681",
            "rhos_709",
            "rhos_754",
            "rhos_865",
            "rhos_884",
            "CI",
            "CI_cyano",
            "NDCI",
            "MPH",
            "FAI",
            "valid_water_mask",
            "generic_flag_exclusion_sensitivity_mask",
            "bloom_rescue_mask",
            "cldice_mask",
            "hisatzen_mask",
            "lake_support_mask",
            "native_source_match_mask",
            "native_match_distance_m",
            "nearest_Hylak_id",
        ],
        "snap_grid": {
            "projection": "EPSG:4326 geographic",
            "grid_method": _m_settings.V267_GRID_METHOD,
            "nominal_resolution_m": float(g["resolution_m"]),
            "reference_latitude": float(g["reference_latitude"]),
            "latitude_step_degrees": float(g["dlat"]),
            "longitude_step_degrees": float(g["dlon"]),
            "nearest_radius_m": float(
                g.get(
                    "nearest_radius_m",
                    max(_m_settings.V267_NEAREST_RADIUS_M, 2.0 * float(g["resolution_m"])),
                )
            ),
            "mapping": "same-lake nearest native observation",
            "interpolation": "none",
            "gap_fill": "none",
            "area_weighting": "not used in production",
        },
        "masking": qa,
        "strict_mask_policy": "Generic LAND+CLDICE+HISATZEN-style exclusion is retained only as a sensitivity mask/count; it is not called CyAN-equivalent because NASA upstream cloud processing differs.",
        "statistics": {
            "source": "native compact OLCI observations",
            "mapped_raster_not_used_for_statistics": True,
        },
        "future_daily_composite": "S3A/S3B daily CIcyano maximum is intentionally not implemented yet; add only after this efficient mapping branch is validated.",
        "important": "The mapped raster is for visualization/GIS/SNAP interoperability. Quantitative lake statistics remain native-observation based.",
    }
    mp = d / f"S3_OLCI_{scene_id}_METADATA.json"
    mp.write_text(json.dumps(meta, indent=2, default=str))
    z = d / f"{str(row.scene_name)}_{cfg}_SHARE_MINIMAL.zip"
    with zipfile.ZipFile(z, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as arc:
        for p in (Path(complete_nc), Path(stats_csv), mp):
            arc.write(p, arcname=p.name)
    return (
        str(z),
        _m_ui_helpers._app_status(
            "Minimal share bundle created",
            f"3 files · nominal {g['resolution_m']:.0f} m EPSG:4326 · same-lake native-nearest · bloom-aware QA · no area-weighting overhead",
        ),
    )


# Reference cell 55, lines 579-583.
def bundle_callback(scene_id, science_resolution_m):
    try:
        return create_scene_share_bundle(scene_id, science_resolution_m)
    except Exception as exc:
        return (None, _m_ui_helpers._app_status("Bundle failed", str(exc), False))


# Reference cell 57, lines 436-450.
def _science_sinusoidal_setup(
    frame, requested_resolution_m=_m_settings.SCIENCE_EXPORT_DEFAULT_RESOLUTION_M
):
    """Exact-resolution globally aligned spherical-sinusoidal grid; never auto-coarsens."""
    from pyproj import CRS, Transformer

    lat = pd.to_numeric(frame["latitude"], errors="coerce").to_numpy(float)
    lon = pd.to_numeric(frame["longitude"], errors="coerce").to_numpy(float)
    finite = np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90) & (np.abs(lon) <= 180)
    if not finite.any():
        raise ValueError("Selected scene contains no finite stored geolocated lake pixels")
    radius = float(_m_settings.CYAN_STYLE_SPHERE_RADIUS_M)
    crs = CRS.from_proj4(f"+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +R={radius:.3f} +units=m +no_defs")
    transformer = Transformer.from_crs(4326, crs, always_xy=True)
    x_all = np.full(lon.shape, np.nan, np.float64)
    y_all = np.full(lat.shape, np.nan, np.float64)
    xx, yy = transformer.transform(lon[finite], lat[finite])
    x_all[finite] = xx
    y_all[finite] = yy
    global_x_min = -math.pi * radius
    global_y_max = 0.5 * math.pi * radius
    res = max(1.0, float(requested_resolution_m))
    cg = np.floor((x_all[finite] - global_x_min) / res).astype(np.int64)
    rg = np.floor((global_y_max - y_all[finite]) / res).astype(np.int64)
    c0, c1 = (int(cg.min()), int(cg.max()))
    r0, r1 = (int(rg.min()), int(rg.max()))
    nx, ny = (c1 - c0 + 1, r1 - r0 + 1)
    finite_idx = np.flatnonzero(finite)
    inside = np.zeros(len(frame), bool)
    inside[finite_idx] = True
    cc = cg - c0
    rr = rg - r0
    flat = rr * nx + cc
    x_left = global_x_min + c0 * res
    y_top = global_y_max - r0 * res
    x = x_left + (np.arange(nx, dtype=np.float64) + 0.5) * res
    y = y_top - (np.arange(ny, dtype=np.float64) + 0.5) * res
    g = {
        "crs": crs,
        "resolution_m": float(res),
        "requested_resolution_m": float(requested_resolution_m),
        "nx": int(nx),
        "ny": int(ny),
        "ncell": int(nx) * int(ny),
        "x": x,
        "y": y,
        "x_left": float(x_left),
        "y_top": float(y_top),
        "global_col0": int(c0),
        "global_row0": int(r0),
        "radius_m": radius,
    }
    return (g, inside, flat)


# Reference cell 55, lines 137-146.
def _binned_mask_grid(values, inside, flat, ncell):
    """Binary-mask aggregation with -1 for bins containing no native observations."""
    vals = np.asarray(values, float)[inside]
    support = np.bincount(flat, minlength=int(ncell)) > 0
    out = np.full(int(ncell), -1, dtype=np.int8)
    out[support] = 0
    positive = np.isfinite(vals) & (vals > 0)
    if positive.any():
        np.maximum.at(out, flat[positive], np.int8(1))
    return out


# Reference cell 59, lines 299-322.
def _share_all_index_stats(scene_id, cfg, row):
    """Lean native-observation statistics; generic CLDICE exclusion is counts-only sensitivity QA."""
    frame = _m_archive._load_scene_pixels(scene_id, "scene", None).copy()
    if frame.empty:
        raise ValueError("No compact lake pixels found for selected scene")
    for c in ("latitude", "longitude", "Hylak_id"):
        frame[c] = pd.to_numeric(frame[c], errors="coerce")
    frame = frame.dropna(subset=["latitude", "longitude", "Hylak_id"]).copy()
    frame["_lat6"] = frame.latitude.round(6)
    frame["_lon6"] = frame.longitude.round(6)
    frame = frame.drop_duplicates(["Hylak_id", "_lat6", "_lon6"]).drop(columns=["_lat6", "_lon6"])
    lakes = _m_state.GLOBAL_STATE.get("lakes")
    lake_meta = {}
    if lakes is not None and (not lakes.empty):
        for _, lr in lakes[
            lakes.Hylak_id.astype(int).isin(frame.Hylak_id.astype(int).unique())
        ].iterrows():
            lake_meta[int(lr.Hylak_id)] = {
                "Lake_name": lr.get("Lake_name", "") or "",
                "Lake_area_km2": float(lr.get("Lake_area", np.nan)),
            }

    def mask(g, n):
        return (
            pd.to_numeric(g.get(n, pd.Series(np.zeros(len(g)), index=g.index)), errors="coerce")
            .fillna(0)
            .to_numpy(float)
            > 0
        )

    def stat(s):
        a = pd.to_numeric(s, errors="coerce").to_numpy(float)
        a = a[np.isfinite(a)]
        if not len(a):
            return {
                "valid_pixels": 0,
                "mean": np.nan,
                "median": np.nan,
                "p10": np.nan,
                "p90": np.nan,
                "max": np.nan,
            }
        return {
            "valid_pixels": int(len(a)),
            "mean": float(np.mean(a)),
            "median": float(np.median(a)),
            "p10": float(np.percentile(a, 10)),
            "p90": float(np.percentile(a, 90)),
            "max": float(np.max(a)),
        }

    recs = []
    for lid, gdf in frame.groupby(frame.Hylak_id.astype(int), sort=True):
        primary = mask(gdf, "valid_water_mask")
        sens = mask(gdf, "cyan_strict_valid_mask")
        rescue = mask(gdf, "bloom_rescue_mask")
        cld = mask(gdf, "cldice_mask")
        cloud_ex = mask(gdf, "cloud_excluded_mask")
        cand = mask(gdf, "ci_candidate_mask")
        det = mask(gdf, "ci_detection_mask")
        ci_valid = mask(gdf, "ci_valid_mask")
        rec = {
            "scene_id": str(scene_id),
            "scene_name": str(row.scene_name),
            "acquisition_start": str(row.acquisition_start),
            "Hylak_id": int(lid),
            **lake_meta.get(int(lid), {}),
            "native_lake_pixels": int(len(gdf)),
            "primary_retained_pixels": int(primary.sum()),
            "generic_flag_exclusion_sensitivity_pixels": int(sens.sum()),
            "bloom_rescued_pixels": int(rescue.sum()),
            "CLDICE_pixels": int(cld.sum()),
            "cloud_excluded_pixels": int(cloud_ex.sum()),
            "bloom_rescue_fraction_of_CLDICE": float(rescue.sum() / cld.sum())
            if cld.sum()
            else np.nan,
            "ci_input_valid_pixels": int(ci_valid.sum()),
            "ci_candidate_pixels": int(cand.sum()),
            "ci_cyano_detection_pixels": int(det.sum()),
        }
        bases = {"CI": cand, "CI_cyano": det, "NDCI": primary, "MPH": primary, "FAI": primary}
        for metric, basis in bases.items():
            s = stat(gdf.loc[basis, metric] if metric in gdf else pd.Series(dtype=float))
            for k, v in s.items():
                rec[f"{metric}_{k}"] = v
        recs.append(rec)
    table = pd.DataFrame(recs).sort_values("Hylak_id").reset_index(drop=True)
    d = _m_settings.EXPORT_DIR / "scene_bundles"
    d.mkdir(parents=True, exist_ok=True)
    p = d / f"S3_OLCI_{scene_id}_ALL_INDEX_STATS.csv"
    table.to_csv(p, index=False)
    return (p, table)


# Reference cell 57, lines 453-454.
def _unique_bin_index(flat, nx):
    unique, inv = np.unique(np.asarray(flat, dtype=np.int64), return_inverse=True)
    rows = (unique // int(nx)).astype(np.int64)
    cols = (unique % int(nx)).astype(np.int64)
    return (unique, inv, rows, cols)


# Reference cell 57, lines 457-463.
def _aggregate_mean_to_unique(values, inside, inv, n_unique, obs_valid=None):
    vals = np.asarray(values, float)[inside]
    good = np.isfinite(vals)
    if obs_valid is not None:
        good &= np.asarray(obs_valid, bool)[inside]
    out = np.full(int(n_unique), np.nan, np.float32)
    if good.any():
        sums = np.bincount(inv[good], weights=vals[good], minlength=int(n_unique))
        counts = np.bincount(inv[good], minlength=int(n_unique))
        nz = counts > 0
        out[nz] = (sums[nz] / counts[nz]).astype(np.float32)
    return out


# Reference cell 57, lines 466-469.
def _aggregate_mask_to_unique(values, inside, inv, n_unique):
    vals = np.asarray(values, float)[inside]
    out = np.zeros(int(n_unique), np.int8)
    pos = np.isfinite(vals) & (vals > 0)
    if pos.any():
        np.maximum.at(out, inv[pos], np.int8(1))
    return out


# Reference cell 57, lines 472-475.
def _chunk_groups(rows, cols, ny, nx, chunk=_m_settings.SCIENCE_EXPORT_CHUNK):
    cy = max(1, min(int(chunk), int(ny)))
    cx = max(1, min(int(chunk), int(nx)))
    nxc = int(math.ceil(nx / cx))
    cid = rows // cy * nxc + cols // cx
    order = np.argsort(cid, kind="mergesort")
    cid_sorted = cid[order]
    cuts = np.flatnonzero(np.r_[True, cid_sorted[1:] != cid_sorted[:-1], True])
    groups = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        groups.append(order[a:b])
    return (groups, cy, cx)


# Reference cell 57, lines 478-483.
def _write_sparse_float(var, values, rows, cols, groups, cy, cx, ny, nx, fill=-9999.0):
    vals = np.asarray(values, np.float32)
    for pos in groups:
        finite = np.isfinite(vals[pos])
        if not finite.any():
            continue
        p = pos[finite]
        r0 = int(rows[p[0]] // cy * cy)
        c0 = int(cols[p[0]] // cx * cx)
        r1 = min(r0 + cy, int(ny))
        c1 = min(c0 + cx, int(nx))
        block = np.full((r1 - r0, c1 - c0), np.float32(fill), np.float32)
        block[rows[p] - r0, cols[p] - c0] = vals[p]
        var[r0:r1, c0:c1] = np.ma.masked_equal(block, np.float32(fill))


# Reference cell 57, lines 486-489.
def _write_sparse_mask(var, values, rows, cols, groups, cy, cx, ny, nx, fill=-1):
    vals = np.asarray(values, np.int8)
    for pos in groups:
        r0 = int(rows[pos[0]] // cy * cy)
        c0 = int(cols[pos[0]] // cx * cx)
        r1 = min(r0 + cy, int(ny))
        c1 = min(c0 + cx, int(nx))
        block = np.full((r1 - r0, c1 - c0), np.int8(fill), np.int8)
        block[rows[pos] - r0, cols[pos] - c0] = vals[pos]
        var[r0:r1, c0:c1] = np.ma.masked_equal(block, np.int8(fill))


# Reference cell 57, lines 493-499.
def _write_sparse_int(var, values, rows, cols, groups, cy, cx, ny, nx, fill=-1):
    vals = np.asarray(values, np.int32)
    for pos in groups:
        r0 = int(rows[pos[0]] // cy * cy)
        c0 = int(cols[pos[0]] // cx * cx)
        r1 = min(r0 + cy, int(ny))
        c1 = min(c0 + cx, int(nx))
        block = np.full((r1 - r0, c1 - c0), np.int32(fill), np.int32)
        block[rows[pos] - r0, cols[pos] - c0] = vals[pos]
        var[r0:r1, c0:c1] = np.ma.masked_equal(block, np.int32(fill))
