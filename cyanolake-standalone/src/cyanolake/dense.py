"""Standalone dense implementation."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Iterable

import geopandas as gpd
import netCDF4
import numpy as np
import pandas as pd
from scipy.ndimage import binary_dilation
from tqdm.auto import tqdm

from . import hydrolakes as _m_hydrolakes
from . import science as _m_science
from . import settings as _m_settings
from . import utils as _m_utils


# Reference cell 30, lines 132-263.
def diagnose_l2_mask(
    l2_path: Path,
    ndci_source: str = None,
    detection_limit: float = _m_settings.CI_DETECTION_LIMIT,
    ss681_max: float = 0.0,
    ss665_min: float = 0.0,
    qa_flags: Iterable[str] = ("NAVFAIL",),
    land_adjacency_pixels: int = _m_settings.LAND_ADJACENCY_PIXELS,
    water_definition: str = "union",
    lakes: gpd.GeoDataFrame | None = None,
    shore_buffer_m: float = None,
    bright_screen: bool = _m_settings.DEFAULT_BRIGHT_SCREEN,
    bright_rhos865: float = _m_settings.DEFAULT_BRIGHT_RHOS865,
) -> dict[str, Any]:
    """Audit exactly where pixels are lost, entirely after L2Gen."""
    if ndci_source is None:
        ndci_source = _m_settings.NDCI_SOURCE
    if shore_buffer_m is None:
        shore_buffer_m = _m_settings.DEFAULT_SHORE_BUFFER_M
    context = _scene_mask_context(
        Path(l2_path),
        qa_flags,
        land_adjacency_pixels,
        water_definition,
        lakes,
        shore_buffer_m,
        bright_screen,
        bright_rhos865,
    )
    with netCDF4.Dataset(l2_path) as src:
        geo = src.groups["geophysical_data"]
        names = {
            wave: _m_science.nearest_product(geo.variables, "rhos", wave)
            for wave in _m_settings.CORE_RHOS_WAVELENGTHS
        }
        core = {wave: _m_science.to_float(geo.variables[name][:]) for wave, name in names.items()}
        finite_core = np.logical_and.reduce([np.isfinite(array) for array in core.values()])
        family = _m_science.ci_family(
            core[620],
            core[665],
            core[681],
            core[709],
            wavelengths=tuple(
                (_m_science.product_wavelength(names[w]) for w in _m_settings.CORE_RHOS_WAVELENGTHS)
            ),
            detection_limit=detection_limit,
            ss681_max=ss681_max,
            ss665_min=ss665_min,
        )
        prefix = "rhos" if ndci_source.lower() == "rhos" else "Rrs"
        red_name = _m_science.nearest_product(geo.variables, prefix, 665)
        edge_name = _m_science.nearest_product(geo.variables, prefix, 709)
        _, ndci_finite = _m_science.ndci_index(
            _m_science.to_float(geo.variables[red_name][:]),
            _m_science.to_float(geo.variables[edge_name][:]),
        )
    water = context["water"]
    retained = context["valid_water"]
    ci_valid = retained & finite_core
    stages = {
        "total_pixels": int(water.size),
        "geolocated_pixels": int(context["geolocated"].sum()),
        "water_pixels": int(water.sum()),
        "after_selected_flags": int((water & ~context["flag_excluded"]).sum()),
        "after_land_adjacency": int(
            (water & ~context["flag_excluded"] & ~context["adjacent"]).sum()
        ),
        "valid_water_pixels": int(retained.sum()),
        "ci_finite_before_qa": int((water & finite_core).sum()),
        "ci_finite_all_geolocated": int((context["geolocated"] & finite_core).sum()),
        "ci_candidates_all_geolocated": int((context["geolocated"] & family["ci_candidate"]).sum()),
        "ci_detections_all_geolocated": int((context["geolocated"] & family["ci_detection"]).sum()),
        "ci_finite_after_qa": int(ci_valid.sum()),
        "ci_candidates_after_qa": int((ci_valid & family["ci_candidate"]).sum()),
        "ci_detections_after_qa": int((ci_valid & family["ci_detection"]).sum()),
        "ndci_valid_before_qa": int((water & ndci_finite).sum()),
        "ndci_valid_all_geolocated": int((context["geolocated"] & ndci_finite).sum()),
        "ndci_valid_after_qa": int((retained & ndci_finite).sum()),
    }
    water_count = max(stages["water_pixels"], 1)
    details = {}
    for name in sorted(context["flag_map"]):
        flagged = water & _m_science.flag_is_set(context["flags"], context["flag_map"], name)
        details[name] = {
            "selected_for_exclusion": name in context["requested_flags"],
            "water_pixel_count": int(flagged.sum()),
            "percent_of_water": 100.0 * float(flagged.sum()) / water_count,
        }
    messages = [
        f"Water definition resolved to: {context['water_definition_resolved']}.",
        "The native L2 NetCDF was not changed; every exclusion shown here is post-L2.",
    ]
    if bright_screen:
        messages.append(
            f"Bright screen removed {int((water & context['bright']).sum()):,} water pixels using {context['bright_product']} > {bright_rhos865:g}. Disable it when testing surface scums."
        )
    if stages["ci_finite_before_qa"] == 0:
        messages.append(
            "No selected-water pixels have four finite CI rhos inputs; inspect L2 rhos coverage."
        )
    elif stages["ci_finite_after_qa"] == 0:
        messages.append(
            "Finite CI inputs exist before QA, so the selected post-mask is the cause of total loss."
        )
    ci_candidates = stages["ci_candidates_after_qa"]
    ci_detections = stages["ci_detections_after_qa"]
    ss665_rejections = int(
        (ci_valid & family["ci_candidate"] & ~(family["SS_665"] > float(ss665_min))).sum()
    )
    threshold_rejections = int(
        (
            ci_valid
            & family["ci_candidate"]
            & (family["SS_665"] > float(ss665_min))
            & ~(family["CI"] > float(detection_limit))
        ).sum()
    )
    if ci_candidates and (not ci_detections):
        messages.append(
            f"The SS681 test found {ci_candidates:,} CI candidate pixels, but none passed both SS665>{float(ss665_min):g} and CI>{float(detection_limit):g}. SS665 rejected {ss665_rejections:,}; the CI detection limit rejected {threshold_rejections:,}. This is an algorithmic non-detection, not a QA/water-mask loss; inspect CI, SS_665 and CI_exclusion_reason."
        )
    elif ci_candidates and ss665_rejections:
        messages.append(
            f"The SS665>{float(ss665_min):g} confirmation gate rejected {ss665_rejections:,} of {ci_candidates:,} retained CI candidates."
        )
    finite_shape = ci_valid & family["finite"]

    def shape_quantiles(array):
        values = array[finite_shape & np.isfinite(array)]
        return (
            {str(q): float(np.nanpercentile(values, q)) for q in (1, 5, 25, 50, 75, 95, 99)}
            if values.size
            else {}
        )

    return {
        "created_utc": _m_utils.utc_now(),
        "source_l2": str(l2_path),
        "water_definition_requested": water_definition,
        "water_definition_resolved": context["water_definition_resolved"],
        "shore_buffer_m": float(shore_buffer_m),
        "bright_screen": bool(bright_screen),
        "bright_rhos865": float(bright_rhos865),
        "qa_flags_requested": list(context["requested_flags"]),
        "qa_flags_missing": context["missing_flags"],
        "land_adjacency_pixels": int(land_adjacency_pixels),
        "spectral_gate_settings": {
            "ss681_max": float(ss681_max),
            "ss665_min": float(ss665_min),
            "ci_detection_limit": float(detection_limit),
        },
        "spectral_shape_audit": {
            "SS681_quantiles_retained_finite": shape_quantiles(family["SS_681"]),
            "SS665_quantiles_retained_finite": shape_quantiles(family["SS_665"]),
            "ss665_rejected_ci_candidates": ss665_rejections,
            "ci_threshold_rejected_candidates": threshold_rejections,
        },
        "stages": stages,
        "flag_details": details,
        "messages": messages,
        "rhos_product_mapping": names,
        "ndci_product_mapping": {665: red_name, 709: edge_name},
        "post_mask_counts": {
            "selected_flag": int((water & context["flag_excluded"]).sum()),
            "land_adjacency": int((water & context["adjacent"]).sum()),
            "bright": int((water & context["bright"]).sum()),
        },
    }


# Reference cell 21, lines 186-223.
def mask_diagnostic_tables(report: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return compact stage and per-flag tables for diagnostic display."""
    labels = {
        "geolocated_pixels": "Geolocated L2 pixels",
        "water_pixels": "Water pixels (LAND not set)",
        "qa_or_adjacency_excluded_water": "Water removed by selected QA/buffer",
        "retained_water_pixels": "Water retained after QA/buffer",
        "ci_finite_before_qa": "Finite 4-band CI inputs before QA",
        "ci_finite_after_qa": "Finite 4-band CI inputs after QA",
        "ci_candidates_after_qa": "CI candidates: SS(681) < 0",
        "ci_detections_after_qa": f"CIcyano detections: SS(665)>0 and CI>{float(report.get('ci_detection_limit', 0)):g}",
        "ndci_valid_before_qa": "Valid NDCI inputs before QA",
        "ndci_valid_after_qa": "Valid NDCI inputs after QA",
    }
    water = report["stages"].get("water_pixels", 0)
    stage_rows = []
    for key, label in labels.items():
        count = int(report["stages"].get(key, 0))
        stage_rows.append(
            {
                "stage": label,
                "pixels": count,
                "percent_of_water": np.nan
                if key == "geolocated_pixels"
                else 100.0 * count / water
                if water
                else np.nan,
            }
        )
    flag_rows = []
    for name, details in report["flag_details"].items():
        if details["selected_for_exclusion"] or details["water_pixel_count"]:
            flag_rows.append(
                {
                    "flag": name,
                    "excluded_by_current_mask": bool(details["selected_for_exclusion"]),
                    "flagged_water_pixels": int(details["water_pixel_count"]),
                    "percent_of_water": details["percent_of_water"],
                }
            )
    return (pd.DataFrame(stage_rows), pd.DataFrame(flag_rows))


# Reference cell 21, lines 226-235.
def _create_float_variable(
    dst: netCDF4.Dataset,
    name: str,
    dimensions: tuple[str, ...],
    chunks: tuple[int, ...],
    attributes: dict[str, Any],
) -> netCDF4.Variable:
    variable = dst.createVariable(
        name,
        "f4",
        dimensions,
        zlib=True,
        complevel=_m_settings.NETCDF_DEFLATE,
        shuffle=True,
        chunksizes=chunks,
        fill_value=_m_settings.DERIVED_FILL,
    )
    variable.setncatts(attributes)
    return variable


# Reference cell 21, lines 238-281.
def validate_derived_file(path: Path) -> dict[str, Any]:
    required = (
        "longitude",
        "latitude",
        "water_mask",
        "valid_water_mask",
        "CI",
        "CI_cyano",
        "NDCI",
        "ci_valid_mask",
        "ci_candidate_mask",
        "ci_detection_mask",
        "ndci_valid_mask",
    )
    with netCDF4.Dataset(path) as dataset:
        missing = [name for name in required if name not in dataset.variables]
        if missing:
            raise RuntimeError(f"Derived file lacks variables: {missing}")
        shape = dataset.variables["latitude"].shape
        if any((dataset.variables[name].shape != shape for name in required)):
            raise RuntimeError("Derived variable shapes differ")
        stride = max(1, int(math.sqrt(max(1, np.prod(shape) / 250000))))
        ci_valid = (
            np.ma.asarray(dataset.variables["ci_valid_mask"][::stride, ::stride]).filled(0) == 1
        )
        ci_candidate = (
            np.ma.asarray(dataset.variables["ci_candidate_mask"][::stride, ::stride]).filled(0) == 1
        )
        ci_detected = (
            np.ma.asarray(dataset.variables["ci_detection_mask"][::stride, ::stride]).filled(0) == 1
        )
        ndci_valid = (
            np.ma.asarray(dataset.variables["ndci_valid_mask"][::stride, ::stride]).filled(0) == 1
        )
        ci = _m_science.to_float(dataset.variables["CI"][::stride, ::stride])
        cyano = _m_science.to_float(dataset.variables["CI_cyano"][::stride, ::stride])
        ndci = _m_science.to_float(dataset.variables["NDCI"][::stride, ::stride])
        ci_finite = np.isfinite(ci)
        cyano_finite = np.isfinite(cyano)
        if not np.array_equal(ci_finite, ci_candidate):
            raise RuntimeError("Finite CI pixels do not match ci_candidate_mask")
        if not np.array_equal(cyano_finite, ci_detected):
            raise RuntimeError("Finite CI_cyano pixels do not match ci_detection_mask")
        if np.any(ci_candidate & ~ci_valid) or np.any(ci_detected & ~ci_candidate):
            raise RuntimeError("CI masks are not nested within CI-input-valid water")
        if ndci_valid.any() and (not np.all(np.isfinite(ndci[ndci_valid]))):
            raise RuntimeError("NDCI is not finite on sampled NDCI-valid water pixels")
        return {
            "path": str(path),
            "shape": tuple(map(int, shape)),
            "size_bytes": Path(path).stat().st_size,
            "sample_ci_valid_water": int(ci_valid.sum()),
            "sample_ci_candidates": int(ci_candidate.sum()),
            "sample_ci_cyano_detections": int(ci_detected.sum()),
            "sample_ndci_valid_water": int(ndci_valid.sum()),
        }


# Reference cell 30, lines 266-514.
def build_indices_netcdf(
    l2_path: Path,
    ndci_source: str = None,
    detection_limit: float = _m_settings.CI_DETECTION_LIMIT,
    ss681_max: float = 0.0,
    ss665_min: float = 0.0,
    qa_flags: Iterable[str] = ("NAVFAIL",),
    land_adjacency_pixels: int = _m_settings.LAND_ADJACENCY_PIXELS,
    mask_profile: str | None = "bloom_retaining",
    destination: Path | None = None,
    rows_per_block: int = None,
    water_definition: str = "union",
    lakes: gpd.GeoDataFrame | None = None,
    shore_buffer_m: float = None,
    bright_screen: bool = _m_settings.DEFAULT_BRIGHT_SCREEN,
    bright_rhos865: float = _m_settings.DEFAULT_BRIGHT_RHOS865,
) -> tuple[Path, dict[str, Any]]:
    """Create a native-swath, selected-water index cube using post-L2 masks."""
    if ndci_source is None:
        ndci_source = _m_settings.NDCI_SOURCE
    if rows_per_block is None:
        rows_per_block = _m_settings.ROW_BLOCK_SIZE
    if shore_buffer_m is None:
        shore_buffer_m = _m_settings.DEFAULT_SHORE_BUFFER_M
    l2_path = Path(l2_path).expanduser().resolve()
    prefix = "rhos" if ndci_source.lower() == "rhos" else "Rrs"
    destination = Path(
        destination or _m_settings.DERIVED_DIR / f"{l2_path.stem}_CI_CIcyano_NDCI_water.nc"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    temporary.unlink(missing_ok=True)
    context = _scene_mask_context(
        l2_path,
        qa_flags,
        land_adjacency_pixels,
        water_definition,
        lakes,
        shore_buffer_m,
        bright_screen,
        bright_rhos865,
    )
    summary: dict[str, Any] = {
        "created_utc": _m_utils.utc_now(),
        "source_l2": str(l2_path),
        "ndci_source": prefix,
        "ci_detection_limit": float(detection_limit),
        "ci_detection_operator": ">",
        "ss681_max": float(ss681_max),
        "ss665_min": float(ss665_min),
        "mask_profile": mask_profile or "custom",
        "water_definition_requested": water_definition,
        "water_definition_resolved": context["water_definition_resolved"],
        "shore_buffer_m": float(shore_buffer_m),
        "qa_flags_requested": list(context["requested_flags"]),
        "qa_flags_missing": context["missing_flags"],
        "land_adjacency_pixels": int(land_adjacency_pixels),
        "bright_screen": bool(bright_screen),
        "bright_rhos865_threshold": float(bright_rhos865),
        "bright_rhos865_product": context["bright_product"],
        "water_pixels": int(context["water"].sum()),
        "valid_water_pixels": int(context["valid_water"].sum()),
        "qa_excluded_pixels": int((context["water"] & context["qa_excluded"]).sum()),
        "bright_excluded_pixels": int((context["water"] & context["bright"]).sum()),
        "ci_valid_pixels": 0,
        "ci_candidates": 0,
        "ci_detections": 0,
        "ndci_valid_pixels": 0,
        "negative_core_rhos_pixels": 0,
    }
    with netCDF4.Dataset(l2_path) as src:
        geo = src.groups["geophysical_data"]
        nav = src.groups["navigation_data"]
        variables = list(geo.variables)
        rhos_names = {
            wave: _m_science.nearest_product(variables, "rhos", wave)
            for wave in _m_settings.CORE_RHOS_WAVELENGTHS
        }
        ndci_names = {
            wave: _m_science.nearest_product(variables, prefix, wave) for wave in (665, 709)
        }
        actual_waves = tuple(
            (
                _m_science.product_wavelength(rhos_names[w])
                for w in _m_settings.CORE_RHOS_WAVELENGTHS
            )
        )
        lat_var = nav.variables["latitude"]
        shape, dimensions = (lat_var.shape, lat_var.dimensions)
        chunks = (min(rows_per_block, shape[0]), min(512, shape[1]))
        with netCDF4.Dataset(temporary, "w", format="NETCDF4") as dst:
            for dimension in dimensions:
                src_dim = src.dimensions[dimension]
                dst.createDimension(dimension, None if src_dim.isunlimited() else len(src_dim))
            dst.setncatts(
                {
                    "Conventions": "CF-1.8 ACDD-1.3",
                    "title": "L2Gen-derived OLCI lake-water CI, CIcyano, and NDCI",
                    "source_l2": l2_path.name,
                    "created_utc": _m_utils.utc_now(),
                    "reflectance_for_CI": "L2Gen rhos; Rayleigh-corrected surface reflectance; dimensionless",
                    "ndci_reflectance_source": prefix,
                    "algorithm_scope": "SS681 plus SS665 gates and strict CI>0.0001 detection limit; research reconstruction, not operational NASA CyAN",
                    "ss681_candidate_max": float(ss681_max),
                    "ss665_confirmation_min": float(ss665_min),
                    "ci_detection_limit": float(detection_limit),
                    "ci_detection_operator": ">",
                    "non_detection_semantics": "CI and CI_cyano valid non-detections are fill/NaN; candidate and detection masks record presence/absence",
                    "water_definition": context["water_definition_resolved"],
                    "shore_buffer_m": float(shore_buffer_m),
                    "mask_profile": mask_profile or "custom",
                    "qa_excluded_flags": " ".join(context["requested_flags"]),
                    "land_adjacency_pixels": int(land_adjacency_pixels),
                    "bright_screen_enabled": int(bool(bright_screen)),
                    "bright_rhos865_threshold": float(bright_rhos865),
                    "l2_internal_output_masks": "disabled during v5 L2Gen call; conditions retained in l2_flags",
                    "core_band_mapping_json": json.dumps(rhos_names),
                    "ndci_band_mapping_json": json.dumps(ndci_names),
                }
            )
            lon_out = _create_float_variable(
                dst,
                "longitude",
                dimensions,
                chunks,
                {
                    "standard_name": "longitude",
                    "units": "degrees_east",
                    "coordinates": "latitude longitude",
                },
            )
            lat_out = _create_float_variable(
                dst,
                "latitude",
                dimensions,
                chunks,
                {
                    "standard_name": "latitude",
                    "units": "degrees_north",
                    "coordinates": "latitude longitude",
                },
            )
            flags_out = dst.createVariable(
                "l2_flags",
                "u4",
                dimensions,
                zlib=True,
                complevel=_m_settings.NETCDF_DEFLATE,
                shuffle=True,
                chunksizes=chunks,
            )
            flags_out.setncatts(
                {
                    "long_name": "source L2Gen flags",
                    "flag_masks": np.asarray(list(context["flag_map"].values()), dtype=np.uint32),
                    "flag_meanings": " ".join(context["flag_map"].keys()),
                    "coordinates": "latitude longitude",
                }
            )

            def mask_var(name: str, long_name: str):
                variable = dst.createVariable(
                    name,
                    "u1",
                    dimensions,
                    zlib=True,
                    complevel=_m_settings.NETCDF_DEFLATE,
                    shuffle=True,
                    chunksizes=chunks,
                    fill_value=np.uint8(255),
                )
                variable.setncatts(
                    {
                        "long_name": long_name,
                        "flag_values": np.array([0, 1], dtype=np.uint8),
                        "flag_meanings": "false true",
                        "coordinates": "latitude longitude",
                    }
                )
                return variable

            water_out = mask_var("water_mask", "water under the selected post-L2 definition")
            l2_water_out = mask_var("l2_water_mask", "geolocated inverse-LAND L2 mask")
            hydro_out = mask_var(
                "hydrolakes_water_mask", "HydroLAKES polygon water after shore buffer"
            )
            valid_out = mask_var("valid_water_mask", "selected water retained after post-L2 QA")
            qa_out = mask_var("qa_excluded_mask", "selected water excluded by post-L2 tests")
            adjacent_out = mask_var(
                "land_adjacency_mask", "non-land pixels in the configured LAND dilation"
            )
            bright_out = mask_var(
                "bright_pixel_mask", "pixels exceeding the optional rhos865 threshold"
            )
            ci_valid_out = mask_var(
                "ci_valid_mask", "retained water with four finite CI rhos inputs"
            )
            ci_candidate_out = mask_var("ci_candidate_mask", "negative-SS681 CI candidate")
            ci_detect_out = mask_var(
                "ci_detection_mask",
                "CIcyano detection after SS681, SS665, and CI detection-limit gates",
            )
            ndci_valid_out = mask_var("ndci_valid_mask", "retained water with finite NDCI")

            def reason_var(name: str, long_name: str, meanings: str):
                variable = dst.createVariable(
                    name,
                    "u1",
                    dimensions,
                    zlib=True,
                    complevel=_m_settings.NETCDF_DEFLATE,
                    shuffle=True,
                    chunksizes=chunks,
                    fill_value=np.uint8(255),
                )
                variable.setncatts(
                    {
                        "long_name": long_name,
                        "flag_values": np.arange(len(meanings.split())),
                        "flag_meanings": meanings,
                        "coordinates": "latitude longitude",
                    }
                )
                return variable

            ci_reason_out = reason_var(
                "CI_exclusion_reason",
                "exclusive reason for CIcyano display outcome",
                "not_geolocated outside_selected_water post_mask_excluded nonfinite_core_rhos SS681_gate_failed SS665_gate_failed below_CI_threshold detected",
            )
            ndci_reason_out = reason_var(
                "NDCI_exclusion_reason",
                "exclusive reason for NDCI display outcome",
                "not_geolocated outside_selected_water post_mask_excluded red_nonfinite red_edge_nonfinite denominator_invalid valid",
            )
            negative_out = dst.createVariable(
                "negative_core_rhos_count",
                "u1",
                dimensions,
                zlib=True,
                complevel=_m_settings.NETCDF_DEFLATE,
                shuffle=True,
                chunksizes=chunks,
                fill_value=np.uint8(255),
            )
            negative_out.setncatts(
                {
                    "long_name": "negative values among rhos 620/665/681/709",
                    "valid_range": np.array([0, 4], dtype=np.uint8),
                }
            )
            ss681_out = _create_float_variable(
                dst,
                "SS_681",
                dimensions,
                chunks,
                {"long_name": "spectral shape near 681 nm", "units": "1"},
            )
            ss665_out = _create_float_variable(
                dst,
                "SS_665",
                dimensions,
                chunks,
                {"long_name": "spectral shape near 665 nm", "units": "1"},
            )
            ci_out = _create_float_variable(
                dst,
                "CI",
                dimensions,
                chunks,
                {
                    "long_name": "detected-only CI before SS665 gate",
                    "units": "1",
                    "algorithm": f"-SS_681 where SS_681<{float(ss681_max):g}; non-detections are fill/NaN",
                },
            )
            cyano_out = _create_float_variable(
                dst,
                "CI_cyano",
                dimensions,
                chunks,
                {
                    "long_name": "detected-only two-stage cyanobacteria index",
                    "units": "1",
                    "algorithm": f"CI where SS_681<{float(ss681_max):g}, SS_665>{float(ss665_min):g}, and CI>threshold; standard threshold=0.0001; non-detections are fill/NaN",
                },
            )
            ci_pre_out = _create_float_variable(
                dst,
                "CI_pre_mask",
                dimensions,
                chunks,
                {
                    "long_name": "CI before the selected water and QA masks",
                    "units": "1",
                    "algorithm": "-SS_681 where finite and SS_681<0; all geolocated pixels",
                },
            )
            cyano_pre_out = _create_float_variable(
                dst,
                "CI_cyano_pre_mask",
                dimensions,
                chunks,
                {
                    "long_name": "CIcyano before the selected water and QA masks",
                    "units": "1",
                    "algorithm": "two-stage spectral gates on all finite geolocated pixels",
                },
            )
            ndci_out = _create_float_variable(
                dst,
                "NDCI",
                dimensions,
                chunks,
                {
                    "long_name": "normalized difference chlorophyll index",
                    "units": "1",
                    "formula": f"({prefix}_709-{prefix}_665)/({prefix}_709+{prefix}_665)",
                },
            )
            ndci_pre_out = _create_float_variable(
                dst,
                "NDCI_pre_mask",
                dimensions,
                chunks,
                {
                    "long_name": "NDCI before the selected water and QA masks",
                    "units": "1",
                    "formula": f"({prefix}_709-{prefix}_665)/({prefix}_709+{prefix}_665)",
                },
            )
            geolocated = context["geolocated"]
            lon_out[:] = np.where(geolocated, context["longitude"], _m_settings.DERIVED_FILL)
            lat_out[:] = np.where(geolocated, context["latitude"], _m_settings.DERIVED_FILL)
            flags_out[:] = context["flags"]
            water_out[:] = context["water"].astype(np.uint8)
            l2_water_out[:] = context["l2_water"].astype(np.uint8)
            if context["hydrolakes_water"] is None:
                hydro_out[:] = np.full(shape, np.uint8(255), dtype=np.uint8)
            else:
                hydro_out[:] = context["hydrolakes_water"].astype(np.uint8)
            valid_out[:] = context["valid_water"].astype(np.uint8)
            qa_out[:] = (context["water"] & context["qa_excluded"]).astype(np.uint8)
            adjacent_out[:] = context["adjacent"].astype(np.uint8)
            bright_out[:] = context["bright"].astype(np.uint8)
            for start in tqdm(range(0, shape[0], rows_per_block), desc="v5 derived row blocks"):
                stop = min(shape[0], start + rows_per_block)
                rs = slice(start, stop)
                arrays = {
                    wave: _m_science.to_float(geo.variables[name][rs, :])
                    for wave, name in rhos_names.items()
                }
                family = _m_science.ci_family(
                    arrays[620],
                    arrays[665],
                    arrays[681],
                    arrays[709],
                    wavelengths=actual_waves,
                    detection_limit=detection_limit,
                    ss681_max=ss681_max,
                    ss665_min=ss665_min,
                )
                ci_valid = context["valid_water"][rs, :] & family["finite"]
                candidate = ci_valid & family["ci_candidate"]
                detected = ci_valid & family["ci_detection"]
                negative_count = sum(((array < 0).astype(np.uint8) for array in arrays.values()))
                red = _m_science.to_float(geo.variables[ndci_names[665]][rs, :])
                edge = _m_science.to_float(geo.variables[ndci_names[709]][rs, :])
                ndci, ndci_finite = _m_science.ndci_index(red, edge)
                ndci_valid = context["valid_water"][rs, :] & ndci_finite
                local_geo = context["geolocated"][rs, :]
                local_water = context["water"][rs, :]
                local_qa = context["qa_excluded"][rs, :]
                ci_reason = np.zeros(ci_valid.shape, dtype=np.uint8)
                ci_reason[local_geo & ~local_water] = 1
                ci_reason[local_water & local_qa] = 2
                eligible_ci = local_water & ~local_qa
                ci_reason[eligible_ci & ~family["finite"]] = 3
                ci_reason[
                    eligible_ci & family["finite"] & ~(family["SS_681"] < float(ss681_max))
                ] = 4
                ci_reason[
                    eligible_ci & family["ci_candidate"] & ~(family["SS_665"] > float(ss665_min))
                ] = 5
                ci_reason[
                    eligible_ci
                    & family["ci_candidate"]
                    & (family["SS_665"] > float(ss665_min))
                    & ~(family["CI"] > float(detection_limit))
                ] = 6
                ci_reason[eligible_ci & family["ci_detection"]] = 7
                ndci_reason = np.zeros(ndci_valid.shape, dtype=np.uint8)
                ndci_reason[local_geo & ~local_water] = 1
                ndci_reason[local_water & local_qa] = 2
                eligible_ndci = local_water & ~local_qa
                ndci_reason[eligible_ndci & ~np.isfinite(red)] = 3
                ndci_reason[eligible_ndci & np.isfinite(red) & ~np.isfinite(edge)] = 4
                ndci_reason[
                    eligible_ndci
                    & np.isfinite(red)
                    & np.isfinite(edge)
                    & ~(np.abs(red + edge) > 1e-08)
                ] = 5
                ndci_reason[eligible_ndci & ndci_finite] = 6
                ci_valid_out[rs, :] = ci_valid.astype(np.uint8)
                ci_candidate_out[rs, :] = np.where(
                    ci_valid, candidate.astype(np.uint8), np.uint8(255)
                )
                ci_detect_out[rs, :] = np.where(ci_valid, detected.astype(np.uint8), np.uint8(255))
                ndci_valid_out[rs, :] = ndci_valid.astype(np.uint8)
                ci_reason_out[rs, :] = ci_reason
                ndci_reason_out[rs, :] = ndci_reason
                negative_out[rs, :] = np.where(ci_valid, negative_count, np.uint8(255))
                ss681_out[rs, :] = np.where(ci_valid, family["SS_681"], _m_settings.DERIVED_FILL)
                ss665_out[rs, :] = np.where(ci_valid, family["SS_665"], _m_settings.DERIVED_FILL)
                ci_out[rs, :] = np.where(candidate, family["CI"], _m_settings.DERIVED_FILL)
                cyano_out[rs, :] = np.where(detected, family["CI_cyano"], _m_settings.DERIVED_FILL)
                ci_pre_out[rs, :] = np.where(
                    local_geo & family["ci_candidate"], family["CI"], _m_settings.DERIVED_FILL
                )
                cyano_pre_out[rs, :] = np.where(
                    local_geo & family["ci_detection"], family["CI_cyano"], _m_settings.DERIVED_FILL
                )
                ndci_out[rs, :] = np.where(ndci_valid, ndci, _m_settings.DERIVED_FILL)
                ndci_pre_out[rs, :] = np.where(
                    local_geo & ndci_finite, ndci, _m_settings.DERIVED_FILL
                )
                summary["ci_valid_pixels"] += int(ci_valid.sum())
                summary["ci_candidates"] += int(candidate.sum())
                summary["ci_detections"] += int(detected.sum())
                summary["ndci_valid_pixels"] += int(ndci_valid.sum())
                summary["negative_core_rhos_pixels"] += int((ci_valid & (negative_count > 0)).sum())
    temporary.replace(destination)
    summary.update(validate_derived_file(destination))
    summary["size_human"] = _m_utils.human_size(destination.stat().st_size)
    summary["ci_detection_fraction_ci_valid"] = (
        summary["ci_detections"] / summary["ci_valid_pixels"]
        if summary["ci_valid_pixels"]
        else None
    )
    _m_utils.write_json(
        _m_settings.EXPORT_DIR / l2_path.stem / "derived_index_summary.json", summary
    )
    print(f"Derived product: {destination} ({summary['size_human']})")
    return (destination, summary)


# Reference cell 30, lines 24-129.
def _scene_mask_context(
    l2_path: Path,
    qa_flags: Iterable[str],
    land_adjacency_pixels: int,
    water_definition: str,
    lakes: gpd.GeoDataFrame | None,
    shore_buffer_m: float,
    bright_screen: bool,
    bright_rhos865: float,
) -> dict[str, Any]:
    """Build masks without changing the native L2Gen file."""
    water_definition = str(water_definition).lower()
    with netCDF4.Dataset(l2_path) as src:
        geo = src.groups["geophysical_data"]
        nav = src.groups["navigation_data"]
        latitude = _m_science.to_float(nav.variables["latitude"][:])
        longitude = _m_science.to_float(nav.variables["longitude"][:])
        flags_var = _m_science.locate_l2_flags(src)
        flags_ma = np.ma.asarray(flags_var[:])
        flags_valid = ~np.ma.getmaskarray(flags_ma)
        flags = flags_ma.filled(0).astype(np.uint32)
        flag_map = _m_science.decode_flag_metadata(flags_var)
        if "LAND" not in flag_map:
            raise RuntimeError("l2_flags metadata lacks LAND")
        geolocated = (
            np.isfinite(latitude)
            & np.isfinite(longitude)
            & (np.abs(latitude) <= 90)
            & (np.abs(longitude) <= 180)
            & flags_valid
        )
        land = _m_science.flag_is_set(flags, flag_map, "LAND")
        l2_water = geolocated & ~land
        hydro_water = _m_hydrolakes.hydrolakes_point_mask(
            latitude, longitude, lakes, shore_buffer_m
        )
        if water_definition == "l2_land":
            water = l2_water
            resolved_water_definition = "inverse L2 LAND flag"
        elif water_definition == "hydrolakes":
            if hydro_water is None:
                raise ValueError("HydroLAKES water was requested but no scene lakes are loaded")
            water = geolocated & hydro_water
            resolved_water_definition = "HydroLAKES polygons"
        elif water_definition == "intersection":
            if hydro_water is None:
                raise ValueError("Intersection water was requested but no scene lakes are loaded")
            water = l2_water & hydro_water
            resolved_water_definition = "HydroLAKES ∩ inverse L2 LAND"
        elif water_definition == "union":
            if hydro_water is None:
                water = l2_water
                resolved_water_definition = "inverse L2 LAND flag (HydroLAKES not loaded)"
            else:
                water = l2_water | geolocated & hydro_water
                resolved_water_definition = "HydroLAKES ∪ inverse L2 LAND"
        elif water_definition == "hydrolakes_if_loaded":
            if hydro_water is None:
                water = l2_water
                resolved_water_definition = "inverse L2 LAND flag (HydroLAKES not loaded)"
            else:
                water = geolocated & hydro_water
                resolved_water_definition = "HydroLAKES polygons"
        else:
            raise ValueError(
                "water_definition must be l2_land, hydrolakes, intersection, union, or hydrolakes_if_loaded"
            )
        requested = tuple(dict.fromkeys((str(name).upper() for name in qa_flags)))
        flag_excluded = np.zeros(latitude.shape, dtype=bool)
        missing_flags = []
        for name in requested:
            if name not in flag_map:
                missing_flags.append(name)
            flag_excluded |= _m_science.flag_is_set(flags, flag_map, name)
        if int(land_adjacency_pixels) > 0:
            adjacent = binary_dilation(land, iterations=int(land_adjacency_pixels)) & ~land
        else:
            adjacent = np.zeros(latitude.shape, dtype=bool)
        bright = np.zeros(latitude.shape, dtype=bool)
        bright_name = None
        if bright_screen:
            try:
                bright_name = _m_science.nearest_product(
                    geo.variables, "rhos", 865, tolerance_nm=25
                )
                rhos865 = _m_science.to_float(geo.variables[bright_name][:])
                bright = np.isfinite(rhos865) & (rhos865 > float(bright_rhos865))
            except KeyError as error:
                raise RuntimeError(
                    "Bright-pixel screen requested, but no rhos band near 865 nm exists"
                ) from error
        qa_excluded = flag_excluded | adjacent | bright
        valid_water = water & ~qa_excluded
        return {
            "latitude": latitude,
            "longitude": longitude,
            "geolocated": geolocated,
            "flags": flags,
            "flag_map": flag_map,
            "land": land,
            "l2_water": l2_water,
            "hydrolakes_water": hydro_water,
            "water": water,
            "valid_water": valid_water,
            "flag_excluded": flag_excluded,
            "adjacent": adjacent,
            "bright": bright,
            "qa_excluded": qa_excluded,
            "requested_flags": requested,
            "missing_flags": missing_flags,
            "water_definition_requested": water_definition,
            "water_definition_resolved": resolved_water_definition,
            "shore_buffer_m": float(shore_buffer_m),
            "bright_screen": bool(bright_screen),
            "bright_rhos865": float(bright_rhos865),
            "bright_product": bright_name,
        }
