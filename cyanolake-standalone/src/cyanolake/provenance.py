"""Standalone provenance implementation."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from . import l2gen as _m_l2gen
from . import settings as _m_settings
from . import utils as _m_utils


# Reference cell 27, lines 1-47.
def write_scene_manifest(
    scene: dict[str, Any] | None,
    raw_zip: Path | None,
    sen3_path: Path,
    l2_path: Path,
    derived_path: Path,
    l2_info: dict[str, Any],
    derived_info: dict[str, Any],
) -> Path:
    scene_name = _m_l2gen.product_basename(sen3_path)
    payload = {
        "created_utc": _m_utils.utc_now(),
        "scene_catalogue_record": scene,
        "files": {
            "raw_l1_zip": str(raw_zip) if raw_zip else None,
            "staged_sen3": str(sen3_path),
            "l2gen_netcdf": str(l2_path),
            "derived_indices_netcdf": str(derived_path),
        },
        "sizes_bytes": {
            "raw_l1_zip": raw_zip.stat().st_size if raw_zip and raw_zip.exists() else None,
            "l2gen_netcdf": l2_path.stat().st_size,
            "derived_indices_netcdf": derived_path.stat().st_size,
        },
        "checksums_sha256": {
            "raw_l1_zip": _m_utils.sha256_file(raw_zip)
            if _m_settings.COMPUTE_SHA256 and raw_zip and raw_zip.exists()
            else None,
            "l2gen_netcdf": _m_utils.sha256_file(l2_path) if _m_settings.COMPUTE_SHA256 else None,
            "derived_indices_netcdf": _m_utils.sha256_file(derived_path)
            if _m_settings.COMPUTE_SHA256
            else None,
        },
        "configuration": {
            "l2_mode": l2_info.get("mode"),
            "ndci_source": derived_info.get("ndci_source"),
            "process_full_scene": l2_info.get("processing_scope") == "full_scene",
            "processing_bbox": l2_info.get("processing_bbox"),
            "full_reflectance_spectra_requested": l2_info.get("full_reflectance_spectra_requested"),
            "get_ancillary": _m_settings.GET_ANCILLARY,
            "allow_climatology_fallback": _m_settings.ALLOW_CLIMATOLOGY_FALLBACK,
            "ci_detection_limit": derived_info.get("ci_detection_limit"),
            "mask_profile": derived_info.get("mask_profile"),
            "qa_exclude_flags": derived_info.get("qa_flags_requested"),
            "land_adjacency_pixels": derived_info.get("land_adjacency_pixels"),
        },
        "l2gen": l2_info,
        "derived_summary": derived_info,
    }
    path = _m_settings.EXPORT_DIR / scene_name / "processing_manifest.json"
    return _m_utils.write_json(path, payload)
