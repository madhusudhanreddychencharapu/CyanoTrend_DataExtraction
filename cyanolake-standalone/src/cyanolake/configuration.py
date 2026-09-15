"""Validated TOML configuration and explicit workspace creation."""

from __future__ import annotations

import copy
import math
import tomllib
from pathlib import Path

from . import settings, state

DEFAULTS = {
    "paths": {"persistent_root": "../data", "scratch_root": "../scratch", "ocssw_root": "~/ocssw"},
    "processing": {
        "l2_mode": "full_swath_rhos",
        "ndci_source": "rhos",
        "get_ancillary": True,
        "allow_climatology_fallback": False,
        "strategy": "adaptive_lake_windows",
        "window_margin_km": 3.0,
        "max_windows": 6,
        "full_scene_fraction_threshold": 0.35,
        "shore_buffer_m": 0.0,
        "rows_per_block": 256,
    },
    "storage": {
        "download_method": "auto",
        "keep_staged_l1": False,
        "keep_raw_zip": False,
        "keep_l2": False,
        "keep_failed_scratch": True,
    },
    "lakes": {"min_area_km2": 10.0, "region": "WORLD"},
}
CURRENT = None


def configure(path: str | Path, *, create: bool = False) -> dict:
    """Use one configured workspace per process. Imports never create directories."""
    global CURRENT
    path = Path(path).expanduser().resolve()
    with path.open("rb") as handle:
        supplied = tomllib.load(handle)
    config = copy.deepcopy(DEFAULTS)
    for section, values in supplied.items():
        if section not in config or not isinstance(values, dict):
            raise ValueError(f"Unknown configuration section: {section}")
        for key, value in values.items():
            if key not in config[section]:
                raise ValueError(f"Unknown configuration key: {section}.{key}")
            expected = config[section][key]
            if isinstance(expected, bool) and not isinstance(value, bool):
                raise ValueError(f"{section}.{key} must be true or false")
            if isinstance(expected, str) and not isinstance(value, str):
                raise ValueError(f"{section}.{key} must be text")
            if isinstance(expected, (int, float)) and not isinstance(expected, bool):
                if (
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not math.isfinite(value)
                ):
                    raise ValueError(f"{section}.{key} must be a finite number")
            config[section][key] = value
    p, s, lakes = config["processing"], config["storage"], config["lakes"]
    if p["l2_mode"] not in ("full_swath_rhos", "standard_dual"):
        raise ValueError("l2_mode must be full_swath_rhos or standard_dual")
    if p["ndci_source"] not in ("rhos", "Rrs") or (
        p["ndci_source"] == "Rrs" and p["l2_mode"] != "standard_dual"
    ):
        raise ValueError("Rrs-based NDCI requires l2_mode = standard_dual")
    if p["strategy"] not in ("adaptive_lake_windows", "full_scene"):
        raise ValueError("strategy must be adaptive_lake_windows or full_scene")
    if s["download_method"] not in ("auto", "s3", "zip"):
        raise ValueError("download_method must be auto, s3 or zip")
    for name in ("max_windows", "rows_per_block"):
        if p[name] < 1 or int(p[name]) != p[name]:
            raise ValueError(f"{name} must be a positive integer")
    if not 0 < p["full_scene_fraction_threshold"] <= 1:
        raise ValueError("full_scene_fraction_threshold must be in (0, 1]")
    if min(p["window_margin_km"], p["shore_buffer_m"], lakes["min_area_km2"]) < 0:
        raise ValueError("Margins and minimum lake area must be nonnegative")
    roots = {}
    for key, value in config["paths"].items():
        root = Path(value).expanduser()
        roots[key] = (path.parent / root).resolve() if not root.is_absolute() else root.resolve()
        if roots[key] in (Path("/"), Path.home().resolve()):
            raise ValueError(f"{key} must be a dedicated directory")
    persistent, scratch, ocssw = (
        roots[k] for k in ("persistent_root", "scratch_root", "ocssw_root")
    )
    for a, b in ((persistent, scratch), (scratch, ocssw), (persistent, ocssw)):
        if a == b or a in b.parents or b in a.parents:
            raise ValueError(
                "Persistent, scratch and OCSSW roots must be separate, non-nested directories"
            )
    settings.PROJECT_ROOT = settings.PERSISTENT_ROOT = persistent
    settings.SCRATCH_ROOT, settings.OCSSWROOT = scratch, ocssw
    settings.OCSSW_INSTALLER_DIR = ocssw.parent / (ocssw.name + "-installer")
    layout = {
        "EXPORT_DIR": persistent / "exports",
        "PAR_DIR": persistent / "par",
        "LOG_DIR": persistent / "logs",
        "HYDRO_ROOT": persistent / "hydrolakes",
        "COMPACT_DIR": persistent / "lake_only",
        "COMPACT_NC_DIR": persistent / "lake_only/netcdf",
        "STATS_SCENE_DIR": persistent / "lake_only/stats_by_scene",
        "CATALOG_DIR": persistent / "catalog",
        "HTML_DIR": persistent / "html_maps",
        "ADMIN_ROOT": persistent / "admin_boundaries",
        "RAW_L1_DIR": scratch / "raw_l1_fallback",
        "STAGING_DIR": scratch / "staging",
        "L2_DIR": scratch / "l2",
        "DERIVED_DIR": scratch / "legacy_dense_derived",
    }
    for key, value in layout.items():
        setattr(settings, key, value)
    settings.REGISTRY_DB = persistent / "scene_registry.sqlite"
    settings.TARGET_LAKES_GPKG = settings.CATALOG_DIR / "target_hydrolakes.gpkg"
    settings.ARCHIVE_INDEX_CSV = settings.CATALOG_DIR / "global_archive_index.csv"
    mapping = {
        "l2_mode": "L2_MODE",
        "ndci_source": "NDCI_SOURCE",
        "get_ancillary": "GET_ANCILLARY",
        "allow_climatology_fallback": "ALLOW_CLIMATOLOGY_FALLBACK",
        "strategy": "PROCESSING_STRATEGY",
        "window_margin_km": "L2_WINDOW_MARGIN_KM",
        "max_windows": "MAX_L2_WINDOWS_PER_SCENE",
        "full_scene_fraction_threshold": "FULL_SCENE_FRACTION_THRESHOLD",
        "shore_buffer_m": "DEFAULT_SHORE_BUFFER_M",
        "rows_per_block": "ROW_BLOCK_SIZE",
    }
    for key, name in mapping.items():
        setattr(settings, name, p[key])
    mapping = {
        "download_method": "DOWNLOAD_METHOD",
        "keep_staged_l1": "KEEP_STAGED_L1",
        "keep_raw_zip": "KEEP_RAW_ZIP_FALLBACK",
        "keep_l2": "KEEP_L2_INTERMEDIATE",
        "keep_failed_scratch": "KEEP_FAILED_SCRATCH",
    }
    for key, name in mapping.items():
        setattr(settings, name, s[key])
    state.GLOBAL_STATE.update(
        lakes=None,
        config_hash=None,
        min_area_km2=lakes["min_area_km2"],
        region=lakes["region"].upper(),
        planned_scene_ids=[],
        planning_shard="",
        planning_geometry=None,
    )
    state.CYAN_STATE.update(source_path=None, source_date=None)
    state.ADM1_STATE.update(country_iso3=None, gdf=None, lookup={}, metadata={})
    settings.TARGET_LAKES_FINGERPRINT = ""
    settings.L2GEN_BIN = settings.GETANC_BIN = None
    settings.OCSSW_SUBPROCESS_ENV = {}
    if create:
        for directory in (persistent, scratch, *layout.values()):
            directory.mkdir(parents=True, exist_ok=True)
    config["paths"] = {k: str(v) for k, v in roots.items()}
    CURRENT = config
    return copy.deepcopy(config)
