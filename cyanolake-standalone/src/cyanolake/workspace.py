"""Canonical lake identity and reproducible configuration selection."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

from . import configuration, settings, state


def lake_fingerprint(lakes) -> str:
    digest = hashlib.sha256()
    for row in lakes.sort_values("Hylak_id").itertuples():
        digest.update(str(int(row.Hylak_id)).encode())
        digest.update(row.geometry.wkb)
    return digest.hexdigest()


def select_configuration(lakes, expected_hash: str | None = None) -> str:
    from .registry import processing_config_hash

    settings.TARGET_LAKES_FINGERPRINT = lake_fingerprint(lakes)
    ch = processing_config_hash(state.GLOBAL_STATE["min_area_km2"], state.GLOBAL_STATE["region"])
    if expected_hash and expected_hash != ch:
        raise ValueError(
            f"Configuration/lake universe differs from {expected_hash}; current hash is {ch}. Restore its saved config and target lakes or create a new plan."
        )
    state.GLOBAL_STATE.update(lakes=lakes, config_hash=ch)
    # Cache products from different science configurations cannot overwrite each other.
    settings.EXPORT_DIR = settings.PERSISTENT_ROOT / "exports" / ch
    settings.HTML_DIR = settings.PERSISTENT_ROOT / "html_maps" / ch
    settings.L2_DIR = settings.SCRATCH_ROOT / "l2" / ch
    settings.PAR_DIR = settings.PERSISTENT_ROOT / "par" / ch
    for folder in (settings.EXPORT_DIR, settings.HTML_DIR, settings.L2_DIR, settings.PAR_DIR):
        folder.mkdir(parents=True, exist_ok=True)
    folder = settings.CATALOG_DIR / "configurations" / ch
    folder.mkdir(parents=True, exist_ok=True)
    snapshot = folder / "configuration.json"
    if not snapshot.exists():
        snapshot.write_text(
            json.dumps(
                {
                    "config_hash": ch,
                    "configuration": configuration.CURRENT,
                    "lake_geometry_sha256": settings.TARGET_LAKES_FINGERPRINT,
                },
                indent=2,
            )
            + "\n"
        )
    target = folder / "target_hydrolakes.gpkg"
    if settings.TARGET_LAKES_GPKG.exists() and not target.exists():
        shutil.copy2(settings.TARGET_LAKES_GPKG, target)
    return ch


def load_lakes(expected_hash: str | None = None):
    import geopandas as gpd

    from .hydrolakes import _normalize_hydrolakes

    if not settings.TARGET_LAKES_GPKG.exists():
        raise FileNotFoundError("Run prepare-lakes before planning or processing")
    lakes = _normalize_hydrolakes(gpd.read_file(settings.TARGET_LAKES_GPKG))
    if lakes.empty:
        raise ValueError("The target HydroLAKES universe is empty")
    ch = select_configuration(lakes, expected_hash)
    return lakes, ch
