"""Standalone state implementation."""

from __future__ import annotations

_S3_CLIENT_CACHE = {}

GLOBAL_STATE = {
    "lakes": None,
    "config_hash": None,
    "min_area_km2": 10.0,
    "region": "WORLD",
    "planned_scene_ids": [],
    "planning_shard": "",
    "planning_geometry": None,
}

ADM1_STATE = {"country_iso3": None, "gdf": None, "lookup": {}, "labels": {}, "metadata": {}}

CYAN_STATE = {"source_path": None, "source_date": None}
