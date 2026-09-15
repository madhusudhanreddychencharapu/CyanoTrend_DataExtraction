"""Standalone windows implementation."""

from __future__ import annotations

import hashlib
import math
from typing import Any

import geopandas as gpd
import numpy as np
from shapely.geometry import box as shapely_box
from shapely.ops import unary_union

from . import catalogue as _m_catalogue
from . import settings as _m_settings


# Reference cell 42, lines 4-10.
def _expanded_lonlat_box(bounds, margin_km: float):
    west, south, east, north = map(float, bounds)
    lat0 = float(np.clip((south + north) / 2.0, -85.0, 85.0))
    dlat = float(margin_km) / 111.32
    dlon = float(margin_km) / max(111.32 * math.cos(math.radians(lat0)), 1.0)
    return shapely_box(
        max(-180, west - dlon), max(-90, south - dlat), min(180, east + dlon), min(90, north + dlat)
    )


# Reference cell 42, lines 13-45.
def plan_scene_windows(
    scene: dict[str, Any],
    scene_lakes: gpd.GeoDataFrame,
    strategy: str = None,
    margin_km: float = None,
    max_windows: int = None,
    full_scene_fraction: float = None,
) -> list[tuple[float, float, float, float] | None]:
    if strategy is None:
        strategy = _m_settings.PROCESSING_STRATEGY
    if margin_km is None:
        margin_km = _m_settings.L2_WINDOW_MARGIN_KM
    if max_windows is None:
        max_windows = _m_settings.MAX_L2_WINDOWS_PER_SCENE
    if full_scene_fraction is None:
        full_scene_fraction = _m_settings.FULL_SCENE_FRACTION_THRESHOLD
    if str(strategy).lower() == "full_scene":
        return [None]
    geom = _m_catalogue.scene_geometry(scene)
    if geom is None or scene_lakes is None or scene_lakes.empty:
        return [None]
    minx, miny, maxx, maxy = geom.bounds
    if maxx - minx > 170:
        return [None]
    boxes = [_expanded_lonlat_box(g.bounds, margin_km) for g in scene_lakes.geometry]
    merged = unary_union(boxes)
    parts = list(merged.geoms) if hasattr(merged, "geoms") else [merged]
    windows = []
    for part in parts:
        clipped = part.intersection(geom)
        if clipped.is_empty:
            continue
        west, south, east, north = clipped.bounds
        if west < east and south < north:
            windows.append((float(west), float(south), float(east), float(north)))
    if not windows or len(windows) > int(max_windows):
        return [None]
    scene_area = max((maxx - minx) * (maxy - miny), 1e-09)
    window_area = sum((max(0, e - w) * max(0, n - s) for w, s, e, n in windows))
    if window_area / scene_area >= float(full_scene_fraction):
        return [None]
    return sorted(windows, key=lambda b: (b[1], b[0]))


# Reference cell 42, lines 48-52.
def window_id(window) -> str:
    if window is None:
        return "FULL"
    text = ",".join((f"{v:.6f}" for v in window))
    return "W_" + hashlib.sha1(text.encode()).hexdigest()[:10]


# Reference cell 42, lines 55-59.
def lakes_for_window(lakes: gpd.GeoDataFrame, window) -> gpd.GeoDataFrame:
    if window is None:
        return lakes.copy().reset_index(drop=True)
    geom = shapely_box(*window)
    return lakes[lakes.intersects(geom)].copy().reset_index(drop=True)
