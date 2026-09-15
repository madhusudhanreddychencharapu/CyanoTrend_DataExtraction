"""Standalone registry implementation."""

from __future__ import annotations

import hashlib
import json
import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd

from . import settings as _m_settings
from . import utils as _m_utils


# Reference cell 36, lines 4-35.
def _registry_connection():
    conn = sqlite3.connect(_m_settings.REGISTRY_DB, timeout=60)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute(
        "\n        CREATE TABLE IF NOT EXISTS scenes (\n            scene_id TEXT NOT NULL,\n            config_hash TEXT NOT NULL,\n            scene_name TEXT,\n            acquisition_start TEXT,\n            acquisition_end TEXT,\n            s3_path TEXT,\n            size_bytes INTEGER,\n            candidate_lakes INTEGER,\n            status TEXT NOT NULL DEFAULT 'queued',\n            attempts INTEGER NOT NULL DEFAULT 0,\n            processing_strategy TEXT,\n            window_count INTEGER,\n            download_method TEXT,\n            compact_bytes INTEGER,\n            stats_rows INTEGER,\n            last_error TEXT,\n            updated_utc TEXT,\n            scene_json TEXT,\n            PRIMARY KEY(scene_id, config_hash)\n        )\n        "
    )
    conn.commit()
    return conn


# Reference cell 57, lines 49-76.
def processing_config_hash(min_area_km2: float, region: str = "WORLD") -> str:
    payload = {
        "l2_mode": _m_settings.L2_MODE,
        "ndci_source": _m_settings.NDCI_SOURCE,
        "analysis_rhos": list(_m_settings.GLOBAL_ANALYSIS_RHOS_WAVELENGTHS),
        "true_color_rhos": list(_m_settings.TRUE_COLOR_RHOS_WAVELENGTHS),
        "indices": ["CI", "CI_cyano", "NDCI", *_m_settings.ADDITIONAL_INDEX_NAMES],
        "full_spectrum": bool(_m_settings.SAVE_FULL_REFLECTANCE_SPECTRA),
        "primary_qa_profile": _m_settings.DEFAULT_MASK_PROFILE,
        "primary_hard_flags": list(_m_settings.BLOOM_AWARE_HARD_FLAGS),
        "cyan_strict_flags": list(_m_settings.CYAN_STRICT_QA_FLAGS),
        "require_l2_nonland": True,
        "bloom_rescue_require_ci_cyano": bool(_m_settings.BLOOM_RESCUE_REQUIRE_CI_CYANO),
        "bloom_rescue_require_positive_afai": bool(_m_settings.BLOOM_RESCUE_REQUIRE_POSITIVE_AFAI),
        "bloom_rescue_afai_min": float(_m_settings.BLOOM_RESCUE_AFAI_MIN),
        "bloom_rescue_require_nir_over_blue": bool(_m_settings.BLOOM_RESCUE_REQUIRE_NIR_OVER_BLUE),
        "detection_limit": float(_m_settings.CI_DETECTION_LIMIT),
        "bright_screen": bool(_m_settings.DEFAULT_BRIGHT_SCREEN),
        "shore_buffer_m": float(_m_settings.DEFAULT_SHORE_BUFFER_M),
        "processing_strategy": _m_settings.PROCESSING_STRATEGY,
        "window_margin_km": float(_m_settings.L2_WINDOW_MARGIN_KM),
        "max_windows": int(_m_settings.MAX_L2_WINDOWS_PER_SCENE),
        "full_scene_fraction": float(_m_settings.FULL_SCENE_FRACTION_THRESHOLD),
        "min_lake_area_km2": float(min_area_km2),
        "region": str(region).upper(),
        "schema": _m_settings.V264_SCHEMA,
    }
    payload.update(
        standalone_schema="cyanolake-2.6.8.1",
        get_ancillary=bool(_m_settings.GET_ANCILLARY),
        allow_climatology_fallback=bool(_m_settings.ALLOW_CLIMATOLOGY_FALLBACK),
        lake_geometry_sha256=getattr(_m_settings, "TARGET_LAKES_FINGERPRINT", ""),
    )
    return hashlib.sha1(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


# Reference cell 36, lines 63-87.
def registry_upsert_scene(scene: dict[str, Any], config_hash: str, candidate_lakes: int) -> None:
    with _registry_connection() as conn:
        conn.execute(
            "\n            INSERT INTO scenes (\n                scene_id, config_hash, scene_name, acquisition_start, acquisition_end,\n                s3_path, size_bytes, candidate_lakes, status, updated_utc, scene_json\n            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, 'queued', ?, ?)\n            ON CONFLICT(scene_id, config_hash) DO UPDATE SET\n                scene_name=excluded.scene_name,\n                acquisition_start=excluded.acquisition_start,\n                acquisition_end=excluded.acquisition_end,\n                s3_path=excluded.s3_path,\n                size_bytes=excluded.size_bytes,\n                candidate_lakes=excluded.candidate_lakes,\n                updated_utc=excluded.updated_utc,\n                scene_json=excluded.scene_json\n            ",
            (
                str(scene["id"]),
                config_hash,
                scene.get("name"),
                scene.get("start"),
                scene.get("end"),
                scene.get("s3_path"),
                int(scene.get("size_bytes") or 0),
                int(candidate_lakes),
                _m_utils.utc_now(),
                json.dumps(scene, default=str),
            ),
        )
        conn.commit()


# Reference cell 36, lines 90-99.
def registry_update(scene_id: str, config_hash: str, **fields) -> None:
    if not fields:
        return
    fields["updated_utc"] = _m_utils.utc_now()
    keys = list(fields)
    sql = ", ".join((f"{k}=?" for k in keys))
    values = [fields[k] for k in keys] + [scene_id, config_hash]
    with _registry_connection() as conn:
        conn.execute(f"UPDATE scenes SET {sql} WHERE scene_id=? AND config_hash=?", values)
        conn.commit()


# Reference cell 36, lines 102-109.
def registry_table(config_hash: str | None = None) -> pd.DataFrame:
    with _registry_connection() as conn:
        if config_hash:
            return pd.read_sql_query(
                "SELECT * FROM scenes WHERE config_hash=? ORDER BY acquisition_start",
                conn,
                params=(config_hash,),
            )
        return pd.read_sql_query("SELECT * FROM scenes ORDER BY acquisition_start", conn)


# Reference cell 36, lines 112-115.
def registry_export_csv() -> Path:
    out = _m_settings.EXPORT_DIR / "scene_registry.csv"
    registry_table().to_csv(out, index=False)
    return out


# Reference cell 48, lines 116-127.
def parse_s3_frame_metadata(name: str) -> dict[str, Any]:
    n = Path(str(name)).name.removesuffix(".SEN3").removesuffix(".zip")
    m = _m_settings.S3_FRAME_RE.match(n)
    if not m:
        raise ValueError(f"Cannot parse Sentinel-3 OLCI frame fields from {n}")
    platform, start, stop, created, duration, cycle, relative_orbit, frame = m.groups()
    return {
        "scene_name": n,
        "platform": platform,
        "sensing_start": start,
        "sensing_stop": stop,
        "creation_time": created,
        "duration_s": int(duration),
        "cycle": int(cycle),
        "relative_orbit": int(relative_orbit),
        "frame": int(frame),
        "native_frame_id": f"{platform}_R{int(relative_orbit):03d}_F{int(frame):04d}",
    }
