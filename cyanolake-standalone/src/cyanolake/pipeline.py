"""Standalone pipeline implementation."""

from __future__ import annotations

import contextlib
import json
import shutil
from pathlib import Path
from typing import Any

import geopandas as gpd
import pandas as pd

from . import catalogue as _m_catalogue
from . import compact as _m_compact
from . import l2gen as _m_l2gen
from . import registry as _m_registry
from . import settings as _m_settings
from . import staging as _m_staging
from . import statistics as _m_statistics
from . import utils as _m_utils
from . import windows as _m_windows


def scene_outputs_valid(scene_id, config_hash):
    path = _m_settings.EXPORT_DIR / scene_id / f"manifest_{config_hash}.json"
    try:
        data = json.loads(path.read_text())
        if data["config_hash"] != config_hash or str(data["scene"]["id"]) != str(scene_id):
            return False
        return all(
            _m_compact.validate_compact_product(p, scene_id=scene_id, config_hash=config_hash)
            for p in data["compact_files"]
        )
    except (OSError, KeyError, ValueError):
        return False


# Reference cell 48, lines 129-139.
def _scene_output_dir(scene: dict[str, Any]) -> Path:
    """Physical persistent layout: platform / relative orbit / frame / date / scene ID."""
    stamp = pd.Timestamp(scene.get("start") or "1970-01-01")
    try:
        meta = _m_registry.parse_s3_frame_metadata(scene.get("name", ""))
        return (
            _m_settings.COMPACT_NC_DIR
            / f"platform={meta['platform']}"
            / f"relative_orbit={meta['relative_orbit']:03d}"
            / f"frame={meta['frame']:04d}"
            / f"year={stamp.year:04d}"
            / f"month={stamp.month:02d}"
            / f"day={stamp.day:02d}"
            / str(scene["id"])
        )
    except Exception:
        return (
            _m_settings.COMPACT_NC_DIR
            / "platform=UNKNOWN"
            / f"year={stamp.year:04d}"
            / f"month={stamp.month:02d}"
            / str(scene["id"])
        )


# Reference cell 46, lines 9-23.
def _cleanup_scene_scratch(
    raw_zip: Path | None, sen3: Path | None, l2_paths: list[Path], success=True
):
    if not success and _m_settings.KEEP_FAILED_SCRATCH:
        return
    if not _m_settings.KEEP_L2_INTERMEDIATE:
        for p in l2_paths:
            Path(p).unlink(missing_ok=True)
    if (
        sen3 is not None
        and (not _m_settings.KEEP_STAGED_L1)
        and Path(sen3).resolve().is_relative_to(_m_settings.STAGING_DIR.resolve())
    ):
        shutil.rmtree(Path(sen3), ignore_errors=True)
        parent = Path(sen3).parent
        if parent != _m_settings.STAGING_DIR and parent.exists():
            with contextlib.suppress(OSError):
                parent.rmdir()
    if (
        raw_zip is not None
        and (not _m_settings.KEEP_RAW_ZIP_FALLBACK)
        and Path(raw_zip).resolve().is_relative_to(_m_settings.RAW_L1_DIR.resolve())
    ):
        Path(raw_zip).unlink(missing_ok=True)


# Reference cell 46, lines 26-127.
def process_scene_once(
    scene: dict[str, Any],
    config_hash: str,
    lakes_all: gpd.GeoDataFrame,
    username: str,
    password: str,
    download_method: str = None,
    strategy: str = None,
) -> pd.DataFrame:
    if download_method is None:
        download_method = _m_settings.DOWNLOAD_METHOD
    if strategy is None:
        strategy = _m_settings.PROCESSING_STRATEGY
    scene_id = str(scene["id"])
    existing = _m_registry.registry_table(config_hash)
    row = existing[existing["scene_id"] == scene_id]
    if not row.empty and row.iloc[0]["status"] == "done":
        stats_path = _m_settings.STATS_SCENE_DIR / f"{scene_id}_{config_hash}.parquet"
        if stats_path.exists() and scene_outputs_valid(scene_id, config_hash):
            print("Already complete; reusing", scene_id)
            return pd.read_parquet(stats_path)
    scene_lakes = _m_catalogue.candidate_lakes_for_scene(scene, lakes_all)
    if scene_lakes.empty:
        _m_registry.registry_update(
            scene_id, config_hash, status="skipped_no_lakes", candidate_lakes=0
        )
        return pd.DataFrame()
    _m_registry.registry_update(
        scene_id,
        config_hash,
        status="running",
        attempts=int(row.iloc[0]["attempts"] if not row.empty else 0) + 1,
    )
    raw_zip = None
    sen3 = None
    l2_paths: list[Path] = []
    success = False
    try:
        raw_zip, sen3, used_method = _m_staging.stage_scene_efficient(
            scene, username, password, download_method
        )
        windows = _m_windows.plan_scene_windows(scene, scene_lakes, strategy=strategy)
        _m_registry.registry_update(
            scene_id,
            config_hash,
            processing_strategy="full_scene" if windows == [None] else "adaptive_lake_windows",
            window_count=len(windows),
            download_method=used_method,
        )
        all_stats = []
        total_compact = 0
        window_provenance = []
        out_dir = _scene_output_dir(scene)
        out_dir.mkdir(parents=True, exist_ok=True)
        for window in windows:
            win_id = _m_windows.window_id(window)
            win_lakes = _m_windows.lakes_for_window(scene_lakes, window)
            compact_path = out_dir / f"{scene_id}_{win_id}_{config_hash}_lakepixels.nc"
            if compact_path.exists() and _m_compact.validate_compact_product(
                compact_path, scene_id=scene_id, config_hash=config_hash
            ):
                print("Reusing completed compact window", compact_path.name)
                stats = _m_statistics.compact_stats(compact_path, win_lakes)
                all_stats.append(stats)
                total_compact += compact_path.stat().st_size
                continue
            l2_path, l2_info = _m_l2gen.run_l2gen(
                sen3,
                mode=_m_settings.L2_MODE,
                process_full_scene=window is None,
                bbox=window,
                full_spectrum=False,
            )
            l2_paths.append(Path(l2_path))
            window_provenance.append(
                {
                    "window_id": win_id,
                    "bbox": None if window is None else list(window),
                    "l2gen": {k: v for k, v in l2_info.items() if k not in {"variables"}},
                }
            )
            compact, info = _m_compact.extract_lake_only_netcdf(
                l2_path,
                scene,
                win_lakes,
                compact_path,
                window=window,
                config_hash=config_hash,
                shore_buffer_m=_m_settings.DEFAULT_SHORE_BUFFER_M,
                rows_per_block=_m_settings.ROW_BLOCK_SIZE,
            )
            if compact is not None:
                stats = _m_statistics.compact_stats(compact, win_lakes)
                all_stats.append(stats)
                total_compact += compact.stat().st_size
            if not _m_settings.KEEP_L2_INTERMEDIATE:
                Path(l2_path).unlink(missing_ok=True)
        stats_table = pd.concat(all_stats, ignore_index=True) if all_stats else pd.DataFrame()
        if not stats_table.empty:
            stats_table = stats_table.drop_duplicates(["scene_id", "window_id", "Hylak_id"])
        _m_settings.STATS_SCENE_DIR.mkdir(parents=True, exist_ok=True)
        stats_path = _m_settings.STATS_SCENE_DIR / f"{scene_id}_{config_hash}.parquet"
        stats_tmp = stats_path.with_suffix(".parquet.tmp")
        stats_table.to_parquet(stats_tmp, index=False)
        stats_tmp.replace(stats_path)
        csv_path = stats_path.with_suffix(".csv")
        stats_table.to_csv(csv_path, index=False)
        manifest = {
            "schema": "global-efficient-v2.5",
            "created_utc": _m_utils.utc_now(),
            "scene": scene,
            "config_hash": config_hash,
            "candidate_lakes": int(len(scene_lakes)),
            "windows": ["FULL" if w is None else list(w) for w in windows],
            "window_provenance": window_provenance,
            "compact_files": sorted(
                (str(p) for p in out_dir.glob(f"*_{config_hash}_lakepixels.nc"))
            ),
            "stats_parquet": str(stats_path),
            "download_method": used_method,
            "cleanup": {
                "keep_staged_l1": _m_settings.KEEP_STAGED_L1,
                "keep_l2_intermediate": _m_settings.KEEP_L2_INTERMEDIATE,
                "keep_raw_zip_fallback": _m_settings.KEEP_RAW_ZIP_FALLBACK,
            },
        }
        _m_utils.write_json(
            _m_settings.EXPORT_DIR / scene_id / f"manifest_{config_hash}.json", manifest
        )
        _m_registry.registry_update(
            scene_id,
            config_hash,
            status="done",
            compact_bytes=int(total_compact),
            stats_rows=int(len(stats_table)),
            last_error=None,
        )
        with contextlib.suppress(Exception):
            _m_registry.registry_export_csv()
        success = True
        return stats_table
    except Exception as exc:
        _m_registry.registry_update(
            scene_id, config_hash, status="failed", last_error=f"{type(exc).__name__}: {exc}"
        )
        with contextlib.suppress(Exception):
            _m_registry.registry_export_csv()
        raise
    finally:
        _cleanup_scene_scratch(raw_zip, sen3, l2_paths, success=success)


# Reference cell 46, lines 130-135.
def load_scene_from_registry(scene_id: str, config_hash: str) -> dict[str, Any]:
    table = _m_registry.registry_table(config_hash)
    row = table[table.scene_id.astype(str) == str(scene_id)]
    if row.empty:
        raise KeyError(scene_id)
    return json.loads(row.iloc[0]["scene_json"])


# Reference cell 46, lines 138-167.
def process_next_queued(
    config_hash: str,
    lakes: gpd.GeoDataFrame,
    username: str,
    password: str,
    n_scenes: int = 1,
    download_method: str = None,
    retry_failed: bool = False,
    scene_ids: list[str] | None = None,
) -> pd.DataFrame:
    if download_method is None:
        download_method = _m_settings.DOWNLOAD_METHOD
    table = _m_registry.registry_table(config_hash)
    if scene_ids is not None:
        wanted = {str(x) for x in scene_ids}
        table = table[table.scene_id.astype(str).isin(wanted)].copy()
    statuses = ["queued", "running"] + (["failed"] if retry_failed else [])
    queue = table[table.status.isin(statuses)].head(int(n_scenes))
    outputs = []
    try:
        for _, row in queue.iterrows():
            scene = json.loads(row.scene_json)
            try:
                result = process_scene_once(
                    scene, config_hash, lakes, username, password, download_method
                )
                if not result.empty:
                    outputs.append(result)
            except Exception as exc:
                print(f"Scene failed {row.scene_name}: {exc}")
    finally:
        with contextlib.suppress(Exception):
            _m_staging.release_temporary_s3_credentials(username, password)
    return pd.concat(outputs, ignore_index=True) if outputs else pd.DataFrame()


# Reference cell 46, lines 181-196.
def storage_summary() -> pd.DataFrame:

    def tree_size(path: Path):
        if not path.exists():
            return 0
        return sum((p.stat().st_size for p in path.rglob("*") if p.is_file()))

    rows = []
    for label, path in (
        ("Runtime lake-only outputs", _m_settings.COMPACT_DIR),
        ("Runtime logs/exports", _m_settings.EXPORT_DIR),
        ("Scratch staged L1", _m_settings.STAGING_DIR),
        ("Scratch temporary L2", _m_settings.L2_DIR),
        ("Scratch fallback ZIP", _m_settings.RAW_L1_DIR),
        ("Scratch compact-file temp", _m_settings.SCRATCH_ROOT / "compact_tmp"),
    ):
        size = tree_size(path)
        rows.append(
            {"category": label, "path": str(path), "bytes": size, "size": _m_utils.human_size(size)}
        )
    return pd.DataFrame(rows)


# Reference cell 46, lines 199-203.
def clear_scratch() -> str:
    for path in (
        _m_settings.STAGING_DIR,
        _m_settings.L2_DIR,
        _m_settings.RAW_L1_DIR,
        _m_settings.DERIVED_DIR,
        _m_settings.SCRATCH_ROOT / "compact_tmp",
    ):
        shutil.rmtree(path, ignore_errors=True)
        path.mkdir(parents=True, exist_ok=True)
    return f"Scratch cleared: {_m_settings.SCRATCH_ROOT}"
