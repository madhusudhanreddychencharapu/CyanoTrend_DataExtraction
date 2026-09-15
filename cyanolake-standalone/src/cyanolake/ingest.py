"""Explicit entry point for existing dense L2Gen products."""

from pathlib import Path

import pandas as pd

from . import catalogue, compact, pipeline, registry, settings, statistics, utils, windows


def ingest_l2(l2_path, scene, lakes, config_hash):
    l2_path = Path(l2_path).expanduser().resolve()
    subset = catalogue.candidate_lakes_for_scene(scene, lakes)
    if subset.empty:
        raise ValueError("Scene footprint intersects no target lakes")
    registry.registry_upsert_scene(scene, config_hash, len(subset))
    out = (
        pipeline._scene_output_dir(scene)
        / f"{scene['id']}_{windows.window_id(None)}_{config_hash}_lakepixels.nc"
    )
    try:
        path, info = compact.extract_lake_only_netcdf(
            l2_path,
            scene,
            subset,
            out,
            config_hash=config_hash,
            shore_buffer_m=settings.DEFAULT_SHORE_BUFFER_M,
            rows_per_block=settings.ROW_BLOCK_SIZE,
        )
        stats = statistics.compact_stats(path, subset) if path else pd.DataFrame()
        stats_path = settings.STATS_SCENE_DIR / f"{scene['id']}_{config_hash}.parquet"
        tmp = stats_path.with_suffix(".parquet.tmp")
        stats.to_parquet(tmp, index=False)
        tmp.replace(stats_path)
        stats.to_csv(stats_path.with_suffix(".csv"), index=False)
        manifest = settings.EXPORT_DIR / scene["id"] / f"manifest_{config_hash}.json"
        utils.write_json(
            manifest,
            {
                "schema": "cyanolake-existing-l2-1",
                "created_utc": utils.utc_now(),
                "scene": scene,
                "config_hash": config_hash,
                "source_l2": str(l2_path),
                "compact_files": [str(path)] if path else [],
                "stats_parquet": str(stats_path),
                "extraction": info,
                "validation_scope": "Existing L2 input; atmospheric correction was not rerun",
            },
        )
        registry.registry_update(
            scene["id"],
            config_hash,
            status="done",
            stats_rows=len(stats),
            compact_bytes=path.stat().st_size if path else 0,
            last_error=None,
            download_method="existing_l2",
        )
        return {
            "compact": path,
            "statistics": stats_path,
            "manifest": manifest,
            "observations": info["nobs"],
        }
    except Exception as exc:
        registry.registry_update(
            scene["id"], config_hash, status="failed", last_error=f"{type(exc).__name__}: {exc}"
        )
        raise
