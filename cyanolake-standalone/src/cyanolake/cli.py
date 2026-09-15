"""Explicit command-line operations; no scene downloads or app launch on import."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path


def positive_int(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def positive_float(value):
    import math

    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be a finite positive number")
    return number


def parser():
    p = argparse.ArgumentParser(description="Standalone OLCI → L2Gen → native HydroLAKES science")
    p.add_argument(
        "--config",
        default="config/default.toml",
        help="TOML settings; relative paths resolve from this file",
    )
    sub = p.add_subparsers(dest="command", required=True)
    sub.add_parser("init", help="Create the configured data and scratch folders")
    q = sub.add_parser("preflight", help="Check local setup without downloading a scene")
    q.add_argument(
        "--without-lakes", action="store_true", help="Check OCSSW before preparing HydroLAKES"
    )
    sub.add_parser(
        "configure-earthdata", help="Prompt for Earthdata credentials and add new .netrc entries"
    )
    q = sub.add_parser(
        "install-ocssw", help="Explicitly download/install OCSSW into an empty configured root"
    )
    q.add_argument(
        "--tag", default="V2026.3", help="Preferred reference tag; installer checks availability"
    )
    q.add_argument(
        "--fallback-tags",
        type=int,
        default=3,
        help="Number of newest operational fallback tags; 0 disables fallback",
    )
    q = sub.add_parser("prepare-lakes", help="Load a local vector or download/cache HydroLAKES")
    q.add_argument("--source", type=Path, help="Existing HydroLAKES GPKG, shapefile or ZIP")
    q = sub.add_parser("list-admin", help="Fetch/cache ADM1 names for a country")
    q.add_argument("--country", required=True, help="Three-letter ISO country code, e.g. USA")
    q = sub.add_parser("plan", help="Discover scenes, intersect lakes and save a reusable plan")
    q.add_argument("--start", required=True, help="First UTC date, YYYY-MM-DD")
    q.add_argument("--end", required=True, help="Last UTC date, inclusive")
    group = q.add_mutually_exclusive_group()
    group.add_argument("--bbox", nargs=4, type=float, metavar=("WEST", "SOUTH", "EAST", "NORTH"))
    group.add_argument("--grid-id", help="For example G05_N40_W085")
    group.add_argument("--admin-key", help="ADM1 key returned by list-admin")
    q.add_argument("--country", help="ISO3 country required with --admin-key")
    q.add_argument(
        "--max-products", type=positive_int, help="Optional catalogue cap for a smoke test"
    )
    q = sub.add_parser("run", help="Process/resume a saved plan in bounded sequential batches")
    q.add_argument("--plan", required=True, type=Path, help="JSON plan written by plan")
    q.add_argument("--max-scenes", type=positive_int, default=1)
    q.add_argument("--retry-failed", action="store_true")
    q = sub.add_parser(
        "process-local", help="Process a local .SEN3 or ZIP using scene metadata JSON"
    )
    q.add_argument("--input", type=Path, required=True)
    q.add_argument("--scene-json", type=Path, required=True)
    q = sub.add_parser(
        "ingest-l2", help="Extract compact lake science from an existing L2Gen NetCDF"
    )
    q.add_argument("--input", type=Path, required=True)
    q.add_argument("--scene-json", type=Path, required=True)
    for name in ("status", "stats", "archive"):
        q = sub.add_parser(name, help=f"Inspect/export {name} for the current configuration")
        q.add_argument("--hash", help="Verify an expected processing configuration hash")
    q = sub.add_parser("export", help="Create the three-file scene share ZIP")
    q.add_argument("--scene-id", required=True)
    q.add_argument("--resolution", type=positive_float, default=300)
    q = sub.add_parser("map", help="Create a portable all-index scene HTML map")
    q.add_argument("--scene-id", required=True)
    q.add_argument("--resolution", type=positive_float, default=1000)
    q = sub.add_parser("compare-cyan", help="Compare a scene/lake with a NASA CyAN GeoTIFF")
    q.add_argument("--scene-id", required=True)
    q.add_argument("--lake-id", type=positive_int, required=True)
    q.add_argument("--date", default="", help="Defaults to the scene UTC date")
    q.add_argument(
        "--source", type=Path, help="Local NASA GeoTIFF; otherwise use Earthdata download"
    )
    q.add_argument("--max-distance", type=positive_float, default=600)
    q.add_argument("--qa", choices=("primary", "strict"), default="primary")
    q = sub.add_parser("diagnose-l2", help="Run the retained dense-swath mask diagnostics")
    q.add_argument("--input", required=True, type=Path)
    q = sub.add_parser("dashboard", help="Launch the complete local Gradio dashboard")
    q.add_argument("--host", default="127.0.0.1")
    q.add_argument("--port", type=positive_int, default=7860)
    return p


def emit(value):
    print(json.dumps(value, indent=2, default=str, allow_nan=False))


def read_scene(path):
    scene = json.loads(Path(path).read_text())
    if "scene" in scene:
        scene = scene["scene"]
    required = ("id", "name", "start", "end", "geofootprint")
    if any(not scene.get(key) for key in required):
        raise ValueError("Scene JSON requires id, name, start, end and GeoJSON geofootprint")
    for key in ("id", "name"):
        value = str(scene[key])
        if value in (".", "..") or "/" in value or "\\" in value:
            raise ValueError(f"Invalid scene {key}")
    return scene


def run_command(args):
    from . import registry, settings, state, workspace

    if args.command == "init":
        registry.registry_table()
        emit({"persistent_root": settings.PERSISTENT_ROOT, "scratch_root": settings.SCRATCH_ROOT})
        return 0
    if args.command == "configure-earthdata":
        from .credentials import configure_earthdata

        emit({"credential_file": configure_earthdata()})
        return 0
    if args.command == "install-ocssw":
        from .installer import install_ocssw_for_olci

        if args.fallback_tags < 0:
            raise ValueError("--fallback-tags must be nonnegative")
        settings.PREFERRED_OCSSW_TAG = args.tag
        settings.FALLBACK_TAG_COUNT = args.fallback_tags
        emit({"installed": install_ocssw_for_olci()})
        return 0
    if args.command == "preflight":
        from .runtime import preflight

        report = preflight(require_lakes=not args.without_lakes)
        emit(report)
        return 0 if report["ok"] else 2
    if args.command == "prepare-lakes":
        from . import catalogue, hydrolakes

        area = state.GLOBAL_STATE["min_area_km2"]
        region = state.GLOBAL_STATE["region"]
        if args.source:
            lakes, source = hydrolakes.load_hydrolakes_file(args.source, area)
        else:
            lakes, source = hydrolakes.prepare_hydrolakes(
                area, "USA" if region == "USA" else "WORLD", delete_global_archive=True
            )
        if region not in ("WORLD", "USA"):
            if "Continent" not in lakes:
                raise ValueError("The supplied HydroLAKES data has no Continent column")
            lakes = lakes[
                lakes.Continent.astype(str).str.lower() == region.replace("_", " ").lower()
            ].copy()
        if lakes.empty:
            raise ValueError("No lakes remain after filtering; target universe was not saved")
        catalogue.save_target_lakes(lakes)
        ch = workspace.select_configuration(lakes)
        emit(
            {
                "lakes": len(lakes),
                "source": source,
                "target": settings.TARGET_LAKES_GPKG,
                "config_hash": ch,
            }
        )
        return 0
    if args.command == "list-admin":
        from .admin import load_adm1_callback

        result = load_adm1_callback(args.country)
        print(result[1])
        emit({"choices": list(state.ADM1_STATE["lookup"])})
        return 0 if state.ADM1_STATE["lookup"] else 2
    if args.command == "diagnose-l2":
        from .dense import diagnose_l2_mask
        from .utils import write_json

        report = diagnose_l2_mask(args.input)
        out = settings.EXPORT_DIR / "l2_mask_diagnostic.json"
        write_json(out, report)
        emit({"diagnostic": out})
        return 0
    if args.command == "dashboard":
        from .runtime import activate_ocssw

        activate_ocssw()
        if settings.TARGET_LAKES_GPKG.exists():
            workspace.load_lakes()
        from .dashboard import build_dashboard

        if not 0 < args.port < 65536:
            raise ValueError("--port must be between 1 and 65535")
        app = build_dashboard()
        app.queue(default_concurrency_limit=1).launch(
            server_name=args.host,
            server_port=args.port,
            share=False,
            inbrowser=False,
            show_error=True,
            allowed_paths=[str(settings.PERSISTENT_ROOT)],
            css=settings.APP_CSS_GLOBAL,
        )
        return 0
    plan = json.loads(args.plan.read_text()) if args.command == "run" else None
    expected = plan["config_hash"] if plan else getattr(args, "hash", None)
    lakes, ch = workspace.load_lakes(expected)
    if args.command == "plan":
        import pandas as pd
        from shapely.geometry import box

        from . import admin, archive, callbacks, catalogue

        if pd.Timestamp(args.end) < pd.Timestamp(args.start):
            raise ValueError("end date must not precede start date")
        bbox = tuple(args.bbox) if args.bbox else None
        geom = None
        label = "WORLD"
        if args.grid_id:
            bbox = archive.grid_id_to_bbox(args.grid_id)
            label = args.grid_id
        elif args.admin_key:
            if not args.country:
                raise ValueError("--country is required with --admin-key")
            admin.load_adm1_callback(args.country)
            geom, label = admin._selected_adm1_geometry(args.country, args.admin_key)
            bbox = tuple(geom.bounds)
        if bbox:
            bbox = catalogue.normalize_bbox(bbox)
            if geom is None:
                geom = box(*bbox)
            if label == "WORLD":
                label = "bbox_" + "_".join(map(str, bbox))
        table, ch, ids, csv, nlakes = callbacks._plan_catalogue_shard(
            args.start,
            args.end,
            lakes,
            state.GLOBAL_STATE["min_area_km2"],
            state.GLOBAL_STATE["region"],
            bbox,
            geom,
            label,
            args.max_products,
        )
        from .utils import write_json

        out = Path(csv).with_suffix(".json")
        scenes = [json.loads(row.scene_json) for row in table.itertuples()]
        write_json(
            out,
            {
                "schema": "cyanolake-plan-1",
                "config_hash": ch,
                "start": args.start,
                "end": args.end,
                "shard": label,
                "bbox": bbox,
                "scene_ids": ids,
                "scenes": scenes,
            },
        )
        for scene in scenes:
            write_json(out.parent / "scenes" / f"{scene['id']}.json", scene)
        emit(
            {
                "plan": out,
                "queue_csv": csv,
                "config_hash": ch,
                "scene_count": len(ids),
                "candidate_lakes": nlakes,
            }
        )
        return 0 if ids else 2
    if args.command in ("status", "stats", "archive"):
        from . import archive, statistics

        if args.command == "status":
            table = registry.registry_table(ch)
            print(
                table[["scene_id", "status", "attempts", "stats_rows", "last_error"]].to_string(
                    index=False
                )
            )
        elif args.command == "stats":
            table = statistics.all_stats(ch)
            out = settings.EXPORT_DIR / "all_lake_statistics.csv"
            table.to_csv(out, index=False)
            emit({"rows": len(table), "csv": out})
        else:
            result = archive.build_archive_index(ch)
            print(result)
        return 0
    if args.command in ("run", "process-local"):
        from . import pipeline, runtime

        report = runtime.preflight()
        if not report["ok"]:
            emit(report)
            return 2
        if args.command == "process-local":
            scene = read_scene(args.scene_json)
            scene["local_input"] = str(args.input.expanduser().resolve())
            registry.registry_upsert_scene(scene, ch, len(lakes))
            table = pipeline.process_scene_once(scene, ch, lakes, "", "")
            emit({"config_hash": ch, "scene_id": scene["id"], "statistics_rows": len(table)})
            return 0
        ids = plan["scene_ids"]
        table = registry.registry_table(ch)
        selected = table[table.scene_id.astype(str).isin(set(ids))]
        if len(selected) != len(set(ids)):
            raise ValueError("Saved plan includes scenes absent from this workspace registry")
        for row in selected.itertuples():
            stats = settings.STATS_SCENE_DIR / f"{row.scene_id}_{ch}.parquet"
            if row.status == "done" and (
                not stats.exists() or not pipeline.scene_outputs_valid(row.scene_id, ch)
            ):
                registry.registry_update(
                    row.scene_id,
                    ch,
                    status="queued",
                    last_error="Saved products are missing or invalid; queued for repair",
                )
        selected = registry.registry_table(ch)
        selected = selected[selected.scene_id.astype(str).isin(set(ids))]
        eligible = selected.status.isin(
            ["queued", "running"] + (["failed"] if args.retry_failed else [])
        )
        if not eligible.any():
            emit(
                {
                    "message": "No eligible scenes in this plan",
                    "statuses": selected.status.value_counts().to_dict(),
                }
            )
            return 2 if (selected.status == "failed").any() else 0
        from .credentials import cdse_credentials

        username, password = cdse_credentials()
        pipeline.process_next_queued(
            ch,
            lakes,
            username,
            password,
            n_scenes=args.max_scenes,
            download_method=settings.DOWNLOAD_METHOD,
            retry_failed=args.retry_failed,
            scene_ids=ids,
        )
        after = registry.registry_table(ch)
        after = after[after.scene_id.astype(str).isin(set(ids))]
        emit({"config_hash": ch, "statuses": after.status.value_counts().to_dict()})
        return 2 if (after.status == "failed").any() else 0
    if args.command == "ingest-l2":
        from .ingest import ingest_l2

        emit(ingest_l2(args.input, read_scene(args.scene_json), lakes, ch))
        return 0
    if args.command == "export":
        from .exports import create_scene_share_bundle

        bundle, status = create_scene_share_bundle(args.scene_id, args.resolution)
        emit({"bundle": bundle, "status": status})
        return 0
    if args.command == "map":
        from .mapping import create_complete_scene_html

        path, grid = create_complete_scene_html(args.scene_id, args.resolution)
        emit({"html": path, "resolution_m": grid["resolution_m"]})
        return 0
    if args.command == "compare-cyan":
        from .comparison import compare_cyan_callback_v265

        result = compare_cyan_callback_v265(
            args.scene_id,
            str(args.lake_id),
            args.date,
            "upload" if args.source else "download",
            str(args.source.resolve()) if args.source else None,
            args.max_distance,
            args.qa,
        )
        print(result[-1])
        # Callback supplies downloadable output files on success, None on failure.
        files = [
            x for x in result if isinstance(x, str) and Path(x).suffix in (".csv", ".nc", ".zip")
        ]
        emit({"files": files, "comparison_directory": settings.EXPORT_DIR / "cyan_comparison"})
        return 0 if files else 2
    raise ValueError(f"Unsupported command: {args.command}")


def main(argv=None):
    args = parser().parse_args(argv)
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("GRADIO_ANALYTICS_ENABLED", "False")
    try:
        from .configuration import configure

        configure(args.config, create=True)
        from .locking import workspace_lock

        with workspace_lock():
            return run_command(args)
    except KeyboardInterrupt:
        print("Interrupted. Valid compact windows can be reused on the next run.", file=sys.stderr)
        return 130
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 2
