"""Standalone dashboard implementation."""

from __future__ import annotations

import gradio as gr
import pandas as pd
import pycountry

from . import admin as _m_admin
from . import archive as _m_archive
from . import callbacks as _m_callbacks
from . import comparison as _m_comparison
from . import exports as _m_exports
from . import mapping as _m_mapping
from . import settings as _m_settings
from . import spectra as _m_spectra

COUNTRY_CHOICES = [
    (f"{c.name} ({c.alpha_3})", c.alpha_3)
    for c in sorted(pycountry.countries, key=lambda c: c.name)
]


def build_dashboard():
    with gr.Blocks(title="Sentinel-3 OLCI Global CyanoLake Processor v2.6.8") as DEMO_GLOBAL:
        gr.Markdown("# Sentinel-3 OLCI → L2Gen → global HydroLAKES processor · v2.6.8")
        gr.Markdown(
            "Production indices are **CI, CIcyano, NDCI, MPH, and OLCI-adapted FAI/AFAI**. The dashboard provides exact-resolution map previews, active State/Province/bbox display filtering, scene-only lake dropdowns, clickable spectral sampling, expanded NASA CyAN bias diagnostics, and a SNAP-friendly EPSG:4326 share raster that uses nearest native OLCI observations rather than center-only binning."
        )
        with gr.Tab("1 · Plan & process scenes"):
            gr.Markdown(
                "### Recommended order\n**A.** Load HydroLAKES. For country/state-by-state production, keep **WORLD** so one canonical lake universe is used.  \n**B.** Choose **State / province (ADM1)**, **Bounding box**, or **5° grid** as the work-tracking shard.  \n**C.** Plan the queue and process a small batch. A shard controls discovery/tracking only; Sentinel-3 scene ID + science configuration prevents duplicate L2Gen processing.  \n**D.** Results persist in your configured data directory. Share ZIPs and registry exports are available after processing."
            )
            with gr.Row():
                region = gr.Dropdown(
                    [
                        "WORLD",
                        "USA",
                        "Africa",
                        "Asia",
                        "Europe",
                        "North America",
                        "South America",
                        "Oceania",
                    ],
                    value="WORLD",
                    label="HydroLAKES universe (WORLD recommended)",
                )
                min_area = gr.Number(value=10.0, label="Minimum lake area (km²)")
                existing_hydro = gr.Textbox(value="", label="Optional existing HydroLAKES path")
                load_lakes = gr.Button("Load / prepare HydroLAKES", variant="primary")
            lake_preview = gr.Dataframe(interactive=False)
            lake_status = gr.HTML()
            load_lakes.click(
                _m_callbacks.load_lakes_callback,
                [region, min_area, existing_hydro],
                [lake_preview, lake_status],
            )
            gr.Markdown("### Choose how to define the next processing shard")
            planning_mode = gr.Radio(
                [
                    ("State / province / region (ADM1)", "adm1"),
                    ("Bounding box", "bbox"),
                    ("Fixed 5° grid", "grid"),
                ],
                value="adm1",
                label="Planning method",
            )
            with gr.Row():
                country_iso3 = gr.Dropdown(
                    COUNTRY_CHOICES, value="USA", label="Country (for ADM1 mode)"
                )
                load_admin = gr.Button("Load states / provinces")
                admin1_choice = gr.Dropdown(label="State / province / first-level region")
            admin_status = gr.HTML()
            load_admin.click(
                _m_admin.load_adm1_callback, country_iso3, [admin1_choice, admin_status]
            )
            gr.Markdown(
                '<span class="small-note">ADM1 boundaries are fetched on demand from geoBoundaries gbOpen. The selected geometry is retained as the **active planning region** and can be used later to filter the map display without changing the full-scene archive/export.</span>'
            )
            with gr.Row():
                start_date = gr.Textbox(value="2024-08-01", label="Start date")
                end_date = gr.Textbox(value="2024-08-07", label="End date")
                bbox_text = gr.Textbox(value="", label="BBox west,south,east,north (bbox mode)")
                grid_id = gr.Textbox(
                    value="", placeholder="G05_N40_W085", label="5° grid ID (grid mode)"
                )
                max_products = gr.Number(value=0, precision=0, label="Catalogue cap (0=all)")
                plan_button = gr.Button("Plan scene queue", variant="primary")
            queue_table = gr.Dataframe(interactive=False)
            plan_snapshot = gr.File(label="Current shard plan CSV")
            plan_status = gr.HTML()
            plan_button.click(
                _m_callbacks.plan_callback,
                [
                    start_date,
                    end_date,
                    planning_mode,
                    country_iso3,
                    admin1_choice,
                    bbox_text,
                    grid_id,
                    max_products,
                ],
                [queue_table, plan_snapshot, plan_status],
            )
            with gr.Row():
                cdse_user = gr.Textbox(label="CDSE username")
                cdse_password = gr.Textbox(label="CDSE password", type="password")
                n_scenes = gr.Slider(1, 20, value=1, step=1, label="Scenes this batch")
                dl_method = gr.Dropdown(["auto", "s3", "zip"], value="auto", label="L1 staging")
                retry_failed = gr.Checkbox(False, label="Retry failed")
                run_button = gr.Button("Run next queued scenes", variant="primary")
            batch_status = gr.HTML()
            storage_after = gr.Dataframe(interactive=False)
            run_button.click(
                _m_callbacks.run_batch_callback,
                [cdse_user, cdse_password, n_scenes, dl_method, retry_failed],
                [queue_table, batch_status, storage_after],
            )
            with gr.Row():
                refresh_queue = gr.Button("Refresh queue")
                registry_file = gr.File(label="Scene registry CSV")
            refresh_queue.click(
                _m_callbacks.refresh_queue_callback, [], [queue_table, registry_file]
            )
        with gr.Tab("2 · Lake time series"):
            with gr.Row():
                refresh_lakes = gr.Button("Refresh processed lakes", variant="primary")
                lake_choice = gr.Dropdown(label="Lake")
                plot_metric = gr.Dropdown(
                    _m_settings.PLOT_METRICS, value=_m_settings.PLOT_METRICS[0], label="Statistic"
                )
                plot_button = gr.Button("Plot lake history")
            lakes_status = gr.HTML()
            lake_plot = gr.Plot()
            lake_table = gr.Dataframe(interactive=False)
            refresh_lakes.click(
                _m_callbacks.refresh_lakes_callback, [], [lake_choice, lakes_status]
            )
            plot_button.click(
                _m_callbacks.lake_timeseries_callback,
                [lake_choice, plot_metric],
                [lake_plot, lake_table],
            )
        with gr.Tab("3 · Scene / multi-lake map"):
            gr.Markdown(
                "**Resolution is now literal.** Choosing 300 m requests an exact 300 m Web-Mercator display; the app no longer silently changes it to ~500 m. If a full-scene 300 m preview is too large for the browser, the app will tell you to use the active planning region, a lake, or explicitly choose a coarser preview. The **Active planning region** option clips the display to the State/Province/bbox/5° shard selected in Tab 1; the stored scene and share ZIP remain full-scene. Mapped cells use the nearest native observation from the **same HydroLAKES lake only**, preventing cross-lake borrowing and avoiding center-binning striping without interpolation."
            )
            with gr.Row():
                refresh_scenes = gr.Button("Refresh completed scenes", variant="primary")
                scene_choice = gr.Dropdown(label="Completed scene")
                map_scope = gr.Radio(
                    [
                        ("Entire processed scene · all lakes", "scene"),
                        ("Active planning region · State/bbox/grid", "planning"),
                        ("Selected lake", "lake"),
                    ],
                    value="planning",
                    label="Map display scope",
                )
                scene_lake = gr.Dropdown(label="Lake within selected scene")
            scene_status = gr.HTML()
            refresh_scenes.click(
                _m_archive.refresh_scene_choices_callback,
                [],
                [scene_choice, scene_lake, scene_status],
            )
            scene_choice.change(_m_archive.scene_changed_callback, scene_choice, scene_lake)
            with gr.Row():
                map_metric = gr.Dropdown(
                    _m_settings.MAP_METRICS, value="CI_cyano", label="Map layer"
                )
                map_resolution = gr.Dropdown(
                    [300, 500, 1000, 2000], value=300, label="Exact display grid (m)"
                )
                basemap = gr.Dropdown(_m_settings.BASEMAPS, value="Esri Dark Gray", label="Basemap")
                labels = gr.Checkbox(True, label="Show Esri gray labels")
                show_lakes = gr.Checkbox(True, label="HydroLAKES boundaries")
            with gr.Row():
                cmap = gr.Dropdown(_m_settings.CMAPS, value="cividis", label="Colormap")
                reverse = gr.Checkbox(False, label="Reverse")
                vmin = gr.Textbox(value="", placeholder="Auto", label="Color minimum")
                vmax = gr.Textbox(value="", placeholder="Auto", label="Color maximum")
                render = gr.Button("Render / export map", variant="primary")
            map_status = gr.HTML()
            map_html = gr.HTML()
            with gr.Row():
                html_file = gr.File(label="Portable interactive HTML · hover values")
                snap_nc = gr.File(label="Displayed-scope CF NetCDF · exact nominal resolution")
                compact_file = gr.File(label="Compact native lake-pixel NetCDF / ZIP")
            render.click(
                _m_mapping.map_callback,
                [
                    scene_choice,
                    map_scope,
                    scene_lake,
                    map_metric,
                    map_resolution,
                    basemap,
                    labels,
                    cmap,
                    reverse,
                    vmin,
                    vmax,
                    show_lakes,
                ],
                [map_html, html_file, snap_nc, compact_file, map_status],
            )
        with gr.Tab("4 · Spectra & lake statistics"):
            gr.Markdown(
                "All spectra use the **native stored OLCI lake pixels**. The coordinate workflow now includes a clickable geographic lake-pixel map; a click is converted to lat/lon and matched to the nearest eligible native OLCI pixel."
            )
            gr.Markdown("### A · Selected-lake mean / median spectrum")
            with gr.Row():
                spec_scene = gr.Dropdown(label="Completed scene")
                spec_lake = gr.Dropdown(label="Lake within scene")
                spec_stat = gr.Radio(["median", "mean"], value="median", label="Spectrum statistic")
                spec_button = gr.Button("Plot lake spectrum", variant="primary")
                spec_refresh = gr.Button("Refresh scenes")
            spec_status = gr.HTML()
            with gr.Row():
                spec_plot = gr.Plot()
                spec_table = gr.Dataframe(interactive=False, label="Lake spectral summary")
            lake_stats_table = gr.Dataframe(interactive=False, label="Stored lake statistics")
            spec_csv = gr.File(label="Lake spectrum CSV")
            spec_refresh.click(
                _m_archive.refresh_scene_choices_callback, [], [spec_scene, spec_lake, spec_status]
            )
            spec_scene.change(_m_archive.scene_changed_callback, spec_scene, spec_lake)
            spec_button.click(
                _m_spectra.scene_spectrum_callback,
                [spec_scene, spec_lake, spec_stat],
                [spec_plot, spec_table, lake_stats_table, spec_csv, spec_status],
            )
            gr.Markdown("### B · Coordinate / map-click native-pixel spectra")
            coordinate_points = gr.Dataframe(
                headers=["label", "latitude", "longitude"],
                datatype=["str", "number", "number"],
                value=pd.DataFrame([{"label": "P1", "latitude": 41.7, "longitude": -82.9}]),
                row_count=(1, "dynamic"),
                column_count=(3, "fixed"),
                type="pandas",
                interactive=True,
                label="Target coordinates · edit manually or click the map",
            )
            with gr.Row():
                coord_label = gr.Textbox(value="P2", label="New point label")
                coord_lat = gr.Number(value=41.8, label="Latitude")
                coord_lon = gr.Number(value=-83.1, label="Longitude")
                coord_add = gr.Button("Add coordinate")
                coord_clear = gr.Button("Clear points")
                spectrum_map_refresh = gr.Button("Refresh clickable map")
            with gr.Row():
                coord_basis = gr.Radio(
                    [
                        ("Any stored HydroLAKES pixel", "stored_lake"),
                        ("Primary bloom-aware valid-water pixel", "valid_water"),
                        ("CI-input-valid pixel", "ci_valid"),
                    ],
                    value="stored_lake",
                    label="Nearest-pixel basis",
                )
                coord_distance = gr.Number(
                    value=2000, minimum=100, maximum=10000, label="Maximum match distance (m)"
                )
                coord_plot_button = gr.Button("Plot coordinate spectra", variant="primary")
            spectrum_map_state = gr.State({})
            with gr.Row():
                spectrum_map = gr.Image(
                    label="Clickable native OLCI lake-pixel map · click to add a spectrum point",
                    type="pil",
                    interactive=True,
                    height=560,
                )
                coord_plot = gr.Plot(label="Native-pixel spectra")
            coord_table = gr.Dataframe(
                interactive=False, label="Matched spectral values and coordinates"
            )
            coord_csv = gr.File(label="Coordinate spectra CSV")
            coord_status = gr.HTML()
            coord_add.click(
                _m_spectra.add_spectrum_point_callback,
                [coord_label, coord_lat, coord_lon, coordinate_points],
                [coordinate_points, coord_status],
            )
            spectrum_map_refresh.click(
                _m_spectra.spectrum_map_refresh_callback,
                [spec_scene, coordinate_points],
                [spectrum_map, spectrum_map_state, coord_status],
            )
            coord_clear.click(
                _m_spectra.clear_spectrum_points_and_map_callback,
                [spec_scene],
                [coordinate_points, spectrum_map, spectrum_map_state, coord_status],
            )
            coord_plot_button.click(
                _m_spectra.coordinate_spectra_with_map_callback,
                [spec_scene, coordinate_points, coord_basis, coord_distance],
                [
                    spectrum_map,
                    spectrum_map_state,
                    coord_plot,
                    coord_table,
                    coord_csv,
                    coord_status,
                ],
            )
            spectrum_map.select(
                _m_spectra.spectrum_map_click_callback,
                [spec_scene, coordinate_points, coord_basis, coord_distance, spectrum_map_state],
                [
                    coordinate_points,
                    spectrum_map,
                    spectrum_map_state,
                    coord_plot,
                    coord_table,
                    coord_csv,
                    coord_status,
                ],
            )
            spec_scene.change(
                _m_spectra.spectrum_map_refresh_callback,
                [spec_scene, coordinate_points],
                [spectrum_map, spectrum_map_state, coord_status],
            )
        with gr.Tab("5 · NASA CyAN comparison"):
            gr.Markdown(
                "Compare one processed scene/lake with same-day NASA **Merged-S3-CYAN**. The NASA daily merged product uses the per-pixel maximum across S3A/S3B, while the local result is one OLCI scene; therefore high-end local-minus-NASA negative bias can be physically expected and should not be removed by empirical scaling. Use **CyAN-strict** QA for the closest product-consistency comparison; bloom-aware QA is provided as a diagnostic sensitivity option."
            )
            with gr.Row():
                cyan_refresh = gr.Button("Refresh completed scenes", variant="primary")
                cyan_scene = gr.Dropdown(label="Completed scene")
                cyan_lake = gr.Dropdown(label="Lake actually stored in selected scene")
            cyan_refresh_status = gr.HTML()
            cyan_refresh.click(
                _m_archive.refresh_scene_choices_callback,
                [],
                [cyan_scene, cyan_lake, cyan_refresh_status],
            )
            cyan_scene.change(_m_archive.scene_changed_callback, cyan_scene, cyan_lake)
            with gr.Row():
                cyan_date = gr.Textbox(
                    value="", placeholder="blank = infer from scene", label="UTC date YYYY-MM-DD"
                )
                cyan_mode = gr.Radio(
                    [
                        ("Download from NASA OBPG using configured Earthdata Login", "download"),
                        ("Upload existing OBPG CIcyano GeoTIFF", "upload"),
                    ],
                    value="download",
                    label="NASA source",
                )
                cyan_upload = gr.File(
                    label="Optional OBPG CIcyano GeoTIFF",
                    type="filepath",
                    file_types=[".tif", ".tiff"],
                )
                cyan_distance = gr.Number(
                    value=450, minimum=50, maximum=1500, label="Max nearest-neighbour distance (m)"
                )
                cyan_qa_basis = gr.Radio(
                    [
                        ("Bloom-aware · recommended science mask", "primary"),
                        ("Generic CLDICE-excluded · sensitivity only", "strict"),
                    ],
                    value="primary",
                    label="Local QA basis",
                )
            with gr.Row():
                cyan_prepare = gr.Button("1 · Download / check NASA GeoTIFF")
                cyan_compare = gr.Button("2 · Compare selected lake", variant="primary")
            cyan_source = gr.File(label="NASA source GeoTIFF")
            cyan_status = gr.HTML()
            cyan_coverage = gr.Dataframe(interactive=False, label="Coverage / pairing audit")
            cyan_metrics = gr.Dataframe(interactive=False, label="Bias / correlation diagnostics")
            cyan_plot = gr.Plot(label="Paired detections + binned median trend")
            with gr.Row():
                cyan_pairs = gr.File(label="Paired cells CSV")
                cyan_metrics_file = gr.File(label="Metrics CSV")
                cyan_nc = gr.File(label="Comparison NetCDF")
                cyan_tif = gr.File(label="NASA source GeoTIFF")
            cyan_prepare.click(
                _m_comparison.prepare_cyan_source_callback,
                [cyan_scene, cyan_date, cyan_mode, cyan_upload],
                [cyan_source, cyan_status],
            )
            cyan_compare.click(
                _m_comparison.compare_cyan_callback_v265,
                [
                    cyan_scene,
                    cyan_lake,
                    cyan_date,
                    cyan_mode,
                    cyan_upload,
                    cyan_distance,
                    cyan_qa_basis,
                ],
                [
                    cyan_plot,
                    cyan_coverage,
                    cyan_metrics,
                    cyan_pairs,
                    cyan_metrics_file,
                    cyan_nc,
                    cyan_tif,
                    cyan_status,
                ],
            )
        with gr.Tab("6 · Global archive & share"):
            gr.Markdown(
                "Use this tab **after processing scenes**. The share ZIP remains deliberately small: **one all-science NetCDF + one all-index per-lake CSV + one metadata JSON**. The all-science NetCDF is now **EPSG:4326 for direct SNAP/GIS geolocation** and uses an exact nominal requested grid (300 m by default). It maps each target lake cell only to the nearest stored native OLCI observation from the same Hylak_id before applying QA. This avoids artificial center-binning striping while preserving cloud/no-data gaps. NASA-style area weighting is not used in the production branch because native observations remain the quantitative source of truth."
            )
            with gr.Row():
                archive_grid = gr.Textbox(
                    value="", placeholder="G05_N40_W085", label="Optional grid filter"
                )
                archive_orbit = gr.Textbox(
                    value="", placeholder="211", label="Optional relative-orbit filter"
                )
                archive_refresh = gr.Button("Refresh archive index", variant="primary")
            archive_table = gr.Dataframe(interactive=False)
            archive_csv = gr.File(label="Global archive index CSV")
            archive_status = gr.HTML()
            archive_refresh.click(
                _m_archive.archive_refresh_callback,
                [archive_grid, archive_orbit],
                [archive_table, archive_csv, archive_status],
            )
            with gr.Row():
                bundle_scene = gr.Dropdown(label="Completed scene")
                bundle_resolution = gr.Dropdown(
                    [300, 500, 1000], value=300, label="Share NetCDF nominal grid resolution (m)"
                )
                bundle_refresh = gr.Button("Refresh completed scenes")
                bundle_button = gr.Button("Create 3-file share ZIP", variant="primary")
            bundle_dummy_lake = gr.Dropdown(visible=False)
            bundle_file = gr.File(label="Shareable scene ZIP · exactly 3 files")
            bundle_status = gr.HTML()
            bundle_refresh.click(
                _m_archive.refresh_scene_choices_callback,
                [],
                [bundle_scene, bundle_dummy_lake, bundle_status],
            )
            bundle_button.click(
                _m_exports.bundle_callback,
                [bundle_scene, bundle_resolution],
                [bundle_file, bundle_status],
            )
        with gr.Tab("7 · Storage & provenance"):
            with gr.Row():
                refresh_storage = gr.Button("Refresh storage")
                clear_button = gr.Button("Clear scratch only")
            storage_table = gr.Dataframe(interactive=False)
            storage_status = gr.HTML()
            registry_download = gr.File(label="Registry CSV")
            refresh_storage.click(
                _m_callbacks.storage_callback, [], [storage_table, registry_download]
            )
            clear_button.click(
                _m_callbacks.clear_scratch_callback, [], [storage_table, storage_status]
            )
    return DEMO_GLOBAL
