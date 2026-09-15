# The data lifecycle

Paths below are relative to the configured roots. The default roots are the package's `data/` and `scratch/` folders, resolved relative to `config/default.toml`.

```text
data/
  scene_registry.sqlite             scene/configuration status ledger
  scene_registry.sqlite-wal         live SQLite WAL, when present
  catalog/
    target_hydrolakes.gpkg           current canonical lake universe
    configurations/<hash>/
      configuration.json            settings + lake-geometry identity
      target_hydrolakes.gpkg         frozen copy for this configuration
    plans/<hash>/
      plan_<region>_<dates>_<time>.json   reusable scene IDs + scene records
      plan_<region>_<dates>_<time>.csv    readable queue snapshot
      scenes/<scene_id>.json         individual scene metadata
    global_archive_index.csv        native orbit/frame and file catalogue
  hydrolakes/                       prepared/cached reference vectors
  admin_boundaries/                 cached ADM1 vectors and metadata
  lake_only/
    netcdf/
      platform=S3A/relative_orbit=012/frame=2520/
        year=2024/month=01/day=01/<scene_id>/
          <scene_id>_<window>_<hash>_lakepixels.nc
    stats_by_scene/
      <scene_id>_<hash>.parquet
      <scene_id>_<hash>.csv
  par/<hash>/                       L2Gen parameters, ancillary records/database
  logs/                             getanc and L2Gen logs
  exports/<hash>/
    <scene_id>/manifest_<hash>.json
    snap_complete/                  complete mapped science NetCDF
    scene_bundles/                  three-file share ZIP + components
    spectra/                        coordinate and lake spectra
    ...                             comparison and other on-demand outputs
  html_maps/<hash>/                 standalone Folium/Leaflet maps

scratch/
  staging/                          extracted or direct-downloaded .SEN3
  raw_l1_fallback/                  fallback ZIP and .part downloads
  l2/<hash>/                        temporary dense L2Gen NetCDF
  legacy_dense_derived/             optional dense diagnostic products
  compact_tmp/                      legacy helper temporary writes, if used
```

## What one compact NetCDF stores

Each row of its one-dimensional `obs` dimension is a stored native OLCI observation: latitude, longitude, `Hylak_id`, source row/column within its L2 processing window, raw `l2_flags`, nine `rhos` bands, spectral shapes, the five production indices, MPH peak wavelength, and QA masks. Non-lake pixels are omitted.

The nominal 885-nm archive band reads L2Gen's actual `rhos_884` product. The compact variable retains the reference name `rhos_885`; the complete GIS export uses `rhos_884` and records the mapping.

QA-excluded science and algorithm non-detections use fill values that decode as NaN. They are not zero concentrations. Boolean masks distinguish a valid non-detection from a QA exclusion.

## What the statistics store

The scene statistics contain one row per lake and processing window, including retained/detected counts, fractions and index summaries. The share statistics consolidate the scene using the reference's native-location deduplication. Do not average window medians to obtain a scene median.

## What is deleted by default

After successful durable output creation, the processor removes its staged L1, fallback ZIP and dense L2. Failed-run scratch is retained by default. Local input files supplied outside the managed scratch folders are preserved.

Retention switches live in `[storage]`. Keeping L2 can help initial validation but increases disk use substantially. Peak space still needs to accommodate the source product plus dense L2; ZIP fallback may temporarily need both ZIP and extracted L1.

## How resume works

Run the same saved plan with the same config and target lake universe. Completed scenes are reused when their statistics and compact files validate. Interrupted scenes and already completed windows can be resumed. Missing outputs return a CLI scene to the queue. Failed scenes require `--retry-failed`.

Discovery is explicit: `run` reads a saved plan and does not rediscover scenes. A new `plan` call creates a new timestamped snapshot.

## Backups and transfers

Keep the entire `data/` tree for a durable archive. Stop writers before copying SQLite, or use a proper SQLite backup tool. A share ZIP contains a mapped science product and summaries; it is not a complete native archive backup.

Manifests and statistics retain absolute paths from the processing machine. The compact loader has a configuration-scoped fallback search under the configured archive root, but moving the tree does not rewrite every provenance path. Preserve the original paths or update deployment paths deliberately.

**Think about this:** deleting `scratch/` means atmospheric correction may need to run again. Deleting `lake_only/` means deleting your persistent measurements. The two folders serve different purposes.
