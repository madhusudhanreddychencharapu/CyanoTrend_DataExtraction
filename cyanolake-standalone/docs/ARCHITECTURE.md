# Architecture, one responsibility at a time

## 1. Configuration decides where and how work runs

`configuration.py` reads TOML, rejects unknown/invalid settings, and creates folders only when requested. `settings.py` contains the final reference defaults. `state.py` holds the current dashboard selection. One process uses one workspace; run separate workspaces for independent jobs.

`workspace.py` fingerprints the canonical lake geometries and selects a science configuration. Configuration snapshots and target lake copies make the processing identity inspectable. Native outputs, scratch L2, and derived exports are separated by configuration where reuse would otherwise be unsafe.

## 2. Discovery determines which scenes are relevant

`catalogue.py` searches the public CDSE catalogue, follows pagination, and intersects scene footprints with HydroLAKES. `admin.py` fetches/caches geoBoundaries ADM1 geometry. `callbacks.py` connects these operations to the dashboard and shard planning.

A country, state, bbox or 5° grid is a discovery region. The expensive unit of processing is a Sentinel-3 scene. A scene crossing two states is stored once per science configuration. Processing uses all eligible lakes from the canonical universe in that scene.

## 3. The registry remembers work

`registry.py` owns SQLite scene status, attempts, scene metadata, and configuration identity. Its primary key is `(scene_id, config_hash)`. Replanning updates metadata without resetting a completed row.

`locking.py` gives one process exclusive write access to a workspace. The OS releases the lock after a crash. The next CLI batch can resume rows left in `running` state. A saved empty plan never expands into the whole queue.

## 4. Staging provides real local L1 files to L2Gen

`auth.py` manages CDSE access tokens. `staging.py` downloads the S3 object tree, falls back to a resumable OData ZIP, validates paths and product structure, and supports existing local inputs. Temporary S3 credentials are released after a batch.

`credentials.py` supplies hidden CLI prompts. The explicit Earthdata setup command adds new entries to `.netrc`, uses mode 0600, and preserves unrelated entries. Existing Earthdata entries must be updated manually. No secret belongs in TOML, scene JSON, Git, or the data archive.

## 5. L2Gen performs atmospheric processing

`runtime.py` finds tools and keeps the OCSSW environment in its subprocesses. Python's GIS libraries use their own PROJ data. `installer.py` exposes explicit OCSSW installation with the reference's preferred tag and optional operational fallbacks. It refuses to erase an existing incomplete installation.

`windows.py` combines buffered lake bounding boxes and switches to a full scene when windows are too numerous or cover too much of the footprint. Dateline-wide footprints use a full-scene processing fallback.

`l2gen.py` preserves product selection, actual band mapping, ancillary queries and validation, parameter files, disabled L2Gen output masks, retries for unavailable optional products, logs, and validated dense L2 reuse.

## 6. Native lake science produces the permanent measurements

`science.py` owns CI-family equations, NDCI, MPH, OLCI-adapted FAI, flag decoding, wavelength lookup, and spherical coordinate conversion. `hydrolakes.py` owns vector preparation, spatial membership and metric shoreline erosion.

`compact.py` walks the L2 swath in row blocks and writes only geolocated observations inside target lakes. Raw bands and flags remain available even where the derived science is QA-excluded. It writes to a temporary NetCDF before renaming it into place.

`pipeline.py` orders staging, windows, L2Gen, compact extraction, statistics, manifests and cleanup. `ingest.py` offers a separate entry point for an existing L2Gen product. `statistics.py` computes native observation summaries. `quality.py` supplies the display/comparison mask selectors.

## 7. Exploration reads saved measurements

`archive.py` loads compact pixels, deduplicates overlapping stored locations at the reference's six-decimal coordinate precision, provides scene/lake selections and builds archive indexes.

`mapping.py` creates exact requested Web Mercator display grids and nominal-resolution geographic grids constrained to each lake. `exports.py` writes CF NetCDFs and three-file share bundles. `spectra.py` handles lake spectra and coordinate/click sampling. `comparison.py` implements v2.6.8 native-first NASA CyAN pairing, QA audit fields, CSV/NetCDF exports, coverage and bias diagnostics.

`dashboard.py` constructs the final UI only when explicitly called. `cli.py` exposes the application operations without a notebook or browser requirement.

## 8. Earlier scientific utilities remain available

`dense.py` preserves dense-swath derived-index creation and mask diagnostics. `legacy_maps.py` preserves native/dense static PNG/PDF maps, GeoTIFFs, Folium maps, lake statistics and pixel spectra utilities. `provenance.py` preserves the dense scene manifest helper. These functions support diagnostics; the global batch uses the compact branch.

`utils.py` supplies atomic JSON writes, sizes, file hashes and subprocess logging.

## Useful reading order

Start with `cli.py`, then `pipeline.py`, then `compact.py`, then `science.py`. Read a whole function by its responsibility. The short source-cell comments and `reference_manifest.json` let you trace its origin without needing the original file at runtime.

**Think about this:** if you change only a map's color scale, which files should change? Only derived visual products. The native observations, atmospheric correction and native statistics remain reusable.
