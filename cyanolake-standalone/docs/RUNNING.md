# Running CyanoLake, step by step

## Before running a real scene

Use Python 3.11 or newer; the local test environment used Python 3.12. The procedure below uses Conda and works as a Linux deployment procedure. Transfer this entire package folder to your server if processing there. Keep your credentials on that machine; do not put them in chat or configuration files.

OCSSW is a separate NASA science installation. Installing this Python package does not install L2Gen. NASA documents the installer, operational `V` tags, and Sentinel-3A/B options in its [official installation instructions](https://www.earthdata.nasa.gov/data/tools/seadas/installers-source-code). The package carries the reference's preferred `V2026.3` tag; availability is checked when you explicitly run the installer. It is not asserted to be the latest tag.

**Two accounts serve different purposes:** CDSE credentials download Sentinel-3 L1 data. Earthdata credentials retrieve NASA ancillary data and optionally NASA CyAN rasters. Public catalogue discovery needs no CDSE password.

## 1. Enter the project folder

`cd` means **change directory**. Quotes keep a path containing spaces together. `/absolute/path/to/cyanolake-standalone` is a placeholder: replace it with the actual folder on your Linux machine.

```bash
cd "/absolute/path/to/cyanolake-standalone"
```

On the Mac where this package was created, the actual folder is:

```text
/Users/madhusudhanreddychencharapu/Documents/ChatGPT/CyanoTrend_DataExtraction/cyanolake-standalone
```

## 2. Create and activate a Python environment

`conda create` creates an environment. `-n cyanolake` gives it the name **cyanolake**. `-c conda-forge` selects the package channel. `python=3.12` requests that Python version; `pip` installs Python's package installer. `conda activate` selects this environment for subsequent commands. If that environment already exists, activate it instead of creating it again.

```bash
conda create -n cyanolake -c conda-forge python=3.12 pip
conda activate cyanolake
```

`python` is the active environment's interpreter. `-m pip` runs that interpreter's pip module. `install` installs a package. `-e` means **editable**: imports use the source files in this folder. `.` means the current folder. `[dev]` adds the testing tools. Quotes prevent the shell from interpreting the square brackets.

```bash
python -m pip install -e ".[dev]"
```

To install a deployment without test tools, use the same command with `.` instead of `".[dev]"`.

## 3. Create your local configuration

`cp` means **copy**. The first path is the source; the second is the destination. The ignored `local.toml` file holds your machine's paths without modifying the supplied defaults.

```bash
cp config/default.toml config/local.toml
```

Edit `config/local.toml` with your preferred editor. The most important entries are:

```toml
[paths]
persistent_root = "../data"
scratch_root = "../scratch"
ocssw_root = "~/ocssw"
```

These defaults are immediately usable when sufficient disk space exists. Relative paths resolve from the **configuration file's directory**, not the terminal's current directory. On a server, replace them with dedicated absolute paths on suitable disks. The three roots must be separate and must not contain one another.

Read the `[storage]` section before processing. Successful runs delete their managed staged L1, fallback ZIP and dense L2 by default. Failed-run scratch remains. For a first real-scene validation, you may set `keep_l2 = true` so you can inspect L2Gen's output afterward.

The final production QA and nine-band compact schema are preserved. This TOML deliberately exposes operational choices and the original shoreline buffer; it does not add a new scientific cloud-mask algorithm.

## 4. Initialize the folders

`cyanolake` is the command installed by pip. `--config` selects the TOML file. Global options go **before** the subcommand. `init` creates the configured directory tree and initializes the SQLite registry.

```bash
cyanolake --config config/local.toml init
```

If the shell cannot find `cyanolake`, `python -m cyanolake` is an equivalent entry point after installation. `-m` means run the named Python module.

## 5. Run the offline tests

`python -m pytest` runs the tests with the active Python interpreter. `-q` means **quiet**: print a short result. These tests use synthetic observations and mocked external processing. They do not download a satellite scene or authenticate to NASA/CDSE.

```bash
python -m pytest -q
```

Passing tests establishes the checked software behavior. It does not establish that your OCSSW installation or a real scene is scientifically correct.

## 6. Install or select NASA OCSSW

If you already have an OLCI-capable installation, point `ocssw_root` to it and proceed to the next section.

`install-ocssw` explicitly downloads NASA's installer and requests SeaDAS processing programs plus OLCI sensor data for S3A and S3B. `--tag V2026.3` supplies the reference's preferred version. `--fallback-tags 0` disables fallback, so a missing preferred tag fails visibly rather than selecting another version.

```bash
cyanolake --config config/local.toml install-ocssw --tag V2026.3 --fallback-tags 0
```

To permit the reference's fallback behavior, change the last number to `3`. That allows the three newest operational tags returned by NASA's installer. Record the actual installed tag and review parity if it differs from the reference run. Detailed installation logs remain in the sibling `ocssw-installer` directory.

The installer options used internally mean:

| Option | Meaning |
|---|---|
| `--list_tags` | Ask NASA which versions are available |
| `--install_dir` | Destination root for OCSSW |
| `--tag` | Requested software/data version |
| `--seadas` | Processing executables and required supporting bundles |
| `--olcis3a` | Sentinel-3A OLCI sensor data |
| `--olcis3b` | Sentinel-3B OLCI sensor data |
| `--verbose` | Detailed installer log messages |

An existing incomplete root is preserved; the command asks you to choose a new root or repair it. Installation is a network operation and can take substantial time. On Apple Silicon, consult NASA's current platform prerequisites rather than assuming that Linux binaries will work.

You do not need to source OCSSW into the shell running this application: the application reads its environment file for OCSSW subprocesses and isolates it from Python's GIS libraries.

## 7. Configure Earthdata

`configure-earthdata` prompts for your NASA username and a hidden password. It adds the NASA machines to `.netrc` with owner-only permissions and preserves unrelated entries. If Earthdata entries already exist, it stops so you can update those entries deliberately.

```bash
cyanolake --config config/local.toml configure-earthdata
```

If you already configured `.netrc`, keep it and run the check below. `preflight` checks local tool/data paths and credential-file structure. `--without-lakes` skips the target-lake requirement because you have not prepared HydroLAKES yet.

```bash
cyanolake --config config/local.toml preflight --without-lakes
```

Proceed only when the report has `"ok": true`. This check verifies credential presence and permissions, not whether the remote server will accept the password.

For a deliberate climatology experiment, `[processing] get_ancillary = false` skips getanc. This changes the processing identity. The default remains time-specific ancillary retrieval with automatic climatology fallback disabled.

## 8. Prepare the canonical HydroLAKES universe

`prepare-lakes` downloads/prepares HydroLAKES and filters by `[lakes] min_area_km2`. The default is WORLD, with a minimum area of 10 km². This may require a large initial reference-data download.

To see the valid lake-universe values from the backend itself, run:

```bash
cyanolake --config config/local.toml list-regions
```

Edit `[lakes] region` in `config/local.toml` before `prepare-lakes` if you want a smaller lake universe. Valid values are `WORLD`, `USA`, `Africa`, `Asia`, `Europe`, `North America`, `South America`, and `Oceania`.

```bash
cyanolake --config config/local.toml prepare-lakes
```

If you already have HydroLAKES, `--source` supplies its local vector path instead of downloading it. Replace the example path below with your real file. A local source is taken as the supplied lake universe after area filtering; choose a WORLD file when you want worldwide scene deduplication.

```bash
cyanolake --config config/local.toml prepare-lakes --source "/absolute/path/to/HydroLAKES.gpkg"
```

Run **one** of those alternatives. Preparing a different lake universe changes its geometry fingerprint and processing hash.

After changing the minimum lake area or region in TOML, rerun `prepare-lakes` before making a new plan. A cached target vector contains the universe you previously prepared; changing TOML cannot restore polygons absent from that vector.

Now `preflight` also requires the saved target-lake file:

```bash
cyanolake --config config/local.toml preflight
```

## 9. Plan a small first batch

`plan` discovers scenes and saves a timestamped JSON plan and CSV queue snapshot. `--start` and `--end` are inclusive UTC dates. `--bbox` is four numbers in **west, south, east, north** order. Negative longitude is west. `--max-products 3` caps catalogue results for a small test; it is not a guarantee of three lake-intersecting scenes.

This is an illustrative one-day Lake Erie-area search. Change the dates and area for your research.

```bash
cyanolake --config config/local.toml plan --start 2024-07-01 --end 2024-07-01 --bbox -84 41 -78 43 --max-products 3
```

Read `scene_count` and the queue CSV. A zero-scene plan exits with code 2; change the search deliberately. The plan command prints the **actual JSON plan path**. Keep this path for processing and resume.

The bbox selects scenes to discover. Each selected scene is processed against the canonical lake universe, so all eligible lakes in that scene can be retained.

For worldwide discovery, omit `--bbox`. Start with a short date range and a bounded batch. For a 5° planning cell, `--grid-id` replaces `--bbox`; the example identifies the cell from 40°N/85°W to 45°N/80°W.

```bash
cyanolake --config config/local.toml plan --start 2024-07-01 --end 2024-07-01 --grid-id G05_N40_W085
```

For an administrative region, `list-admin` prints state/province names and their exact internal keys. `--country` is an ISO3 code. For normal use, plan with `--admin-name` and the readable name. For exact reproducible reruns, use `--admin-key` and copy the key printed beside that name.

```bash
cyanolake --config config/local.toml list-admin --country USA
cyanolake --config config/local.toml plan --start 2024-07-01 --end 2024-07-01 --country USA --admin-name "Georgia" --max-products 3
```

The key-based form is still available when you need to pin one exact boundary row:

```bash
cyanolake --config config/local.toml plan --start 2024-07-01 --end 2024-07-01 --country USA --admin-key "KEY_FROM_LIST_ADMIN"
```

## 10. Process one scene

`PLAN_JSON` is a shell variable you choose for convenience. The value below is a placeholder; replace it with the exact path printed by `plan`. No spaces surround `=`. `"$PLAN_JSON"` expands that variable while keeping its path together.

`run` processes a saved plan. `--plan` selects that JSON file. `--max-scenes 1` limits this invocation to one eligible scene. CDSE credentials are requested in the terminal; the password is hidden. Existing `CDSE_USERNAME`/`CDSE_PASSWORD` environment variables are also supported, but do not type passwords into shell history.

```bash
PLAN_JSON="/actual/path/printed/by/plan.json"
cyanolake --config config/local.toml run --plan "$PLAN_JSON" --max-scenes 1
```

Keep the process running until it completes. The data lifecycle is staging → getanc → L2Gen → compact extraction → native statistics → manifest → successful-run cleanup.

`status` shows the current configuration's scene IDs, state, attempts and errors. `stats` writes a consolidated CSV from the per-window native statistics.

```bash
cyanolake --config config/local.toml status
cyanolake --config config/local.toml stats
```

## 11. Continue or resume

The same `run` command continues the frozen plan. `--max-scenes 5` attempts at most five eligible scenes this time. It does not mean five workers. Processing remains sequential; S3 staging can use several object-download threads.

```bash
cyanolake --config config/local.toml run --plan "$PLAN_JSON" --max-scenes 5
```

`--retry-failed` includes failed scenes after you resolve their cause. Valid completed windows can be reused. Rows left `running` after a stopped process are eligible for recovery under the workspace lock.

```bash
cyanolake --config config/local.toml run --plan "$PLAN_JSON" --max-scenes 1 --retry-failed
```

If the configuration or lake universe changed, the saved plan is rejected with a hash mismatch. Restore the original configuration and target lakes from `catalog/configurations/<hash>/`, or create a new plan for the new setup. This prevents accidental reuse under different processing settings.

## 12. Export a completed scene

`SCENE_ID` stores an actual ID from `status`. Replace the example value; it is a placeholder. `export` creates the three-file share ZIP. `--scene-id` identifies the saved scene. `--resolution 300` requests the reference's nominal 300 m geographic GIS grid. There is no silent coarsening.

```bash
SCENE_ID="ACTUAL_SCENE_UUID_FROM_STATUS"
cyanolake --config config/local.toml export --scene-id "$SCENE_ID" --resolution 300
```

`map` writes a portable all-index HTML map. A 1000 m Web Mercator display grid is lighter for browsing; this changes only the visualization. `archive` writes an orbit/frame/file index for the current configuration.

```bash
cyanolake --config config/local.toml map --scene-id "$SCENE_ID" --resolution 1000
cyanolake --config config/local.toml archive
```

Open the generated HTML in a browser. Basemap tiles require internet access. The share NetCDF can be opened in SNAP/GIS. Its statistics still come from native compact observations, not from counting mapped cells.

## 13. NASA CyAN comparison

`compare-cyan` compares one stored scene and lake. `--lake-id` is a real `Hylak_id`, not a scene ID. `--source` uses an existing NASA OBPG CIcyano GeoTIFF; omitting it requests the same-day NASA download through Earthdata. `--max-distance 600` sets the maximum native match distance in metres. `--qa primary` selects bloom-aware QA. `strict` is only the generic flag-exclusion sensitivity case.

Replace `1234` and the path with your real lake and NASA file.

```bash
cyanolake --config config/local.toml compare-cyan --scene-id "$SCENE_ID" --lake-id 1234 --source "/absolute/path/to/NASA_CYAN.tif" --max-distance 600 --qa primary
```

The nearest same-lake native observation is chosen first; its QA is then checked. Pair tables record its coordinates, distance, selected QA, flags and bloom rescue. NASA's merged daily product and one local scene have different temporal support: correlation is a comparison diagnostic, not proof of algorithm equivalence.

## 14. Use an existing local input

`process-local` uses a local `.SEN3` or ZIP. `--input` gives that path. `--scene-json` supplies metadata saved under a plan's `scenes/` directory. The processor preserves original inputs outside its managed scratch tree.

```bash
cyanolake --config config/local.toml process-local --input "/absolute/path/to/PRODUCT.SEN3" --scene-json "/absolute/path/to/scene.json"
```

`ingest-l2` starts from an existing L2Gen NetCDF and performs compact extraction/statistics without rerunning atmospheric correction. Use matching scene metadata. Only use this route when you know how that L2 was generated; its existence does not prove reference-compatible parameters.

```bash
cyanolake --config config/local.toml ingest-l2 --input "/absolute/path/to/L2GEN.nc" --scene-json "/absolute/path/to/scene.json"
```

Scene JSON requires `id`, `name`, `start`, `end`, and a GeoJSON `geofootprint`. `plan` saves these records for you. `diagnose-l2` runs the preserved dense mask diagnostic on a local L2 file:

```bash
cyanolake --config config/local.toml diagnose-l2 --input "/absolute/path/to/L2GEN.nc"
```

## 15. Launch the dashboard

`dashboard` constructs and starts the full Gradio application. `--host 127.0.0.1` listens only on this machine. `--port 7860` selects the local TCP port. Open the printed local URL. Use Ctrl+C in the terminal to stop it. The dashboard owns the workspace lock while running; stop it before running CLI writers against the same workspace.

```bash
cyanolake --config config/local.toml dashboard --host 127.0.0.1 --port 7860
```

The seven tabs cover planning/processing, time series, maps, spectra/statistics, NASA comparison, global archive/share, and storage/provenance. The dashboard has the reference's additional interactive controls and downloads. It does not create a public Gradio share tunnel.

For a Linux server, leave the dashboard listening on loopback and use SSH forwarding from your own computer. `ssh` opens an encrypted connection. `-N` means do not run a remote shell command. `-L 7860:127.0.0.1:7860` forwards your local port to the server's loopback port. `YOUR_USER@YOUR_SERVER` is a placeholder. Keep this connection open and browse `http://127.0.0.1:7860` locally.

```bash
ssh -N -L 7860:127.0.0.1:7860 YOUR_USER@YOUR_SERVER
```

## Common failures

| Message / symptom | Next step |
|---|---|
| `l2gen is missing` | Correct `ocssw_root` or install the science processors |
| Missing S3A/S3B sensor data | Install both OLCI bundles |
| Earthdata 401/403 | Check the local account, permissions and NASA authorization; do not paste credentials into logs/chat |
| No target lakes | Prepare the canonical vector and inspect its area filter |
| Zero discovered scenes | Inspect dates, bbox/ADM1 and catalogue limit |
| Required L2 product unavailable | Inspect installed tag and L2Gen log; the code does not silently discard required science bands |
| Exact map grid exceeds its cap | Select a smaller scope or deliberately request a coarser display resolution |
| Another process owns workspace | Stop the other dashboard/CLI writer or use a separate workspace |
| Config hash differs | Restore the saved setup or plan a new configuration |

Exit code 0 means the requested command completed. Exit code 2 means an actionable setup, planning or processing failure. Ctrl+C returns 130. `run` can complete a bounded batch while more queued scenes remain; inspect `status`.

For built-in help, `--help` prints usage and exits:

```bash
cyanolake --help
cyanolake run --help
```

**Think before scaling up:** after the first real scene, can you explain why a particular pixel is retained, cloud-excluded, bloom-rescued, or a valid non-detection? Inspect its native flags, bands and masks before increasing the batch size.
