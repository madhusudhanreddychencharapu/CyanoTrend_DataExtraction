# CyanoLake — standalone OLCI lake processing

This package implements the **effective final v2.6.8 workflow** from the supplied scientific reference as ordinary Python modules. It does not contain, open, import, or execute a notebook. No Jupyter runtime is required.

## Start here

1. [Run guide](docs/RUNNING.md): Linux/Conda setup, every command and flag explained before use, first scene, resume, exports, and dashboard.
2. [Architecture](docs/ARCHITECTURE.md): what each module is responsible for.
3. [Data lifecycle](docs/DATA_LAYOUT.md): where files go, which files remain, and what can be reproduced.
4. [Science and reference parity](docs/SCIENCE_PARITY.md): equations, final patch selection, scientific caveats, and operational changes.
5. [Validation](docs/VALIDATION.md): executed checks and the remaining real-scene validation boundary.

## Workflow

```text
CDSE catalogue + canonical HydroLAKES
          ↓
saved plan + SQLite scene registry
          ↓
S3 object staging / ZIP fallback / local .SEN3
          ↓
getanc → L2Gen over adaptive lake windows or the full scene
          ↓
native lake observations → bloom-aware QA → CI / CIcyano / NDCI / MPH / FAI
          ↓
compact NetCDF + lake statistics + provenance
          ↓
maps, spectra, same-lake native-nearest GIS export, NASA CyAN comparison
```

Native lake observations are the numerical source for statistics. Mapped grids are derived products. The nearest same-lake source is selected **before QA**: an invalid nearest source stays no-data.

## Interfaces

The command line supports setup, planning, bounded sequential batches, resume, local inputs, existing L2 ingestion, diagnostics, archive inspection, share bundles, HTML maps, and NASA comparison. The complete Gradio dashboard retains ADM1/bbox/grid planning, map controls, spectra sampling, archive downloads, and comparison diagnostics.

For backend-only use, start with the discoverability commands before processing:

```bash
python -m cyanolake --config config/local.toml guide
python -m cyanolake --config config/local.toml list-regions
python -m cyanolake --config config/local.toml list-admin --country USA
```

Use readable state/province names in planning when possible:

```bash
python -m cyanolake --config config/local.toml plan --start 2024-07-01 --end 2024-07-01 --country USA --admin-name Georgia --max-products 3
```

The three-file share ZIP contains a complete science NetCDF, all-index native statistics CSV, and metadata JSON. The package preserves the nine-band production archive and its final bloom-aware QA behavior.

**Validation boundary:** local synthetic tests verify the implemented behavior. An authenticated real-scene L2Gen run and numerical comparison against the reference's real outputs are still required before calling a deployment scientifically validated.
