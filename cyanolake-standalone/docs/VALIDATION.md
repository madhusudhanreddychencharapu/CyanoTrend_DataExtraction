# Validation record

## Executed locally

The standalone test suite currently contains **22 passing tests**, using Python 3.12 on macOS ARM64. Tests run without downloading real scenes or providing credentials.

Coverage includes:

- CI gates, strict detection threshold, finite negative inputs, NDCI denominator handling, MPH/FAI and the 885→884 L2 band mapping.
- Synthetic L2 extraction into the native compact schema, dual QA masks, cloud/bloom rescue, valid non-detection fill semantics, lake IDs and statistics.
- Complete NetCDF export, signed-byte classic-format masks, primary-QA fill behavior and the three-file share ZIP.
- Same-lake native-first assignment, including a nearest invalid source and a closer observation from another lake.
- Portable HTML maps, coordinate spectra and synthetic NASA GeoTIFF comparison under primary and strict QA.
- Configuration errors, geometry/configuration fingerprints, registry deduplication, catalogue pagination and inclusive end dates.
- ZIP path traversal rejection, managed scratch cleanup, preservation of supplied local inputs, empty-plan isolation and recovery of interrupted scenes.
- L2Gen parameter generation with a mocked executable; successful scene orchestration, cleanup and reuse.
- Separate CLI processes performing existing-L2 ingestion, status inspection and share export.
- Construction of the final Gradio dashboard: 181 components and 28 event handlers.

The test warnings observed are dependency deprecations from Rasterio/Affine and Gradio's legacy row-count form; they do not change the checked numerical results.

## Direct reference comparison

The one-time development audit verified 225 source-definition hashes, compared 10,000 synthetic spectra, compared all 36 compact variables/masks for a 12-observation L2 scene, and compared the checked native-first geographic assignment. Exact values and masks matched. See [reference_parity_result.json](reference_parity_result.json).

The shipped tests do not need the original notebook. The notebook audit/conversion tools are not distributed with the package.

## Additional checks

- Static undefined-name/import checks and source formatting.
- Package import and CLI entry-point checks.
- Built and installed a Python wheel; its `cyanolake --help` entry point worked outside the source directory. The wheel contains 37 Python modules and no notebook files.
- Started the dashboard on localhost, successfully fetched its HTTP configuration (181 components, 28 handlers), and shut down the temporary server.
- No notebook file, notebook loader, `exec`/`eval` runtime, IPython or Colab imports in application code.

The exact selected dependency versions used for validation are recorded in `validation_environment.json`. They describe this test environment, not a claim that every Linux/HPC environment has been validated.

## Not executed here

- Installing NASA OCSSW on the user's processing system.
- Authenticating CDSE or Earthdata, live S3/ZIP scene downloads, or ancillary retrieval.
- A real Sentinel-3 L1 → L2Gen → compact archive run.
- Numerical parity against a real reference scene, real NASA CyAN rasters or field observations.
- Full-world scalability, peak-disk/memory benchmarking, polar/dateline science validation or multi-user dashboard operation.

The implementation is runnable standalone code with tested synthetic behavior. Scientific deployment acceptance still requires the real-scene checks described in [SCIENCE_PARITY.md](SCIENCE_PARITY.md).

## Reproduce the shipped tests

`python -m pytest` runs the tests with the active interpreter; `-q` requests concise output. Run this from the installed package's root directory:

```bash
python -m pytest -q
```
