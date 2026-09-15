# Science, reference selection and parity

## Reference identity

The supplied reference was `S3_OLCI_L2Gen_CyanoLake_GLOBAL_EFFICIENT_v2_6_8_NATIVE_FIRST_CYAN_COMPARISON.ipynb`.

SHA-256: `044cbc3395645362788f0da2669c8f23d1ad9deff77d92738de21fe00fdf6c67`.

The notebook is not distributed in this package. No source loader, notebook executor, Jupyter dependency, or conversion script is part of the application. The Python modules are standalone source files. `reference_manifest.json` records where the 225 selected definitions came from and the hashes of their original source blocks.

## Which version of each function was used?

The reference contains successive redefinitions. Executing it top to bottom gives the final production behavior below. Earlier descriptions of different default masks or sinusoidal production exports do not override later code.

| Responsibility | Effective source |
|---|---|
| Ancillary retrieval and L2Gen parameter logic | Cell 17, with the extended L2 validator from cell 29 |
| CI spectral gates | Cell 30 |
| Core/extra reflectance products and MPH/FAI | Cell 48 |
| Compact extraction, bloom-aware QA and per-window statistics | Cell 57 |
| Exact display grid and geographic grid setup | Cell 58 |
| Same-lake native-first maps, complete GIS export and three-file share | Cell 59 |
| NASA CyAN native-first comparison | Cell 60 |
| Complete final dashboard | Cell 61 |

Cells are numbered from zero, as in `reference_manifest.json`. Shared helpers from earlier cells remain ordinary functions. The earlier base L2 validator is retained explicitly because the final validator calls it. Old overridden versions are not replayed at runtime.

## L2Gen settings

The default `full_swath_rhos` mode preserves `proc_ocean=2`, `proc_land=1`, `aer_opt=-99`, `brdf_opt=0`, and `gas_opt=15`. L2Gen's land, cloud, glint, solar/sensor zenith, high-radiance and straylight output masks are disabled in the same parameter block as the reference. Their flags remain available for post-L2 QA.

The code requests the required nine-band production set at nominal 490, 560, 620, 665, 681, 709, 754, 865 and 885 nm. Nominal 885 maps to the actual L2Gen name `rhos_884`. The full-spectrum/true-color helper capabilities remain, while the compact production branch keeps the final nine-band schema. Optional unavailable L2 products may be retried without that product; a required product cannot be silently dropped.

Ancillary retrieval keeps the explicit-time query, OLCI-manifest fallback query, parameter validation, logs and optional climatology policy. Time-specific retrieval is enabled and climatology fallback is disabled by default. L2 and ancillary data are external inputs whose versions affect scientific reproducibility.

## Index definitions

The equations use the actual wavelengths identified in the L2 variable names. With nominal wavelengths shown for readability:

```text
SS681 = r681 - r665 - (r709 - r665) × (681 - 665)/(709 - 665)
SS665 = r665 - r620 - (r681 - r620) × (665 - 620)/(681 - 620)

CI candidate: finite inputs and SS681 < 0
CI magnitude: -SS681 for candidates
CIcyano detection: candidate and SS665 > 0 and CI > 0.0001

NDCI = (r709 - r665)/(r709 + r665), with the reference's finite/denominator test

FAI/AFAI = r754 - the straight spectral baseline between r665 and r865 at 754 nm
MPH = maximum line height at 681, 709 or 754 nm above the 665-to-885 baseline
```

FAI is the **OLCI-adapted AFAI form**, not an assertion that OLCI has a SWIR band. MPH's actual long-wave endpoint is taken from the mapped L2 product. FLH and OLH are excluded from the final production set.

The detection operator is strictly greater than 0.0001. Finite negative `rhos` values are retained as inputs rather than removed by an added positivity threshold. Invalid inputs, QA exclusions and detection-only non-detections are encoded using the reference's NaN/fill semantics and accompanying masks.

## Primary QA and sensitivity QA

Native extraction first establishes HydroLAKES membership, finite geolocation and valid raw flags. LAND, HISATZEN and NAVFAIL are excluded from primary quantitative water. CLDICE is normally excluded, but a CLDICE pixel can be rescued when it also satisfies all three independent bloom-evidence conditions:

1. Standard CIcyano spectral detection.
2. Positive OLCI AFAI.
3. Reflectance at 754 nm greater than reflectance at 490 nm.

The source's additional shoreline/bright-screen controls and reason masks are retained. The production defaults have no shoreline erosion, no adjacency dilation and no rhos865 bright screen. The separate `cyan_strict_valid_mask` excludes CLDICE without bloom recovery; its exported label identifies it as **generic flag-exclusion sensitivity**, not NASA-equivalent cloud processing.

The earlier dense diagnostic utilities keep their earlier bloom-retaining/default union-water semantics. They are diagnostic helpers, not a substitute for inspecting the final compact product's bloom-aware masks.

The reference records missing QA flags but treats an unavailable flag as unset. Inspect `missing_flags` and flag metadata when using unusual L2 inputs. The standalone package preserves this scientific behavior rather than silently inventing a new exclusion policy.

## Native-first mapping and comparison

Each geographic target lake cell selects the nearest finite stored native observation **from the same Hylak_id**, before inspecting QA. If the selected native observation is invalid, the mapped science stays fill/NaN; it is not replaced by a farther clear pixel. Values are not averaged or interpolated.

The complete mapped product uses an exact requested **nominal** resolution in EPSG:4326. Degree spacing corresponds to the requested metre spacing at the reference latitude; geographic cells do not have equal physical width at every latitude. Web Mercator previews use exact requested projection-coordinate spacing and the reference's approximate local ground-distance handling. These are not NASA Integerized Sinusoidal products or an area-weighted l2bin/l3mapgen reconstruction.

The final NASA comparison follows the same native-first order. CSV/NetCDF outputs include chosen native coordinates, match distance, QA validity, CLDICE/HISATZEN/NAVFAIL and bloom rescue, along with coverage and detection pairing diagnostics. The default basis is primary bloom-aware QA.

NASA Merged-S3-CYAN is interpreted by the supplied reference as a same-day S3A/S3B maximum, while a local compact product represents one scene. This temporal difference and upstream processing differences limit the interpretation of correlations and biases. This code does not implement a new daily local maximum composite.

## Operational changes made for a standalone application

These changes do not introduce a new index or QA equation:

- Explicit imports and one responsibility per module replace notebook execution order. Configuration-sensitive operational defaults resolve at call time.
- TOML selects persistent/scratch/OCSSW paths. Dependency installation, credentials, OCSSW setup and dashboard launch are explicit commands.
- The application uses a single workspace writer lock and a sequential scene batch. S3 object downloading retains its thread pool.
- Configuration identity additionally records ancillary policy and canonical lake geometry; exports and dense L2 caches are separated by configuration. The original hash payload remains part of the calculation, but standalone hash values are not expected to equal old notebook-run hashes.
- Timestamped plans persist the exact scene list. Empty plans cannot accidentally process other queued scenes. CLI resume recovers interrupted `running` rows and validates saved products.
- Stats writes are atomic; complete manifests list only the current configuration's compact files. Supplied local L1 inputs outside managed scratch survive cleanup.
- OCSSW setup refuses to erase an existing incomplete installation. Failed fallback attempts are preserved when moving to another tag. `.netrc` setup preserves unrelated entries.
- Colab-only UI wording and setup side effects are replaced with local instructions. The Gradio application is private to loopback by default.
- Existing-L2 ingestion, preflight, command-line operations and configuration snapshots are added around the copied scientific routines.

Some scientific metadata schema strings intentionally retain their source patch versions (for example compact v2.6.4 and share v2.6.7). Those identify the unchanged product schema within the final v2.6.8 workflow.

## Direct comparison completed

`reference_parity_result.json` records a one-time development comparison against definitions read directly from the original reference. That development script is outside the delivered package and is not needed to run it.

- All 225 recorded original definition hashes were verified against the supplied file.
- 10,000 deterministic synthetic spectra produced exact matching CI-family, NDCI, MPH and FAI arrays/masks.
- Independent reference and standalone extraction of a miniature L2 scene produced exactly matching values and masks for **36 compact variables across 12 observations**.
- The checked same-lake native-first geographic assignment case matched exactly.

This is stronger evidence than compilation alone. It is still **synthetic parity**, not a real-scene atmospheric correction comparison, independent scientific validation, or proof of equivalence to operational NASA CyAN.

## Before scientific production

Use a real matched L1 scene, the same OCSSW tag/sensor tables/ancillary data, and the same lake geometries. Compare L2 bands, compact observations, masks, detections, statistics, native-first assignments and exported variables. Inspect known clear water, cloud, bright bloom and shoreline locations. Record differences before increasing the batch size.

The reference's dateline, very large lake, metric-erosion and global display limitations have not been expanded into a new algorithm here. The bounded synthetic tests do not cover every world geometry.
