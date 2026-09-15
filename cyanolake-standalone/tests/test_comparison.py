from pathlib import Path

import netCDF4
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import from_origin

from cyanolake.comparison import compare_cyan_callback_v265


def test_nasa_comparison_keeps_rejected_nearest_pixel(ingested, tmp_path):
    _, _, scene, _ = ingested
    tif = tmp_path / "synthetic_nasa_dn.tif"
    with rasterio.open(
        tif,
        "w",
        driver="GTiff",
        width=11,
        height=6,
        count=1,
        dtype="uint8",
        crs="EPSG:4326",
        transform=from_origin(-85.012, 40.010, 0.002, 0.002),
        nodata=255,
    ) as ds:
        ds.write(np.full((1, 6, 11), 180, dtype=np.uint8))
    outputs = compare_cyan_callback_v265(scene["id"], "1", "", "upload", str(tif), 600, "primary")
    assert outputs[3] is not None, outputs[-1]
    pairs = pd.read_csv(outputs[3])
    invalid = pairs.local_match_found & ~pairs.nearest_selected_QA_valid
    assert invalid.any()
    assert not pairs.loc[invalid, "local_ci_cyano_detected"].any()
    assert not pairs.loc[invalid, "paired_above_0_0001_for_correlation"].any()
    with netCDF4.Dataset(outputs[5]) as ds:
        assert len(ds.variables) > 10
    strict = compare_cyan_callback_v265(scene["id"], "1", "", "upload", str(tif), 600, "strict")
    assert strict[3] is not None, strict[-1]
    strict_pairs = pd.read_csv(strict[3])
    assert np.allclose(
        pairs.nearest_native_latitude, strict_pairs.nearest_native_latitude, equal_nan=True
    )
    assert np.allclose(
        pairs.nearest_native_longitude, strict_pairs.nearest_native_longitude, equal_nan=True
    )
    assert strict_pairs.nearest_selected_QA_valid.sum() <= pairs.nearest_selected_QA_valid.sum()
