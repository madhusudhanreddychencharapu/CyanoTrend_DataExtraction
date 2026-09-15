import json
import zipfile

import netCDF4
import numpy as np
import pandas as pd

from cyanolake import compact, exports, mapping, pipeline


def test_compact_qa_fill_and_native_storage(ingested):
    result, lakes, scene, ch = ingested
    assert result["observations"] == 12
    assert compact.validate_compact_product(result["compact"], scene_id=scene["id"], config_hash=ch)
    with netCDF4.Dataset(result["compact"]) as ds:
        assert set(ds.dimensions) == {"obs"}
        assert ds["valid_water_mask"][:6].tolist() == [1, 1, 0, 0, 0, 0]
        assert ds["bloom_rescue_mask"][:6].tolist() == [0, 1, 0, 0, 0, 0]
        assert ds["cyan_strict_valid_mask"][:6].tolist() == [1, 0, 0, 0, 0, 0]
        assert np.ma.getmaskarray(ds["CI_cyano"][:])[2:6].all()
        assert np.ma.getmaskarray(ds["CI_cyano"][:])[6]
        assert ds["valid_water_mask"][6] == 1
        assert ds["rhos_620"][7] < 0
        assert ds["ci_valid_mask"][7] == 1
        assert ds["rhos_885"].shape == (12,)
    assert pipeline.scene_outputs_valid(scene["id"], ch)
    stats = pd.read_parquet(result["statistics"])
    assert len(stats) == 2 and set(stats.Hylak_id) == {1, 2}


def test_all_science_export_and_share_bundle(ingested):
    _, _, scene, _ = ingested
    path, g = exports.export_complete_scene_science_netcdf(scene["id"], 300)
    assert g["resolution_m"] == 300
    with netCDF4.Dataset(path) as ds:
        assert ds.data_model == "NETCDF4_CLASSIC"
        assert "rhos_884" in ds.variables and "rhos_885" not in ds.variables
        assert {
            "CI",
            "CI_cyano",
            "NDCI",
            "MPH",
            "FAI",
            "nearest_Hylak_id",
            "native_match_distance_m",
        } <= set(ds.variables)
        assert ds["valid_water_mask"].dtype == np.dtype("int8")
        assert ds["valid_water_mask"]._FillValue == -1
        invalid = np.ma.filled(ds["valid_water_mask"][:], 0) == 0
        assert np.ma.getmaskarray(ds["CI_cyano"][:])[invalid].all()
    bundle, _ = exports.create_scene_share_bundle(scene["id"], 300)
    with zipfile.ZipFile(bundle) as z:
        assert len(z.namelist()) == 3
        assert sorted(n.rsplit(".", 1)[1] for n in z.namelist()) == ["csv", "json", "nc"]
        meta = json.loads(z.read(next(n for n in z.namelist() if n.endswith(".json"))))
        assert meta["statistics"]["mapped_raster_not_used_for_statistics"] is True


def test_same_lake_native_first_assignment(ingested):
    _, lakes, scene, _ = ingested
    frame = pd.DataFrame(
        {
            "latitude": [40.001, 40.001, 40.001],
            "longitude": [-85.005, -85.003, -85.00501],
            "Hylak_id": [1, 1, 2],
            "valid_water_mask": [0, 1, 1],
            "CI_cyano": [np.nan, 0.002, 0.5],
        }
    )
    g = {
        "latitude": np.array([40.001]),
        "longitude": np.array([-85.005]),
        "ny": 1,
        "nx": 1,
        "resolution_m": 300,
    }
    rows, cols, src, dist = mapping._native_assignments_geographic_v267(frame, g, scene["id"])
    assert src.tolist() == [0]
    assert np.isnan(frame.iloc[src[0]].CI_cyano) and dist[0] < 1


def test_portable_map_and_spectra(ingested):
    _, _, scene, _ = ingested
    path, g = mapping.create_complete_scene_html(scene["id"], 300)
    assert "CI_cyano" in open(path).read()
    from cyanolake.spectra import coordinate_spectra_callback

    result = coordinate_spectra_callback(
        scene["id"], pd.DataFrame({"label": ["A"], "latitude": [40.001], "longitude": [-85.009]})
    )
    assert result[2] is not None, result[-1]
    assert len(result[1]) > 0
