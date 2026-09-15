import geopandas as gpd
import netCDF4
import numpy as np
import pytest
from shapely.geometry import box, mapping

from cyanolake import catalogue, configuration, workspace


@pytest.fixture
def configured(tmp_path):
    cfg = tmp_path / "settings.toml"
    cfg.write_text(f'''[paths]
persistent_root = "{tmp_path}/data"
scratch_root = "{tmp_path}/scratch"
ocssw_root = "{tmp_path}/ocssw"
[processing]
get_ancillary = false
rows_per_block = 1
''')
    configuration.configure(cfg, create=True)
    return cfg


@pytest.fixture
def lake_scene(configured):
    lakes = gpd.GeoDataFrame(
        {"Hylak_id": [1, 2], "Lake_name": ["Test A", "Test B"], "Lake_area": [20.0, 12.0]},
        geometry=[box(-85.012, 39.999, -85.000, 40.010), box(-85.000, 39.999, -84.990, 40.010)],
        crs=4326,
    )
    catalogue.save_target_lakes(lakes)
    ch = workspace.select_configuration(lakes)
    scene = {
        "id": "00000000-0000-4000-8000-000000000001",
        "name": "S3A_OL_1_EFR____20240101T120000_20240101T120300_20240102T120000_0180_108_012_2520_LN1_O_NT_002",
        "start": "2024-01-01T12:00:00Z",
        "end": "2024-01-01T12:03:00Z",
        "geofootprint": mapping(box(-85.02, 39.99, -84.98, 40.02)),
        "size_bytes": 0,
    }
    return lakes, scene, ch


@pytest.fixture
def synthetic_l2(tmp_path, lake_scene):
    path = tmp_path / "synthetic_l2.nc"
    with netCDF4.Dataset(path, "w") as ds:
        ds.createDimension("number_of_lines", 2)
        ds.createDimension("pixels_per_line", 6)
        dims = ("number_of_lines", "pixels_per_line")
        nav = ds.createGroup("navigation_data")
        geo = ds.createGroup("geophysical_data")
        nav.createVariable("latitude", "f4", dims)[:] = np.tile([[40.001], [40.004]], (1, 6))
        nav.createVariable("longitude", "f4", dims)[:] = np.tile(
            [-85.009, -85.007, -85.005, -85.003, -84.997, -84.995], (2, 1)
        )
        fv = geo.createVariable("l2_flags", "u4", dims)
        fv.flag_meanings = "LAND CLDICE HISATZEN NAVFAIL"
        fv.flag_masks = np.array([1, 2, 4, 8], dtype=np.uint32)
        fv[:] = np.array([[0, 2, 2, 4, 1, 8], [0, 0, 0, 0, 0, 0]], dtype=np.uint32)
        vals = {
            490: 0.01,
            560: 0.02,
            620: 0.01,
            665: 0.03,
            681: 0.02,
            709: 0.04,
            754: 0.05,
            865: 0.02,
            884: 0.01,
        }
        for w, value in vals.items():
            data = np.full((2, 6), value, np.float32)
            if w == 754:
                data[0, 2] = 0.001
            if w == 681:
                data[1, 0] = 0.05
            if w == 620:
                data[1, 1] = -0.01
            geo.createVariable(f"rhos_{w}", "f4", dims, fill_value=-32767)[:] = data
    return path


@pytest.fixture
def ingested(synthetic_l2, lake_scene):
    from cyanolake.ingest import ingest_l2

    lakes, scene, ch = lake_scene
    result = ingest_l2(synthetic_l2, scene, lakes, ch)
    return result, lakes, scene, ch
