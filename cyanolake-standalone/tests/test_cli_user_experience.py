import json
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely.geometry import box, mapping

from cyanolake import catalogue, cli, settings, workspace


def _write_adm1_cache(country="USA"):
    """Create a tiny cached ADM1 file so CLI tests do not use the network."""
    settings.ADMIN_ROOT.mkdir(parents=True, exist_ok=True)
    metadata = {
        "boundaryName": "United States of America",
        "boundaryYearRepresented": "2018",
        "boundarySource": "Synthetic ADM1 test data",
        "boundaryLicense": "Test fixture",
    }
    (settings.ADMIN_ROOT / f"geoBoundaries_{country}_ADM1_metadata.json").write_text(
        json.dumps(metadata)
    )
    gdf = gpd.GeoDataFrame(
        {
            "shapeName": ["Georgia", "Ohio"],
            "shapeISO": ["US-GA", "US-OH"],
            "shapeID": ["test-georgia", "test-ohio"],
        },
        geometry=[box(-85, 30, -80, 35), box(-84, 39, -80, 42)],
        crs=4326,
    )
    gdf.to_file(settings.ADMIN_ROOT / f"geoBoundaries_{country}_ADM1_simplified.geojson")


def test_list_regions_is_human_readable(configured, capsys):
    code = cli.main(["--config", str(configured), "list-regions"])
    output = capsys.readouterr().out
    assert code == 0
    assert "Supported [lakes].region values" in output
    assert "North America" in output
    assert 'region = "North America"' in output


def test_list_admin_prints_names_and_keys(configured, capsys):
    _write_adm1_cache()
    code = cli.main(["--config", str(configured), "list-admin", "--country", "USA"])
    output = capsys.readouterr().out
    assert code == 0
    assert "Georgia (US-GA)" in output
    assert "test-georgia" in output
    assert "--admin-name" in output
    assert "State/province list for United States of America" in output


def test_plan_accepts_readable_admin_name(configured, monkeypatch, capsys):
    _write_adm1_cache()
    lakes = gpd.GeoDataFrame(
        {"Hylak_id": [1], "Lake_name": ["Readable Lake"], "Lake_area": [20.0]},
        geometry=[box(-84, 33, -83, 34)],
        crs=4326,
    )
    catalogue.save_target_lakes(lakes)
    workspace.select_configuration(lakes)

    scene = {
        "id": "scene-georgia",
        "name": "S3A_SYNTHETIC_GEORGIA",
        "start": "2024-07-01T00:00:00Z",
        "end": "2024-07-01T00:03:00Z",
        "geofootprint": mapping(box(-85, 30, -80, 35)),
        "size_bytes": 0,
    }

    def fake_plan(*args, **kwargs):
        table = pd.DataFrame(
            [{"scene_id": scene["id"], "status": "queued", "scene_json": json.dumps(scene)}]
        )
        return table, "testhash", [scene["id"]], settings.CATALOG_DIR / "plan.csv", 1

    monkeypatch.setattr(cli, "emit", lambda value: print(json.dumps(value, default=str)))
    from cyanolake import callbacks

    monkeypatch.setattr(callbacks, "_plan_catalogue_shard", fake_plan)
    code = cli.main(
        [
            "--config",
            str(configured),
            "plan",
            "--start",
            "2024-07-01",
            "--end",
            "2024-07-01",
            "--country",
            "USA",
            "--admin-name",
            "Georgia",
            "--max-products",
            "1",
        ]
    )
    output = capsys.readouterr().out
    assert code == 0
    payload = json.loads(output)
    saved_plan = json.loads(Path(payload["plan"]).read_text())
    assert saved_plan["scene_ids"] == ["scene-georgia"]
    assert saved_plan["shard"] == "Georgia"
