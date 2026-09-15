import json
import zipfile
from pathlib import Path

import pandas as pd
import pytest

from cyanolake import catalogue, configuration, pipeline, registry, settings, staging, workspace


def test_invalid_configuration_stops_before_directory_creation(tmp_path):
    p = tmp_path / "bad.toml"
    p.write_text('[paths]\npersistent_root="/"\n')
    with pytest.raises(ValueError, match="dedicated"):
        configuration.configure(p, create=True)
    p.write_text('[storage]\nkeep_l2="false"\n')
    with pytest.raises(ValueError, match="true or false"):
        configuration.configure(p)
    p.write_text("[processing]\nmax_windows=0\n")
    with pytest.raises(ValueError, match="positive integer"):
        configuration.configure(p)


def test_registry_replanning_does_not_reset_completed_scene(lake_scene):
    lakes, scene, ch = lake_scene
    registry.registry_upsert_scene(scene, ch, len(lakes))
    registry.registry_update(scene["id"], ch, status="done", attempts=2)
    registry.registry_upsert_scene(scene, ch, len(lakes))
    table = registry.registry_table(ch)
    assert len(table) == 1 and table.iloc[0].status == "done" and table.iloc[0].attempts == 2


def test_hash_changes_for_science_or_lake_geometry(lake_scene, monkeypatch):
    lakes, scene, ch = lake_scene
    assert registry.processing_config_hash(10, "WORLD") == ch
    monkeypatch.setattr(settings, "GET_ANCILLARY", True)
    assert registry.processing_config_hash(10, "WORLD") != ch
    shifted = lakes.copy()
    shifted.geometry = shifted.geometry.translate(xoff=0.001)
    assert workspace.lake_fingerprint(lakes) != workspace.lake_fingerprint(shifted)


def test_zip_rejects_path_traversal(tmp_path):
    p = tmp_path / "bad.zip"
    with zipfile.ZipFile(p, "w") as z:
        z.writestr("../escaped.txt", "bad")
    with pytest.raises(RuntimeError, match="Unsafe ZIP"):
        staging.safe_extract_zip(p, tmp_path / "staged")
    assert not (tmp_path / "escaped.txt").exists()


def test_local_inputs_survive_cleanup(configured, tmp_path):
    source = tmp_path / "original.SEN3"
    source.mkdir()
    (source / "valuable.txt").write_text("keep")
    zipped = tmp_path / "original.zip"
    zipped.write_bytes(b"keep")
    pipeline._cleanup_scene_scratch(zipped, source, [], success=True)
    assert source.exists() and zipped.exists()
    owned = settings.STAGING_DIR / "staged.SEN3"
    owned.mkdir()
    pipeline._cleanup_scene_scratch(None, owned, [], success=False)
    assert owned.exists()
    pipeline._cleanup_scene_scratch(None, owned, [], success=True)
    assert not owned.exists()


def test_empty_plan_cannot_run_other_queued_scenes(lake_scene, monkeypatch):
    lakes, scene, ch = lake_scene
    registry.registry_upsert_scene(scene, ch, len(lakes))

    def fail(*a, **k):
        pytest.fail("An empty plan must not stage any scene")

    monkeypatch.setattr(pipeline, "process_scene_once", fail)
    monkeypatch.setattr(staging, "release_temporary_s3_credentials", lambda *a: None)
    result = pipeline.process_next_queued(ch, lakes, "", "", scene_ids=[])
    assert result.empty


def test_running_scene_is_resumable_after_process_exit(lake_scene, monkeypatch):
    lakes, scene, ch = lake_scene
    registry.registry_upsert_scene(scene, ch, len(lakes))
    registry.registry_update(scene["id"], ch, status="running")
    seen = []

    def process(s, *args, **kwargs):
        seen.append(s["id"])
        return pd.DataFrame()

    monkeypatch.setattr(pipeline, "process_scene_once", process)
    monkeypatch.setattr(staging, "release_temporary_s3_credentials", lambda *a: None)
    pipeline.process_next_queued(ch, lakes, "", "", scene_ids=[scene["id"]])
    assert seen == [scene["id"]]


def test_missing_compact_invalidates_done_manifest(ingested):
    result, lakes, scene, ch = ingested
    assert pipeline.scene_outputs_valid(scene["id"], ch)
    result["compact"].unlink()
    assert not pipeline.scene_outputs_valid(scene["id"], ch)


def test_catalogue_pagination_preserves_last_day(monkeypatch):
    calls = []

    class Response:
        def __init__(self, data):
            self.data = data

        def raise_for_status(self):
            pass

        def json(self):
            return self.data

    def get(url, params=None, **kwargs):
        calls.append((url, params))
        if len(calls) == 1:
            return Response(
                {
                    "value": [{"Id": "a", "Name": "S3A_test.SEN3"}],
                    "@odata.nextLink": "/odata/v1/Products?next=2",
                }
            )
        return Response({"value": [{"Id": "b", "Name": "S3B_test.SEN3"}]})

    monkeypatch.setattr(catalogue.requests, "get", get)
    result = catalogue.search_olci_l1_catalog("2024-01-01", "2024-01-02")
    assert [r["id"] for r in result] == ["a", "b"]
    assert "2024-01-03T00:00:00.000Z" in calls[0][1]["$filter"]
    assert calls[1][1] is None
