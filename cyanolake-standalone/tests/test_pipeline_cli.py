import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from cyanolake import l2gen, pipeline, registry, runtime, settings, staging, windows


def test_l2gen_parameters_preserve_reference_science(
    synthetic_l2, lake_scene, tmp_path, monkeypatch
):
    lakes, scene, ch = lake_scene
    sen3 = tmp_path / (scene["name"] + ".SEN3")
    sen3.mkdir()
    (sen3 / "xfdumanifest.xml").write_text("<manifest/>")
    for n in range(1, 22):
        (sen3 / f"Oa{n:02}_radiance.nc").write_bytes(b"placeholder")
    monkeypatch.setattr(settings, "L2GEN_BIN", Path("/synthetic/l2gen"))
    monkeypatch.setattr(runtime, "command_version", lambda path: "synthetic test executable")
    captured = {}

    def run(command, log_path, **kw):
        par = Path(command[-1].split("=", 1)[1])
        values = dict(
            line.split("=", 1)
            for line in par.read_text().splitlines()
            if "=" in line and not line.startswith("#")
        )
        captured.update(values)
        shutil.copy2(synthetic_l2, values["ofile"])

    monkeypatch.setattr(l2gen._m_utils, "run_logged", run)
    output, info = l2gen.run_l2gen(sen3, process_full_scene=True)
    assert output.exists()
    assert (
        captured["gas_opt"] == "15" and captured["aer_opt"] == "-99" and captured["brdf_opt"] == "0"
    )
    assert captured["proc_ocean"] == "2" and captured["proc_land"] == "1"
    for key in (
        "maskland",
        "maskcloud",
        "maskglint",
        "masksunzen",
        "masksatzen",
        "maskhilt",
        "maskstlight",
    ):
        assert captured[key] == "0"
    assert info["internal_output_masks_disabled"] is True
    assert "rhos_884" in captured["l2prod"]


def test_pipeline_completes_and_reuses_valid_scene(synthetic_l2, lake_scene, monkeypatch):
    lakes, scene, ch = lake_scene
    registry.registry_upsert_scene(scene, ch, len(lakes))
    calls = []

    def stage(*args):
        calls.append("stage")
        target = settings.STAGING_DIR / "managed.SEN3"
        target.mkdir()
        return None, target, "synthetic"

    def l2(*args, **kwargs):
        calls.append("l2gen")
        target = settings.L2_DIR / "dense.nc"
        shutil.copy2(synthetic_l2, target)
        return target, {"synthetic": True}

    monkeypatch.setattr(staging, "stage_scene_efficient", stage)
    monkeypatch.setattr(l2gen, "run_l2gen", l2)
    monkeypatch.setattr(windows, "plan_scene_windows", lambda *a, **k: [None])
    stats = pipeline.process_scene_once(scene, ch, lakes, "", "")
    assert len(stats) == 2
    assert registry.registry_table(ch).iloc[0].status == "done"
    assert not (settings.STAGING_DIR / "managed.SEN3").exists()
    assert not (settings.L2_DIR / "dense.nc").exists()
    pipeline.process_scene_once(scene, ch, lakes, "", "")
    assert calls == ["stage", "l2gen"]


def test_cli_existing_l2_to_share_zip(configured, synthetic_l2, lake_scene, tmp_path):
    _, scene, _ = lake_scene
    scene_json = tmp_path / "scene.json"
    scene_json.write_text(json.dumps(scene))
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(__file__).resolve().parents[1] / "src")

    def command(*args):
        result = subprocess.run(
            [sys.executable, "-m", "cyanolake", "--config", str(configured), *args],
            env=env,
            capture_output=True,
            text=True,
            timeout=60,
        )
        assert result.returncode == 0, result.stdout + result.stderr
        return result.stdout

    output = command("ingest-l2", "--input", str(synthetic_l2), "--scene-json", str(scene_json))
    assert "observations" in output
    assert "done" in command("status")
    output = command("export", "--scene-id", scene["id"], "--resolution", "300")
    assert "SHARE_MINIMAL.zip" in output


def test_dashboard_constructs_without_launch(configured):
    from cyanolake.dashboard import build_dashboard

    app = build_dashboard()
    assert len(app.config["components"]) > 100
    assert len(app.config["dependencies"]) >= 20
    app.close()
