"""Standalone runtime implementation."""

from __future__ import annotations

import os
import shlex
import shutil
import sqlite3 as _sqlite3
import subprocess
from pathlib import Path

from . import settings as _m_settings


# Reference cell 11, lines 38-54.
def read_bash_environment(script: Path) -> dict[str, str]:
    """Return a sourced Bash environment without mutating this Python process."""
    environment = os.environ.copy()
    environment["OCSSWROOT"] = str(_m_settings.OCSSWROOT)
    script = Path(script)
    if not script.is_file():
        return environment
    command = f"source {shlex.quote(str(script))} >/dev/null 2>&1 && env -0"
    result = subprocess.run(
        ["bash", "-lc", command], capture_output=True, check=True, env=environment
    )
    for item in result.stdout.split(b"\x00"):
        if b"=" not in item:
            continue
        key, value = item.split(b"=", 1)
        environment[key.decode(errors="ignore")] = value.decode(errors="ignore")
    return environment


# Reference cell 11, lines 86-93.
def _proj_database_layout(database: Path) -> tuple[int, int]:
    with _sqlite3.connect(database) as connection:
        rows = dict(
            connection.execute(
                "SELECT key, value FROM metadata WHERE key IN ('DATABASE.LAYOUT.VERSION.MAJOR', 'DATABASE.LAYOUT.VERSION.MINOR')"
            ).fetchall()
        )
    return (int(rows["DATABASE.LAYOUT.VERSION.MAJOR"]), int(rows["DATABASE.LAYOUT.VERSION.MINOR"]))


# Reference cell 11, lines 116-125.
def find_ocssw_tool(name: str) -> Path | None:
    candidates = [
        _m_settings.OCSSWROOT / "bin" / name,
        _m_settings.OCSSWROOT / "scripts" / name,
        _m_settings.OCSSWROOT / "scripts" / f"{name}.py",
    ]
    found = shutil.which(name) or shutil.which(f"{name}.py")
    if found:
        candidates.append(Path(found))
    return next((path for path in candidates if path.is_file()), None)


# Reference cell 11, lines 132-143.
def command_version(path: Path | None) -> str:
    if path is None:
        return "not found"
    for option in ("--version", "version"):
        result = subprocess.run(
            [str(path), option],
            capture_output=True,
            text=True,
            env=_m_settings.OCSSW_SUBPROCESS_ENV,
            timeout=30,
        )
        text = (result.stdout + "\n" + result.stderr).strip()
        if text:
            return text.splitlines()[0]
    return str(path)


# Reference cell 12, lines 67-71.
def require_ocssw(require_getanc: bool = False) -> None:
    if _m_settings.L2GEN_BIN is None:
        raise RuntimeError("l2gen was not found. Complete the local OCSSW setup first.")
    if require_getanc and _m_settings.GETANC_BIN is None:
        raise RuntimeError("getanc was not found in the active OCSSW installation.")


def activate_ocssw() -> dict:
    """Keep OCSSW libraries in subprocesses, separate from Python's GIS stack."""
    environment = read_bash_environment(_m_settings.OCSSWROOT / "OCSSW_bash.env")
    environment["OCSSWROOT"] = str(_m_settings.OCSSWROOT)
    _m_settings.OCSSW_SUBPROCESS_ENV = environment
    for name in (
        "PROJ_LIB",
        "PROJ_DATA",
        "GDAL_DATA",
        "LD_LIBRARY_PATH",
        "PYTHONPATH",
        "GDAL_DRIVER_PATH",
    ):
        value = os.environ.get(name, "")
        if str(_m_settings.OCSSWROOT) in value:
            retained = [
                part for part in value.split(os.pathsep) if str(_m_settings.OCSSWROOT) not in part
            ]
            if retained:
                os.environ[name] = os.pathsep.join(retained)
            else:
                os.environ.pop(name, None)
    # Unlike Colab wheel-only paths, this also works with Conda's PROJ data.
    import pyproj
    import rasterio

    pyproj.CRS.from_epsg(4326)
    rasterio.crs.CRS.from_epsg(4326)
    _m_settings.PYTHON_PROJ_DATA_DIR = Path(pyproj.datadir.get_data_dir())
    _m_settings.PYTHON_PROJ_DATABASE = _m_settings.PYTHON_PROJ_DATA_DIR / "proj.db"
    _m_settings.L2GEN_BIN = find_ocssw_tool("l2gen")
    _m_settings.GETANC_BIN = find_ocssw_tool("getanc")
    return {
        "ocssw_root": str(_m_settings.OCSSWROOT),
        "l2gen": str(_m_settings.L2GEN_BIN) if _m_settings.L2GEN_BIN else None,
        "getanc": str(_m_settings.GETANC_BIN) if _m_settings.GETANC_BIN else None,
        "python_proj_data": str(_m_settings.PYTHON_PROJ_DATA_DIR),
    }


def preflight(*, require_lakes: bool = True) -> dict:
    """Report actionable setup failures without downloading a scene."""
    import netrc
    import platform

    from . import configuration

    report = activate_ocssw()
    errors = []
    if not (_m_settings.OCSSWROOT / "OCSSW_bash.env").is_file():
        errors.append("OCSSW_bash.env is missing from the configured OCSSW root")
    if not _m_settings.L2GEN_BIN:
        errors.append("l2gen is missing: install OCSSW or correct paths.ocssw_root")
    elif not os.access(_m_settings.L2GEN_BIN, os.X_OK):
        errors.append("The selected l2gen file is not executable")
    if _m_settings.GET_ANCILLARY:
        if not _m_settings.GETANC_BIN:
            errors.append("getanc is missing")
        auth_path = Path.home() / ".netrc"
        try:
            if auth_path.stat().st_mode & 0o077:
                errors.append("Earthdata .netrc permissions must be 0600")
            credentials = netrc.netrc(str(auth_path)).authenticators("urs.earthdata.nasa.gov")
            if not credentials or not credentials[0] or not credentials[2]:
                errors.append("Earthdata login is absent from .netrc")
        except (OSError, netrc.NetrcParseError):
            errors.append("Create a valid Earthdata .netrc using configure-earthdata")
    for sensor in ("s3a", "s3b"):
        if not (_m_settings.OCSSWROOT / "share" / "olci" / sensor).is_dir():
            errors.append(f"OCSSW OLCI {sensor} sensor data is missing")
    if require_lakes and not _m_settings.TARGET_LAKES_GPKG.is_file():
        errors.append("Target HydroLAKES file is missing: run prepare-lakes")
    root = _m_settings.SCRATCH_ROOT
    existing = next((p for p in (root, *root.parents) if p.exists()), Path("/"))
    report.update(
        ok=not errors,
        errors=errors,
        platform=platform.platform(),
        scratch_free_gib=round(shutil.disk_usage(existing).free / 1024**3, 2),
        configuration=configuration.CURRENT,
        validation_scope="Local installation and credentials presence; no remote login or science validation",
    )
    return report
