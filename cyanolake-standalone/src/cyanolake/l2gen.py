"""Standalone l2gen implementation."""

from __future__ import annotations

import datetime as dt
import hashlib
import math
import os
import re
from pathlib import Path
from typing import Any, Iterable

import netCDF4
import numpy as np

from . import catalogue as _m_catalogue
from . import runtime as _m_runtime
from . import science as _m_science
from . import settings as _m_settings
from . import staging as _m_staging
from . import utils as _m_utils


# Reference cell 48, lines 69-82.
def l2_products(mode: str, full_spectrum: bool = False) -> tuple[list[str], set[str]]:
    mode = mode.lower()
    if mode in {"full_swath_rhos", "rayleigh_ci"}:
        mode = "full_swath_rhos"
    if mode not in {"full_swath_rhos", "standard_dual"}:
        raise ValueError("mode must be 'full_swath_rhos' or 'standard_dual'")
    waves = (
        _m_settings.SPECTRAL_RHOS_WAVELENGTHS
        if full_spectrum
        else _m_settings.GLOBAL_ANALYSIS_RHOS_WAVELENGTHS
    )
    products = [_l2gen_rhos_product_name(wave) for wave in waves]
    required = {
        _l2gen_rhos_product_name(wave) for wave in _m_settings.GLOBAL_ANALYSIS_RHOS_WAVELENGTHS
    }
    if mode == "standard_dual":
        rrs_waves = _m_settings.AQUATIC_RRS_WAVELENGTHS if full_spectrum else (665, 709)
        products.extend((f"Rrs_{wave}" for wave in rrs_waves))
        required.update({"Rrs_665", "Rrs_709"})
    return (list(dict.fromkeys(products)), required)


# Reference cell 17, lines 27-38.
def write_par_file(path: Path, parameters: dict[str, Any]) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for key, value in parameters.items():
        if value is None:
            continue
        if isinstance(value, bool):
            value = int(value)
        lines.append(f"{key}={value}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# Reference cell 17, lines 41-42.
def product_basename(sen3_path: Path) -> str:
    return Path(sen3_path).name.removesuffix(".SEN3")


# Reference cell 17, lines 45-64.
def olci_acquisition_interval(sen3_path: Path) -> tuple[str, str]:
    """Return ISO start/stop times from a standard Sentinel-3 product name."""
    name = product_basename(Path(sen3_path))
    stamps = re.findall("\\d{8}T\\d{6}", name)
    if len(stamps) < 2:
        raise ValueError(
            f"Could not parse the acquisition start/stop times from the OLCI product name: {name}"
        )
    start_dt = dt.datetime.strptime(stamps[0], "%Y%m%dT%H%M%S")
    stop_dt = dt.datetime.strptime(stamps[1], "%Y%m%dT%H%M%S")
    duration = (stop_dt - start_dt).total_seconds()
    if duration <= 0 or duration > 2 * 3600:
        raise ValueError(f"Implausible OLCI acquisition interval ({duration:g} seconds) in {name}")
    return (start_dt.strftime("%Y-%m-%dT%H:%M:%S"), stop_dt.strftime("%Y-%m-%dT%H:%M:%S"))


# Reference cell 17, lines 67-100.
def validate_ancillary_par(path: Path) -> dict[str, Any]:
    """Require a non-empty getanc par file with usable MET/OZONE inputs."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError("getanc returned no usable ancillary parameter file")
    parameters: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        parameters[key.strip()] = value.strip()
    required = ("met1", "met2", "met3", "ozone1", "ozone2", "ozone3")
    missing_keys = [key for key in required if not parameters.get(key)]
    if missing_keys:
        raise RuntimeError(
            "ancillary parameter file lacks required entries: " + ", ".join(missing_keys)
        )
    missing_files = []
    for key, value in parameters.items():
        expanded = Path(os.path.expandvars(value)).expanduser()
        if not expanded.is_file():
            missing_files.append(f"{key}={expanded}")
    if missing_files:
        raise RuntimeError(
            "ancillary parameter file references missing downloads: "
            + "; ".join(missing_files[:12])
        )
    return {
        "path": str(path),
        "size_bytes": path.stat().st_size,
        "parameter_keys": sorted(parameters),
        "referenced_file_count": len(parameters),
    }


# Reference cell 17, lines 103-183.
def run_getanc(sen3_path: Path, allow_fallback: bool = False) -> Path | None:
    """Download time-specific ancillary data and return its .anc parameter file."""
    if not _m_settings.GET_ANCILLARY:
        return None
    _m_runtime.require_ocssw(require_getanc=True)
    sen3_path = _m_staging.validate_sen3(sen3_path)
    name = product_basename(sen3_path)
    anc_path = _m_settings.PAR_DIR / f"{name}.anc"
    if anc_path.is_file() and anc_path.stat().st_size:
        try:
            validate_ancillary_par(anc_path)
            return anc_path
        except RuntimeError:
            anc_path.unlink(missing_ok=True)
    start_time, stop_time = olci_acquisition_interval(sen3_path)
    ancdb_path = _m_settings.PAR_DIR / "ancillary_data.db"
    attempts = [
        (
            "mission-neutral explicit-time query",
            [
                str(_m_settings.GETANC_BIN),
                "--start",
                start_time,
                "--stop",
                stop_time,
                "--ofile",
                str(anc_path),
                "--ancdb",
                str(ancdb_path),
                "--noprint",
                "--verbose",
            ],
        ),
        (
            "standard OLCI manifest query",
            [
                str(_m_settings.GETANC_BIN),
                str(sen3_path / "xfdumanifest.xml"),
                "--ofile",
                str(anc_path),
                "--ancdb",
                str(ancdb_path),
                "--noprint",
                "--verbose",
                "--refreshDB",
            ],
        ),
    ]
    errors: list[str] = []
    for number, (label, command) in enumerate(attempts, start=1):
        log_path = _m_settings.LOG_DIR / f"{name}_getanc_attempt{number}.log"
        try:
            result = _m_utils.run_logged(command, log_path, timeout=1800, cwd=_m_settings.PAR_DIR)
            try:
                validation = validate_ancillary_par(anc_path)
            except RuntimeError as validation_error:
                output_tail = (
                    ((result.stdout or "") + "\n" + (result.stderr or "")).strip()
                    or "getanc produced no console output"
                )[-3000:]
                raise RuntimeError(
                    f"{label}: getanc exited 0 but its result was unusable: {validation_error}. Log: {log_path}\n{output_tail}"
                ) from validation_error
            _m_utils.write_json(
                _m_settings.EXPORT_DIR / name / "ancillary_summary.json",
                {
                    "created_utc": _m_utils.utc_now(),
                    "query_method": label,
                    "acquisition_start": start_time,
                    "acquisition_stop": stop_time,
                    **validation,
                },
            )
            print(f"Ancillary data ready via {label}: {anc_path}")
            return anc_path
        except Exception as error:
            errors.append(f"Attempt {number} ({label})\n{error}")
            anc_path.unlink(missing_ok=True)
    if allow_fallback:
        print(
            "WARNING: getanc failed after both OLCI strategies; L2Gen will use its climatological ancillary defaults because fallback was explicitly enabled."
        )
        for error in errors:
            print("\n" + error[-2500:])
        return None
    raise RuntimeError(
        f"No usable time-specific ancillary file was produced, and climatology fallback is disabled. The detailed attempt logs are under {_m_settings.LOG_DIR}. If you deliberately accept climatological defaults for a test run, set ALLOW_CLIMATOLOGY_FALLBACK=True and rerun the configuration and dashboard cells.\n\n"
        + "\n\n".join(errors)
    )


# Reference cell 29, lines 52-78.
def validate_l2_file(path: Path, required_products: Iterable[str]) -> dict[str, Any]:
    info = _BASE_VALIDATE_L2_FILE_V5(path, required_products)
    with netCDF4.Dataset(path) as root:
        geo = root.groups["geophysical_data"]
        nav = root.groups["navigation_data"]
        rhos_name = _m_science.nearest_product(geo.variables, "rhos", 681)
        shape = nav.variables["latitude"].shape
        stride = max(1, int(math.sqrt(max(1, np.prod(shape) / 350000))))
        rhos = _m_science.to_float(geo.variables[rhos_name][::stride, ::stride])
        flags_var = _m_science.locate_l2_flags(root)
        flags = np.ma.asarray(flags_var[::stride, ::stride]).filled(0).astype(np.uint32)
        flag_map = _m_science.decode_flag_metadata(flags_var)
        land = _m_science.flag_is_set(flags, flag_map, "LAND")
        finite = np.isfinite(rhos)
        info["sampled_rhos_coverage"] = {
            "product": rhos_name,
            "stride": stride,
            "finite_all_pixels": int(finite.sum()),
            "sampled_all_pixels": int(finite.size),
            "finite_fraction_all": float(finite.mean()),
            "finite_land_pixels": int((finite & land).sum()),
            "sampled_land_pixels": int(land.sum()),
            "finite_fraction_land": float((finite & land).sum() / land.sum())
            if land.any()
            else None,
            "finite_nonland_pixels": int((finite & ~land).sum()),
            "sampled_nonland_pixels": int((~land).sum()),
            "finite_fraction_nonland": float((finite & ~land).sum() / (~land).sum())
            if (~land).any()
            else None,
        }
    return info


# Reference cell 17, lines 228-344.
def run_l2gen(
    sen3_path: Path,
    mode: str = None,
    process_full_scene: bool = None,
    bbox: tuple[float, float, float, float] | None = None,
    full_spectrum: bool = None,
) -> tuple[Path, dict[str, Any]]:
    """Run one OLCI scene through L2Gen, with optional-product retry."""
    if mode is None:
        mode = _m_settings.L2_MODE
    if process_full_scene is None:
        process_full_scene = _m_settings.PROCESS_FULL_SCENE
    if full_spectrum is None:
        full_spectrum = _m_settings.SAVE_FULL_REFLECTANCE_SPECTRA
    _m_runtime.require_ocssw(require_getanc=_m_settings.GET_ANCILLARY)
    sen3_path = _m_staging.validate_sen3(sen3_path)
    mode = mode.lower()
    name = product_basename(sen3_path)
    suffix = (
        "FULLSWATH_RHOS" if mode in {"full_swath_rhos", "rayleigh_ci"} else "FULLSWATH_RHOS_RRS"
    )
    processing_bbox = None
    if process_full_scene:
        scope_tag = "FULL"
    else:
        if bbox is None:
            raise ValueError("bbox is required when process_full_scene=False")
        processing_bbox = _m_catalogue.normalize_bbox(bbox)
        bbox_text = ",".join((f"{value:.6f}" for value in processing_bbox))
        scope_tag = "BBOX_" + hashlib.sha1(bbox_text.encode("utf-8")).hexdigest()[:10]
    spectrum_tag = "SPEC" if full_spectrum else "CORE"
    output = _m_settings.L2_DIR / f"{name}_L2GEN_{suffix}_{scope_tag}_{spectrum_tag}.nc"
    products, required = l2_products(mode, full_spectrum)
    if output.exists():
        try:
            info = validate_l2_file(output, required)
            info.update(
                {
                    "mode": mode,
                    "processing_scope": "full_scene" if process_full_scene else "bbox",
                    "processing_bbox": processing_bbox,
                    "full_reflectance_spectra_requested": bool(full_spectrum),
                    "current_l2gen_version": _m_runtime.command_version(_m_settings.L2GEN_BIN),
                    "reused_existing": True,
                }
            )
            print(f"Using validated existing L2 file: {output}")
            return (output, info)
        except Exception:
            stamp = dt.datetime.now().strftime("%Y%m%dT%H%M%S")
            incomplete = output.with_name(f"{output.stem}.incompatible_{stamp}{output.suffix}")
            output.replace(incomplete)
            print(f"Moved incompatible prior output to {incomplete}")
    anc_path = run_getanc(sen3_path, allow_fallback=_m_settings.ALLOW_CLIMATOLOGY_FALLBACK)
    run_tag = f"{mode}_{scope_tag}_{spectrum_tag}".lower()
    par_path = _m_settings.PAR_DIR / f"{name}_l2gen_{run_tag}.par"
    log_path = _m_settings.LOG_DIR / f"{name}_l2gen_{run_tag}.log"
    while True:
        parameters: dict[str, Any] = {
            "ifile": sen3_path / "xfdumanifest.xml",
            "ofile": output,
            "l2prod": " ".join(products),
            "proc_ocean": 2 if mode in {"full_swath_rhos", "rayleigh_ci"} else 1,
            "proc_land": 1,
            "maskland": 0,
            "maskcloud": 0,
            "maskglint": 0,
            "masksunzen": 0,
            "masksatzen": 0,
            "maskhilt": 0,
            "maskstlight": 0,
            "proc_uncertainty": 0,
            "gas_opt": 15,
            "deflate": int(_m_settings.NETCDF_DEFLATE),
        }
        if mode in {"full_swath_rhos", "rayleigh_ci"}:
            parameters.update({"aer_opt": -99, "brdf_opt": 0})
        if not process_full_scene:
            west, south, east, north = processing_bbox
            parameters.update({"west": west, "south": south, "east": east, "north": north})
        write_par_file(par_path, parameters)
        command = [str(_m_settings.L2GEN_BIN)]
        if anc_path is not None:
            command.append(f"par={anc_path}")
        command.append(f"par={par_path}")
        try:
            _m_utils.run_logged(command, log_path, timeout=4 * 3600)
            break
        except RuntimeError as error:
            text = str(error)
            match = re.search("product\\s+([A-Za-z0-9_]+)\\s+not found", text, re.IGNORECASE)
            if not match:
                output.unlink(missing_ok=True)
                raise
            unavailable = match.group(1)
            if unavailable in required or unavailable not in products:
                output.unlink(missing_ok=True)
                raise RuntimeError(f"Required L2Gen product unavailable: {unavailable}") from error
            products.remove(unavailable)
            output.unlink(missing_ok=True)
            print(f"Optional product {unavailable} is unavailable; retrying without it.")
    info = validate_l2_file(output, required)
    info.update(
        {
            "mode": mode,
            "l2gen_version": _m_runtime.command_version(_m_settings.L2GEN_BIN),
            "processing_scope": "full_scene" if process_full_scene else "bbox",
            "processing_bbox": processing_bbox,
            "full_reflectance_spectra_requested": bool(full_spectrum),
            "reused_existing": False,
            "internal_output_masks_disabled": True,
            "rhos_scope": "complete geolocated swath where L2Gen computes finite values; LAND retained in l2_flags",
            "products_requested_final": products,
            "parameter_file": str(par_path),
            "ancillary_file": str(anc_path) if anc_path else None,
            "log": str(log_path),
        }
    )
    print(f"L2 complete: {output} ({_m_utils.human_size(output.stat().st_size)})")
    return (output, info)


# Reference cell 48, lines 66-67.
def _l2gen_rhos_product_name(wave: int) -> str:
    return f"rhos_{_m_settings.OLCI_L2GEN_RHOS_NAME_OVERRIDES.get(int(wave), int(wave))}"


# Reference cell 17, lines 186-225.
def _BASE_VALIDATE_L2_FILE_V5(path: Path, required_products: Iterable[str]) -> dict[str, Any]:
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing/empty L2 output: {path}")
    with netCDF4.Dataset(path) as root:
        if "geophysical_data" not in root.groups or "navigation_data" not in root.groups:
            raise RuntimeError("L2 output lacks geophysical_data/navigation_data groups")
        geo = root.groups["geophysical_data"]
        nav = root.groups["navigation_data"]
        flags_present = (
            "l2_flags" in geo.variables
            or "l2_flags" in nav.variables
            or "l2_flags" in root.variables
        )
        if not flags_present:
            raise RuntimeError("L2 output lacks l2_flags required for water/QA masking")
        missing = [name for name in required_products if name not in geo.variables]
        if missing:
            raise RuntimeError(f"L2 output lacks required products: {missing}")
        for coordinate in ("latitude", "longitude"):
            if coordinate not in nav.variables:
                raise RuntimeError(f"L2 output lacks navigation variable {coordinate}")
        shape = nav.variables["latitude"].shape
        if nav.variables["longitude"].shape != shape:
            raise RuntimeError("L2 latitude/longitude shapes differ")
        provenance_attributes = {}
        for attribute in (
            "product_name",
            "instrument",
            "platform",
            "processing_version",
            "software_name",
            "software_version",
            "processing_time",
            "processing_control",
            "history",
        ):
            if attribute in root.ncattrs():
                provenance_attributes[attribute] = str(root.getncattr(attribute))
        return {
            "path": str(path),
            "shape": tuple(map(int, shape)),
            "products": sorted(geo.variables),
            "size_bytes": path.stat().st_size,
            "source_global_attributes": provenance_attributes,
        }
