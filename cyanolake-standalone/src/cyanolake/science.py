"""Standalone science implementation."""

from __future__ import annotations

from typing import Any, Iterable

import netCDF4
import numpy as np

from . import settings as _m_settings


# Reference cell 19, lines 5-6.
def to_float(values: Any) -> np.ndarray:
    return np.ma.asarray(values).filled(np.nan).astype(np.float32, copy=False)


# Reference cell 19, lines 9-13.
def product_wavelength(name: str) -> float:
    match = _m_settings.PRODUCT_PATTERN.match(name)
    if not match:
        raise ValueError(f"Cannot parse wavelength from product name: {name}")
    return float(match.group("wavelength"))


# Reference cell 19, lines 16-30.
def nearest_product(
    variables: Iterable[str], prefix: str, target: float, tolerance_nm: float = 12.0
) -> str:
    candidates = []
    for name in variables:
        if not name.startswith(prefix + "_"):
            continue
        try:
            wavelength = product_wavelength(name)
        except ValueError:
            continue
        candidates.append((abs(wavelength - target), name))
    if not candidates or min(candidates)[0] > tolerance_nm:
        raise KeyError(f"No {prefix} product within {tolerance_nm} nm of {target}")
    return min(candidates)[1]


# Reference cell 19, lines 33-38.
def decode_flag_metadata(variable: netCDF4.Variable) -> dict[str, np.uint32]:
    meanings = str(getattr(variable, "flag_meanings", "")).replace(",", " ").split()
    masks = np.asarray(getattr(variable, "flag_masks", []), dtype=np.uint32).ravel()
    if len(meanings) != len(masks):
        raise RuntimeError("l2_flags flag_meanings and flag_masks metadata do not align")
    return {meaning.upper(): np.uint32(mask) for meaning, mask in zip(meanings, masks)}


# Reference cell 19, lines 41-45.
def flag_is_set(flags: np.ndarray, flag_map: dict[str, np.uint32], name: str) -> np.ndarray:
    mask = flag_map.get(name.upper())
    if mask is None:
        return np.zeros(flags.shape, dtype=bool)
    return flags.astype(np.uint32, copy=False) & mask != 0


# Reference cell 19, lines 48-55.
def locate_l2_flags(root: netCDF4.Dataset) -> netCDF4.Variable:
    for group_name in ("geophysical_data", "navigation_data"):
        group = root.groups.get(group_name)
        if group is not None and "l2_flags" in group.variables:
            return group.variables["l2_flags"]
    if "l2_flags" in root.variables:
        return root.variables["l2_flags"]
    raise KeyError("l2_flags is absent from the L2 file")


# Reference cell 30, lines 1-21.
def ci_family(
    r620: np.ndarray,
    r665: np.ndarray,
    r681: np.ndarray,
    r709: np.ndarray,
    wavelengths: tuple[float, float, float, float] = (620.0, 665.0, 681.0, 709.0),
    detection_limit: float = _m_settings.CI_DETECTION_LIMIT,
    ss681_max: float = 0.0,
    ss665_min: float = 0.0,
) -> dict[str, np.ndarray]:
    """Two-stage CI family with explicit, tunable spectral gates."""
    r620, r665, r681, r709 = [
        np.asarray(array, dtype=np.float32) for array in (r620, r665, r681, r709)
    ]
    w620, w665, w681, w709 = map(float, wavelengths)
    finite = np.isfinite(r620) & np.isfinite(r665) & np.isfinite(r681) & np.isfinite(r709)
    ss681 = r681 - r665 - (r709 - r665) * ((w681 - w665) / (w709 - w665))
    ss665 = r665 - r620 - (r681 - r620) * ((w665 - w620) / (w681 - w620))
    ss681 = np.where(finite, ss681, np.nan).astype(np.float32)
    ss665 = np.where(finite, ss665, np.nan).astype(np.float32)
    candidate = finite & (ss681 < float(ss681_max))
    ci = np.where(candidate, -ss681, np.nan).astype(np.float32)
    detection = candidate & (ss665 > float(ss665_min)) & (ci > float(detection_limit))
    return {
        "SS_681": ss681,
        "SS_665": ss665,
        "CI": ci,
        "CI_cyano": np.where(detection, ci, np.nan).astype(np.float32),
        "finite": finite,
        "ci_candidate": candidate,
        "ci_detection": detection,
    }


# Reference cell 20, lines 39-46.
def ndci_index(red: np.ndarray, red_edge: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    red = np.asarray(red, dtype=np.float32)
    red_edge = np.asarray(red_edge, dtype=np.float32)
    denominator = red_edge + red
    valid = np.isfinite(red) & np.isfinite(red_edge) & (np.abs(denominator) > 1e-08)
    result = np.full(red.shape, np.nan, dtype=np.float32)
    np.divide(red_edge - red, denominator, out=result, where=valid)
    return (result, valid)


# Reference cell 48, lines 11-16.
def _spectral_line_height(center, left, right, wl_center, wl_left, wl_right):
    center = np.asarray(center, dtype=float)
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    baseline = left + (right - left) * (
        (float(wl_center) - float(wl_left)) / (float(wl_right) - float(wl_left))
    )
    return center - baseline


# Reference cell 48, lines 19-56.
def additional_olci_indices(
    rhos: dict[int, np.ndarray], wavelengths: dict[int, float] | None = None
):
    """Compute production OLCI bloom indices from Rayleigh-corrected rhos.

    MPH: maximum peak height among ~681/~709/~754 nm above a common 665--~885 nm baseline.
    FAI: OLCI-adapted floating-algae index (AFAI form), 754 nm above the 665--865 nm baseline.
         This is not the original Hu SWIR-based FAI because OLCI has no suitable SWIR band.
    """
    need = (665, 681, 709, 754, 865, 885)
    missing = [w for w in need if w not in rhos]
    if missing:
        raise KeyError(f"Missing rhos bands required for MPH/FAI: {missing}")
    wl = {w: float((wavelengths or {}).get(w, w)) for w in need}
    finite = {w: np.isfinite(rhos[w]) for w in need}
    fai = _spectral_line_height(rhos[754], rhos[665], rhos[865], wl[754], wl[665], wl[865])
    fai_valid = finite[754] & finite[665] & finite[865]
    mph_candidates = []
    mph_valid_components = []
    for peak in (681, 709, 754):
        h = _spectral_line_height(rhos[peak], rhos[665], rhos[885], wl[peak], wl[665], wl[885])
        v = finite[peak] & finite[665] & finite[885]
        mph_candidates.append(np.where(v, h, np.nan))
        mph_valid_components.append(v)
    stack = np.stack(mph_candidates, axis=0)
    with np.errstate(all="ignore"):
        mph = np.nanmax(stack, axis=0)
    mph_valid = np.logical_or.reduce(mph_valid_components) & np.isfinite(mph)
    peak_waves = np.asarray([wl[681], wl[709], wl[754]], dtype=np.float32)
    safe = np.where(np.isfinite(stack), stack, -np.inf)
    peak_idx = np.argmax(safe, axis=0)
    mph_peak_nm = peak_waves[peak_idx]
    mph_peak_nm[~mph_valid] = np.nan
    return {
        "MPH": mph,
        "MPH_valid": mph_valid,
        "MPH_peak_nm": mph_peak_nm,
        "FAI": fai,
        "FAI_valid": fai_valid,
    }


# Reference cell 50, lines 1011-1013.
def _unit_sphere_xyz(lat, lon):
    lat = np.deg2rad(np.asarray(lat, float))
    lon = np.deg2rad(np.asarray(lon, float))
    c = np.cos(lat)
    return np.column_stack((c * np.cos(lon), c * np.sin(lon), np.sin(lat)))
