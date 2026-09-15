"""Standalone legacy maps implementation."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import folium
import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.patheffects as path_effects
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import rasterio
from branca.colormap import LinearColormap
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from folium import plugins as folium_plugins
from matplotlib.patches import Patch
from pyproj import CRS, Transformer
from pyresample import geometry as pr_geometry
from pyresample import kd_tree as pr_kd_tree
from rasterio.transform import from_bounds
from shapely import contains_xy
from tqdm.auto import tqdm

from . import catalogue as _m_catalogue
from . import hydrolakes as _m_hydrolakes
from . import science as _m_science
from . import settings as _m_settings
from . import utils as _m_utils

LAYER_STYLE = {
    "CI": ("magma", (0.0, 0.003), "CI (dimensionless)"),
    "CI_cyano": ("viridis", (0.0001, 0.003), "CI$_{cyano}$ detection (>0.0001)"),
    "NDCI": ("RdYlGn", (-0.25, 0.5), "NDCI (dimensionless)"),
    "SS_681": ("RdBu_r", (-0.003, 0.003), "SS(681)"),
    "SS_665": ("RdBu_r", (-0.003, 0.003), "SS(665)"),
    "negative_core_rhos_count": ("plasma", (0, 4), "Negative core rhos band count"),
}

MASK_LAYER_STYLE = {
    "water_mask": ("Blues", "Water mask (0/1)"),
    "valid_water_mask": ("Greens", "Retained water after QA (0/1)"),
    "qa_excluded_mask": ("Reds", "QA/buffer-excluded water (0/1)"),
    "land_adjacency_mask": ("Oranges", "Land-adjacency mask (0/1)"),
    "ci_valid_mask": ("Greens", "Finite CI inputs after QA (0/1)"),
    "ci_candidate_mask": ("Purples", "SS(681) CI candidate mask (0/1)"),
    "ci_detection_mask": ("viridis", "CIcyano detection mask (0/1)"),
    "ndci_valid_mask": ("Greens", "Valid NDCI mask (0/1)"),
}

LAYER_STYLE.update(
    {
        "CI": ("magma", (0.0, 0.003), "CI (dimensionless)"),
        "CI_cyano": ("cividis", (0.0001, 0.003), "CI$_{cyano}$ detection (>0.0001)"),
        "NDCI": ("BrBG", (-0.25, 0.5), "NDCI (dimensionless)"),
    }
)

LAYER_STYLE.update(
    {
        "CI_pre_mask": ("magma", (0.0, 0.003), "CI before water/QA mask"),
        "CI_cyano_pre_mask": (
            "cividis",
            (0.0001, 0.003),
            "CI$_{cyano}$ detection (>0.0001) before water/QA mask",
        ),
        "NDCI_pre_mask": ("BrBG", (-0.25, 0.5), "NDCI before water/QA mask"),
        "rhos_620": ("viridis", (-0.02, 0.12), "$\\rho_s$(620)"),
        "rhos_665": ("viridis", (-0.02, 0.12), "$\\rho_s$(665)"),
        "rhos_681": ("viridis", (-0.02, 0.12), "$\\rho_s$(681)"),
        "rhos_709": ("viridis", (-0.02, 0.12), "$\\rho_s$(709)"),
        "rhos_865": ("viridis", (-0.02, 0.12), "$\\rho_s$(865)"),
        "CI_exclusion_reason": ("turbo", (0, 7), "CIcyano outcome reason (code)"),
        "NDCI_exclusion_reason": ("turbo", (0, 6), "NDCI outcome reason (code)"),
    }
)


# Reference cell 25, lines 1-20.
def bbox_slices(
    latitude: np.ndarray,
    longitude: np.ndarray,
    bbox: tuple[float, float, float, float] | None,
    margin: int = 4,
) -> tuple[slice, slice]:
    if bbox is None:
        return (slice(None), slice(None))
    west, south, east, north = _m_catalogue.normalize_bbox(bbox)
    inside = (
        np.isfinite(latitude)
        & np.isfinite(longitude)
        & (longitude >= west)
        & (longitude <= east)
        & (latitude >= south)
        & (latitude <= north)
    )
    rows, columns = np.where(inside)
    if not rows.size:
        raise ValueError("No geolocated pixels intersect the requested bbox")
    return (
        slice(max(0, rows.min() - margin), min(latitude.shape[0], rows.max() + margin + 1)),
        slice(max(0, columns.min() - margin), min(latitude.shape[1], columns.max() + margin + 1)),
    )


# Reference cell 25, lines 23-24.
def preview_stride(shape: tuple[int, int], max_points: int = 350000) -> int:
    return max(1, int(math.ceil(math.sqrt(np.prod(shape) / max_points))))


# Reference cell 25, lines 27-34.
def robust_limits(values: np.ndarray, default: tuple[float, float]) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if not finite.size:
        return default
    low, high = np.nanpercentile(finite, [2, 98])
    if not np.isfinite(low) or not np.isfinite(high) or low == high:
        return default
    return (float(low), float(high))


# Reference cell 25, lines 37-65.
def load_true_color(l2_path: Path, row_slice: slice, col_slice: slice, stride: int):
    with netCDF4.Dataset(l2_path) as root:
        geo = root.groups["geophysical_data"]
        nav = root.groups["navigation_data"]
        try:
            names = [
                _m_science.nearest_product(geo.variables, "rhot", target, 15)
                for target in (665, 560, 443)
            ]
            title = "L2Gen input TOA true color"
        except KeyError:
            names = [
                _m_science.nearest_product(geo.variables, "rhos", target, 15)
                for target in (665, 560, 443)
            ]
            title = "L2Gen Rayleigh-corrected pseudo true color"
        rs = slice(row_slice.start, row_slice.stop, stride)
        cs = slice(col_slice.start, col_slice.stop, stride)
        latitude = _m_science.to_float(nav.variables["latitude"][rs, cs])
        longitude = _m_science.to_float(nav.variables["longitude"][rs, cs])
        rgb = np.dstack([_m_science.to_float(geo.variables[name][rs, cs]) for name in names])
    finite = np.all(np.isfinite(rgb), axis=2)
    stretched = np.zeros_like(rgb, dtype=np.float32)
    for channel in range(3):
        values = rgb[..., channel]
        sample = values[finite]
        low, high = (0.0, 0.12) if not sample.size else np.nanpercentile(sample, [1, 99])
        if not np.isfinite(high) or high <= low:
            low, high = (0.0, 0.12)
        stretched[..., channel] = np.clip((values - low) / (high - low), 0, 1)
    stretched = np.power(stretched, 0.8)
    return (latitude, longitude, stretched, finite, title)


# Reference cell 31, lines 519-598.
def render_layer(
    ax,
    derived_path,
    l2_path,
    layer,
    bbox=None,
    max_points=350000,
    vmin=None,
    vmax=None,
    resolution_m=_m_settings.DEFAULT_MAP_RESOLUTION_M,
    lakes=None,
    show_hydrolakes=True,
    show_context=False,
):
    """Geographic publication raster without a second Cartopy reprojection."""
    grid = regular_geographic_raster(
        derived_path, l2_path, layer, bbox, resolution_m, max_cells=_m_settings.MAP_MAX_CELLS
    )
    data = grid["data"]
    west, south, east, north = grid["extent"]
    kwargs = dict(
        extent=(west, east, south, north),
        origin="upper",
        transform=ccrs.PlateCarree(),
        interpolation="nearest",
        rasterized=True,
    )
    title = layer.replace("_", " ")
    if layer == "true_color":
        artist = ax.imshow(np.nan_to_num(data, nan=0.92), **kwargs)
        title = (
            "Sentinel-3 OLCI top-of-atmosphere true colour"
            if grid["source_products"]["source"] == "rhot"
            else "Sentinel-3 OLCI Rayleigh-corrected pseudo true colour"
        )
    else:
        if layer in MASK_LAYER_STYLE or layer in {
            "l2_water_mask",
            "hydrolakes_water_mask",
            "bright_pixel_mask",
        }:
            cmap, label = MASK_LAYER_STYLE.get(layer, ("viridis", layer.replace("_", " ")))
            lower, upper = (0 if vmin is None else float(vmin), 1 if vmax is None else float(vmax))
        else:
            cmap, defaults, label = LAYER_STYLE[layer]
            lower, upper = robust_limits(data, defaults)
            lower = lower if vmin is None else float(vmin)
            upper = upper if vmax is None else float(vmax)
        cmap_obj = plt.get_cmap(cmap).copy()
        cmap_obj.set_bad((0.96, 0.96, 0.96, 0))
        if layer in {"CI", "CI_cyano", "NDCI"}:

            def companion(name):
                lat, lon, values, _ = _native_scene_layer(derived_path, l2_path, name, bbox)
                return _resample_nearest(
                    lon, lat, values, grid["area"], grid["effective_resolution_m"]
                )

            water = companion("water_mask") >= 0.5
            retained = companion("valid_water_mask") >= 0.5
            finite = np.isfinite(data)
            rgba = np.zeros(data.shape + (4,), dtype=np.float32)
            rgba[water & ~retained] = mcolors.to_rgba("#252525", 0.92)
            rgba[retained & ~finite] = mcolors.to_rgba("#c7c7c7", 0.92)
            rgba[finite] = cmap_obj(mcolors.Normalize(lower, upper, clip=True)(data[finite]))
            artist = ax.imshow(rgba, **kwargs)
            scalar = plt.cm.ScalarMappable(norm=mcolors.Normalize(lower, upper), cmap=cmap_obj)
            scalar.set_array([])
            plt.colorbar(scalar, ax=ax, fraction=0.043, pad=0.04, label=label)
            ax.legend(
                handles=[
                    Patch(facecolor="#c7c7c7", label="Valid water: algorithm non-detection"),
                    Patch(facecolor="#252525", label="Excluded by selected post-mask"),
                ],
                loc="lower left",
                fontsize=7,
                framealpha=0.88,
            )
            title = f"{title} — detections coloured; non-detections grey"
        else:
            artist = ax.imshow(data, cmap=cmap_obj, vmin=lower, vmax=upper, **kwargs)
            plt.colorbar(artist, ax=ax, fraction=0.043, pad=0.04, label=label)
    ax.set_extent((west, east, south, north), crs=ccrs.PlateCarree())
    if show_context:
        try:
            ax.add_feature(
                cfeature.LAKES.with_scale("10m"),
                facecolor="none",
                edgecolor="#607d8b",
                linewidth=0.45,
            )
            ax.add_feature(
                cfeature.COASTLINE.with_scale("10m"), edgecolor="#37474f", linewidth=0.55
            )
            ax.add_feature(cfeature.BORDERS.with_scale("10m"), edgecolor="#546e7a", linewidth=0.45)
            ax.add_feature(cfeature.STATES.with_scale("10m"), edgecolor="#78909c", linewidth=0.35)
        except Exception as error:
            print(f"Context boundaries unavailable: {error}")
    if show_hydrolakes and lakes is not None and (not lakes.empty):
        for geometry in lakes.to_crs(4326).geometry:
            ax.add_geometries(
                [geometry],
                ccrs.PlateCarree(),
                facecolor="none",
                edgecolor="#00bcd4",
                linewidth=0.7,
                zorder=6,
            )
    gridlines = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.45,
        color="#455a64",
        alpha=0.5,
        linestyle="--",
    )
    gridlines.top_labels = False
    gridlines.right_labels = False
    gridlines.xformatter = LONGITUDE_FORMATTER
    gridlines.yformatter = LATITUDE_FORMATTER
    gridlines.xlabel_style = {"size": 8}
    gridlines.ylabel_style = {"size": 8}
    ax.set_title(title, weight="bold", pad=10)
    ax.annotate(
        "N",
        xy=(0.94, 0.93),
        xytext=(0.94, 0.81),
        xycoords="axes fraction",
        ha="center",
        weight="bold",
        arrowprops=dict(facecolor="white", edgecolor="black", width=2),
    )
    ax.text(
        0.995,
        0.006,
        f"nearest-neighbour • {grid['width']}×{grid['height']} cells • {grid['effective_resolution_m']:.0f} m display grid",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color="white",
        path_effects=[path_effects.withStroke(linewidth=2, foreground="black")],
    )
    ax._olci_grid = grid
    return artist


# Reference cell 31, lines 601-609.
def plot_layer(
    derived_path,
    l2_path,
    layer="CI_cyano",
    bbox=None,
    vmin=None,
    vmax=None,
    figsize=(10, 8),
    resolution_m=_m_settings.DEFAULT_MAP_RESOLUTION_M,
    lakes=None,
    show_context=False,
):
    figure = plt.figure(figsize=figsize, constrained_layout=False)
    axis = figure.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    render_layer(
        axis,
        derived_path,
        l2_path,
        layer,
        bbox,
        vmin=vmin,
        vmax=vmax,
        resolution_m=resolution_m,
        lakes=lakes,
        show_context=show_context,
    )
    figure.subplots_adjust(left=0.1, right=0.9, bottom=0.09, top=0.91)
    return figure


# Reference cell 31, lines 333-380.
def export_map_bundle(
    derived_path: Path,
    l2_path: Path,
    bbox: tuple[float, float, float, float] | None = None,
    dpi: int = 400,
    resolution_m: float = _m_settings.DEFAULT_MAP_RESOLUTION_M,
    lakes: gpd.GeoDataFrame | None = None,
    export_cogs: bool = True,
) -> dict[str, Path]:
    scene_dir = _m_settings.EXPORT_DIR / Path(l2_path).stem / "publication_maps"
    scene_dir.mkdir(parents=True, exist_ok=True)
    outputs: dict[str, Path] = {}
    layers = ("true_color", "CI", "CI_cyano", "NDCI")
    grids = {}
    for layer in layers:
        figure = plot_layer(
            derived_path, l2_path, layer, bbox, resolution_m=resolution_m, lakes=lakes
        )
        png = scene_dir / f"{Path(l2_path).stem}_{layer}_{int(resolution_m)}m.png"
        pdf = png.with_suffix(".pdf")
        figure.savefig(png, dpi=int(dpi), facecolor="white")
        figure.savefig(pdf, dpi=int(dpi), facecolor="white")
        plt.close(figure)
        outputs[f"{layer}_png"] = png
        outputs[f"{layer}_pdf"] = pdf
        if export_cogs:
            grid = regular_scene_raster(derived_path, l2_path, layer, bbox, resolution_m)
            grids[layer] = grid
            tif = write_regular_geotiff(grid, png.with_suffix(".tif"))
            outputs[f"{layer}_geotiff"] = tif
    combined, axes = plt.subplots(
        2,
        2,
        figsize=(15.5, 12.5),
        constrained_layout=False,
        subplot_kw={"projection": ccrs.PlateCarree()},
    )
    for axis, layer in zip(axes.ravel(), layers):
        render_layer(
            axis, derived_path, l2_path, layer, bbox, resolution_m=resolution_m, lakes=lakes
        )
    combined.subplots_adjust(left=0.07, right=0.94, bottom=0.06, top=0.95, wspace=0.22, hspace=0.2)
    combined_png = scene_dir / f"{Path(l2_path).stem}_four_panel_{int(resolution_m)}m.png"
    combined_pdf = combined_png.with_suffix(".pdf")
    combined.savefig(combined_png, dpi=int(dpi), facecolor="white")
    combined.savefig(combined_pdf, dpi=int(dpi), facecolor="white")
    plt.close(combined)
    outputs.update({"combined_png": combined_png, "combined_pdf": combined_pdf})
    metadata = {
        "created_utc": _m_utils.utc_now(),
        "source_l2": str(l2_path),
        "source_derived": str(derived_path),
        "requested_resolution_m": float(resolution_m),
        "resampling": "nearest-neighbour",
        "bbox": bbox,
        "dpi": int(dpi),
        "hydrolakes_overlay": bool(lakes is not None and (not lakes.empty)),
        "native_swath_statistics": True,
        "regular_grid_used_only_for_maps_and_GeoTIFF": True,
    }
    metadata_path = _m_utils.write_json(scene_dir / "publication_map_metadata.json", metadata)
    outputs["metadata_json"] = metadata_path
    print(f"Saved {len(outputs)} publication files under {scene_dir}")
    return outputs


# Reference cell 26, lines 34-47.
def nearest_valid_pixel(
    derived_path: Path, target_lat: float, target_lon: float
) -> tuple[int, int, float, float]:
    with netCDF4.Dataset(derived_path) as derived:
        latitude = _m_science.to_float(derived.variables["latitude"][:])
        longitude = _m_science.to_float(derived.variables["longitude"][:])
        valid = np.ma.asarray(derived.variables["ci_valid_mask"][:]).filled(0) == 1
    scale = max(math.cos(math.radians(float(target_lat))), 0.1)
    distance = (latitude - target_lat) ** 2 + ((longitude - target_lon) * scale) ** 2
    distance[~valid | ~np.isfinite(distance)] = np.inf
    if not np.isfinite(distance).any():
        raise RuntimeError("No CI-valid water pixel is available")
    row, column = np.unravel_index(np.argmin(distance), distance.shape)
    return (int(row), int(column), float(latitude[row, column]), float(longitude[row, column]))


# Reference cell 26, lines 50-111.
def export_pixel_spectrum(
    l2_path: Path, derived_path: Path, target_lat: float, target_lon: float
) -> tuple[Path, Path, pd.DataFrame]:
    row, column, pixel_lat, pixel_lon = nearest_valid_pixel(derived_path, target_lat, target_lon)
    with netCDF4.Dataset(l2_path) as l2:
        geo = l2.groups["geophysical_data"]
        records_by_wavelength: dict[float, dict[str, Any]] = {}
        for prefix in ("rhos", "Rrs"):
            names = []
            for name in geo.variables:
                if not name.startswith(prefix + "_"):
                    continue
                try:
                    wavelength = _m_science.product_wavelength(name)
                except ValueError:
                    continue
                names.append((wavelength, name))
            for wavelength, name in sorted(names):
                record = records_by_wavelength.setdefault(
                    wavelength,
                    {
                        "target_latitude": float(target_lat),
                        "target_longitude": float(target_lon),
                        "pixel_latitude": pixel_lat,
                        "pixel_longitude": pixel_lon,
                        "row": row,
                        "column": column,
                        "wavelength_nm": wavelength,
                    },
                )
                record[f"{prefix}_product"] = name
                record[prefix] = float(_m_science.to_float(geo.variables[name][row, column]))
    table = pd.DataFrame([records_by_wavelength[key] for key in sorted(records_by_wavelength)])
    scene_dir = _m_settings.EXPORT_DIR / Path(l2_path).stem / "spectra"
    scene_dir.mkdir(parents=True, exist_ok=True)
    stem = f"spectrum_{pixel_lat:.5f}_{pixel_lon:.5f}".replace("-", "m")
    csv_path = scene_dir / f"{stem}.csv"
    png_path = scene_dir / f"{stem}.png"
    table.to_csv(csv_path, index=False)
    panels = [name for name in ("rhos", "Rrs") if name in table and table[name].notna().any()]
    figure, axes = plt.subplots(
        len(panels),
        1,
        figsize=(9, 4.5 * len(panels)),
        sharex=True,
        constrained_layout=True,
        squeeze=False,
    )
    for axis, prefix in zip(axes.ravel(), panels):
        axis.plot(table["wavelength_nm"], table[prefix], "o-", lw=1.5, ms=4)
        for wavelength in (620, 665, 681, 709):
            axis.axvline(wavelength, color="0.65", lw=0.7, ls="--")
        axis.axhline(0, color="0.35", lw=0.8)
        axis.set_ylabel(
            "L2Gen rhos (dimensionless)" if prefix == "rhos" else "L2Gen Rrs (sr$^{-1}$)"
        )
        axis.grid(alpha=0.2)
    axes[-1, 0].set_xlabel("Wavelength (nm)")
    figure.suptitle(f"Nearest CI-valid water pixel: {pixel_lat:.5f}, {pixel_lon:.5f}")
    figure.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(figure)
    print(f"Spectrum CSV: {csv_path}")
    return (csv_path, png_path, table)


# Reference cell 31, lines 9-17.
def _target_crs(longitude: np.ndarray, latitude: np.ndarray) -> CRS:
    valid = np.isfinite(longitude) & np.isfinite(latitude)
    if not valid.any():
        raise RuntimeError("No finite geolocation is available for map projection")
    center_lon = float(np.nanmedian(longitude[valid]))
    center_lat = float(np.nanmedian(latitude[valid]))
    zone = int(np.clip(math.floor((center_lon + 180.0) / 6.0) + 1, 1, 60))
    epsg = (32600 if center_lat >= 0 else 32700) + zone
    return CRS.from_epsg(epsg)


# Reference cell 31, lines 20-64.
def regular_area_definition(
    longitude: np.ndarray,
    latitude: np.ndarray,
    resolution_m: float = _m_settings.DEFAULT_MAP_RESOLUTION_M,
    max_cells: int = _m_settings.MAP_MAX_CELLS,
) -> tuple[pr_geometry.AreaDefinition, dict[str, Any]]:
    """Create a north-up square-pixel UTM grid and coarsen only for safety."""
    resolution_m = float(resolution_m)
    if resolution_m <= 0:
        raise ValueError("Map resolution must be positive")
    valid = (
        np.isfinite(longitude)
        & np.isfinite(latitude)
        & (np.abs(longitude) <= 180)
        & (np.abs(latitude) <= 90)
    )
    if not valid.any():
        raise RuntimeError("No geolocated source pixels are available")
    crs = _target_crs(longitude[valid], latitude[valid])
    transformer = Transformer.from_crs(4326, crs, always_xy=True)
    x, y = transformer.transform(longitude[valid], latitude[valid])
    xmin, xmax = (float(np.nanmin(x)), float(np.nanmax(x)))
    ymin, ymax = (float(np.nanmin(y)), float(np.nanmax(y)))
    pad = max(resolution_m, 0.002 * max(xmax - xmin, ymax - ymin))
    xmin, xmax, ymin, ymax = (xmin - pad, xmax + pad, ymin - pad, ymax + pad)
    width = max(1, int(math.ceil((xmax - xmin) / resolution_m)))
    height = max(1, int(math.ceil((ymax - ymin) / resolution_m)))
    effective = resolution_m
    if width * height > int(max_cells):
        factor = math.sqrt(width * height / int(max_cells))
        effective *= factor
        width = max(1, int(math.ceil((xmax - xmin) / effective)))
        height = max(1, int(math.ceil((ymax - ymin) / effective)))
        print(
            f"Requested map grid exceeded {max_cells:,} cells; effective resolution was increased to {effective:.1f} m for this display/export."
        )
    extent = (xmin, ymin, xmax, ymax)
    area = pr_geometry.AreaDefinition(
        "olci_regular_grid", "OLCI regular projected grid", "olci_utm", crs, width, height, extent
    )
    return (
        area,
        {
            "crs": crs,
            "crs_wkt": crs.to_wkt(),
            "epsg": crs.to_epsg(),
            "extent": extent,
            "width": width,
            "height": height,
            "requested_resolution_m": resolution_m,
            "effective_resolution_m": effective,
        },
    )


# Reference cell 31, lines 67-95.
def _resample_nearest(
    longitude: np.ndarray,
    latitude: np.ndarray,
    values: np.ndarray,
    area: pr_geometry.AreaDefinition,
    resolution_m: float,
) -> np.ndarray:
    geolocated = (
        np.isfinite(longitude)
        & np.isfinite(latitude)
        & (np.abs(longitude) <= 180)
        & (np.abs(latitude) <= 90)
    )
    swath = pr_geometry.SwathDefinition(
        lons=np.ma.array(longitude, mask=~geolocated), lats=np.ma.array(latitude, mask=~geolocated)
    )
    data = np.asarray(values)
    if data.shape[:2] != longitude.shape:
        raise ValueError(f"Data shape {data.shape} does not match geolocation {longitude.shape}")
    radius = max(600.0, 2.25 * float(resolution_m))
    if data.ndim == 2:
        return np.asarray(
            pr_kd_tree.resample_nearest(
                swath, data, area, radius_of_influence=radius, fill_value=np.nan, epsilon=0.2
            ),
            dtype=np.float32,
        )
    channels = [
        pr_kd_tree.resample_nearest(
            swath, data[..., band], area, radius_of_influence=radius, fill_value=np.nan, epsilon=0.2
        )
        for band in range(data.shape[2])
    ]
    return np.dstack(channels).astype(np.float32)


# Reference cell 31, lines 98-113.
def _rgb_stretch(rgb: np.ndarray) -> np.ndarray:
    output = np.zeros_like(rgb, dtype=np.float32)
    valid = np.all(np.isfinite(rgb), axis=2)
    for band in range(3):
        sample = rgb[..., band][valid]
        if sample.size:
            low, high = np.nanpercentile(sample, [1.0, 99.5])
        else:
            low, high = (0.0, 0.15)
        if not np.isfinite(high) or high <= low:
            low, high = (0.0, 0.15)
        output[..., band] = np.clip((rgb[..., band] - low) / (high - low), 0, 1)
    output = np.power(output, 0.82)
    output[~valid] = np.nan
    return output


# Reference cell 31, lines 496-508.
def regular_scene_raster(
    derived_path,
    l2_path,
    layer,
    bbox=None,
    resolution_m=_m_settings.DEFAULT_MAP_RESOLUTION_M,
    max_cells=_m_settings.MAP_MAX_CELLS,
):
    latitude, longitude, values, source = _native_scene_layer(derived_path, l2_path, layer, bbox)
    area, metadata = regular_area_definition(longitude, latitude, resolution_m, max_cells=max_cells)
    raster = _resample_nearest(
        longitude, latitude, values, area, metadata["effective_resolution_m"]
    )
    if layer == "true_color":
        raster = _rgb_stretch(raster)
        title = (
            "Sentinel-3 OLCI top-of-atmosphere true colour"
            if source["source"] == "rhot"
            else "Sentinel-3 OLCI Rayleigh-corrected pseudo true colour"
        )
    else:
        title = layer.replace("_", " ")
    return {
        "layer": layer,
        "data": raster,
        "area": area,
        **metadata,
        "title": title,
        "source_products": source,
        "bbox": bbox,
    }


# Reference cell 31, lines 162-165.
def _project_lakes(lakes: gpd.GeoDataFrame | None, crs: CRS):
    if lakes is None or lakes.empty:
        return None
    return lakes.to_crs(crs)


# Reference cell 31, lines 168-183.
def _add_scale_and_north(axis: plt.Axes, grid: dict[str, Any]) -> None:
    xmin, ymin, xmax, ymax = grid["extent"]
    width = xmax - xmin
    candidate = width * 0.2
    power = 10 ** math.floor(math.log10(max(candidate, 1)))
    scale = max(power, round(candidate / power) * power)
    x0, y0 = (xmin + width * 0.06, ymin + (ymax - ymin) * 0.065)
    axis.plot(
        [x0, x0 + scale],
        [y0, y0],
        color="white",
        lw=5,
        solid_capstyle="butt",
        path_effects=[path_effects.Stroke(linewidth=7, foreground="black"), path_effects.Normal()],
    )
    axis.text(
        x0 + scale / 2,
        y0 + (ymax - ymin) * 0.025,
        f"{scale / 1000:g} km",
        ha="center",
        va="bottom",
        color="white",
        fontsize=9,
        weight="bold",
        path_effects=[path_effects.withStroke(linewidth=2.5, foreground="black")],
    )
    axis.annotate(
        "N",
        xy=(0.94, 0.91),
        xytext=(0.94, 0.79),
        xycoords="axes fraction",
        ha="center",
        va="center",
        color="white",
        fontsize=11,
        weight="bold",
        arrowprops=dict(facecolor="white", edgecolor="black", width=2.5, headwidth=9),
        path_effects=[path_effects.withStroke(linewidth=2, foreground="black")],
    )


# Reference cell 31, lines 666-708.
def interactive_regular_map(
    derived_path,
    l2_path,
    layer="true_color",
    bbox=None,
    resolution_m=_m_settings.DEFAULT_MAP_RESOLUTION_M,
    lakes=None,
    vmin=None,
    vmax=None,
    destination=None,
):
    """Return a standalone Leaflet map with true longitude/latitude interaction."""
    grid = regular_geographic_raster(
        derived_path,
        l2_path,
        layer,
        bbox,
        resolution_m,
        max_cells=_m_settings.INTERACTIVE_MAX_CELLS,
    )
    west, south, east, north = grid["extent"]
    water = retained = None
    if layer in {"CI", "CI_cyano", "NDCI"}:

        def companion(name):
            lat, lon, values, _ = _native_scene_layer(derived_path, l2_path, name, bbox)
            return _resample_nearest(lon, lat, values, grid["area"], grid["effective_resolution_m"])

        water = companion("water_mask") >= 0.5
        retained = companion("valid_water_mask") >= 0.5
    rgba, legend = _folium_rgba(grid["data"], layer, vmin, vmax, water, retained)
    map_object = folium.Map(
        location=[(south + north) / 2, (west + east) / 2],
        zoom_start=8,
        tiles=None,
        control_scale=True,
        prefer_canvas=True,
    )
    folium.TileLayer("CartoDB positron", name="Light basemap", control=True).add_to(map_object)
    folium.TileLayer("OpenStreetMap", name="OpenStreetMap", control=True).add_to(map_object)
    folium.raster_layers.ImageOverlay(
        rgba,
        bounds=[[south, west], [north, east]],
        name=f"OLCI {layer}",
        opacity=1.0,
        interactive=True,
        cross_origin=False,
        zindex=3,
    ).add_to(map_object)
    if lakes is not None and (not lakes.empty):
        folium.GeoJson(
            lakes.to_crs(4326).__geo_interface__,
            name="HydroLAKES",
            style_function=lambda feature: {"color": "#00bcd4", "weight": 1.2, "fillOpacity": 0},
        ).add_to(map_object)
    if legend is not None:
        legend.add_to(map_object)
    folium_plugins.Fullscreen(position="topright").add_to(map_object)
    folium_plugins.MousePosition(
        position="bottomright", separator=" | ", prefix="Lat / Lon", num_digits=5
    ).add_to(map_object)
    folium.LatLngPopup().add_to(map_object)
    folium_plugins.MeasureControl(position="topleft", primary_length_unit="kilometers").add_to(
        map_object
    )
    folium.LayerControl(collapsed=False).add_to(map_object)
    map_object.fit_bounds([[south, west], [north, east]])
    if destination is None:
        destination = (
            _m_settings.EXPORT_DIR
            / Path(l2_path).stem
            / "interactive_maps"
            / f"{layer}_leaflet.html"
        )
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    map_object.save(str(destination))
    map_object._olci_html_path = destination
    print(f"Interactive latitude/longitude map: {destination}")
    return map_object


# Reference cell 31, lines 299-330.
def write_regular_geotiff(grid: dict[str, Any], path: Path) -> Path:
    """Write a cloud-optimized GeoTIFF when the rasterio COG driver is available."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = grid["data"]
    if data.ndim == 3:
        array = np.moveaxis(np.uint8(np.clip(np.nan_to_num(data, nan=0), 0, 1) * 255), -1, 0)
        dtype, nodata, count = ("uint8", 0, array.shape[0])
    else:
        array = np.where(np.isfinite(data), data, _m_settings.DERIVED_FILL).astype(np.float32)[
            None, ...
        ]
        dtype, nodata, count = ("float32", float(_m_settings.DERIVED_FILL), 1)
    transform = from_bounds(*grid["extent"], grid["width"], grid["height"])
    profile = dict(
        driver="COG",
        height=grid["height"],
        width=grid["width"],
        count=count,
        dtype=dtype,
        crs=grid["crs"],
        transform=transform,
        nodata=nodata,
        compress="DEFLATE",
        blocksize=512,
        overview_resampling="nearest",
    )
    try:
        with rasterio.open(path, "w", **profile) as destination:
            destination.write(array)
            destination.update_tags(
                source_layer=grid["layer"],
                resampling="nearest",
                effective_resolution_m=grid["effective_resolution_m"],
            )
    except Exception:
        path.unlink(missing_ok=True)
        profile.update(driver="GTiff", tiled=True, blockxsize=512, blockysize=512)
        profile.pop("blocksize", None)
        profile.pop("overview_resampling", None)
        with rasterio.open(path, "w", **profile) as destination:
            destination.write(array)
            levels = [
                level
                for level in (2, 4, 8, 16)
                if min(grid["width"], grid["height"]) // level >= 128
            ]
            if levels:
                destination.build_overviews(levels, rasterio.enums.Resampling.nearest)
    return path


# Reference cell 31, lines 383-442.
def analyze_lakes(
    derived_path: Path, lakes: gpd.GeoDataFrame, shore_buffer_m: float = 0.0
) -> pd.DataFrame:
    """Native-swath per-lake statistics; map resampling never enters these values."""
    if lakes is None or lakes.empty:
        raise ValueError("No scene HydroLAKES polygons are loaded")
    analysis_lakes = _m_hydrolakes.erode_lake_geometries(lakes, shore_buffer_m)
    with netCDF4.Dataset(derived_path) as ds:
        lat = _m_science.to_float(ds.variables["latitude"][:])
        lon = _m_science.to_float(ds.variables["longitude"][:])
        retained = np.ma.asarray(ds.variables["valid_water_mask"][:]).filled(0) == 1
        ci_valid = np.ma.asarray(ds.variables["ci_valid_mask"][:]).filled(0) == 1
        detected = np.ma.asarray(ds.variables["ci_detection_mask"][:]).filled(0) == 1
        bright = (
            np.ma.asarray(ds.variables["bright_pixel_mask"][:]).filled(0) == 1
            if "bright_pixel_mask" in ds.variables
            else np.zeros(lat.shape, bool)
        )
        ci = _m_science.to_float(ds.variables["CI"][:])
        cyano = _m_science.to_float(ds.variables["CI_cyano"][:])
        ndci = _m_science.to_float(ds.variables["NDCI"][:])
    records = []
    for _, lake in tqdm(
        analysis_lakes.iterrows(), total=len(analysis_lakes), desc="Lake statistics"
    ):
        west, south, east, north = lake.geometry.bounds
        window = (lon >= west) & (lon <= east) & (lat >= south) & (lat <= north)
        rows, columns = np.where(window)
        if not rows.size:
            polygon = np.zeros(lat.shape, dtype=bool)
        else:
            rs = slice(rows.min(), rows.max() + 1)
            cs = slice(columns.min(), columns.max() + 1)
            local = contains_xy(lake.geometry, lon[rs, cs], lat[rs, cs])
            polygon = np.zeros(lat.shape, dtype=bool)
            polygon[rs, cs] = local
        kept = polygon & retained
        valid_ci = polygon & ci_valid
        hits = polygon & detected
        ci_values = ci[polygon & np.isfinite(ci)]
        cyano_values = cyano[hits & np.isfinite(cyano)]
        ndci_values = ndci[kept & np.isfinite(ndci)]

        def percentile(values, q):
            return float(np.nanpercentile(values, q)) if values.size else np.nan

        record = {
            "Hylak_id": int(lake["Hylak_id"]),
            "Lake_name": lake.get("Lake_name", "") or "",
            "Lake_area_km2": float(lake.get("Lake_area", np.nan)),
            "native_polygon_pixels": int(polygon.sum()),
            "retained_pixels": int(kept.sum()),
            "retained_percent": 100.0 * kept.sum() / polygon.sum() if polygon.any() else np.nan,
            "bright_screen_pixels": int((polygon & bright).sum()),
            "ci_input_valid_pixels": int(valid_ci.sum()),
            "ci_candidate_pixels": int((polygon & np.isfinite(ci)).sum()),
            "ci_cyano_detection_pixels": int(hits.sum()),
            "ci_cyano_detection_fraction_of_ci_valid": hits.sum() / valid_ci.sum()
            if valid_ci.any()
            else np.nan,
            "ci_cyano_detection_fraction_of_retained": hits.sum() / kept.sum()
            if kept.any()
            else np.nan,
            "CI_median_candidates": percentile(ci_values, 50),
            "CI_p10_candidates": percentile(ci_values, 10),
            "CI_p90_candidates": percentile(ci_values, 90),
            "CIcyano_median_detections": percentile(cyano_values, 50),
            "CIcyano_p10_detections": percentile(cyano_values, 10),
            "CIcyano_p90_detections": percentile(cyano_values, 90),
            "NDCI_median_retained": percentile(ndci_values, 50),
            "NDCI_p10_retained": percentile(ndci_values, 10),
            "NDCI_p90_retained": percentile(ndci_values, 90),
        }
        records.append(record)
    return pd.DataFrame(records).sort_values(["Lake_area_km2", "Hylak_id"], ascending=[False, True])


# Reference cell 31, lines 445-450.
def export_lake_statistics(derived_path: Path, table: pd.DataFrame) -> Path:
    destination = (
        _m_settings.EXPORT_DIR / Path(derived_path).stem / "hydrolakes_native_swath_statistics.csv"
    )
    destination.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(destination, index=False)
    print(f"Lake statistics: {destination}")
    return destination


# Reference cell 31, lines 468-493.
def _native_scene_layer(derived_path, l2_path, layer, bbox=None):
    with netCDF4.Dataset(derived_path) as derived:
        lat_full = _m_science.to_float(derived.variables["latitude"][:])
        lon_full = _m_science.to_float(derived.variables["longitude"][:])
        row_slice, col_slice = bbox_slices(lat_full, lon_full, bbox, margin=8)
        latitude = lat_full[row_slice, col_slice]
        longitude = lon_full[row_slice, col_slice]
        if layer not in {"true_color"} and (not layer.startswith("rhos_")):
            if layer not in derived.variables:
                raise KeyError(f"Layer not in derived NetCDF: {layer}")
            return (
                latitude,
                longitude,
                _m_science.to_float(derived.variables[layer][row_slice, col_slice]),
                None,
            )
    with netCDF4.Dataset(l2_path) as root:
        geo = root.groups["geophysical_data"]
        if layer.startswith("rhos_"):
            wave = int(layer.split("_")[1])
            name = _m_science.nearest_product(geo.variables, "rhos", wave, 25)
            values = _m_science.to_float(geo.variables[name][row_slice, col_slice])
            return (latitude, longitude, values, [name])
        try:
            names = [
                _m_science.nearest_product(geo.variables, "rhot", wave, 20)
                for wave in (665, 560, 443)
            ]
            source = "rhot"
        except KeyError:
            names = [
                _m_science.nearest_product(geo.variables, "rhos", wave, 20)
                for wave in (665, 560, 443)
            ]
            source = "rhos"
        values = np.dstack(
            [_m_science.to_float(geo.variables[name][row_slice, col_slice]) for name in names]
        )
        return (latitude, longitude, values, {"products": names, "source": source})


# Reference cell 31, lines 511-516.
def _grid_lonlat_bounds(grid):
    xmin, ymin, xmax, ymax = grid["extent"]
    transformer = Transformer.from_crs(grid["crs"], 4326, always_xy=True)
    xs = np.array([xmin, xmax, xmax, xmin])
    ys = np.array([ymin, ymin, ymax, ymax])
    lon, lat = transformer.transform(xs, ys)
    return (float(np.min(lon)), float(np.min(lat)), float(np.max(lon)), float(np.max(lat)))


# Reference cell 31, lines 612-640.
def regular_geographic_raster(
    derived_path,
    l2_path,
    layer="true_color",
    bbox=None,
    resolution_m=_m_settings.DEFAULT_MAP_RESOLUTION_M,
    max_cells=_m_settings.MAP_MAX_CELLS,
):
    """Regular EPSG:4326 preview for browser maps; native science arrays are unchanged."""
    latitude, longitude, values, source = _native_scene_layer(derived_path, l2_path, layer, bbox)
    valid = np.isfinite(latitude) & np.isfinite(longitude)
    west, east = (float(np.nanmin(longitude[valid])), float(np.nanmax(longitude[valid])))
    south, north = (float(np.nanmin(latitude[valid])), float(np.nanmax(latitude[valid])))
    center_lat = 0.5 * (south + north)
    dy = float(resolution_m) / 111320.0
    dx = float(resolution_m) / max(111320.0 * math.cos(math.radians(center_lat)), 1.0)
    width = max(1, int(math.ceil((east - west) / dx)))
    height = max(1, int(math.ceil((north - south) / dy)))
    effective = float(resolution_m)
    if width * height > int(max_cells):
        factor = math.sqrt(width * height / int(max_cells))
        effective *= factor
        dx *= factor
        dy *= factor
        width = max(1, int(math.ceil((east - west) / dx)))
        height = max(1, int(math.ceil((north - south) / dy)))
        print(
            f"Display grid was coarsened to {effective:.1f} m ({width * height:,} cells); native NetCDF values are unchanged."
        )
    area = pr_geometry.AreaDefinition(
        "olci_leaflet",
        "OLCI geographic preview",
        "latlon",
        CRS.from_epsg(4326),
        width,
        height,
        (west, south, east, north),
    )
    data = _resample_nearest(longitude, latitude, values, area, effective)
    if layer == "true_color":
        data = _rgb_stretch(data)
    return {
        "layer": layer,
        "data": data,
        "area": area,
        "crs": CRS.from_epsg(4326),
        "epsg": 4326,
        "extent": (west, south, east, north),
        "width": width,
        "height": height,
        "requested_resolution_m": float(resolution_m),
        "effective_resolution_m": effective,
        "bbox": bbox,
        "source_products": source,
    }


# Reference cell 31, lines 643-663.
def _folium_rgba(data, layer, vmin=None, vmax=None, water=None, retained=None):
    if layer == "true_color":
        rgb = np.uint8(np.clip(np.nan_to_num(data, nan=0), 0, 1) * 255)
        alpha = np.uint8(np.all(np.isfinite(data), axis=2) * 255)
        return (np.dstack([rgb, alpha]), None)
    if layer in MASK_LAYER_STYLE or layer in {
        "l2_water_mask",
        "hydrolakes_water_mask",
        "bright_pixel_mask",
    }:
        cmap_name, label, low, high = ("viridis", layer.replace("_", " "), 0, 1)
    else:
        cmap_name, defaults, label = LAYER_STYLE[layer]
        low, high = robust_limits(data, defaults)
    low = low if vmin is None else float(vmin)
    high = high if vmax is None else float(vmax)
    norm = mcolors.Normalize(low, high, clip=True)
    rgba = np.uint8(plt.get_cmap(cmap_name)(norm(np.nan_to_num(data, nan=low))) * 255)
    rgba[..., 3] = np.uint8(np.isfinite(data) * 220)
    if layer in {"CI", "CI_cyano", "NDCI"} and water is not None and (retained is not None):
        rgba[water & ~retained] = np.array([37, 37, 37, 235], dtype=np.uint8)
        rgba[retained & ~np.isfinite(data)] = np.array([199, 199, 199, 235], dtype=np.uint8)
    colors = [mcolors.to_hex(plt.get_cmap(cmap_name)(value)) for value in np.linspace(0, 1, 7)]
    legend = LinearColormap(colors, vmin=low, vmax=high, caption=label)
    return (rgba, legend)
