"""Standalone spectra implementation."""

from __future__ import annotations

import contextlib

import gradio as gr
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree

from . import archive as _m_archive
from . import callbacks as _m_callbacks
from . import quality as _m_quality
from . import science as _m_science
from . import settings as _m_settings
from . import ui_helpers as _m_ui_helpers


# Reference cell 50, lines 965-977.
def scene_spectrum_callback(scene_id, lake_id, statistic="median"):
    try:
        if not scene_id or not lake_id:
            raise ValueError("Choose a completed scene and lake")
        f = _m_archive._load_scene_pixels(scene_id, "lake", lake_id)
        waves = list(_m_settings.GLOBAL_ANALYSIS_RHOS_WAVELENGTHS)
        rows = []
        for w in waves:
            a = pd.to_numeric(f[f"rhos_{w}"], errors="coerce").to_numpy(float)
            a = a[np.isfinite(a)]
            rows.append(
                {
                    "wavelength_nm": w,
                    "mean_rhos": float(np.nanmean(a)) if a.size else np.nan,
                    "median_rhos": float(np.nanmedian(a)) if a.size else np.nan,
                    "std_rhos": float(np.nanstd(a)) if a.size else np.nan,
                    "n": int(a.size),
                }
            )
        tab = pd.DataFrame(rows)
        y = "median_rhos" if statistic == "median" else "mean_rhos"
        fig = px.line(
            tab,
            x="wavelength_nm",
            y=y,
            markers=True,
            title=f"Hylak_id {lake_id} · {statistic} Rayleigh-corrected rhos",
        )
        stats = _m_callbacks._stats()
        stats = (
            stats[
                (stats.scene_id.astype(str) == str(scene_id))
                & (pd.to_numeric(stats.Hylak_id, errors="coerce") == int(lake_id))
            ]
            if not stats.empty
            else pd.DataFrame()
        )
        d = _m_settings.EXPORT_DIR / "spectra"
        d.mkdir(parents=True, exist_ok=True)
        csv = d / f"{scene_id}_Hylak_{int(lake_id)}_spectrum.csv"
        tab.to_csv(csv, index=False)
        return (
            fig,
            tab,
            stats,
            str(csv),
            _m_ui_helpers._app_status(
                "Spectrum ready", f"{len(f):,} retained/native lake-pixel records"
            ),
        )
    except Exception as exc:
        return (
            None,
            pd.DataFrame(),
            pd.DataFrame(),
            None,
            _m_ui_helpers._app_status("Spectrum failed", str(exc), False),
        )


# Reference cell 50, lines 982-996.
def _spectrum_points_frame(points):
    if points is None:
        frame = pd.DataFrame(columns=["label", "latitude", "longitude"])
    elif isinstance(points, pd.DataFrame):
        frame = points.copy()
    else:
        frame = pd.DataFrame(points, columns=["label", "latitude", "longitude"])
    frame = frame.rename(columns={c: str(c).strip().lower() for c in frame.columns})
    req = ["label", "latitude", "longitude"]
    if any((c not in frame.columns for c in req)):
        raise ValueError("Point table requires label, latitude, longitude")
    frame = frame[req].copy()
    frame["label"] = frame["label"].fillna("").astype(str).str.strip()
    frame["latitude"] = pd.to_numeric(frame["latitude"], errors="coerce")
    frame["longitude"] = pd.to_numeric(frame["longitude"], errors="coerce")
    frame = frame.dropna(subset=["latitude", "longitude"]).reset_index(drop=True)
    if frame.empty:
        raise ValueError("Add at least one latitude/longitude point")
    if len(frame) > 20:
        raise ValueError("Maximum 20 coordinate spectrum points")
    for i in frame.index:
        if not frame.at[i, "label"]:
            frame.at[i, "label"] = f"P{i + 1}"
    if ((frame.latitude < -90) | (frame.latitude > 90)).any() or (
        (frame.longitude < -180) | (frame.longitude > 180)
    ).any():
        raise ValueError("Invalid latitude/longitude")
    return frame


# Reference cell 50, lines 999-1004.
def add_spectrum_point_callback(label, latitude, longitude, points):
    try:
        f = (
            pd.DataFrame(columns=["label", "latitude", "longitude"])
            if points is None
            else points.copy()
            if isinstance(points, pd.DataFrame)
            else pd.DataFrame(points, columns=["label", "latitude", "longitude"])
        )
        row = pd.DataFrame(
            [
                {
                    "label": str(label or f"P{len(f) + 1}"),
                    "latitude": float(latitude),
                    "longitude": float(longitude),
                }
            ]
        )
        return (
            pd.concat([f, row], ignore_index=True),
            _m_ui_helpers._app_status("Spectrum point added"),
        )
    except Exception as exc:
        return (points, _m_ui_helpers._app_status("Could not add point", str(exc), False))


# Reference cell 50, lines 1007-1008.
def clear_spectrum_points_callback():
    return (
        pd.DataFrame(columns=["label", "latitude", "longitude"]),
        _m_ui_helpers._app_status("Spectrum points cleared"),
    )


# Reference cell 50, lines 1016-1040.
def coordinate_spectra_callback(scene_id, points, pixel_basis="stored_lake", max_distance_m=2000):
    try:
        if not scene_id:
            raise ValueError("Choose a completed scene")
        targets = _spectrum_points_frame(points)
        f = _m_archive._load_scene_pixels(scene_id, "scene", None)
        mask = np.isfinite(f.latitude) & np.isfinite(f.longitude)
        if pixel_basis == "valid_water" and "valid_water_mask" in f:
            mask &= pd.to_numeric(f.valid_water_mask, errors="coerce").fillna(0).to_numpy(float) > 0
        if pixel_basis == "ci_valid" and "ci_valid_mask" in f:
            mask &= pd.to_numeric(f.ci_valid_mask, errors="coerce").fillna(0).to_numpy(float) > 0
        source = f.loc[mask].reset_index(drop=True)
        if source.empty:
            raise ValueError("No stored pixels satisfy the selected nearest-pixel basis")
        tree = cKDTree(_m_science._unit_sphere_xyz(source.latitude, source.longitude))
        records = []
        for _, p in targets.iterrows():
            chord, idx = tree.query(
                _m_science._unit_sphere_xyz([p.latitude], [p.longitude])[0], k=1
            )
            pix = source.iloc[int(idx)]
            angle = 2 * np.arcsin(min(1.0, float(chord) / 2))
            dist = 6371000.0 * angle
            if dist > float(max_distance_m):
                raise ValueError(
                    f"{p.label}: nearest stored lake pixel is {dist:.0f} m away, beyond {float(max_distance_m):.0f} m limit"
                )
            for w in _m_settings.GLOBAL_ANALYSIS_RHOS_WAVELENGTHS:
                col = f"rhos_{w}"
                value = float(pix[col]) if col in pix and np.isfinite(pix[col]) else np.nan
                records.append(
                    {
                        "point": p.label,
                        "target_latitude": float(p.latitude),
                        "target_longitude": float(p.longitude),
                        "pixel_latitude": float(pix.latitude),
                        "pixel_longitude": float(pix.longitude),
                        "Hylak_id": int(pix.Hylak_id),
                        "match_distance_m": float(dist),
                        "pixel_basis": pixel_basis,
                        "product": "rhos_884" if w == 885 else col,
                        "wavelength_nm": 884.0 if w == 885 else float(w),
                        "value": value,
                    }
                )
        tab = pd.DataFrame(records)
        fig = go.Figure()
        symbols = [
            "circle",
            "square",
            "diamond",
            "cross",
            "triangle-up",
            "triangle-down",
            "star",
            "hexagon",
            "x",
        ]
        for i, (label, gp) in enumerate(tab.groupby("point", sort=False)):
            finite = gp[np.isfinite(gp.value)]
            fig.add_trace(
                go.Scatter(
                    x=finite.wavelength_nm,
                    y=finite.value,
                    mode="lines+markers",
                    name=str(label),
                    marker={"symbol": symbols[i % len(symbols)], "size": 7},
                    customdata=finite[
                        ["pixel_latitude", "pixel_longitude", "match_distance_m", "Hylak_id"]
                    ].to_numpy(),
                    hovertemplate="%{fullData.name}<br>%{x:.1f} nm<br>%{y:.7g}<br>pixel %{customdata[0]:.5f}, %{customdata[1]:.5f}<br>distance %{customdata[2]:.0f} m · Hylak %{customdata[3]}<extra></extra>",
                )
            )
        fig.update_layout(
            title=f"Sentinel-3 OLCI coordinate spectra · {len(targets)} point(s)",
            xaxis_title="Wavelength (nm)",
            yaxis_title="Rayleigh-corrected rhos",
            template="plotly_white",
            height=560,
            hovermode="closest",
        )
        d = _m_settings.EXPORT_DIR / "spectra"
        d.mkdir(parents=True, exist_ok=True)
        csv = d / f"{scene_id}_coordinate_spectra.csv"
        tab.to_csv(csv, index=False)
        return (
            fig,
            tab.round(8),
            str(csv),
            _m_ui_helpers._app_status(
                "Coordinate spectra ready",
                f"{len(targets)} requested point(s); nearest stored HydroLAKES lake pixels; non-lake pixels are not retained in this efficient archive",
            ),
        )
    except Exception as exc:
        return (
            None,
            pd.DataFrame(),
            None,
            _m_ui_helpers._app_status("Coordinate spectra failed", str(exc), False),
        )


# Reference cell 58, lines 516-565.
def _spectrum_map_image(scene_id, points=None, matched_table=None, width=1000, height=560):
    if not scene_id:
        return (None, {})
    f = _m_archive._load_scene_pixels(scene_id, "scene", None)
    finite = np.isfinite(pd.to_numeric(f.latitude, errors="coerce")) & np.isfinite(
        pd.to_numeric(f.longitude, errors="coerce")
    )
    f = f.loc[finite].copy()
    if f.empty:
        return (None, {})
    if len(f) > 180000:
        take = np.linspace(0, len(f) - 1, 180000, dtype=int)
        draw_f = f.iloc[take]
    else:
        draw_f = f
    west, east = (float(f.longitude.min()), float(f.longitude.max()))
    south, north = (float(f.latitude.min()), float(f.latitude.max()))
    dx = max(east - west, 0.0001)
    dy = max(north - south, 0.0001)
    west -= 0.03 * dx
    east += 0.03 * dx
    south -= 0.03 * dy
    north += 0.03 * dy
    img = Image.new("RGB", (int(width), int(height)), (244, 247, 250))
    dr = ImageDraw.Draw(img)
    dr.rectangle([0, 0, width - 1, height - 1], outline=(80, 90, 100), width=1)

    def xy(lon, lat):
        x = (float(lon) - west) / (east - west) * (width - 1)
        y = (north - float(lat)) / (north - south) * (height - 1)
        return (int(round(x)), int(round(y)))

    with contextlib.suppress(Exception):
        b = _m_archive._lake_boundaries_for_ids(draw_f.Hylak_id.unique())
        if b is not None and (not b.empty):
            b = b.to_crs(4326)

            def draw_geom(geom):
                if geom is None or geom.is_empty:
                    return
                if geom.geom_type == "Polygon":
                    coords = list(geom.exterior.coords)
                    pts = [xy(xx, yy) for xx, yy in coords]
                    if len(pts) > 1:
                        dr.line(pts, fill=(0, 165, 195), width=1)
                elif geom.geom_type in {"MultiPolygon", "GeometryCollection"}:
                    for gg in geom.geoms:
                        draw_geom(gg)

            for geom in b.geometry:
                draw_geom(geom)
    valid = _m_quality._primary_selector(draw_f)
    lonv = draw_f.longitude.to_numpy(float)
    latv = draw_f.latitude.to_numpy(float)
    for x0, y0, v in zip(lonv, latv, valid):
        x, y = xy(x0, y0)
        col = (17, 94, 140) if v else (175, 182, 190)
        dr.point((x, y), fill=col)
    if points is not None:
        with contextlib.suppress(Exception):
            pf = _spectrum_points_frame(points)
            for _, r in pf.iterrows():
                x, y = xy(r.longitude, r.latitude)
                dr.ellipse([x - 5, y - 5, x + 5, y + 5], outline=(220, 38, 38), width=2)
                dr.line([x - 7, y, x + 7, y], fill=(220, 38, 38), width=1)
                dr.line([x, y - 7, x, y + 7], fill=(220, 38, 38), width=1)
    if (
        matched_table is not None
        and isinstance(matched_table, pd.DataFrame)
        and (not matched_table.empty)
    ):
        mt = matched_table.drop_duplicates("point")
        for _, r in mt.iterrows():
            if np.isfinite(r.get("pixel_longitude", np.nan)) and np.isfinite(
                r.get("pixel_latitude", np.nan)
            ):
                x, y = xy(r.pixel_longitude, r.pixel_latitude)
                dr.rectangle([x - 4, y - 4, x + 4, y + 4], outline=(0, 150, 70), width=2)
    dr.rectangle([6, 6, 365, 45], fill=(255, 255, 255), outline=(180, 185, 190))
    dr.text((12, 10), "Click map to add spectrum point", fill=(20, 30, 40))
    dr.text(
        (12, 26),
        "blue=valid native pixels · gray=QA excluded · red=target · green=matched",
        fill=(55, 65, 75),
    )
    state = {
        "west": west,
        "east": east,
        "south": south,
        "north": north,
        "width": int(width),
        "height": int(height),
        "scene_id": str(scene_id),
    }
    return (img, state)


# Reference cell 58, lines 568-573.
def spectrum_map_refresh_callback(scene_id, points):
    try:
        img, state = _spectrum_map_image(scene_id, points)
        return (
            img,
            state,
            _m_ui_helpers._app_status(
                "Spectrum map ready",
                "Click a location to append a coordinate and immediately plot the nearest native-pixel spectrum.",
            ),
        )
    except Exception as exc:
        return (None, {}, _m_ui_helpers._app_status("Spectrum map failed", str(exc), False))


# Reference cell 58, lines 576-582.
def coordinate_spectra_with_map_callback(
    scene_id, points, pixel_basis="stored_lake", max_distance_m=2000
):
    fig, tab, csv, status = coordinate_spectra_callback(
        scene_id, points, pixel_basis, max_distance_m
    )
    try:
        img, state = _spectrum_map_image(scene_id, points, tab)
    except Exception:
        img, state = (None, {})
    return (img, state, fig, tab, csv, status)


# Reference cell 58, lines 585-591.
def clear_spectrum_points_and_map_callback(scene_id):
    empty = pd.DataFrame(columns=["label", "latitude", "longitude"])
    try:
        img, state = _spectrum_map_image(scene_id, empty)
    except Exception:
        img, state = (None, {})
    return (empty, img, state, _m_ui_helpers._app_status("Spectrum points cleared"))


# Reference cell 58, lines 594-614.
def spectrum_map_click_callback(
    scene_id, points, pixel_basis, max_distance_m, map_state, evt: gr.SelectData
):
    try:
        if not scene_id:
            raise ValueError("Choose a completed scene")
        st = dict(map_state or {})
        if str(st.get("scene_id")) != str(scene_id):
            _, st = _spectrum_map_image(scene_id, points)
        idx = evt.index
        if not isinstance(idx, (tuple, list)) or len(idx) < 2:
            raise ValueError("Could not read map click coordinates")
        x = float(idx[0])
        y = float(idx[1])
        width = float(st["width"])
        height = float(st["height"])
        lon = float(st["west"] + x / max(width - 1, 1) * (st["east"] - st["west"]))
        lat = float(st["north"] - y / max(height - 1, 1) * (st["north"] - st["south"]))
        current = (
            pd.DataFrame(columns=["label", "latitude", "longitude"])
            if points is None
            else points.copy()
            if isinstance(points, pd.DataFrame)
            else pd.DataFrame(points, columns=["label", "latitude", "longitude"])
        )
        label = f"Click{len(current) + 1}"
        current = pd.concat(
            [current, pd.DataFrame([{"label": label, "latitude": lat, "longitude": lon}])],
            ignore_index=True,
        )
        fig, tab, csv, status = coordinate_spectra_callback(
            scene_id, current, pixel_basis, max_distance_m
        )
        img, new_state = _spectrum_map_image(scene_id, current, tab)
        return (
            current,
            img,
            new_state,
            fig,
            tab,
            csv,
            _m_ui_helpers._app_status(
                "Clicked spectrum ready",
                f"{label}: target {lat:.5f}, {lon:.5f}. The spectrum uses the nearest eligible stored native OLCI pixel.",
            ),
        )
    except Exception as exc:
        return (
            points,
            None,
            map_state,
            None,
            pd.DataFrame(),
            None,
            _m_ui_helpers._app_status("Map-click spectrum failed", str(exc), False),
        )
