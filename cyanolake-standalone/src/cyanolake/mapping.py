"""Standalone mapping implementation."""

from __future__ import annotations

import base64
import html as html_lib
import io
import json
import math
from urllib.parse import quote

import folium
import geopandas as gpd
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from folium import plugins as folium_plugins
from PIL import Image
from pyproj import Transformer
from scipy.spatial import cKDTree
from shapely import contains_xy

from . import archive as _m_archive
from . import hydrolakes as _m_hydrolakes
from . import quality as _m_quality
from . import science as _m_science
from . import settings as _m_settings
from . import state as _m_state
from . import ui_helpers as _m_ui_helpers

_TO_WEBMERC = Transformer.from_crs(4326, 3857, always_xy=True)

_FROM_WEBMERC = Transformer.from_crs(3857, 4326, always_xy=True)


# Reference cell 58, lines 69-101.
def _grid_geometry_mercator(lat, lon, resolution_m, max_cells=_m_settings.V265_DISPLAY_MAX_CELLS):
    """Exact EPSG:3857 display grid; never silently changes the user-requested resolution."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    finite = np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) < 85.0)
    if not finite.any():
        raise ValueError("No finite Web-Mercator-compatible geolocation")
    x, y = _TO_WEBMERC.transform(lon[finite], lat[finite])
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    res = float(resolution_m)
    if res <= 0:
        raise ValueError("Display resolution must be positive")
    xmin = math.floor(float(np.nanmin(x)) / res) * res - res
    xmax = math.ceil(float(np.nanmax(x)) / res) * res + res
    ymin = math.floor(float(np.nanmin(y)) / res) * res - res
    ymax = math.ceil(float(np.nanmax(y)) / res) * res + res
    nx = max(2, int(round((xmax - xmin) / res)))
    ny = max(2, int(round((ymax - ymin) / res)))
    xmax = xmin + nx * res
    ymax = ymin + ny * res
    ncell = int(nx) * int(ny)
    if ncell > int(max_cells):
        raise RuntimeError(
            f"Exact {res:.0f} m display would require {ncell / 1000000.0:.1f} million cells, above the {int(max_cells) / 1000000.0:.0f}-million browser safety limit. Choose 'Active planning region' or 'Selected lake', or explicitly choose a coarser preview. The app will not silently change resolution."
        )
    west, south = _FROM_WEBMERC.transform(xmin, ymin)
    east, north = _FROM_WEBMERC.transform(xmax, ymax)
    return {
        "xmin": xmin,
        "xmax": xmax,
        "ymin": ymin,
        "ymax": ymax,
        "nx": nx,
        "ny": ny,
        "ncell": ncell,
        "resolution_m": res,
        "west": float(west),
        "east": float(east),
        "south": float(south),
        "north": float(north),
        "crs": "EPSG:3857",
    }


# Reference cell 59, lines 177-191.
def _regularize_dataframe_mercator(frame, columns, resolution_m, scene_id=None, clip_geom=None):
    """Exact requested display grid with same-lake nearest-native assignment."""
    if not scene_id:
        raise ValueError("scene_id is required for same-lake Web-Mercator remapping")
    lat = pd.to_numeric(frame.latitude, errors="coerce").to_numpy(float)
    lon = pd.to_numeric(frame.longitude, errors="coerce").to_numpy(float)
    g = _grid_geometry_mercator(
        lat, lon, float(resolution_m), max_cells=_m_settings.V265_DISPLAY_MAX_CELLS
    )
    rows, cols, src, dist = _native_assignments_mercator_v267(
        frame, g, scene_id, clip_geom=clip_geom
    )
    out = {c: np.full((g["ny"], g["nx"]), np.nan, np.float32) for c in columns}
    matched = src >= 0
    for c in columns:
        if c not in frame.columns or not matched.any():
            continue
        vals = pd.to_numeric(frame[c], errors="coerce").to_numpy(float)
        good = matched & np.isfinite(vals[np.maximum(src, 0)])
        if good.any():
            out[c][rows[good], cols[good]] = vals[src[good]].astype(np.float32)
    g.update(
        {
            "mapping_method": _m_settings.V267_DISPLAY_METHOD,
            "target_lake_cells": int(len(rows)),
            "matched_target_cells": int(matched.sum()),
            "nearest_radius_m": max(
                _m_settings.V267_NEAREST_RADIUS_M, 2.0 * float(g["resolution_m"])
            ),
        }
    )
    return (out, g)


# Reference cell 50, lines 448-464.
def _grid_geometry_geographic(lat, lon, resolution_m, max_cells=_m_settings.DISPLAY_MAX_CELLS):
    """Consistent EPSG:4326 grid used only for SNAP/GIS NetCDF export."""
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    finite = np.isfinite(lat) & np.isfinite(lon)
    if not finite.any():
        raise ValueError("No finite geolocation")
    mean_lat = float(np.nanmean(lat[finite]))
    res = float(resolution_m)
    for _ in range(10):
        dlat = res / 111320.0
        dlon = res / max(111320.0 * np.cos(np.deg2rad(mean_lat)), 1000.0)
        west = math.floor(float(np.nanmin(lon[finite])) / dlon) * dlon - dlon
        east = math.ceil(float(np.nanmax(lon[finite])) / dlon) * dlon + dlon
        south = math.floor(float(np.nanmin(lat[finite])) / dlat) * dlat - dlat
        north = math.ceil(float(np.nanmax(lat[finite])) / dlat) * dlat + dlat
        nx = max(2, int(round((east - west) / dlon)))
        ny = max(2, int(round((north - south) / dlat)))
        east = west + nx * dlon
        north = south + ny * dlat
        if nx * ny <= max_cells:
            break
        res *= math.sqrt(nx * ny / max_cells) * 1.05
    return {
        "west": west,
        "east": east,
        "south": south,
        "north": north,
        "dlat": dlat,
        "dlon": dlon,
        "nx": nx,
        "ny": ny,
        "resolution_m": res,
        "crs": "EPSG:4326",
    }


# Reference cell 59, lines 112-130.
def _regularize_dataframe_geographic(frame, columns, resolution_m, scene_id=None, clip_geom=None):
    """Exact nominal-resolution EPSG:4326 remap constrained by the same HydroLAKES polygon."""
    if scene_id is None:
        scene_id = (
            str(frame.get("scene_id", pd.Series([""])).iloc[0])
            if isinstance(frame, pd.DataFrame) and "scene_id" in frame.columns
            else None
        )
    if not scene_id:
        raise ValueError("scene_id is required for same-lake geographic remapping")
    g = _science_geographic_setup_v265(frame, float(resolution_m))
    rows, cols, src, dist = _native_assignments_geographic_v267(
        frame, g, scene_id, clip_geom=clip_geom
    )
    out = {c: np.full((g["ny"], g["nx"]), np.nan, np.float32) for c in columns}
    matched = src >= 0
    for c in columns:
        if c not in frame.columns or not matched.any():
            continue
        vals = pd.to_numeric(frame[c], errors="coerce").to_numpy(float)
        good = matched & np.isfinite(vals[np.maximum(src, 0)])
        if good.any():
            out[c][rows[good], cols[good]] = vals[src[good]].astype(np.float32)
    g.update(
        {
            "mapping_method": _m_settings.V267_GRID_METHOD,
            "target_lake_cells": int(len(rows)),
            "matched_target_cells": int(matched.sum()),
            "nearest_radius_m": max(
                _m_settings.V267_NEAREST_RADIUS_M, 2.0 * float(g["resolution_m"])
            ),
        }
    )
    return (out, g)


# Reference cell 50, lines 480-481.
def _regularize_dataframe(frame, columns, resolution_m):
    return _regularize_dataframe_geographic(frame, columns, resolution_m)


# Reference cell 50, lines 484-492.
def _auto_limits(values, metric):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if not v.size:
        return (0.0, 1.0)
    if metric == "CI_cyano":
        v = v[v > 0]
    if not v.size:
        return (0.0, 1.0)
    lo, hi = np.nanpercentile(v, [2, 98]) if v.size > 20 else (np.nanmin(v), np.nanmax(v))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo = float(np.nanmin(v))
        hi = float(np.nanmax(v))
    if hi <= lo:
        hi = lo + 1e-06
    return (float(lo), float(hi))


# Reference cell 50, lines 495-496.
def _optional_float(v):
    t = str(v or "").strip().lower()
    return None if not t or t == "auto" else float(t)


# Reference cell 50, lines 499-503.
def _rgba_data_url(grid, metric, cmap_name, reverse, vmin=None, vmax=None):
    vmin0, vmax0 = _auto_limits(grid, metric)
    vmin = vmin0 if vmin is None else float(vmin)
    vmax = vmax0 if vmax is None else float(vmax)
    if vmax <= vmin:
        raise ValueError("Color maximum must exceed minimum")
    cm = plt.get_cmap(cmap_name + ("_r" if reverse and (not cmap_name.endswith("_r")) else ""))
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    rgba = (cm(norm(np.nan_to_num(grid, nan=vmin))) * 255).astype(np.uint8)
    rgba[..., 3] = np.where(np.isfinite(grid), 230, 0).astype(np.uint8)
    im = Image.fromarray(rgba, "RGBA")
    b = io.BytesIO()
    im.save(b, format="PNG", optimize=True)
    return ("data:image/png;base64," + base64.b64encode(b.getvalue()).decode(), vmin, vmax, cm.name)


# Reference cell 50, lines 506-512.
def _truecolor_data_url(grids):
    chans = []
    for n in ("rhos_665", "rhos_560", "rhos_490"):
        a = np.asarray(grids[n], float)
        v = a[np.isfinite(a)]
        if not v.size:
            chans.append(np.zeros(a.shape, np.float32))
            continue
        lo, hi = np.nanpercentile(v, [2, 98])
        x = np.clip((a - lo) / max(hi - lo, 1e-08), 0, 1)
        x = np.power(x, 0.85)
        chans.append(x)
    rgb = np.stack(chans, -1)
    alpha = np.where(
        np.all(
            np.isfinite(
                np.stack([grids["rhos_665"], grids["rhos_560"], grids["rhos_490"]], axis=0)
            ),
            axis=0,
        ),
        255,
        0,
    ).astype(np.uint8)
    rgba = np.dstack([(np.nan_to_num(rgb) * 255).astype(np.uint8), alpha])
    im = Image.fromarray(rgba, "RGBA")
    b = io.BytesIO()
    im.save(b, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(b.getvalue()).decode()


# Reference cell 50, lines 515-519.
def _cmap_gradient_css(cmap_name, reverse=False, n=12):
    cm = plt.get_cmap(cmap_name + ("_r" if reverse and (not cmap_name.endswith("_r")) else ""))
    stops = []
    for i in range(n):
        r, g, b, _ = cm(i / (n - 1))
        stops.append(f"rgb({int(r * 255)},{int(g * 255)},{int(b * 255)}) {100 * i / (n - 1):.1f}%")
    return ",".join(stops)


# Reference cell 50, lines 522-528.
def _add_basemap(m, basemap, labels=True):
    if basemap == "None":
        return
    if basemap == "OpenStreetMap":
        folium.TileLayer(
            "https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png",
            attr="© OpenStreetMap contributors",
            name="OpenStreetMap",
            show=True,
        ).add_to(m)
        return
    if basemap == "Esri World Imagery":
        folium.TileLayer(
            "https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}",
            attr="Esri",
            name="Esri World Imagery",
            show=True,
        ).add_to(m)
        return
    dark = basemap == "Esri Dark Gray"
    stem = "Canvas/World_Dark_Gray_Base" if dark else "Canvas/World_Light_Gray_Base"
    ref = "Canvas/World_Dark_Gray_Reference" if dark else "Canvas/World_Light_Gray_Reference"
    folium.TileLayer(
        f"https://server.arcgisonline.com/ArcGIS/rest/services/{stem}/MapServer/tile/{{z}}/{{y}}/{{x}}",
        attr="Esri",
        name=basemap,
        show=True,
    ).add_to(m)
    if labels:
        folium.TileLayer(
            f"https://server.arcgisonline.com/ArcGIS/rest/services/{ref}/MapServer/tile/{{z}}/{{y}}/{{x}}",
            attr="Esri",
            name="Labels",
            overlay=True,
            control=True,
            show=True,
        ).add_to(m)


# Reference cell 50, lines 537-538.
def _add_scalar_legend(m, metric, cmap_name, reverse, vmin, vmax):
    grad = _cmap_gradient_css(cmap_name, reverse)
    legend = f"<div style='position:fixed;bottom:28px;left:24px;z-index:9999;background:rgba(255,255,255,.94);padding:10px 12px;border:1px solid #777;border-radius:5px;font:12px sans-serif;min-width:220px'><b>{html_lib.escape(metric)}</b><div style='height:14px;margin:5px 0 2px;background:linear-gradient(to right,{grad})'></div><div style='display:flex;justify-content:space-between'><span>{vmin:.5g}</span><span>{vmax:.5g}</span></div></div>"
    m.get_root().html.add_child(folium.Element(legend))


# Reference cell 50, lines 541-571.
def _hover_payload(frame, metric, resolution_m, max_points=180000):
    """Build a compact spatial-hash payload for portable HTML hover.

    The older implementation scanned every retained lake point on every mousemove,
    which became slow for a full OLCI scene.  This version buckets points in a
    lon/lat grid; JavaScript checks only the cursor bucket and its eight neighbors.
    """
    cols = ["latitude", "longitude", "Hylak_id"] + (
        [metric] if metric != "true_color" else ["rhos_665", "rhos_560", "rhos_490"]
    )
    f = frame[cols].copy()
    f = f[
        np.isfinite(pd.to_numeric(f["latitude"], errors="coerce"))
        & np.isfinite(pd.to_numeric(f["longitude"], errors="coerce"))
    ]
    if len(f) > max_points:
        take = np.linspace(0, len(f) - 1, max_points, dtype=int)
        f = f.iloc[take]
    cell = 0.02
    buckets = {}
    for r in f.itertuples(index=False):
        lat = float(r.latitude)
        lon = float(r.longitude)
        key = f"{math.floor(lon / cell)}:{math.floor(lat / cell)}"
        if metric == "true_color":
            rec = [
                round(lat, 6),
                round(lon, 6),
                int(r.Hylak_id),
                float(getattr(r, "rhos_665")) if np.isfinite(getattr(r, "rhos_665")) else None,
                float(getattr(r, "rhos_560")) if np.isfinite(getattr(r, "rhos_560")) else None,
                float(getattr(r, "rhos_490")) if np.isfinite(getattr(r, "rhos_490")) else None,
            ]
        else:
            val = getattr(r, metric)
            rec = [
                round(lat, 6),
                round(lon, 6),
                int(r.Hylak_id),
                float(val) if np.isfinite(val) else None,
            ]
        buckets.setdefault(key, []).append(rec)
    return {"cell": cell, "buckets": buckets, "metric": metric}


# Reference cell 50, lines 575-599.
def _inject_hover_js(m, frame, metric, resolution_m):
    payload = json.dumps(_hover_payload(frame, metric, resolution_m), separators=(",", ":"))
    map_name = m.get_name()
    radius = max(float(resolution_m) * 1.8, 450.0)
    js = f"<script>(function(){{\n    var P={payload}; var cell=P.cell, buckets=P.buckets, metric=P.metric; var map={map_name};\n    var info=L.control({{position:'topright'}});\n    info.onAdd=function(){{this._div=L.DomUtil.create('div','hover-readout');\n      this._div.style.cssText='background:rgba(255,255,255,.95);padding:8px 10px;border:1px solid #777;border-radius:4px;font:12px sans-serif;min-width:220px';\n      this._div.innerHTML='<b>Hover lake pixels</b><br>Move over a retained lake pixel'; return this._div;}}; info.addTo(map);\n    function hav(a,b,c,d){{var R=6371000,p=Math.PI/180,x=(d-b)*p*Math.cos((a+c)*p/2),y=(c-a)*p;return Math.sqrt(x*x+y*y)*R;}}\n    map.on('mousemove',function(e){{\n      var ix=Math.floor(e.latlng.lng/cell), iy=Math.floor(e.latlng.lat/cell), best=null, bd=1e99;\n      for(var dx=-1;dx<=1;dx++) for(var dy=-1;dy<=1;dy++){{\n        var arr=buckets[(ix+dx)+':'+(iy+dy)]; if(!arr) continue;\n        for(var j=0;j<arr.length;j++){{var p=arr[j],d=hav(e.latlng.lat,e.latlng.lng,p[0],p[1]);if(d<bd){{bd=d;best=p;}}}}\n      }}\n      if(!best||bd>{radius:.1f}){{info._div.innerHTML='<b>Hover lake pixels</b><br>No nearby retained pixel';return;}}\n      var txt='<b>Hylak_id '+best[2]+'</b><br>'+best[0].toFixed(5)+', '+best[1].toFixed(5)+'<br>';\n      if(metric==='true_color') txt+='rhos R/G/B: '+best[3]+', '+best[4]+', '+best[5];\n      else txt+=metric+': '+(best[3]===null?'NaN':Number(best[3]).toPrecision(6));\n      txt+='<br><small>nearest stored native OLCI lake pixel; distance '+Math.round(bd)+' m</small>'; info._div.innerHTML=txt;\n    }});\n    }})();</script>"
    m.get_root().html.add_child(folium.Element(js))


# Reference cell 59, lines 194-223.
def _build_map(
    scene_id,
    scope,
    lake_id,
    metric,
    resolution_m,
    basemap,
    labels,
    cmap,
    reverse,
    vmin_text,
    vmax_text,
    show_lakes=True,
):
    if not scene_id:
        raise ValueError("Choose a completed scene")
    scope = str(scope or "scene")
    planning_geom = None
    if scope == "lake":
        raw = _m_archive._load_scene_pixels(scene_id, "lake", lake_id)
        scope_label = f"lake {lake_id}"
    else:
        raw = _m_archive._load_scene_pixels(scene_id, "scene", None)
        if scope == "planning":
            planning_geom = _m_state.GLOBAL_STATE.get("planning_geometry")
            raw = _planning_geometry_frame(raw)
            scope_label = _m_state.GLOBAL_STATE.get("planning_shard") or "active planning region"
        else:
            scope_label = "entire processed scene"
    cols = ["rhos_665", "rhos_560", "rhos_490"] if metric == "true_color" else [metric]
    frame = _m_quality._masked_display_frame(raw, cols)
    grids, g = _regularize_dataframe_mercator(
        frame, cols, float(resolution_m), scene_id=scene_id, clip_geom=planning_geom
    )
    bounds = [[g["south"], g["west"]], [g["north"], g["east"]]]
    m = folium.Map(tiles=None, control_scale=True, prefer_canvas=True)
    _add_basemap(m, basemap, bool(labels))
    if metric == "true_color":
        folium.raster_layers.ImageOverlay(
            _truecolor_data_url(grids),
            bounds=bounds,
            opacity=0.95,
            name="True color · bloom-aware QA",
        ).add_to(m)
    else:
        img, vmin, vmax, cm = _rgba_data_url(
            grids[metric],
            metric,
            cmap,
            bool(reverse),
            _optional_float(vmin_text),
            _optional_float(vmax_text),
        )
        folium.raster_layers.ImageOverlay(
            img, bounds=bounds, opacity=0.88, name=f"{metric} · bloom-aware QA"
        ).add_to(m)
        _add_scalar_legend(m, metric, cmap, bool(reverse), vmin, vmax)
    if show_lakes:
        b = _m_archive._lake_boundaries_for_ids(raw.Hylak_id.unique())
        if not b.empty:
            folium.GeoJson(
                json.loads(b.to_crs(4326).to_json()),
                name="HydroLAKES boundaries",
                style_function=lambda f: {"color": "#00e5ff", "weight": 1.2, "fillOpacity": 0},
            ).add_to(m)
    if planning_geom is not None and (not planning_geom.is_empty):
        reg = gpd.GeoDataFrame({"name": [scope_label]}, geometry=[planning_geom], crs=4326)
        folium.GeoJson(
            json.loads(reg.to_json()),
            name="Active planning region",
            style_function=lambda f: {"color": "#ffcc00", "weight": 2.2, "fillOpacity": 0},
        ).add_to(m)
        west, south, east, north = planning_geom.bounds
        fit_bounds = [[south, west], [north, east]]
    else:
        fit_bounds = bounds
    folium_plugins.Fullscreen(position="topleft").add_to(m)
    folium_plugins.MousePosition(
        position="bottomright", separator=" | ", prefix="Lat / Lon", num_digits=5
    ).add_to(m)
    folium_plugins.MeasureControl(position="topleft", primary_length_unit="kilometers").add_to(m)
    folium.LayerControl(collapsed=False).add_to(m)
    m.fit_bounds(fit_bounds)
    hover_frame = frame.loc[_m_quality._primary_selector(frame)].copy()
    _inject_hover_js(m, hover_frame, metric, g["resolution_m"])
    out = (
        _m_settings.HTML_DIR
        / f"scene_{scene_id}_{scope}_{lake_id or 'ALL'}_{metric}_{int(round(g['resolution_m']))}m.html"
    )
    m.save(out)
    map_src = f"/gradio_api/file={quote(str(out), safe='/')}?v={out.stat().st_mtime_ns}"
    iframe = (
        "<iframe style='width:100%;height:740px;border:1px solid #ccd5df;border-radius:8px;background:white' loading='eager' src='"
        + map_src
        + "'></iframe>"
    )
    return (iframe, str(out), frame, grids, g)


# Reference cell 50, lines 791-810.
def _inject_all_science_hover_js(m, frame, resolution_m):
    cols = ["latitude", "longitude", "Hylak_id", "CI", "CI_cyano", "NDCI", "MPH", "FAI"]
    f = frame[[c for c in cols if c in frame.columns]].copy()
    f = f[
        np.isfinite(pd.to_numeric(f.latitude, errors="coerce"))
        & np.isfinite(pd.to_numeric(f.longitude, errors="coerce"))
    ]
    if len(f) > 180000:
        f = f.iloc[np.linspace(0, len(f) - 1, 180000, dtype=int)]
    cell = 0.02
    buckets = {}
    for r in f.itertuples(index=False):
        key = f"{math.floor(float(r.longitude) / cell)}:{math.floor(float(r.latitude) / cell)}"
        vals = []
        for name in ("CI", "CI_cyano", "NDCI", "MPH", "FAI"):
            v = getattr(r, name, np.nan)
            vals.append(float(v) if np.isfinite(v) else None)
        buckets.setdefault(key, []).append(
            [round(float(r.latitude), 6), round(float(r.longitude), 6), int(r.Hylak_id), *vals]
        )
    payload = json.dumps({"cell": cell, "buckets": buckets}, separators=(",", ":"))
    map_name = m.get_name()
    radius = max(float(resolution_m) * 2.0, 650.0)
    js = f"<script>(function(){{var P={payload},cell=P.cell,buckets=P.buckets,map={map_name};\n    var info=L.control({{position:'topright'}});info.onAdd=function(){{this._div=L.DomUtil.create('div','all-science-hover');this._div.style.cssText='background:rgba(255,255,255,.96);padding:8px 10px;border:1px solid #777;border-radius:4px;font:12px sans-serif;min-width:250px';this._div.innerHTML='<b>Hover native lake pixels</b>';return this._div;}};info.addTo(map);\n    function hav(a,b,c,d){{var R=6371000,p=Math.PI/180,x=(d-b)*p*Math.cos((a+c)*p/2),y=(c-a)*p;return Math.sqrt(x*x+y*y)*R;}}\n    map.on('mousemove',function(e){{var ix=Math.floor(e.latlng.lng/cell),iy=Math.floor(e.latlng.lat/cell),best=null,bd=1e99;for(var dx=-1;dx<=1;dx++)for(var dy=-1;dy<=1;dy++){{var a=buckets[(ix+dx)+':'+(iy+dy)];if(!a)continue;for(var j=0;j<a.length;j++){{var p=a[j],d=hav(e.latlng.lat,e.latlng.lng,p[0],p[1]);if(d<bd){{bd=d;best=p;}}}}}}if(!best||bd>{radius:.1f}){{info._div.innerHTML='<b>Hover native lake pixels</b><br>No nearby retained pixel';return;}}var n=['CI','CIcyano','NDCI','MPH','FAI/AFAI'],txt='<b>Hylak_id '+best[2]+'</b><br>'+best[0].toFixed(5)+', '+best[1].toFixed(5);for(var k=0;k<5;k++)txt+='<br>'+n[k]+': '+(best[3+k]===null?'NaN':Number(best[3+k]).toPrecision(6));txt+='<br><small>nearest stored native OLCI lake pixel · '+Math.round(bd)+' m</small>';info._div.innerHTML=txt;}});}})();</script>"
    m.get_root().html.add_child(folium.Element(js))


# Reference cell 59, lines 338-354.
def create_complete_scene_html(scene_id, resolution_m=1000, basemap="Esri Dark Gray", force=False):
    """Portable multi-layer HTML using the efficient same-lake mapper."""
    if not scene_id:
        raise ValueError("Choose a completed scene")
    frame = _m_archive._load_scene_pixels(scene_id, "scene", None)
    cols = ["rhos_665", "rhos_560", "rhos_490", "CI", "CI_cyano", "NDCI", "MPH", "FAI"]
    display = _m_quality._masked_display_frame(frame, cols)
    grids, g = _regularize_dataframe_mercator(
        display, cols, float(resolution_m), scene_id=scene_id, clip_geom=None
    )
    out = (
        _m_settings.HTML_DIR / f"scene_{scene_id}_ALL_SCIENCE_{int(round(g['resolution_m']))}m.html"
    )
    if out.exists() and (not force):
        return (str(out), g)
    bounds = [[g["south"], g["west"]], [g["north"], g["east"]]]
    m = folium.Map(tiles=None, control_scale=True, prefer_canvas=True)
    _add_basemap(m, basemap, True)
    tc = _truecolor_data_url({k: grids[k] for k in ("rhos_665", "rhos_560", "rhos_490")})
    folium.raster_layers.ImageOverlay(
        tc, bounds=bounds, opacity=0.95, name="True color rhos 665/560/490", show=False
    ).add_to(m)
    legend_rows = []
    palette = {
        "CI": "cividis",
        "CI_cyano": "cividis",
        "NDCI": "BrBG",
        "MPH": "plasma",
        "FAI": "viridis",
    }
    for metric in ("CI", "CI_cyano", "NDCI", "MPH", "FAI"):
        img, vmin, vmax, _ = _rgba_data_url(
            grids[metric], metric, palette[metric], False, None, None
        )
        folium.raster_layers.ImageOverlay(
            img, bounds=bounds, opacity=0.88, name=metric, show=metric == "CI_cyano"
        ).add_to(m)
        legend_rows.append(
            f"<tr><td><b>{html_lib.escape(metric)}</b></td><td>{vmin:.4g}</td><td>{vmax:.4g}</td></tr>"
        )
    b = _m_archive._lake_boundaries_for_ids(frame.Hylak_id.unique())
    if not b.empty:
        folium.GeoJson(
            json.loads(b.to_crs(4326).to_json()),
            name="HydroLAKES boundaries",
            style_function=lambda f: {"color": "#00e5ff", "weight": 1.1, "fillOpacity": 0},
        ).add_to(m)
    panel = (
        "<div style='position:fixed;bottom:28px;left:24px;z-index:9999;background:rgba(255,255,255,.95);padding:8px 10px;border:1px solid #777;border-radius:5px;font:11px sans-serif'><b>Auto display ranges</b><table><tr><th>Layer</th><th>min</th><th>max</th></tr>"
        + "".join(legend_rows)
        + "</table><small>Use layer control at upper right.</small></div>"
    )
    m.get_root().html.add_child(folium.Element(panel))
    folium_plugins.Fullscreen(position="topleft").add_to(m)
    folium_plugins.MousePosition(
        position="bottomright", separator=" | ", prefix="Lat / Lon", num_digits=5
    ).add_to(m)
    folium_plugins.MeasureControl(position="topleft", primary_length_unit="kilometers").add_to(m)
    folium.LayerControl(collapsed=False).add_to(m)
    m.fit_bounds(bounds)
    _inject_all_science_hover_js(m, display, g["resolution_m"])
    m.save(out)
    return (str(out), g)


# Reference cell 58, lines 281-293.
def map_callback(
    scene_id,
    scope,
    lake_id,
    metric,
    resolution_m,
    basemap,
    labels,
    cmap,
    reverse,
    vmin,
    vmax,
    show_lakes,
):
    from . import exports as _m_exports

    try:
        iframe, html_path, frame, grids, g = _build_map(
            scene_id,
            scope,
            lake_id,
            metric,
            resolution_m,
            basemap,
            labels,
            cmap,
            reverse,
            vmin,
            vmax,
            show_lakes,
        )
        quick_scope = scope if scope in {"scene", "planning", "lake"} else "scene"
        nc = _m_exports.export_snap_gridded_netcdf(
            scene_id, quick_scope, lake_id, metric, frame, resolution_m
        )
        compact = _m_archive._compact_download(scene_id)
        scope_text = (
            _m_state.GLOBAL_STATE.get("planning_shard", "active planning region")
            if scope == "planning"
            else f"lake {lake_id}"
            if scope == "lake"
            else "full scene"
        )
        return (
            iframe,
            html_path,
            nc,
            compact,
            _m_ui_helpers._app_status(
                "Map ready",
                f"{len(frame):,} stored lake pixels; {frame.Hylak_id.nunique():,} lakes; {scope_text}; exact Web-Mercator display grid = {g['resolution_m']:.0f} m.",
            ),
        )
    except Exception as exc:
        return (
            _m_ui_helpers._app_status("Map failed", str(exc), False),
            None,
            None,
            None,
            _m_ui_helpers._app_status("Map failed", str(exc), False),
        )


# Reference cell 58, lines 53-66.
def _planning_geometry_frame(frame):
    geom = _m_state.GLOBAL_STATE.get("planning_geometry")
    if geom is None or getattr(geom, "is_empty", True):
        raise ValueError(
            "No active planning region. Plan a State/Province, bbox, or 5° grid in Tab 1 first."
        )
    lon = pd.to_numeric(frame.longitude, errors="coerce").to_numpy(float)
    lat = pd.to_numeric(frame.latitude, errors="coerce").to_numpy(float)
    finite = np.isfinite(lon) & np.isfinite(lat)
    keep = np.zeros(len(frame), dtype=bool)
    if finite.any():
        keep[finite] = contains_xy(geom, lon[finite], lat[finite])
    out = frame.loc[keep].reset_index(drop=True)
    if out.empty:
        raise ValueError(
            "The selected processed scene contains no stored lake pixels inside the active planning region."
        )
    return out


# Reference cell 58, lines 296-321.
def _science_geographic_setup_v265(
    frame, requested_resolution_m=_m_settings.SCIENCE_EXPORT_DEFAULT_RESOLUTION_M
):
    """SNAP-friendly regular EPSG:4326 grid with nominal metric spacing at scene median latitude."""
    lat = pd.to_numeric(frame.latitude, errors="coerce").to_numpy(float)
    lon = pd.to_numeric(frame.longitude, errors="coerce").to_numpy(float)
    finite = np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90) & (np.abs(lon) <= 180)
    if not finite.any():
        raise ValueError("No finite stored scene geolocation")
    res = float(requested_resolution_m)
    ref_lat = float(np.nanmedian(lat[finite]))
    dlat = res / 111320.0
    dlon = res / max(111320.0 * math.cos(math.radians(ref_lat)), 1000.0)
    west = math.floor(float(np.nanmin(lon[finite])) / dlon) * dlon - dlon
    east = math.ceil(float(np.nanmax(lon[finite])) / dlon) * dlon + dlon
    south = math.floor(float(np.nanmin(lat[finite])) / dlat) * dlat - dlat
    north = math.ceil(float(np.nanmax(lat[finite])) / dlat) * dlat + dlat
    nx = max(2, int(round((east - west) / dlon)))
    ny = max(2, int(round((north - south) / dlat)))
    east = west + nx * dlon
    north = south + ny * dlat
    longitude = west + (np.arange(nx, dtype=np.float64) + 0.5) * dlon
    latitude = south + (np.arange(ny, dtype=np.float64) + 0.5) * dlat
    return {
        "resolution_m": res,
        "reference_latitude": ref_lat,
        "dlat": dlat,
        "dlon": dlon,
        "west": west,
        "east": east,
        "south": south,
        "north": north,
        "nx": nx,
        "ny": ny,
        "ncell": int(nx) * int(ny),
        "longitude": longitude,
        "latitude": latitude,
        "crs": "EPSG:4326",
    }


# Reference cell 58, lines 324-332.
def _scene_lake_union(scene_id, frame):
    try:
        b = _m_archive._lake_boundaries_for_ids(
            pd.to_numeric(frame.Hylak_id, errors="coerce").dropna().astype(np.int64).unique()
        )
        if b is None or b.empty:
            return None
        b = b.to_crs(4326)
        return b.geometry.union_all()
    except Exception:
        return None


# Reference cell 58, lines 335-353.
def _candidate_geographic_chunks(frame, g, chunk=_m_settings.V265_EXPORT_CHUNK):
    lat = pd.to_numeric(frame.latitude, errors="coerce").to_numpy(float)
    lon = pd.to_numeric(frame.longitude, errors="coerce").to_numpy(float)
    finite = np.isfinite(lat) & np.isfinite(lon)
    rr = np.floor((lat[finite] - g["south"]) / g["dlat"]).astype(np.int64)
    cc = np.floor((lon[finite] - g["west"]) / g["dlon"]).astype(np.int64)
    good = (rr >= 0) & (rr < g["ny"]) & (cc >= 0) & (cc < g["nx"])
    rr = rr[good]
    cc = cc[good]
    nxc = int(math.ceil(g["nx"] / chunk))
    nyc = int(math.ceil(g["ny"] / chunk))
    base = np.unique(rr // chunk * nxc + cc // chunk)
    keys = set()
    for cid in base.tolist():
        cr = int(cid // nxc)
        ccol = int(cid % nxc)
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                ar, ac = (cr + dr, ccol + dc)
                if 0 <= ar < nyc and 0 <= ac < nxc:
                    keys.add((ar, ac))
    return sorted(keys)


# Reference cell 59, lines 15-43.
def _scene_analysis_lakes_v267(scene_id, frame, clip_geom=None):
    """Return the exact scene lake polygons used as mapped support.

    The shoreline erosion stored in the compact-product metadata is reapplied so the
    mapped product does not expand back into shoreline pixels intentionally removed
    during native-pixel extraction.
    """
    ids = (
        pd.to_numeric(frame.get("Hylak_id", pd.Series(dtype=float)), errors="coerce")
        .dropna()
        .astype(np.int64)
        .unique()
    )
    if not len(ids):
        raise ValueError("No Hylak_id values are available in the selected scene/scope")
    lakes = _m_archive._lake_boundaries_for_ids(ids)
    if lakes is None or lakes.empty:
        raise RuntimeError(
            "HydroLAKES polygons for the stored scene lakes are unavailable; reload the target HydroLAKES universe before mapping/export."
        )
    lakes = lakes.to_crs(4326).copy()
    qa = _m_archive._compact_qa_metadata(scene_id)
    shore = float(qa.get("shore_buffer_m", _m_settings.DEFAULT_SHORE_BUFFER_M) or 0.0)
    if shore > 0:
        lakes = _m_hydrolakes.erode_lake_geometries(lakes, shore)
    if clip_geom is not None and (not getattr(clip_geom, "is_empty", True)):
        cg = clip_geom
        try:
            lakes["geometry"] = lakes.geometry.intersection(cg)
            lakes = lakes[~lakes.geometry.is_empty].copy()
        except Exception:
            pass
    if lakes.empty:
        raise ValueError("No mapped lake polygons remain in the selected display/export scope")
    return lakes


# Reference cell 59, lines 46-59.
def _dedupe_target_assignments_v267(rows, cols, src, dist, nx):
    """Keep one source per target cell, preferring a real match and then shortest distance."""
    rows = np.asarray(rows, np.int64)
    cols = np.asarray(cols, np.int64)
    src = np.asarray(src, np.int64)
    dist = np.asarray(dist, float)
    if not len(rows):
        return (rows, cols, src, dist)
    flat = rows * int(nx) + cols
    rank = np.where(src >= 0, dist, np.inf)
    order = np.lexsort((rank, flat))
    flat_o = flat[order]
    keep = np.r_[True, flat_o[1:] != flat_o[:-1]]
    take = order[keep]
    return (rows[take], cols[take], src[take], dist[take])


# Reference cell 59, lines 62-109.
def _native_assignments_geographic_v267(frame, g, scene_id, clip_geom=None, radius_m=None):
    """Generate target lake cells and same-lake nearest-native source assignments.

    A target cell is never allowed to borrow an observation from a different Hylak_id.
    The nearest source is chosen before QA is inspected. Thus a cloud/invalid native
    observation remains no-data rather than being replaced by a farther clear pixel.
    """
    radius_m = max(
        float(radius_m or _m_settings.V267_NEAREST_RADIUS_M), 2.0 * float(g["resolution_m"])
    )
    lat = pd.to_numeric(frame.latitude, errors="coerce").to_numpy(float)
    lon = pd.to_numeric(frame.longitude, errors="coerce").to_numpy(float)
    hid = pd.to_numeric(frame.Hylak_id, errors="coerce").fillna(-1).to_numpy(np.int64)
    finite = np.isfinite(lat) & np.isfinite(lon)
    lakes = _scene_analysis_lakes_v267(scene_id, frame, clip_geom=clip_geom)
    all_r = []
    all_c = []
    all_s = []
    all_d = []
    lats_grid = np.asarray(g["latitude"], float)
    lons_grid = np.asarray(g["longitude"], float)
    for _, lr in lakes.iterrows():
        lid = int(lr.Hylak_id)
        geom = lr.geometry
        src_idx = np.flatnonzero(finite & (hid == lid))
        if not len(src_idx) or geom is None or geom.is_empty:
            continue
        tree = cKDTree(_m_science._unit_sphere_xyz(lat[src_idx], lon[src_idx]))
        west, south, east, north = map(float, geom.bounds)
        c0 = max(0, int(np.searchsorted(lons_grid, west, side="left") - 1))
        c1 = min(g["nx"], int(np.searchsorted(lons_grid, east, side="right") + 1))
        r0 = max(0, int(np.searchsorted(lats_grid, south, side="left") - 1))
        r1 = min(g["ny"], int(np.searchsorted(lats_grid, north, side="right") + 1))
        if c1 <= c0 or r1 <= r0:
            continue
        width = max(1, c1 - c0)
        block_rows = max(1, min(r1 - r0, int(max(1, 350000 // width))))
        for rb0 in range(r0, r1, block_rows):
            rb1 = min(rb0 + block_rows, r1)
            lon2, lat2 = np.meshgrid(lons_grid[c0:c1], lats_grid[rb0:rb1])
            support = contains_xy(geom, lon2, lat2)
            pos = np.flatnonzero(support.ravel())
            if not len(pos):
                continue
            rr_local, cc_local = np.unravel_index(pos, support.shape)
            rr = rr_local.astype(np.int64) + rb0
            cc = cc_local.astype(np.int64) + c0
            xyz = _m_science._unit_sphere_xyz(lat2.ravel()[pos], lon2.ravel()[pos])
            chord, ii = tree.query(xyz, k=1, workers=-1)
            dist = 6371000.0 * (2.0 * np.arcsin(np.minimum(1.0, np.asarray(chord, float) / 2.0)))
            matched = np.isfinite(dist) & (dist <= radius_m) & (np.asarray(ii) < len(src_idx))
            src = np.full(len(pos), -1, np.int64)
            dd = np.full(len(pos), np.nan, float)
            if matched.any():
                src[matched] = src_idx[np.asarray(ii[matched], np.int64)]
                dd[matched] = dist[matched]
            all_r.append(rr)
            all_c.append(cc)
            all_s.append(src)
            all_d.append(dd)
    if not all_r:
        return (np.empty(0, np.int64),) * 3 + (np.empty(0, float),)
    rows = np.concatenate(all_r)
    cols = np.concatenate(all_c)
    src = np.concatenate(all_s)
    dist = np.concatenate(all_d)
    return _dedupe_target_assignments_v267(rows, cols, src, dist, g["nx"])


# Reference cell 59, lines 133-174.
def _native_assignments_mercator_v267(frame, g, scene_id, clip_geom=None, radius_m=None):
    """Same-lake nearest-native assignments for the Leaflet/Web-Mercator display grid."""
    radius_ground = max(
        float(radius_m or _m_settings.V267_NEAREST_RADIUS_M), 2.0 * float(g["resolution_m"])
    )
    lat = pd.to_numeric(frame.latitude, errors="coerce").to_numpy(float)
    lon = pd.to_numeric(frame.longitude, errors="coerce").to_numpy(float)
    hid = pd.to_numeric(frame.Hylak_id, errors="coerce").fillna(-1).to_numpy(np.int64)
    finite = np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) < 85)
    x = np.full(len(frame), np.nan, float)
    y = np.full(len(frame), np.nan, float)
    if finite.any():
        xx, yy = _TO_WEBMERC.transform(lon[finite], lat[finite])
        x[finite] = xx
        y[finite] = yy
    lakes = _scene_analysis_lakes_v267(scene_id, frame, clip_geom=clip_geom).to_crs(3857)
    all_r = []
    all_c = []
    all_s = []
    all_d = []
    xs = g["xmin"] + (np.arange(g["nx"], dtype=float) + 0.5) * g["resolution_m"]
    for _, lr in lakes.iterrows():
        lid = int(lr.Hylak_id)
        geom = lr.geometry
        src_idx = np.flatnonzero(finite & (hid == lid))
        if not len(src_idx) or geom is None or geom.is_empty:
            continue
        tree = cKDTree(np.column_stack([x[src_idx], y[src_idx]]))
        lat_med = float(np.nanmedian(lat[src_idx]))
        scale = 1.0 / max(math.cos(math.radians(lat_med)), 0.2)
        radius_web = radius_ground * scale
        minx, miny, maxx, maxy = map(float, geom.bounds)
        c0 = max(0, int(math.floor((minx - g["xmin"]) / g["resolution_m"])) - 1)
        c1 = min(g["nx"], int(math.ceil((maxx - g["xmin"]) / g["resolution_m"])) + 1)
        r0 = max(0, int(math.floor((g["ymax"] - maxy) / g["resolution_m"])) - 1)
        r1 = min(g["ny"], int(math.ceil((g["ymax"] - miny) / g["resolution_m"])) + 1)
        if c1 <= c0 or r1 <= r0:
            continue
        width = max(1, c1 - c0)
        block_rows = max(1, min(r1 - r0, int(max(1, 350000 // width))))
        for rb0 in range(r0, r1, block_rows):
            rb1 = min(rb0 + block_rows, r1)
            ys = g["ymax"] - (np.arange(rb0, rb1, dtype=float) + 0.5) * g["resolution_m"]
            xx2, yy2 = np.meshgrid(xs[c0:c1], ys)
            support = contains_xy(geom, xx2, yy2)
            pos = np.flatnonzero(support.ravel())
            if not len(pos):
                continue
            rr_local, cc_local = np.unravel_index(pos, support.shape)
            rr = rr_local.astype(np.int64) + rb0
            cc = cc_local.astype(np.int64) + c0
            dist, ii = tree.query(
                np.column_stack([xx2.ravel()[pos], yy2.ravel()[pos]]), k=1, workers=-1
            )
            matched = np.isfinite(dist) & (dist <= radius_web) & (np.asarray(ii) < len(src_idx))
            src = np.full(len(pos), -1, np.int64)
            dd = np.full(len(pos), np.nan, float)
            if matched.any():
                src[matched] = src_idx[np.asarray(ii[matched], np.int64)]
                dd[matched] = np.asarray(dist, float)[matched] / scale
            all_r.append(rr)
            all_c.append(cc)
            all_s.append(src)
            all_d.append(dd)
    if not all_r:
        return (np.empty(0, np.int64),) * 3 + (np.empty(0, float),)
    return _dedupe_target_assignments_v267(
        np.concatenate(all_r),
        np.concatenate(all_c),
        np.concatenate(all_s),
        np.concatenate(all_d),
        g["nx"],
    )
