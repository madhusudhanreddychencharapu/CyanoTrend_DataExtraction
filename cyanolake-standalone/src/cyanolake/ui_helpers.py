"""Standalone ui helpers implementation."""

from __future__ import annotations

import html as html_lib

from . import catalogue as _m_catalogue


# Reference cell 50, lines 60-62.
def _app_status(title, detail="", ok=True):
    color = "#157347" if ok else "#b42318"
    return f"<div style='padding:10px 12px;border-left:5px solid {color};background:#f8fafc;border-radius:6px'><b style='color:{color}'>{html_lib.escape(str(title))}</b><br>{html_lib.escape(str(detail)).replace(chr(10), '<br>')}</div>"


# Reference cell 50, lines 220-225.
def _parse_bbox(text):
    text = str(text or "").strip()
    if not text:
        return None
    parts = [float(x.strip()) for x in text.split(",")]
    if len(parts) != 4:
        raise ValueError("BBox must be west,south,east,north")
    return _m_catalogue.normalize_bbox(parts)
