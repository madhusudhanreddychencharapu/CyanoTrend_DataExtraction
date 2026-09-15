"""Standalone quality implementation."""

from __future__ import annotations

import numpy as np
import pandas as pd


# Reference cell 57, lines 355-358.
def _primary_selector(frame):
    if "valid_water_mask" not in frame.columns:
        return np.ones(len(frame), dtype=bool)
    return pd.to_numeric(frame["valid_water_mask"], errors="coerce").fillna(0).to_numpy(float) > 0


# Reference cell 57, lines 361-364.
def _cyan_comparison_selector(frame):
    if (
        "cyan_strict_valid_mask" in frame.columns
        and pd.to_numeric(frame["cyan_strict_valid_mask"], errors="coerce").notna().any()
    ):
        return (
            pd.to_numeric(frame["cyan_strict_valid_mask"], errors="coerce")
            .fillna(0)
            .to_numpy(float)
            > 0
        )
    return _primary_selector(frame)


# Reference cell 57, lines 367-373.
def _masked_display_frame(frame, columns):
    out = frame.copy()
    keep = _primary_selector(out)
    for c in columns:
        if c in out.columns:
            a = pd.to_numeric(out[c], errors="coerce").to_numpy(float)
            a[~keep] = np.nan
            out[c] = a
    return out
