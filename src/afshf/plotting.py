"""Plotting helpers (spec section 14).

Skeleton: figure generation lives behind a feature flag so the test suite can
run without matplotlib being installed yet.
"""

from __future__ import annotations

from typing import Iterable

import pandas as pd


def quadrant_data(df: pd.DataFrame,
                  x_col: str,
                  y_col: str,
                  target_col: str = "target_id") -> pd.DataFrame:
    """Project per-target effects onto a 2D quadrant plane (e.g., AF vs HF)."""
    out = df[[target_col, x_col, y_col]].copy()
    out["quadrant"] = [
        ("benefit_benefit" if x < 0 and y < 0 else
         "benefit_harm" if x < 0 and y > 0 else
         "harm_benefit" if x > 0 and y < 0 else
         "harm_harm")
        for x, y in zip(out[x_col], out[y_col])
    ]
    return out


def lollipop_table(df: pd.DataFrame,
                   target_col: str = "target_id",
                   value_col: str = "theta",
                   ci_lower_col: str | None = None,
                   ci_upper_col: str | None = None) -> pd.DataFrame:
    cols = [target_col, value_col]
    if ci_lower_col:
        cols.append(ci_lower_col)
    if ci_upper_col:
        cols.append(ci_upper_col)
    return df[cols].copy()


def heatmap_table(df: pd.DataFrame,
                  index_col: str,
                  column_col: str,
                  value_col: str) -> pd.DataFrame:
    return df.pivot_table(index=index_col, columns=column_col, values=value_col, aggfunc="first")


def render_figures_if_available(figure_specs: Iterable[dict]) -> list[str]:
    """Render figures only if matplotlib is importable. Returns paths produced."""
    try:
        import matplotlib  # noqa: F401
    except ImportError:
        return []
    return []  # actual rendering implemented in Phase 6
