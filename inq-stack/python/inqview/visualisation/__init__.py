"""inqview.visualisation — the rendering layer (ADR 0003).

The ONLY sub-package permitted to import plotting/VTK libraries. Holds the
canonical theme (ADR 0004) and, as the restructure proceeds, all renderers that
consume ``inqview.analysis`` result dataclasses.
"""
from __future__ import annotations

from .style import (
    ONE_COL_IN,
    TWO_COL_W_IN,
    ROLE_CMAP,
    STYLE_CONFIG,
    apply_theme,
    cmap_for,
    figure_one_col,
    figure_two_col,
)
from .energy_components import (
    render_breakdown_gif,
    render_flow_lines,
    render_initial_vs_final_bars,
)
from .field_io import VtiField, load_vti
from .density_gifs import make_density_gif_battery, make_twin_density_matrix

__all__ = [
    "VtiField",
    "load_vti",
    "make_density_gif_battery",
    "make_twin_density_matrix",
    "ONE_COL_IN",
    "TWO_COL_W_IN",
    "ROLE_CMAP",
    "STYLE_CONFIG",
    "apply_theme",
    "cmap_for",
    "figure_one_col",
    "figure_two_col",
    "render_initial_vs_final_bars",
    "render_flow_lines",
    "render_breakdown_gif",
]
