"""Report 1 (MSc dissertation) figure-generation scripts.

Each fig0N_<slug>.py module is a runnable script that loads simulation
data, applies the shared matplotlib style via apply_style(), and writes
docs/reports/report1/figures/fig0N_<slug>.png.

Style and palette definitions live in _shared_style; the human-readable
specification is docs/reports/report1/figures/global_style.md.
"""

from ._shared_style import (
    apply_style,
    palette_sweep5,
    palette_sweep3,
    palette_regime3,
    regime_tints,
    references,
    column_widths_in,
    panel_label,
    TufteCritic,
)

__all__ = [
    "apply_style",
    "palette_sweep5",
    "palette_sweep3",
    "palette_regime3",
    "regime_tints",
    "references",
    "column_widths_in",
    "panel_label",
    "TufteCritic",
]
