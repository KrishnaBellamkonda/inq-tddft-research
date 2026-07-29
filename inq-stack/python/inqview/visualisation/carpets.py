"""Carpet (waterfall/contour) renderers (IV-M05, Runfeng's idea).

Renders a (time × bin) array as a filled contour — e.g. the WP momentum
distribution P(|k|, t) as a (k, t) carpet. Renderers are not unit-tested
(IV-M10); they consume already-validated arrays.
"""
from __future__ import annotations

import numpy as np

from .style import apply_theme, cmap_for, figure_two_col


def momentum_carpet(times_au, k_vals, grid, *, log: bool = False, levels: int = 40):
    """Filled-contour carpet of P(|k|, t).

    times_au : (n_t,) · k_vals : (n_k,) · grid : (n_t, n_k) intensity.
    Returns (fig, ax). x = |k| (Bohr⁻¹), y = time (a.u.), colour = intensity.
    """
    apply_theme()
    g = np.asarray(grid, dtype=float)
    z = np.log1p(g) if log else g
    fig, ax = figure_two_col(height_in=4.0)
    cs = ax.contourf(np.asarray(k_vals), np.asarray(times_au), z,
                     levels=levels, cmap=cmap_for("sequential"))
    fig.colorbar(cs, ax=ax, label=("log(1+P)" if log else "P(|k|)"))
    ax.set_xlabel(r"$|k|$ (bohr$^{-1}$)")
    ax.set_ylabel("time (a.u.)")
    return fig, ax
