"""Canonical inqview plotting theme (ADR 0004).

Promoted from ``report1/_shared_style.py`` + the ``report-figures`` skill +
``docs/reports/report1/figures/global_style.md`` into the library. The single
visual standard for ALL inqview output:

- **semantic cmap roles** — phases ask for a ROLE, never a literal cmap:
  ``sequential → inferno``, ``diverging → RdBu_r`` (zero-centred),
  ``phase → twilight_shifted``.
- **fixed-dimension figure factory** — ``figure_one_col()`` = 3.5×3.0 in with a
  FIXED axes rectangle (every one-column panel shares an identical data box),
  ``figure_two_col()`` = 7.0 in wide. Individual plots only; panel composition
  is a downstream LaTeX concern.

Known pitfall (memory ``reference_fixed_dimension_plot_pitfalls``): tight-bbox /
constrained-layout silently break the fixed width. The factory therefore uses an
explicit ``add_axes`` rectangle and NO constrained layout; save with
``bbox_inches=None`` to preserve the on-page size.

Report vs presentation saving (memory ``feedback`` FB-001, 2026-06-25):
- **Report** figures have NO on-canvas title and the colorbar fits inside the
  fixed rect — save with ``bbox_inches=None`` (the fixed-dimension idiom).
- **Presentation** figures DO carry an on-canvas ``set_title`` and often an
  external (``make_axes_locatable``) colorbar whose label spills past the fixed
  rect. The fixed rect leaves ~3.5 % headroom above the axes, so a bare
  ``savefig`` clips the title/colorbar off-canvas. Such figures MUST be saved
  via :func:`save_presentation` (``bbox_inches="tight"``), which sacrifices the
  fixed on-page size — acceptable because a slide deck contain-fits each image.
  Never save a titled / external-colorbar figure with the bare report save.
"""
from __future__ import annotations

import math

import matplotlib as mpl
import matplotlib.pyplot as plt

# --- semantic cmap roles (the designed standard) ---------------------------
ROLE_CMAP: dict[str, str] = {
    "sequential": "inferno",          # positive-magnitude maps, diffraction intensity
    "diverging": "RdBu_r",            # signed Δn / difference maps (zero-centred)
    "phase": "twilight_shifted",      # cyclic data
}


def cmap_for(role: str) -> str:
    """Return the canonical cmap name for a semantic role (never a literal)."""
    try:
        return ROLE_CMAP[role]
    except KeyError:
        raise ValueError(
            f"unknown cmap role {role!r}; valid: {sorted(ROLE_CMAP)}"
        ) from None


# --- canonical axis units (the single units standard; ADR 0004, Cluster R) --
# Promoted from the tufte skill's "project annotation rule 5". Every inqview
# axis uses these — NEVER atomic units on a public time axis.
UNITS: dict[str, str] = {
    "energy": "eV",
    "length": "Bohr",
    "time": "fs",                 # femtoseconds, not atomic units
    "momentum": r"Bohr$^{-1}$",
    "stopping_power": "eV/Bohr",
}


def axis_label(quantity: str, symbol: str | None = None) -> str:
    """Canonical axis label ``"<symbol or quantity> (<unit>)"`` for a known
    quantity (e.g. ``axis_label("time")`` -> ``"time (fs)"``)."""
    try:
        unit = UNITS[quantity]
    except KeyError:
        raise ValueError(
            f"unknown quantity {quantity!r}; valid: {sorted(UNITS)}"
        ) from None
    name = symbol if symbol is not None else quantity.replace("_", " ")
    return f"{name} ({unit})"


# --- fixed-dimension scheme (report1 STANDARD FIXED-COLUMN, 2026-05-29) -----
ONE_COL_IN: tuple[float, float] = (3.5, 3.0)   # (width, height) one-column plot
TWO_COL_W_IN: float = 7.0                      # width of a two-column plot

# FIXED axes rectangle (figure fractions) so every one-column plot shares an
# identical data box and panels align. Margins sized for the worst-case y-label.
_ONE_COL_AXES_RECT: tuple[float, float, float, float] = (0.180, 0.160, 0.785, 0.805)


# --- rcParams standard ------------------------------------------------------
STYLE_CONFIG: dict[str, object] = {
    "font_family": "serif",
    "font_size": 10,
    "axes_label": 10,
    "axes_title": 10,
    "tick_label": 9,
    "legend": 9,
    "axes_linewidth": 0.8,
    "line_width": 1.2,
    "tick_direction": "in",
    "tick_top": True,
    "tick_right": True,
    "tick_major_size": 3.5,
    "tick_major_width": 0.6,
    "cmap_sequential": ROLE_CMAP["sequential"],
    "cmap_diverging": ROLE_CMAP["diverging"],
    "save_dpi": 600,
    "preview_dpi": 150,
}


def apply_theme() -> None:
    """Install the canonical rcParams. Call once at the top of a plotting run."""
    c = STYLE_CONFIG
    mpl.rcParams.update({
        "font.family": c["font_family"],
        "font.size": c["font_size"],
        "axes.labelsize": c["axes_label"],
        "axes.titlesize": c["axes_title"],
        "xtick.labelsize": c["tick_label"],
        "ytick.labelsize": c["tick_label"],
        "legend.fontsize": c["legend"],
        "axes.linewidth": c["axes_linewidth"],
        "lines.linewidth": c["line_width"],
        "xtick.direction": c["tick_direction"],
        "ytick.direction": c["tick_direction"],
        "xtick.top": c["tick_top"],
        "ytick.right": c["tick_right"],
        "xtick.major.size": c["tick_major_size"],
        "ytick.major.size": c["tick_major_size"],
        "xtick.major.width": c["tick_major_width"],
        "ytick.major.width": c["tick_major_width"],
        "image.cmap": c["cmap_sequential"],
        "savefig.dpi": c["save_dpi"],
        "figure.dpi": c["preview_dpi"],
        # Scientific notation must render as ×10^n (mathtext superscript), NEVER
        # as "1e-03" on any axis/colorbar offset or tick (user rule FB-003).
        "axes.formatter.use_mathtext": True,
    })


def figure_one_col(*, with_colorbar: bool = False):
    """A one-column figure (3.5×3.0 in) with the FIXED axes rectangle.

    Returns ``(fig, ax)``. With ``with_colorbar`` the axes width shrinks to leave
    room for a colourbar while the figure size stays fixed.
    """
    fig = plt.figure(figsize=ONE_COL_IN)
    left, bottom, width, height = _ONE_COL_AXES_RECT
    if with_colorbar:
        width *= 0.86
    ax = fig.add_axes((left, bottom, width, height))
    return fig, ax


def figure_two_col(*, height_in: float = 3.0):
    """A two-column figure (7.0 in wide) at the given height."""
    fig = plt.figure(figsize=(TWO_COL_W_IN, height_in))
    ax = fig.add_subplot(1, 1, 1)
    return fig, ax


def save_presentation(fig, path, *, dpi: int | None = None,
                      pad_inches: float = 0.02, transparent: bool = False):
    """Save a PRESENTATION figure without clipping its title / external colorbar.

    Uses ``bbox_inches="tight"`` so an on-canvas ``set_title`` and a
    ``make_axes_locatable`` colorbar (whose label spills past the fixed
    ``figure_one_col`` axes rectangle) are captured in full. This deliberately
    abandons the fixed on-page size — correct ONLY for deck/slide figures, which
    are contain-fit into their slide box. For REPORT figures (no title, colorbar
    inside the rect) keep the bare ``savefig(bbox_inches=None)`` so the
    fixed-dimension idiom holds. See FB-001 / module docstring.
    """
    fig.savefig(path, dpi=dpi if dpi is not None else STYLE_CONFIG["save_dpi"],
                bbox_inches="tight", pad_inches=pad_inches,
                transparent=transparent)
    plt.close(fig)


def sci_notation(x: float, sig: int = 2) -> str:
    """Format ``x`` as a mathtext fragment ``m\\times10^{e}`` for embedding.

    Returns the body WITHOUT ``$`` delimiters so it drops straight into an
    existing mathtext string, e.g. ``rf"$n_0={sci_notation(n)}$/Bohr$^3$"``.
    Enforces the house rule (FB-003): scientific notation in plot text is
    ALWAYS rendered as ×10^n, never the bare-Python ``"1.30e-03"`` form.
    """
    if x == 0 or not math.isfinite(x):
        return "0"
    exp = int(math.floor(math.log10(abs(x))))
    mant = x / 10.0 ** exp
    return rf"{mant:.{sig}f}\times10^{{{exp}}}"
