"""Shared matplotlib style for Report 1 figures.

Executable counterpart to docs/reports/report1/figures/global_style.md.
Every fig0N_<slug>.py script must call apply_style() at the top of main()
so the rcParams stay in lock-step with the specification.

Tufte critique is enabled by default. To disable for a specific figure:

    from inqview.report1 import TufteCritic
    TufteCritic.enabled = False

Or per-save:

    critic = TufteCritic()
    critic.save(fig, "path/to/fig.png")           # critique + save
    critic.save(fig, "path/to/fig.png", skip=True) # save without critique

Usage
-----
    from inqview.report1 import apply_style, palette_sweep5, column_widths_in
    apply_style()
    fig, ax = plt.subplots(figsize=(column_widths_in["single"], 2.2))
    ax.plot(x, y, color=palette_sweep5[0])
    fig.savefig("docs/reports/report1/figures/fig01_<slug>.png",
                dpi=600, bbox_inches="tight", pad_inches=0.02)
"""

from __future__ import annotations

import warnings
from contextlib import contextmanager
from typing import Iterator

import matplotlib as mpl


# ═══════════════════════════════════════════════════════════════════════════════
# TUNEABLE STYLE CONFIG — change here to restyle all figures at once
# ═══════════════════════════════════════════════════════════════════════════════

STYLE_CONFIG = {
    # --- LaTeX font packages (swap this one line to change all fonts) ---
    # Current: scalable Computer Modern (lmodern + mlmodern)
    # Alternative: Times — r"\usepackage{newtxtext}\usepackage{newtxmath}"
    "latex_packages": r"\usepackage{lmodern}\usepackage{mlmodern}",
    "latex_extra":    r"\usepackage{amsmath}\usepackage{siunitx}",
    "font_family":    "serif",
    "font_serif":     ["Computer Modern Roman"],

    # --- Font sizes (pt) — SAVE-AT-FINAL-WIDTH: plots are saved at their true
    #     on-page width so these render at their literal pt size (8-12 band).
    #     Do NOT bump to compensate for LaTeX downscaling (ad-hoc fix rejected
    #     2026-05-28). Tune once here; confirm adequacy in LaTeX. ---
    "font_size":      10,
    "axes_label":     10,
    "axes_title":     10,
    "tick_label":     9,
    "legend":         9,

    # --- Line widths ---
    "axes_linewidth":  0.8,
    "line_width":      1.2,
    "patch_linewidth": 1.0,

    # --- Tick geometry ---
    "tick_direction":      "in",
    "tick_top":            True,
    "tick_right":          True,
    "tick_major_size":     3.5,
    "tick_major_width":    0.6,
    "tick_minor_size":     2.0,
    "tick_minor_width":    0.5,

    # --- Default colormaps (resolved 2026-05-28) ---
    "cmap_diverging":  "RdBu_r",   # signed / difference maps (zero-centred)
    "cmap_sequential": "inferno",  # positive-only magnitude maps

    # --- Legend ---
    "legend_frameon":      True,
    "legend_framealpha":   1.0,
    "legend_edgecolor":    "#404040",
    "legend_fancybox":     False,

    # --- Figure / save ---
    "preview_dpi":   150,
    "save_dpi":      600,
    "save_pad":      0.02,
}


column_widths_in = {
    "single": 5.00,    # 127 mm  one-column moderate (for single-panel plots)
    "1.5col": 5.50,    # 140 mm  intermediate
    "full":   6.50,    # 165 mm  one-column full text width (elsarticle 3p)
}

# draft5 = single-column `article`, A4, 27 mm side margins → textwidth 156 mm.
DRAFT5_TEXTWIDTH_IN = 156.0 / 25.4   # 6.142 in  (\textwidth, \figwidthwide)
draft5_widths_in = {
    "full":  DRAFT5_TEXTWIDTH_IN,            # 6.14  \figwidthwide / minipage \textwidth
    "large": 0.85 * DRAFT5_TEXTWIDTH_IN,     # 5.22  \figwidthlarge
    "med":   0.65 * DRAFT5_TEXTWIDTH_IN,     # 3.99  \figwidthmed
    "small": 0.45 * DRAFT5_TEXTWIDTH_IN,     # 2.76  \figwidthsmall
    "half":  0.48 * DRAFT5_TEXTWIDTH_IN,     # 2.95  0.48\textwidth (2-per-row minipage)
}

# STANDARD FIXED-COLUMN SIZE SCHEME (set 2026-05-29). Plot width is fixed by
# how many columns it occupies, NOT by the current draft's textwidth:
#   - one-column plot : 3.5 × 3.0 in
#   - two-column plot : 7.0 in wide
# A sub-plot in a 2-up row is a one-column plot (3.5 in); the pair spans 7 in.
ONE_COL_IN = (3.5, 3.0)    # (width, height) for a single-column plot
TWO_COL_W_IN = 7.0         # width for a plot spanning both columns

# FIXED axes rectangle (figure fractions) for one-column plots, so EVERY
# one-column plot has an identical data box and the panels align. Margins are
# sized to fit the worst-case y-label; plots with short labels get extra
# whitespace rather than a bigger axes. Use instead of layout="constrained".
# Resulting axes ≈ (0.775×3.5) × (0.805×3.0) = 2.71 × 2.42 in for all.
# (Plots with a colourbar override `right` per-plot — handled case by case.)
ONE_COL_MARGINS = dict(left=0.200, right=0.975, bottom=0.160, top=0.965)


def fix_one_col_axes(fig):
    """Pin a one-column figure's axes to ONE_COL_MARGINS (no constrained layout)."""
    fig.subplots_adjust(**ONE_COL_MARGINS)


# Square one-column heatmap (3.5 × 3.5 in). Pinned so every square heatmap has
# an identical data box. `cbar=True` reserves room on the right for a colourbar.
SQ_MARGINS          = dict(left=0.150, right=0.965, bottom=0.135, top=0.965)
SQ_MARGINS_CBAR     = dict(left=0.150, right=0.820, bottom=0.135, top=0.965)


def fix_square_axes(fig, *, cbar=False):
    """Pin a 3.5×3.5 heatmap's axes so all square heatmaps share one data box."""
    fig.subplots_adjust(**(SQ_MARGINS_CBAR if cbar else SQ_MARGINS))


# Pixel-sampled from QuantumKickExtension/iter 5/*.png
palette_sweep5 = [
    "#881818",   # wine red    — low
    "#C03828",   # brick red   — low-mid
    "#783898",   # purple      — mid
    "#2070A0",   # steel blue  — mid-high
    "#185070",   # deep navy   — high
]

palette_sweep3 = [
    "#881818",   # wine red    — low
    "#783898",   # purple      — mid
    "#185070",   # deep navy   — high
]

palette_regime3 = [
    "#881818",   # wine red    — Low-v / Characteristic
    "#185070",   # deep navy   — Mid-v / Softening
    "#188048",   # forest green— High-v / Hardening
]

regime_tints = {
    "characteristic": "#FCEAEA",   # palette_regime3[0] @ alpha 0.15 on white
    "softening":      "#E6ECF1",   # palette_regime3[1] @ alpha 0.15 on white
    "hardening":      "#E6F2EB",   # palette_regime3[2] @ alpha 0.15 on white
}

references = {
    "asymptote":    {"color": "#808080", "linestyle": "--", "linewidth": 0.9},
    "theory":       {"color": "#000000", "linestyle": "--", "linewidth": 1.0},
    "fit_overlay":  {"color": "#881818", "linestyle": "-",  "linewidth": 1.0},
    "annotation":   {"color": "#404040"},
}


def apply_style() -> None:
    """Set rcParams to the report-1 standard. Idempotent.

    All values are read from STYLE_CONFIG so a single edit restyles every
    figure on the next run.
    """
    c = STYLE_CONFIG
    preamble = c["latex_packages"] + c["latex_extra"]

    mpl.rcParams.update({
        # text rendering
        "text.usetex": True,
        "font.family": c["font_family"],
        "font.serif":  c["font_serif"],
        "text.latex.preamble": preamble,

        # font sizes
        "axes.labelsize":  c["axes_label"],
        "axes.titlesize":  c["axes_title"],
        "xtick.labelsize": c["tick_label"],
        "ytick.labelsize": c["tick_label"],
        "legend.fontsize": c["legend"],
        "font.size":       c["font_size"],

        # line widths
        "axes.linewidth":  c["axes_linewidth"],
        "lines.linewidth": c["line_width"],
        "patch.linewidth": c["patch_linewidth"],

        # ticks
        "xtick.direction":    c["tick_direction"],
        "ytick.direction":    c["tick_direction"],
        "xtick.top":          c["tick_top"],
        "ytick.right":        c["tick_right"],
        "xtick.major.size":   c["tick_major_size"],
        "ytick.major.size":   c["tick_major_size"],
        "xtick.major.width":  c["tick_major_width"],
        "ytick.major.width":  c["tick_major_width"],
        "xtick.minor.size":   c["tick_minor_size"],
        "ytick.minor.size":   c["tick_minor_size"],
        "xtick.minor.width":  c["tick_minor_width"],
        "ytick.minor.width":  c["tick_minor_width"],

        # legend
        "legend.frameon":     c["legend_frameon"],
        "legend.framealpha":  c["legend_framealpha"],
        "legend.edgecolor":   c["legend_edgecolor"],
        "legend.fancybox":    c["legend_fancybox"],

        # figure
        "figure.dpi":         c["preview_dpi"],
        "savefig.dpi":        c["save_dpi"],
        # "standard" = save the full fixed canvas so the PNG's physical size
        # == figsize exactly (save-at-final-width). Do NOT use "tight" here:
        # it re-trims the canvas and breaks the width guarantee. Legacy scripts
        # that want trimming still pass bbox_inches="tight" explicitly.
        "savefig.bbox":       "standard",
        "savefig.pad_inches": c["save_pad"],
        "savefig.facecolor":  "white",
        "figure.facecolor":   "white",
        "axes.facecolor":     "white",

        # default colour cycle
        "axes.prop_cycle": mpl.cycler(color=palette_sweep5),
    })


def panel_label(ax, label: str, *, x: float = 0.02, y: float = 0.95) -> None:
    """Place a panel label '(a)', '(b)', ... at the top-left inside `ax`.

    The default position is governed by global_style.md section 7.
    """
    ax.text(
        x, y, label,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=mpl.rcParams["axes.labelsize"],
    )


# ---------------------------------------------------------------------------
# Tufte Critic — enabled by default, removable
# ---------------------------------------------------------------------------

class TufteCritic:
    """Automated Tufte-principles checker for matplotlib figures.

    Enabled by default. Runs a battery of heuristic checks against Tufte's
    core principles (data-ink ratio, chartjunk, graphical integrity, labeling,
    data density). Emits warnings for violations.

    Disable globally:
        TufteCritic.enabled = False

    Disable per-context:
        with TufteCritic.disabled():
            fig.savefig(...)
    """

    enabled: bool = True

    @classmethod
    @contextmanager
    def disabled(cls) -> Iterator[None]:
        """Context manager to temporarily suppress Tufte checks."""
        prev = cls.enabled
        cls.enabled = False
        try:
            yield
        finally:
            cls.enabled = prev

    def critique(self, fig) -> list[str]:
        """Run Tufte checks on a matplotlib Figure. Returns a list of warnings."""
        if not self.enabled:
            return []

        issues: list[str] = []
        axes_list = fig.get_axes()

        for i, ax in enumerate(axes_list):
            tag = f"axes[{i}]"

            # --- 1. Chartjunk: figure-level title is discouraged ---
            if ax.get_title():
                issues.append(
                    f"[chartjunk] {tag} has a title '{ax.get_title()[:40]}…' — "
                    f"captions belong in the LaTeX body, not on the figure."
                )

            # --- 2. Data-ink: check for visible grid ---
            if ax.xaxis.get_gridlines() and any(
                gl.get_visible() for gl in ax.xaxis.get_gridlines()
            ):
                issues.append(
                    f"[data-ink] {tag} has visible x-gridlines. "
                    f"Remove or mute to gray ≤0.3 alpha."
                )
            if ax.yaxis.get_gridlines() and any(
                gl.get_visible() for gl in ax.yaxis.get_gridlines()
            ):
                issues.append(
                    f"[data-ink] {tag} has visible y-gridlines. "
                    f"Remove or mute to gray ≤0.3 alpha."
                )

            # --- 3. Labeling: axes without labels ---
            if not ax.get_xlabel():
                issues.append(
                    f"[labeling] {tag} has no x-axis label. "
                    f"Every axis with units should be labeled."
                )
            if not ax.get_ylabel():
                issues.append(
                    f"[labeling] {tag} has no y-axis label."
                )

            # --- 4. Integrity: check for 3D projection ---
            if hasattr(ax, "get_proj") or ax.name == "3d":
                issues.append(
                    f"[integrity] {tag} uses 3D projection. "
                    f"3D charts distort proportions — prefer 2D with "
                    f"layering/color for the extra dimension."
                )

            # --- 5. Data density: empty axes ---
            has_data = bool(ax.lines or ax.collections or ax.patches
                           or ax.images or ax.texts)
            if not has_data:
                issues.append(
                    f"[density] {tag} appears to contain no data elements."
                )

            # --- 6. Legend placement: outside the axes box ---
            legend = ax.get_legend()
            if legend is not None:
                bbox = legend.get_window_extent(
                    fig.canvas.get_renderer()
                    if hasattr(fig.canvas, "get_renderer")
                    else None
                )
                if bbox is not None:
                    ax_bbox = ax.get_window_extent()
                    if (bbox.x0 < ax_bbox.x0 - 20
                            or bbox.x1 > ax_bbox.x1 + 20):
                        issues.append(
                            f"[data-ink] {tag} legend extends beyond the "
                            f"axes box. Prefer labels near the data."
                        )

        # --- 7. Figure-level: suptitle ---
        if fig._suptitle is not None:
            issues.append(
                "[chartjunk] Figure has a suptitle — remove it; "
                "captions go in the LaTeX body."
            )

        return issues

    def warn(self, fig) -> list[str]:
        """Run critique and emit Python warnings for each issue."""
        issues = self.critique(fig)
        for issue in issues:
            warnings.warn(f"TufteCritic: {issue}", stacklevel=2)
        return issues

    def save(self, fig, path: str, *, skip: bool = False, **kwargs) -> list[str]:
        """Critique (unless skip=True), then save the figure.

        Keyword arguments are forwarded to fig.savefig().
        Returns the list of Tufte issues found (empty if skip=True).
        """
        issues = [] if skip else self.warn(fig)
        defaults = {"dpi": 600, "bbox_inches": "tight", "pad_inches": 0.02}
        defaults.update(kwargs)
        fig.savefig(path, **defaults)
        return issues
