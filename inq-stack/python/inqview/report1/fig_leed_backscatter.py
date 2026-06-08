"""fig_leed_backscatter — Backscattering LEED: simulation vs Tsubonoya.

2x2 layout:
  (a) Coronene target positions (center + C-C bond)
  (b) Tsubonoya extracted LEED (center impact, Fig 2a)
  (c) This work: screen 14, step 330, center target
  (d) This work: screen 14, C-C bond target

Screen 14 is at z = +14.4 Bohr, comparable to Tsubonoya's observation
distance D = 6.35 Å ≈ 12 Bohr.

Run:
    python -m inqview.report1.fig_leed_backscatter
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path
from PIL import Image

from inqview.screens import load_leed_pattern
from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    panel_label,
    TufteCritic,
)

SCREEN_IDX = 14
CENTER_STEP = 330
CCBOND_STEP = 330
DT = 0.020

CENTER_SNAP = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens_snapshots"
)
CCBOND_SNAP = Path(
    "ResearchProject/systems/coronene/run_cc_bond/"
    "results/raw/screens/instantaneous"
)
TSUB_IMG = Path(
    "docs/reports/report1/figures/tsubonoya_extracted/tsubonoya-comparison.png"
)
TARGET_IMG = Path(
    "docs/reports/report1/figures/fig_coronene_target.png"
)
OUT = Path("docs/reports/report1/figures/fig_leed_backscatter.png")


def fix_snapshot_dx(pattern, ref_screen_path: Path):
    """Override dx/dy from a time-integrated screen that has the full header."""
    ref = load_leed_pattern(ref_screen_path)
    pattern.dx_bohr = ref.dx_bohr
    pattern.dy_bohr = ref.dy_bohr
    pattern.origin_x_bohr = -0.5 * pattern.nx * ref.dx_bohr
    pattern.origin_y_bohr = -0.5 * pattern.ny * ref.dy_bohr
    return pattern


def plot_screen(ax, pattern, *, cmap="inferno"):
    """Plot a LeedPattern with correct physical coordinates (linear scale)."""
    im = ax.imshow(
        pattern.data, origin="lower",
        extent=pattern.extent_bohr,
        aspect="equal",
        cmap=cmap,
        vmin=0,
        interpolation="nearest",
    )
    ax.set_xlabel(r"$x$ (Bohr)")
    ax.set_ylabel(r"$y$ (Bohr)")
    return im


def main() -> None:
    apply_style()

    # Load screen data using the canonical loader (handles fftshift + origin)
    # Snapshot files lack dx metadata — get it from the time-integrated screens
    center_ref = Path("ResearchProject/systems/coronene/run_propagate_paper_replica/"
                      f"results/screens/screen_{SCREEN_IDX:02d}.dat")
    ccbond_ref = Path("ResearchProject/systems/coronene/run_cc_bond/"
                      f"results/raw/screens/total/screen_{SCREEN_IDX:02d}.dat")

    center_path = CENTER_SNAP / f"step_{CENTER_STEP:06d}" / f"screen_{SCREEN_IDX:02d}.dat"
    center_pat = fix_snapshot_dx(load_leed_pattern(center_path), center_ref)
    print(f"Center: screen {SCREEN_IDX}, step {CENTER_STEP} (t={CENTER_STEP*DT:.1f} a.u.), "
          f"z={center_pat.z_bohr:.2f} Bohr, dx={center_pat.dx_bohr:.4f}, extent={center_pat.extent_bohr}")

    ccbond_path = CCBOND_SNAP / f"screen_{SCREEN_IDX:02d}_t{CCBOND_STEP:06d}.dat"
    ccbond_pat = fix_snapshot_dx(load_leed_pattern(ccbond_path), ccbond_ref)
    print(f"C-C bond: screen {SCREEN_IDX}, step {CCBOND_STEP} (t={CCBOND_STEP*DT:.1f} a.u.), "
          f"z={ccbond_pat.z_bohr:.2f} Bohr, dx={ccbond_pat.dx_bohr:.4f}, extent={ccbond_pat.extent_bohr}")

    # Load reference images
    target_img = np.array(Image.open(TARGET_IMG))
    tsub_img = np.array(Image.open(TSUB_IMG))

    # Build 2x2 figure — all panels identical physical size via gridspec
    W = column_widths_in["full"]
    panel_h = W * 0.42  # height of each panel row
    fig = plt.figure(figsize=(W, 2 * panel_h + 0.55))  # extra for colorbar
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1, 1], height_ratios=[1, 1],
        wspace=0.08, hspace=0.25,
        left=0.08, right=0.92, top=0.97, bottom=0.10,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])

    # (a) Target positions — image fills the subplot, no axes frame
    ax_a.imshow(target_img, aspect="auto")
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    panel_label(ax_a, "(a)", x=0.03, y=0.97)

    # (b) Tsubonoya reference — image fills the subplot, no axes frame
    ax_b.imshow(tsub_img, aspect="auto")
    ax_b.set_xticks([])
    ax_b.set_yticks([])
    for spine in ax_b.spines.values():
        spine.set_visible(False)
    panel_label(ax_b, "(b)", x=0.03, y=0.97)

    # (c) This work: center target
    im_c = plot_screen(ax_c, center_pat)
    panel_label(ax_c, "(c)", x=0.03, y=0.97)

    # (d) This work: C-C bond target
    im_d = plot_screen(ax_d, ccbond_pat)
    ax_d.set_ylabel("")
    ax_d.set_yticklabels([])
    panel_label(ax_d, "(d)", x=0.03, y=0.97)

    # Shared colorbar for (c) and (d) — positioned below bottom row
    cbar_ax = fig.add_axes([0.15, 0.04, 0.70, 0.018])
    cbar = fig.colorbar(im_c, cax=cbar_ax, orientation="horizontal")
    cbar.set_label(r"Screen density (arb.\ units)", fontsize=9)
    cbar.ax.tick_params(labelsize=9)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(str(OUT), dpi=600, bbox_inches="tight", pad_inches=0.03)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
