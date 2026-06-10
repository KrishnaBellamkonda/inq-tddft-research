"""fig_leed_transmission — Transmission LEED pattern vs GS density FFT.

2x2 layout:
  (a) Coronene target positions
  (b) |FFT[n_GS(x,y)]|² — analytical diffraction prediction
  (c) Transmission screen 5, step 450, center target
  (d) Transmission screen 5, step 450, C-C bond target

Screen 5 is at z = -14.46 Bohr (transmission side, symmetric to
backscatter screen 14 at z = +14.40 Bohr).

Run:
    python -m inqview.report1.fig_leed_transmission
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from PIL import Image

from inqview.io.leed import load_leed_pattern
from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    panel_label,
    TufteCritic,
)

SCREEN_IDX = 5
STEP = 450

CENTER_SNAP = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens_snapshots"
)
CCBOND_SNAP = Path(
    "ResearchProject/systems/coronene/run_cc_bond/"
    "results/raw/screens/instantaneous"
)
GS_VTI = Path(
    "ResearchProject/systems/coronene/run_save_gs_paper_replica/"
    "results/density_gs/density_t000000.vti"
)
TARGET_IMG = Path(
    "docs/reports/report1/figures/fig_coronene_target.png"
)
OUT = Path("docs/reports/report1/figures/fig_leed_transmission.png")

# Reference screen with full header (for dx recovery)
CENTER_REF = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    f"results/screens/screen_{SCREEN_IDX:02d}.dat"
)
CCBOND_REF = Path(
    "ResearchProject/systems/coronene/run_cc_bond/"
    f"results/raw/screens/total/screen_{SCREEN_IDX:02d}.dat"
)


def fix_snapshot_dx(pattern, ref_screen_path: Path):
    """Override dx/dy from a time-integrated screen that has the full header."""
    ref = load_leed_pattern(ref_screen_path)
    pattern.dx_bohr = ref.dx_bohr
    pattern.dy_bohr = ref.dy_bohr
    pattern.origin_x_bohr = -0.5 * pattern.nx * ref.dx_bohr
    pattern.origin_y_bohr = -0.5 * pattern.ny * ref.dy_bohr
    return pattern


def load_gs_density_xy(vti_path: Path) -> tuple[np.ndarray, float]:
    """Load GS density VTI, integrate over z, return n_2D(x,y) and dx."""
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(vti_path))
    reader.Update()
    data = reader.GetOutput()

    dims = data.GetDimensions()
    spacing = data.GetSpacing()

    arr = data.GetPointData().GetArray(0)
    n_pts = arr.GetNumberOfTuples()
    flat = np.array([arr.GetValue(i) for i in range(n_pts)])

    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    dx = spacing[0]
    dz = spacing[2]
    n2d = rho.sum(axis=2) * dz
    return n2d, dx


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

    # Load transmission screen data
    center_path = CENTER_SNAP / f"step_{STEP:06d}" / f"screen_{SCREEN_IDX:02d}.dat"
    center_pat = fix_snapshot_dx(load_leed_pattern(center_path), CENTER_REF)
    print(f"Center: screen {SCREEN_IDX}, step {STEP}, "
          f"z={center_pat.z_bohr:.2f} Bohr, dx={center_pat.dx_bohr:.4f}")

    ccbond_path = CCBOND_SNAP / f"screen_{SCREEN_IDX:02d}_t{STEP:06d}.dat"
    ccbond_pat = fix_snapshot_dx(load_leed_pattern(ccbond_path), CCBOND_REF)
    print(f"C-C bond: screen {SCREEN_IDX}, step {STEP}, "
          f"z={ccbond_pat.z_bohr:.2f} Bohr, dx={ccbond_pat.dx_bohr:.4f}")

    # Load GS density n_2D(x,y) for comparison
    print("Loading GS density...")
    n2d, gs_dx = load_gs_density_xy(GS_VTI)
    nx_gs = n2d.shape[0]
    gs_extent = (-0.5 * nx_gs * gs_dx, 0.5 * nx_gs * gs_dx,
                 -0.5 * nx_gs * gs_dx, 0.5 * nx_gs * gs_dx)
    print(f"GS density: {nx_gs}x{nx_gs}, dx={gs_dx:.4f}, extent={gs_extent}")

    # Load target image
    target_img = np.array(Image.open(TARGET_IMG))

    # Build 2x2 figure
    W = column_widths_in["full"]
    fig, axes = plt.subplots(2, 2, figsize=(W, W * 0.92),
                              gridspec_kw={"wspace": 0.05, "hspace": 0.20})

    # (a) Target positions — no outline
    axes[0, 0].imshow(target_img, aspect="equal")
    axes[0, 0].set_xticks([])
    axes[0, 0].set_yticks([])
    for spine in axes[0, 0].spines.values():
        spine.set_visible(False)
    panel_label(axes[0, 0], "(a)", x=0.03, y=0.97)

    # (b) n_GS(x,y) — ground-state density for comparison
    im_b = axes[0, 1].imshow(
        n2d, origin="lower",
        extent=gs_extent,
        aspect="equal",
        cmap="inferno",
        vmin=0,
        interpolation="nearest",
    )
    screen_lim = abs(center_pat.extent_bohr[0])
    axes[0, 1].set_xlim(-screen_lim, screen_lim)
    axes[0, 1].set_ylim(-screen_lim, screen_lim)
    axes[0, 1].set_xlabel(r"$x$ (Bohr)")
    axes[0, 1].set_ylabel(r"$y$ (Bohr)")
    panel_label(axes[0, 1], "(b)", x=0.03, y=0.97)

    # (c) Transmission screen: center target
    im_c = plot_screen(axes[1, 0], center_pat)
    panel_label(axes[1, 0], "(c)", x=0.03, y=0.97)

    # (d) Transmission screen: C-C bond target
    im_d = plot_screen(axes[1, 1], ccbond_pat)
    axes[1, 1].set_ylabel("")
    axes[1, 1].set_yticklabels([])
    panel_label(axes[1, 1], "(d)", x=0.03, y=0.97)

    # Shared colorbar for (c) and (d)
    fig.subplots_adjust(bottom=0.10)
    cbar_ax = fig.add_axes([0.15, 0.03, 0.70, 0.015])
    cbar = fig.colorbar(im_c, cax=cbar_ax, orientation="horizontal")
    cbar.set_label(r"Screen density (arb.\ units)", fontsize=8)
    cbar.ax.tick_params(labelsize=6)

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
