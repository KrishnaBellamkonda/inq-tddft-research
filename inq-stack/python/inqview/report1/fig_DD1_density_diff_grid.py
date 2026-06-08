"""fig_DD1_density_diff_grid — Bath response: delta-density xz-slices.

Shows Δn(r,t) = n(r,t) - n_GS(r) at 4 time snapshots as xz-slices
through the y=0 plane, using the density_delta VTI output from the
WP-in-jellium run (E=100 eV, sigma=1 Bohr, N=162, L=50 Bohr).

Run:
    python -m inqview.report1.fig_DD1_density_diff_grid
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from pathlib import Path

from inqview.report1._shared_style import (
    apply_style, column_widths_in, panel_label, TufteCritic,
)

# ---------------------------------------------------------------------------
# Data source
# ---------------------------------------------------------------------------
DELTA_DIR = Path(
    "ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1/"
    "results/raw/vti/density_delta"
)
DT_AU = 0.02          # time step (a.u.)
L_BOHR = 50.0
N_PANELS = 6
TARGET_TIMES_AU = [0.4, 2.2, 4.0, 5.8, 7.6, 9.5]

OUT = "docs/reports/report1/figures/fig_DD1_density_diff_grid.png"


def load_vti(path: str):
    """Load a VTI file, return (rho_3d, origin, spacing).

    rho_3d is shaped (nz, ny, nx) following VTK point-data convention.
    """
    import vtk
    from vtk.util.numpy_support import vtk_to_numpy

    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    img = reader.GetOutput()
    dims = img.GetDimensions()      # (nx, ny, nz)
    spacing = img.GetSpacing()
    origin = img.GetOrigin()
    arr = vtk_to_numpy(img.GetPointData().GetArray(0))
    rho = arr.reshape(dims[2], dims[1], dims[0])   # (nz, ny, nx)
    return rho, origin, spacing


def main() -> None:
    apply_style()

    # Discover available frames, skip t=0 (all zeros)
    vti_files = sorted(DELTA_DIR.glob("density_delta_t*.vti"))
    steps = []
    for f in vti_files:
        s = int(f.stem.split("_t")[1])
        if s > 0:
            steps.append(s)

    # Pick frames closest to the target times
    all_times = np.array([s * DT_AU for s in steps])
    selected = []
    for t_target in TARGET_TIMES_AU:
        idx = np.argmin(np.abs(all_times - t_target))
        selected.append(steps[idx])
    print(f"Selected steps: {selected}")

    # Load slices
    panels = []
    for step in selected:
        path = DELTA_DIR / f"density_delta_t{step:06d}.vti"
        rho, origin, spacing = load_vti(str(path))
        nz, ny, nx = rho.shape
        slice_xz = rho[:, ny // 2, :]       # (nz, nx)
        t_au = step * DT_AU
        panels.append((step, t_au, slice_xz))
        print(f"  Step {step}: t = {t_au:.2f} a.u., "
              f"|Δn|_max = {np.abs(slice_xz).max():.4e} e/Bohr^3")

    # Coordinate axes (Bohr)
    nz, ny, nx = rho.shape
    x = origin[0] + np.arange(nx + 1) * spacing[0]   # cell edges for pcolormesh
    z = origin[2] + np.arange(nz + 1) * spacing[2]

    # Shared colour scale
    all_max = max(np.abs(d).max() for _, _, d in panels)
    linthresh = 1e-3 * all_max
    norm = SymLogNorm(linthresh=linthresh, linscale=1.0,
                      vmin=-all_max, vmax=all_max)

    # Figure: 1 row x 4 panels, full width
    W = column_widths_in["full"]
    panel_w = W / N_PANELS
    fig_h = panel_w * (nz / nx) * 1.15 + 0.45   # aspect + colourbar room
    fig, axes = plt.subplots(1, N_PANELS, figsize=(W, fig_h),
                              sharey=True,
                              gridspec_kw={"wspace": 0.08})

    labels = ["(a)", "(b)", "(c)", "(d)", "(e)", "(f)"]
    im = None
    for i, (step, t_au, slice_xz) in enumerate(panels):
        ax = axes[i]
        im = ax.pcolormesh(x, z, slice_xz, cmap="RdBu_r", norm=norm,
                           shading="flat", rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlabel(r"$x$ (Bohr)")
        if i == 0:
            ax.set_ylabel(r"$z$ (Bohr)")

        # Time annotation in whitespace
        ax.text(0.97, 0.97, rf"$t = {t_au:.1f}$ a.u.",
                transform=ax.transAxes, fontsize=7, ha="right", va="top",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
        panel_label(ax, labels[i], x=0.03, y=0.97)

    # Shared colourbar
    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.90, 0.18, 0.015, 0.65])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(r"$\Delta n$ (e/Bohr$^3$)", fontsize=8)
    cbar.ax.tick_params(labelsize=6)

    # Tufte critique
    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    Path(OUT).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
