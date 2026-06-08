"""fig_delta_density_multiples — Small multiples of delta-density xz slices.

Shows Dn(x,z) = n(t) - n(0) integrated over y at several time steps for the
WP-in-jellium run (E=100 eV, sigma=1 Bohr, N=162, L=50).

Run:
    python -m inqview.report1.fig_delta_density_multiples
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from pathlib import Path

from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    panel_label,
    TufteCritic,
    palette_sweep5,
)

# --- Paths -------------------------------------------------------------------
VTI_DIR = Path(
    "ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1/"
    "results/raw/vti/density_delta"
)
OUT = Path("docs/reports/report1/figures/fig_delta_density_multiples.png")

# --- Simulation parameters ---------------------------------------------------
DT_AU = 0.020
L_BOHR = 50.0
AU_TO_FS = 1.0 / 41.341374575751


def load_vti_density(path: str) -> tuple[np.ndarray, tuple, tuple]:
    """Load a VTI density file. Returns (3D array [ix,iy,iz], origin, spacing)."""
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()

    dims = data.GetDimensions()
    spacing = data.GetSpacing()
    origin = data.GetOrigin()

    arr = data.GetPointData().GetArray(0)
    n_pts = arr.GetNumberOfTuples()
    flat = np.array([arr.GetValue(i) for i in range(n_pts)])

    # VTK: Fortran ordering (x fastest)
    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    return rho, origin, spacing


def main() -> None:
    apply_style()

    # --- Select time steps ---------------------------------------------------
    # Available steps: 0 to 474 in steps of 2
    # Choose 6 evenly spaced frames covering the WP transit
    all_files = sorted(VTI_DIR.glob("density_delta_t*.vti"))
    if not all_files:
        print(f"ERROR: No VTI files found in {VTI_DIR}")
        return

    # Extract step numbers
    steps_available = []
    for f in all_files:
        step_str = f.stem.split("_t")[-1]
        steps_available.append(int(step_str))
    steps_available = sorted(steps_available)

    # Choose 6 frames: skip t=0 (delta=0 by construction), then evenly space
    # through the propagation
    n_panels = 6
    # Skip early transient (first ~20 steps) and take from step 40 onward
    useful_steps = [s for s in steps_available if s >= 20]
    indices = np.linspace(0, len(useful_steps) - 1, n_panels, dtype=int)
    selected_steps = [useful_steps[i] for i in indices]
    print(f"Selected steps: {selected_steps}")

    # --- Load data -----------------------------------------------------------
    panels_data = []
    origin = spacing = None
    for step in selected_steps:
        fname = VTI_DIR / f"density_delta_t{step:06d}.vti"
        print(f"  Loading {fname.name}...")
        rho, origin, spacing = load_vti_density(str(fname))
        # Integrate over y to get Dn(x, z)
        dy = spacing[1]
        dn_xz = rho.sum(axis=1) * dy  # shape (nx, nz)
        panels_data.append(dn_xz)

    dx, dy_sp, dz = spacing
    nx, nz = panels_data[0].shape
    x = origin[0] + np.arange(nx) * dx
    z = origin[2] + np.arange(nz) * dz

    # --- Generate both log and linear versions --------------------------------
    _make_figure(panels_data, selected_steps, x, z,
                 use_log=True, out=OUT)
    _make_figure(panels_data, selected_steps, x, z,
                 use_log=False, out=OUT.with_stem(OUT.stem + "_linear"))


def _make_figure(panels_data, selected_steps, x, z, *, use_log: bool, out: Path):
    W = column_widths_in["full"]
    n_cols = 3
    n_rows = 2
    panel_h = W / n_cols * 0.9
    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(W, panel_h * n_rows * 1.15),
        sharex=True, sharey=True,
    )
    axes_flat = axes.flatten()

    all_max = max(np.abs(dn).max() for dn in panels_data)
    if use_log:
        linthresh = 1e-4 * all_max
        norm = SymLogNorm(linthresh=linthresh, linscale=1.0,
                          vmin=-all_max, vmax=all_max)
    else:
        norm = plt.Normalize(vmin=-all_max, vmax=all_max)

    labels = "(a) (b) (c) (d) (e) (f)".split()

    for i, (dn_xz, step) in enumerate(zip(panels_data, selected_steps)):
        ax = axes_flat[i]
        t_fs = step * DT_AU * AU_TO_FS

        im = ax.pcolormesh(
            z, x, dn_xz,
            cmap="RdBu_r",
            norm=norm,
            shading="auto",
            rasterized=True,
        )
        ax.set_aspect("equal")

        ax.text(
            0.97, 0.95, rf"$t = {t_fs:.3f}$ fs",
            transform=ax.transAxes,
            fontsize=7, ha="right", va="top",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1),
        )
        panel_label(ax, labels[i], x=0.02, y=0.95)

        if i >= n_cols:
            ax.set_xlabel(r"$z$ (Bohr)")
        if i % n_cols == 0:
            ax.set_ylabel(r"$x$ (Bohr)")

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(r"$\Delta n(x,z)$ (e/Bohr$^2$)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic ({out.name}): {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(str(out), dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
