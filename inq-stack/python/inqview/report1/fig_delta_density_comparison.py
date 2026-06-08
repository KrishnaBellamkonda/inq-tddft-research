"""fig_delta_density_comparison_standard — Δn(jellium - free) small multiples.

Shows the bath response by subtracting free-WP propagation from the
WP-in-jellium run. Isolates wake formation and screening response.

Run:
    python -m inqview.report1.fig_delta_density_comparison_standard
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

RUN_JELL = Path("ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1")
RUN_FREE = Path("ResearchProject/systems/jellium/run_free_wp_L50_E100_sigma1")
L_BOHR = 50.0
AU_TO_FS = 0.02418884
DT_JELL = 0.02

OUT = "docs/reports/report1/figures/fig_delta_density_comparison.png"


def load_vti(path):
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    data = reader.GetOutput()
    dims = data.GetDimensions()
    spacing = data.GetSpacing()
    origin = data.GetOrigin()
    arr = data.GetPointData().GetArray(0)
    flat = np.array([arr.GetValue(i) for i in range(arr.GetNumberOfTuples())])
    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    return rho, origin, spacing


def main() -> None:
    apply_style()

    # Select 6 time steps evenly across the trajectory
    jell_vti_dir = RUN_JELL / "results/raw/vti/density_total"
    free_vti_dir = RUN_FREE / "results/raw/vti/density_rt_total"

    jell_files = sorted(jell_vti_dir.glob("density_t*.vti"))
    # Extract step numbers
    jell_steps = []
    for f in jell_files:
        s = int(f.stem.split("_t")[1])
        jell_steps.append(s)

    # Select 6 evenly spaced (skip t=0 which has no difference)
    useful = [s for s in jell_steps if s >= 20]
    indices = np.linspace(0, len(useful) - 1, 6, dtype=int)
    selected = [useful[i] for i in indices]
    print(f"Selected steps: {selected}")

    panels_data = []
    for step in selected:
        jell_path = jell_vti_dir / f"density_t{step:06d}.vti"
        free_path = free_vti_dir / f"density_t{step:06d}.vti"

        if not free_path.exists():
            # Free run may have different step numbering (every step vs every 2)
            # The jellium run writes every 2 steps, free run every step
            # So jellium step N corresponds to free step N
            free_path = free_vti_dir / f"density_t{step:06d}.vti"

        if not jell_path.exists() or not free_path.exists():
            print(f"  SKIP step {step}: jell={jell_path.exists()}, free={free_path.exists()}")
            continue

        rho_j, origin, spacing = load_vti(str(jell_path))
        rho_f, _, _ = load_vti(str(free_path))

        # Subtract: Δn = n_jellium - n_free
        # Both are total densities; jellium includes bath electrons
        # We want the bath response, so: Δn = n_jell_total - n_free_total
        delta = rho_j - rho_f

        # Take xz slice at y = mid
        ny = delta.shape[1]
        delta_xz = delta[:, ny // 2, :]

        panels_data.append((step, delta_xz))
        print(f"  Step {step}: t={step*DT_JELL:.2f} a.u. = {step*DT_JELL*AU_TO_FS:.3f} fs, "
              f"|Δn|_max={np.abs(delta_xz).max():.4e}")

    if len(panels_data) < 6:
        print(f"WARNING: only {len(panels_data)} panels (expected 6)")

    # Coordinate axes
    dx = spacing[0]
    dz = spacing[2]
    nx, nz = panels_data[0][1].shape
    x = origin[0] + np.arange(nx) * dx
    z = origin[2] + np.arange(nz) * dz

    # Figure
    W = column_widths_in["full"]
    n_cols, n_rows = 3, 2
    panel_h = W / n_cols * 0.9
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(W, panel_h * n_rows * 1.15),
                              sharex=True, sharey=True)
    axes_flat = axes.flatten()

    all_max = max(np.abs(d).max() for _, d in panels_data)
    linthresh = 1e-4 * all_max
    norm = SymLogNorm(linthresh=linthresh, linscale=1.0, vmin=-all_max, vmax=all_max)
    labels = "(a) (b) (c) (d) (e) (f)".split()

    for i, (step, delta_xz) in enumerate(panels_data):
        ax = axes_flat[i]
        t_fs = step * DT_JELL * AU_TO_FS

        im = ax.pcolormesh(z, x, delta_xz, cmap="RdBu_r", norm=norm,
                           shading="auto", rasterized=True)
        ax.set_aspect("equal")
        ax.text(0.97, 0.95, rf"$t = {t_fs:.3f}$ fs", transform=ax.transAxes,
                fontsize=7, ha="right", va="top",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.7, pad=1))
        panel_label(ax, labels[i], x=0.02, y=0.95)

        if i >= n_cols:
            ax.set_xlabel(r"$z$ (Bohr)")
        if i % n_cols == 0:
            ax.set_ylabel(r"$x$ (Bohr)")

    fig.subplots_adjust(right=0.88, hspace=0.15, wspace=0.08)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label(r"$\Delta n_{\mathrm{jell}} - \Delta n_{\mathrm{free}}$ (e/Bohr$^3$)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
