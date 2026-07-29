"""fig_density_profile_comparison_standard — z-profile: jellium vs free WP.

Small multiples comparing n_WP(z,t) in jellium against free propagation.
Shows how the bath modifies the WP shape over time.

Run:
    python -m applications.report1.fig_density_profile_comparison_standard
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from applications.report1._shared_style import (
    apply_style, palette_sweep5, column_widths_in, panel_label, TufteCritic,
)

RUN_JELL = Path("ResearchProject/systems/jellium/run_wp_n162_L50_E100_sigma1")
RUN_FREE = Path("ResearchProject/systems/jellium/run_free_wp_L50_E100_sigma1")
L_BOHR = 50.0
DT_JELL = 0.02
AU_TO_FS = 0.02418884

OUT = "docs/reports/report1/figures/fig_density_profile_comparison.png"


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

    jell_vti_dir = RUN_JELL / "results/raw/vti/density_total"
    jell_gs_dir = RUN_JELL / "results/raw/vti/density_gs_system"
    free_vti_dir = RUN_FREE / "results/raw/vti/density_rt_total"

    # Load GS density (static background to subtract)
    gs_files = sorted(jell_gs_dir.glob("*.vti"))
    rho_gs, origin, spacing = load_vti(str(gs_files[0]))
    dx, dy, dz = spacing
    nz = rho_gs.shape[2]
    z = origin[2] + np.arange(nz) * dz
    gs_z_profile = rho_gs.sum(axis=(0, 1)) * dx * dy

    # Select 6 timesteps
    jell_steps = sorted([int(f.stem.split("_t")[1])
                         for f in jell_vti_dir.glob("density_t*.vti")])
    useful = [s for s in jell_steps if s >= 20]
    indices = np.linspace(0, len(useful) - 1, 6, dtype=int)
    selected = [useful[i] for i in indices]

    W = column_widths_in["full"]
    fig, axes = plt.subplots(2, 3, figsize=(W, W * 0.55), sharex=True, sharey=True)
    axes_flat = axes.flatten()
    labels = "(a) (b) (c) (d) (e) (f)".split()

    for i, step in enumerate(selected):
        jell_path = jell_vti_dir / f"density_t{step:06d}.vti"
        free_path = free_vti_dir / f"density_t{step:06d}.vti"

        if not jell_path.exists() or not free_path.exists():
            print(f"  SKIP step {step}")
            continue

        rho_j, _, _ = load_vti(str(jell_path))
        rho_f, _, _ = load_vti(str(free_path))

        # Jellium WP contribution: n_total - n_GS
        jell_wp_z = rho_j.sum(axis=(0, 1)) * dx * dy - gs_z_profile
        # Free WP: just the total (single WP, no bath)
        free_z = rho_f.sum(axis=(0, 1)) * dx * dy

        t_fs = step * DT_JELL * AU_TO_FS
        ax = axes_flat[i]

        ax.plot(z, jell_wp_z, "-", color=palette_sweep5[0], linewidth=0.9,
                label="Jellium" if i == 0 else None)
        ax.plot(z, free_z, "--", color=palette_sweep5[4], linewidth=0.9,
                label="Free" if i == 0 else None)

        ax.text(0.97, 0.95, rf"$t = {t_fs:.3f}$ fs", transform=ax.transAxes,
                fontsize=6, ha="right", va="top")
        panel_label(ax, labels[i], x=0.02, y=0.95)

        if i >= 3:
            ax.set_xlabel(r"$z$ (Bohr)")
        if i % 3 == 0:
            ax.set_ylabel(r"$n(z)$ (e/Bohr)")

    if len(selected) > 0:
        axes_flat[0].legend(fontsize=5.5, loc="upper left", frameon=False)

    fig.tight_layout(h_pad=0.5, w_pad=0.3)

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
