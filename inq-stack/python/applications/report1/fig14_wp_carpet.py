"""fig14 — Wave-packet z-vs-t carpet plot.

2D heatmap of the WP projected density ρ_z(z,t) = ∫∫|ψ_wp|² dx dy,
showing propagation (centroid moves diagonally) and spreading (band widens).

Uses the coronene WP density VTI series (61 frames).

Run:
    python -m applications.report1.fig14_wp_carpet
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    column_widths_in,
    panel_label,
    TufteCritic,
)

WP_DIR = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/density_rt_wp"
)


def load_vti_density(path: str) -> tuple[np.ndarray, tuple, tuple]:
    """Load VTI density. Returns (rho[ix,iy,iz], origin, spacing)."""
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

    # discover VTI files
    vti_files = sorted(WP_DIR.glob("density_t*.vti"))
    if not vti_files:
        print("No VTI files found in", WP_DIR)
        return

    print(f"Loading {len(vti_files)} WP density frames...")

    # load first to get grid
    rho0, origin, spacing = load_vti_density(vti_files[0])
    dx, dy, dz = spacing
    nx, ny, nz = rho0.shape
    z = origin[2] + np.arange(nz) * dz

    # extract time steps from filenames
    steps = []
    for f in vti_files:
        s = int(f.stem.split("_t")[1])
        steps.append(s)
    steps = np.array(steps)

    # assume dt = 0.02 a.u. per step (from run_summary); adjust if needed
    dt_au = 0.02
    times = steps * dt_au

    # build carpet: ρ_z(z, t) = ∫∫ ρ dx dy
    carpet = np.zeros((len(vti_files), nz))
    for i, f in enumerate(vti_files):
        if i == 0:
            rho = rho0
        else:
            rho, _, _ = load_vti_density(f)
        carpet[i, :] = rho.sum(axis=(0, 1)) * dx * dy

    # compute centroid and sigma
    z_centroid = np.zeros(len(vti_files))
    z_sigma = np.zeros(len(vti_files))
    for i in range(len(vti_files)):
        rho_z = carpet[i, :]
        norm = rho_z.sum() * dz
        if norm > 1e-12:
            z_centroid[i] = (rho_z * z).sum() * dz / norm
            z_sigma[i] = np.sqrt((rho_z * (z - z_centroid[i])**2).sum() * dz / norm)

    # plot
    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 1.1))

    vmax = carpet.max()
    vmin = max(1e-4 * vmax, carpet[carpet > 0].min()) if (carpet > 0).any() else 1e-10

    im = ax.pcolormesh(z, times, carpet, cmap="viridis",
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading="auto", rasterized=True)

    # centroid overlay
    ax.plot(z_centroid, times, color="white", linewidth=0.8, alpha=0.9)
    # ±σ band
    ax.plot(z_centroid - z_sigma, times, color="white", linewidth=0.5,
            linestyle="--", alpha=0.6)
    ax.plot(z_centroid + z_sigma, times, color="white", linewidth=0.5,
            linestyle="--", alpha=0.6)

    # free-propagation reference
    if len(z_centroid) > 1:
        v_est = (z_centroid[-1] - z_centroid[0]) / (times[-1] - times[0])
        z_free = z_centroid[0] + v_est * (times - times[0])
        ax.plot(z_free, times, color="#FFD700", linewidth=0.5, linestyle=":",
                alpha=0.7)

    ax.set_xlabel(r"$z$ (Bohr)")
    ax.set_ylabel(r"$t$ (a.u.)")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
    cbar.set_label(r"$\rho_z(z,t)$ (Bohr$^{-1}$)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig14_wp_carpet.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
