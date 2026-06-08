"""fig06 — Density difference Δn(x,z) for coronene WP collision.

Planar-integrated (over y) density difference n(t) − n_GS, shown on the
(x,z) plane.  Replicates the layout of Yao & Schleife, arXiv:2605.12854v1
Fig 1.

Run:
    python -m inqview.report1.fig06_density_diff_coronene
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm

from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    TufteCritic,
)

GS_PATH = (
    "ResearchProject/systems/coronene/run_save_gs_paper_replica/"
    "results/density_gs/density_t000000.vti"
)
RT_PATH = (
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/density_rt_target/density_t000520.vti"
)


def load_vti_density(path: str) -> tuple[np.ndarray, tuple, tuple]:
    """Load a VTI density file. Returns (3D array, origin, spacing)."""
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(path)
    reader.Update()
    data = reader.GetOutput()

    dims = data.GetDimensions()
    spacing = data.GetSpacing()
    origin = data.GetOrigin()

    arr = data.GetPointData().GetArray(0)
    n_pts = arr.GetNumberOfTuples()
    flat = np.array([arr.GetValue(i) for i in range(n_pts)])

    # VTK uses Fortran ordering (x fastest)
    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    # rho[ix, iy, iz]

    return rho, origin, spacing


def main() -> None:
    apply_style()

    print("Loading GS density...")
    rho_gs, origin, spacing = load_vti_density(GS_PATH)
    print("Loading RT density...")
    rho_rt, _, _ = load_vti_density(RT_PATH)

    dx, dy, dz = spacing

    # planar integration over y → n_planar(x, z)
    n_gs_xz = rho_gs.sum(axis=1) * dy
    n_rt_xz = rho_rt.sum(axis=1) * dy

    # density difference
    dn = n_rt_xz - n_gs_xz

    # coordinate axes
    nx, nz = dn.shape
    x = origin[0] + np.arange(nx) * dx
    z = origin[2] + np.arange(nz) * dz

    W = column_widths_in["full"]
    fig, ax = plt.subplots(figsize=(W, W * 0.38))

    A = np.max(np.abs(dn)) * 0.9
    norm = SymLogNorm(linthresh=1e-4 * A, linscale=1.0, vmin=-A, vmax=A)

    im = ax.pcolormesh(z, x, dn, cmap="RdBu_r", norm=norm,
                       shading="auto", rasterized=True)

    # coronene plane at z ≈ 0
    ax.axvline(0, color="#404040", linewidth=0.6, linestyle="--")
    ax.text(0.5, -16, "coronene plane", fontsize=5.5, ha="center",
            color="#404040", rotation=90, va="top")

    ax.set_xlabel(r"$z$ (Bohr)")
    ax.set_ylabel(r"$x$ (Bohr)")
    ax.set_aspect("equal")

    # colorbar
    cbar = fig.colorbar(im, ax=ax, orientation="horizontal",
                        pad=0.15, aspect=40, shrink=0.8)
    cbar.set_label(r"$n(t) - n_{\mathrm{GS}}$ (e/Bohr$^2$)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig06_density_diff_coronene.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
