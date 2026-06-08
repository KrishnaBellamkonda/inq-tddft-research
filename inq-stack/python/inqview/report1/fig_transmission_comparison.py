"""fig_transmission_comparison — Transmission: GS density FFT vs simulation screen.

Two side-by-side panels for the coronene center-target run:
  (a) |FFT[n_GS(x,y)]|² — analytical prediction
  (b) Transmission screen from simulation (screen 5, step 450)

Run:
    python -m inqview.report1.fig_transmission_comparison
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from inqview.screens import load_leed_pattern
from inqview.report1._shared_style import (
    apply_style, column_widths_in, panel_label, TufteCritic,
)

SCREEN_IDX = 5
STEP = 450

CENTER_SNAP = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens_snapshots"
)
CENTER_REF = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    f"results/screens/screen_{SCREEN_IDX:02d}.dat"
)
GS_VTI = Path(
    "ResearchProject/systems/coronene/run_save_gs_paper_replica/"
    "results/density_gs/density_t000000.vti"
)
OUT = "docs/reports/report1/figures/fig_transmission_comparison.png"


def fix_snapshot_dx(pattern, ref_path):
    ref = load_leed_pattern(ref_path)
    pattern.dx_bohr = ref.dx_bohr
    pattern.dy_bohr = ref.dy_bohr
    pattern.origin_x_bohr = -0.5 * pattern.nx * ref.dx_bohr
    pattern.origin_y_bohr = -0.5 * pattern.ny * ref.dy_bohr
    return pattern


def load_gs_density_xy(vti_path):
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(str(vti_path))
    reader.Update()
    data = reader.GetOutput()
    dims = data.GetDimensions()
    spacing = data.GetSpacing()
    arr = data.GetPointData().GetArray(0)
    flat = np.array([arr.GetValue(i) for i in range(arr.GetNumberOfTuples())])
    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    n2d = rho.sum(axis=2) * spacing[2]
    return n2d, spacing[0]


def main() -> None:
    apply_style()

    # Load transmission screen
    snap_path = CENTER_SNAP / f"step_{STEP:06d}" / f"screen_{SCREEN_IDX:02d}.dat"
    pat = fix_snapshot_dx(load_leed_pattern(snap_path), CENTER_REF)

    # Load GS density and FFT
    n2d, gs_dx = load_gs_density_xy(GS_VTI)
    nx_gs = n2d.shape[0]
    fft_gs = np.fft.fftshift(np.fft.fft2(n2d))
    power_gs = np.abs(fft_gs) ** 2
    cx = nx_gs // 2
    power_gs[cx - 1:cx + 2, cx - 1:cx + 2] = 0
    dk = 2 * np.pi / (nx_gs * gs_dx)
    k_extent = (-nx_gs / 2 * dk, nx_gs / 2 * dk, -nx_gs / 2 * dk, nx_gs / 2 * dk)

    W = column_widths_in["full"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(W, W * 0.45),
                                      gridspec_kw={"wspace": 0.25})

    # (a) |FFT[n_GS]|²
    vmax_a = np.percentile(power_gs[power_gs > 0], 99)
    ax_a.imshow(power_gs, origin="lower", extent=k_extent, aspect="equal",
                cmap="inferno", vmin=0, vmax=vmax_a,
                interpolation="nearest")
    k_lim = 5.5
    ax_a.set_xlim(-k_lim, k_lim)
    ax_a.set_ylim(-k_lim, k_lim)
    ax_a.set_xlabel(r"$k_x$ (Bohr$^{-1}$)")
    ax_a.set_ylabel(r"$k_y$ (Bohr$^{-1}$)")
    panel_label(ax_a, "(a)", x=0.03, y=0.97)

    # (b) Transmission screen (real space, linear)
    ax_b.imshow(pat.data, origin="lower", extent=pat.extent_bohr,
                aspect="equal", cmap="inferno", vmin=0,
                interpolation="nearest")
    ax_b.set_xlabel(r"$x$ (Bohr)")
    ax_b.set_ylabel(r"$y$ (Bohr)")
    panel_label(ax_b, "(b)", x=0.03, y=0.97)

    critic = TufteCritic()
    for iss in critic.critique(fig):
        print(f"  TufteCritic: {iss}")

    fig.savefig(OUT, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
