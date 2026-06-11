"""fig09 — Transmission LEED vs analytic 2D-FT of GS density.

Two-panel: (a) far-field transmission screen, (b) |FFT[n_GS(x,y)]|²,
both with Bragg peak overlay.

Run:
    python -m applications.report1.fig09_transmission_vs_fft
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
    column_widths_in,
    panel_label,
    TufteCritic,
)

SCREEN_PATH = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens/screen_00.dat"
)
GS_VTI = (
    "ResearchProject/systems/coronene/run_save_gs_paper_replica/"
    "results/density_gs/density_t000000.vti"
)


def load_screen(path):
    lines = open(path).readlines()
    z_pos = float(lines[0].split("z=")[1].split()[0])
    parts = lines[1].lstrip("# ").split()
    nx = int(parts[0].split("=")[1])
    dx = float(parts[2].split("=")[1])
    data = np.loadtxt(path, skiprows=2)
    return np.fft.fftshift(data.reshape(nx, nx)), nx, dx


def load_gs_density_2d():
    """Load GS density, integrate over z → n_2D(x,y), compute |FFT|²."""
    import vtk
    reader = vtk.vtkXMLImageDataReader()
    reader.SetFileName(GS_VTI)
    reader.Update()
    data = reader.GetOutput()
    dims = data.GetDimensions()
    spacing = data.GetSpacing()
    arr = data.GetPointData().GetArray(0)
    flat = np.array([arr.GetValue(i) for i in range(arr.GetNumberOfTuples())])
    rho = flat.reshape(dims[2], dims[1], dims[0]).transpose(2, 1, 0)
    # integrate over z → n_2D(x,y)
    n2d = rho.sum(axis=2) * spacing[2]
    # 2D FFT
    F = np.fft.fft2(n2d)
    F_shifted = np.fft.fftshift(F)
    power = np.abs(F_shifted)**2
    # normalise to peak = 1
    power /= power.max()
    return power, dims[0], spacing[0]


def main() -> None:
    apply_style()

    screen, nx_s, dx_s = load_screen(SCREEN_PATH)
    print("Loading GS density for FFT...")
    fft_power, nx_g, dx_g = load_gs_density_2d()

    # k-space for screen
    dk_s = 2 * np.pi / (nx_s * dx_s)
    kx_s = (np.arange(nx_s) - nx_s // 2) * dk_s

    # k-space for FFT
    dk_g = 2 * np.pi / (nx_g * dx_g)
    kx_g = (np.arange(nx_g) - nx_g // 2) * dk_g

    # normalise screen to peak=1
    screen_norm = screen / screen.max()

    W = column_widths_in["full"]
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(W, W * 0.45))

    # Bragg peaks
    a_cc = 2.68
    k_inner = (2 * np.pi / a_cc) / np.sqrt(3)
    k_outer = 2 * (2 * np.pi / a_cc) / np.sqrt(3)

    for ax, data, kx, label in [
        (ax_a, screen_norm, kx_s, "(a)"),
        (ax_b, fft_power, kx_g, "(b)"),
    ]:
        vmax = data.max()
        vmin = max(1e-3 * vmax, data[data > 0].min()) if (data > 0).any() else 1e-10
        im = ax.pcolormesh(kx, kx, data, cmap="inferno",
                           norm=LogNorm(vmin=vmin, vmax=vmax),
                           shading="auto", rasterized=True)
        ax.set_aspect("equal")

        for k_bragg, angles in [
            (k_inner, [30, 90, 150, 210, 270, 330]),
            (k_outer, [0, 60, 120, 180, 240, 300]),
        ]:
            for a_deg in angles:
                a_rad = np.radians(a_deg)
                ax.plot(k_bragg * np.cos(a_rad), k_bragg * np.sin(a_rad),
                        "o", markersize=4, markerfacecolor="none",
                        markeredgecolor="cyan", markeredgewidth=0.6)

        k_lim = min(k_outer * 2.5, kx.max())
        ax.set_xlim(-k_lim, k_lim)
        ax.set_ylim(-k_lim, k_lim)
        ax.set_xlabel(r"$k_x$ (Bohr$^{-1}$)")
        panel_label(ax, label)

    ax_a.set_ylabel(r"$k_y$ (Bohr$^{-1}$)")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig09_transmission_vs_fft.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
