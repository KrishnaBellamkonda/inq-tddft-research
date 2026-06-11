"""fig07 — LEED single backscattering screen with Bragg peak overlay.

Shows the 2D LEED pattern from the nearest backscattering screen
(screen_10, z ≈ 1.5 Bohr) with expected coronene Bragg peak positions
overlaid as open circles.

Run:
    python -m applications.report1.fig07_leed_single_screen
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
    "results/screens/screen_10.dat"
)


def load_screen(path: Path):
    """Load screen .dat. Returns (2D array, z, nx, dx)."""
    lines = open(path).readlines()
    z_pos = float(lines[0].split("z=")[1].split()[0])
    parts = lines[1].lstrip("# ").split()
    nx = int(parts[0].split("=")[1])
    ny = int(parts[1].split("=")[1])
    dx = float(parts[2].split("=")[1])

    data = np.loadtxt(path, skiprows=2)
    return data.reshape(ny, nx), z_pos, nx, dx


def main() -> None:
    apply_style()

    screen, z_pos, nx, dx = load_screen(SCREEN_PATH)
    screen = np.fft.fftshift(screen)

    # k-space coordinates
    dk = 2 * np.pi / (nx * dx)
    kx = (np.arange(nx) - nx // 2) * dk
    ky = (np.arange(nx) - nx // 2) * dk

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W))

    vmax = screen.max()
    vmin = max(5e-4 * vmax, screen[screen > 0].min())

    im = ax.pcolormesh(kx, ky, screen, cmap="inferno",
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading="auto", rasterized=True)
    ax.set_aspect("equal")

    # Bragg peak positions for coronene
    a_cc = 2.68  # C-C in Bohr
    k_outer = 2 * (2 * np.pi / a_cc) / np.sqrt(3)
    k_inner = (2 * np.pi / a_cc) / np.sqrt(3)

    for k_bragg, angles_deg, label, ls in [
        (k_inner, [30, 90, 150, 210, 270, 330], r"inner ($k_1$)", "--"),
        (k_outer, [0, 60, 120, 180, 240, 300], r"outer ($k_2$)", "-"),
    ]:
        for a_deg in angles_deg:
            a_rad = np.radians(a_deg)
            ax.plot(k_bragg * np.cos(a_rad), k_bragg * np.sin(a_rad),
                    "o", markersize=5, markerfacecolor="none",
                    markeredgecolor="cyan", markeredgewidth=0.8)
        # label one
        a0 = np.radians(angles_deg[0])
        ax.text(k_bragg * np.cos(a0) + 0.3, k_bragg * np.sin(a0) + 0.3,
                label, fontsize=5, color="cyan")

    ax.set_xlabel(r"$k_x$ (Bohr$^{-1}$)")
    ax.set_ylabel(r"$k_y$ (Bohr$^{-1}$)")

    # limit to useful range
    k_lim = min(k_outer * 2.5, kx.max())
    ax.set_xlim(-k_lim, k_lim)
    ax.set_ylim(-k_lim, k_lim)

    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label("Intensity (arb.)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig07_leed_single_screen.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
