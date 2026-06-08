"""fig10 — LEED waterfall across transmission (forward) screens.

Radially-integrated LEED intensity as a 2D heatmap: (k_r, screen_z).
Uses screens 0-9 (z < 0 = transmission side).

Run:
    python -m inqview.report1.fig10_leed_waterfall_fwd
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pathlib import Path

from inqview.report1._shared_style import (
    apply_style,
    column_widths_in,
    panel_label,
    TufteCritic,
)

SCREEN_DIR = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens"
)

# forward/transmission screens: indices 0-9, z < 0
BACK_INDICES = list(range(0, 10))


def load_screen(path: Path) -> tuple[np.ndarray, float, int, float]:
    """Load a screen .dat file. Returns (2D array, z_pos, nx, dx)."""
    with open(path) as f:
        header = f.readline()

    # parse header: # label=screen_NN z=ZZ total_time=TT n_accum=NN
    z_pos = float(header.split("z=")[1].split()[0])

    header2_line = open(path).readlines()[1].lstrip("# ")
    parts = header2_line.split()
    nx = int(parts[0].split("=")[1])
    ny = int(parts[1].split("=")[1])
    dx = float(parts[2].split("=")[1])

    data = np.loadtxt(path, skiprows=2)
    screen = data.reshape(ny, nx)

    return screen, z_pos, nx, dx


def radial_profile(screen: np.ndarray, nx: int, dx: float,
                   n_bins: int = 60) -> tuple[np.ndarray, np.ndarray]:
    """Compute radially-integrated intensity I(k_r) from a 2D screen."""
    # k-space coordinates
    dk = 2 * np.pi / (nx * dx)
    kx = (np.arange(nx) - nx // 2) * dk
    ky = (np.arange(nx) - nx // 2) * dk
    KX, KY = np.meshgrid(kx, ky)
    KR = np.sqrt(KX**2 + KY**2)

    kr_max = KR.max() * 0.9
    kr_bins = np.linspace(0, kr_max, n_bins + 1)
    kr_centers = 0.5 * (kr_bins[:-1] + kr_bins[1:])

    # shift screen to center (fftshift-like)
    screen_shifted = np.fft.fftshift(screen)

    profile = np.zeros(n_bins)
    for i in range(n_bins):
        mask = (KR >= kr_bins[i]) & (KR < kr_bins[i + 1])
        if mask.any():
            profile[i] = screen_shifted[mask].sum()

    return kr_centers, profile


def main() -> None:
    apply_style()

    # load all backscattering screens
    z_positions = []
    profiles = []
    kr_ref = None

    for idx in BACK_INDICES:
        path = SCREEN_DIR / f"screen_{idx:02d}.dat"
        if not path.exists():
            continue
        screen, z_pos, nx, dx = load_screen(path)
        kr, prof = radial_profile(screen, nx, dx)
        z_positions.append(z_pos)
        profiles.append(prof)
        if kr_ref is None:
            kr_ref = kr

    if not profiles:
        print("No backscattering screen data found")
        return

    # build 2D array
    carpet = np.array(profiles)  # shape: (n_screens, n_bins)
    z_arr = np.array(z_positions)

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 0.85))

    vmax = carpet.max()
    vmin = max(1e-3 * vmax, carpet[carpet > 0].min()) if (carpet > 0).any() else 1e-10

    im = ax.pcolormesh(kr_ref, z_arr, carpet, cmap="magma",
                       norm=LogNorm(vmin=vmin, vmax=vmax),
                       shading="auto", rasterized=True)

    ax.set_xlabel(r"$k_r$ (Bohr$^{-1}$)")
    ax.set_ylabel(r"$z_{\mathrm{screen}}$ (Bohr)")

    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label("Radial intensity (arb.)", fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    # expected Bragg peak positions for coronene
    a_cc = 2.68  # C-C distance in Bohr (1.42 Å)
    k_outer = 2 * (2 * np.pi / a_cc) / np.sqrt(3)
    k_inner = (2 * np.pi / a_cc) / np.sqrt(3)
    for k_bragg, label in [(k_inner, r"$k_1$"), (k_outer, r"$k_2$")]:
        ax.axvline(k_bragg, color="cyan", linewidth=0.5, linestyle="--",
                   alpha=0.7)
        ax.text(k_bragg + 0.05, z_arr[-1] + 0.3, label, fontsize=5,
                color="cyan", va="bottom")

    ax.text(0.97, 0.03, r"\textit{transmission}",
            transform=ax.transAxes, fontsize=5, color="#d0d0d0",
            ha="right", va="bottom")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig10_leed_waterfall_fwd.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
