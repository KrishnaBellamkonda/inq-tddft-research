"""fig_leed_comparison — Backscattering LEED pattern comparison with Tsubonoya.

Two-panel figure:
(a) LEED-window time-integrated backscattering pattern (FFT of real-space
    screen density -> k-space diffraction)
(b) Best single-snapshot diffraction pattern

Both panels show coronene hexagonal Bragg peak positions overlaid.

Run:
    python -m applications.report1.fig_leed_comparison
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

# --- Paths -------------------------------------------------------------------
SCREEN_DIR = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens"
)
SNAP_DIR = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens_snapshots"
)
LEED_WINDOW_DIR = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens_leed_window"
)
OUT = Path("docs/reports/report1/figures/fig_leed_comparison.png")


def load_screen(path: Path):
    """Load screen .dat. Returns (2D array, metadata dict).

    Handles two header formats:
    - Time-integrated: 2 header lines (line 1: label/z/total_time, line 2: nx/ny/dx/dy)
    - Instantaneous snapshot: 1 header line (label/z/t), grid inferred from data
    """
    lines = open(path).readlines()
    hdr1 = lines[0]
    z_pos = float(hdr1.split("z=")[1].split()[0])

    if lines[1].startswith("#"):
        hdr2 = lines[1].lstrip("# ")
        parts = hdr2.split()
        meta = {}
        for p in parts:
            if "=" in p:
                k, v = p.split("=")
                meta[k] = v
        nx = int(meta["nx"])
        dx = float(meta["dx"])
        skip = 2
    else:
        first_row = np.fromstring(lines[1], sep=" ")
        nx = len(first_row)
        dx = 0.289758
        skip = 1

    data = np.loadtxt(path, skiprows=skip)
    return data.reshape(nx, nx), {"z": z_pos, "nx": nx, "dx": dx}


def screen_to_leed(screen: np.ndarray, dx: float):
    """FFT real-space screen density to get k-space LEED pattern.

    Returns (power_spectrum, k-axis).
    """
    nx = screen.shape[0]
    # 2D FFT -> shift zero-frequency to center
    fft = np.fft.fft2(screen)
    fft_shifted = np.fft.fftshift(fft)
    power = np.abs(fft_shifted)**2

    # k-space coordinates
    dk = 2 * np.pi / (nx * dx)
    k = (np.arange(nx) - nx // 2) * dk
    return power, k


def plot_leed_panel(ax, power, k, *, vmin=None, vmax=None, suppress_dc=True):
    """Plot LEED diffraction pattern with Bragg peak overlay."""
    nx = len(k)

    if suppress_dc:
        # Suppress the central DC peak for better contrast
        cx = nx // 2
        power = power.copy()
        power[cx-1:cx+2, cx-1:cx+2] = 0

    if vmax is None:
        vmax = power.max()
    if vmin is None:
        vmin = max(1e-3 * vmax, power[power > 0].min())

    im = ax.pcolormesh(
        k, k, power,
        cmap="inferno",
        norm=LogNorm(vmin=vmin, vmax=vmax),
        shading="auto",
        rasterized=True,
    )
    ax.set_aspect("equal")

    # Bragg peak positions for coronene (honeycomb lattice)
    # C-C bond length ~1.42 Ang = 2.68 Bohr
    a_cc = 2.68  # Bohr
    # Honeycomb reciprocal lattice vectors magnitude
    k_outer = 2 * (2 * np.pi / a_cc) / np.sqrt(3)
    k_inner = (2 * np.pi / a_cc) / np.sqrt(3)

    for k_bragg, angles_deg in [
        (k_inner, [30, 90, 150, 210, 270, 330]),
        (k_outer, [0, 60, 120, 180, 240, 300]),
    ]:
        for a_deg in angles_deg:
            a_rad = np.radians(a_deg)
            ax.plot(
                k_bragg * np.cos(a_rad), k_bragg * np.sin(a_rad),
                "o", markersize=4.5, markerfacecolor="none",
                markeredgecolor="white", markeredgewidth=0.8,
            )

    k_lim = k_outer * 1.8
    ax.set_xlim(-k_lim, k_lim)
    ax.set_ylim(-k_lim, k_lim)
    ax.set_xlabel(r"$k_x$ (Bohr$^{-1}$)")
    ax.set_ylabel(r"$k_y$ (Bohr$^{-1}$)")

    return im


def main() -> None:
    apply_style()

    # --- Load LEED-windowed screen -------------------------------------------
    screen_win, meta_win = load_screen(LEED_WINDOW_DIR / "screen_10.dat")
    power_win, k_win = screen_to_leed(screen_win, meta_win["dx"])
    print(f"LEED-window screen: z={meta_win['z']:.3f}")

    # --- Load best instantaneous snapshot ------------------------------------
    best_power = None
    best_max = 0
    best_step = None
    best_k = None
    for step in [150, 180, 210, 240, 270, 300]:
        snap_path = SNAP_DIR / f"step_{step:06d}" / "screen_10.dat"
        if snap_path.exists():
            s, m = load_screen(snap_path)
            p, kk = screen_to_leed(s, m["dx"])
            # Suppress DC for comparison
            cx = len(kk) // 2
            p_test = p.copy()
            p_test[cx-1:cx+2, cx-1:cx+2] = 0
            pmax = p_test.max()
            if pmax > best_max:
                best_max = pmax
                best_power = p
                best_step = step
                best_k = kk
    if best_power is not None:
        print(f"Best snapshot: step {best_step}")

    # --- Figure: 2-panel -----------------------------------------------------
    W = column_widths_in["full"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(W, W * 0.48))

    # Panel (a): LEED-windowed
    im1 = plot_leed_panel(ax1, power_win, k_win)
    ax1.text(0.03, 0.97, "(a)", transform=ax1.transAxes, fontsize=10,
             ha="left", va="top", color="white",
             bbox=dict(facecolor="black", edgecolor="none", alpha=0.5, pad=1))

    # Panel (b): instantaneous snapshot
    if best_power is not None:
        im2 = plot_leed_panel(ax2, best_power, best_k)
        ax2.text(0.03, 0.97, "(b)", transform=ax2.transAxes, fontsize=10,
                 ha="left", va="top", color="white",
                 bbox=dict(facecolor="black", edgecolor="none", alpha=0.5, pad=1))
    else:
        screen_int, meta_int = load_screen(SCREEN_DIR / "screen_10.dat")
        power_int, k_int = screen_to_leed(screen_int, meta_int["dx"])
        im2 = plot_leed_panel(ax2, power_int, k_int)
        ax2.text(0.03, 0.97, "(b)", transform=ax2.transAxes, fontsize=10,
                 ha="left", va="top", color="white",
                 bbox=dict(facecolor="black", edgecolor="none", alpha=0.5, pad=1))

    # Shared colorbar
    fig.subplots_adjust(right=0.88, wspace=0.35)
    cbar_ax = fig.add_axes([0.90, 0.15, 0.015, 0.7])
    cbar = fig.colorbar(im1, cax=cbar_ax)
    cbar.set_label(r"$|I(\mathbf{k})|^2$ (arb.\ units)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    # --- Save ----------------------------------------------------------------
    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    fig.savefig(str(OUT), dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved -> {OUT}")
    plt.close(fig)


if __name__ == "__main__":
    main()
