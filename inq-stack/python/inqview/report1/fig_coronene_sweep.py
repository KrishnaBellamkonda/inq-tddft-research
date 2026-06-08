"""fig_coronene_sweep — LEED backscattering patterns at multiple screen distances.

Small-multiples display: 2x5 grid of I(k_x, k_y) patterns from
10 backscattering screens at different z-positions. Uses the windowed
LEED screens from the coronene paper-replica run.

Run:
    python -m inqview.report1.fig_coronene_sweep
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
    "results/screens_leed_window"
)

# backscattering screens: indices 10-19, z > 0
BACK_INDICES = list(range(10, 20))


def load_screen(path: Path) -> tuple[np.ndarray, float, int, float]:
    """Load a screen .dat file. Returns (2D array fftshifted, z_pos, nx, dx)."""
    lines = open(path).readlines()

    # parse header: # label=screen_NN z=ZZ total_time=TT n_accum=NN
    z_pos = float(lines[0].split("z=")[1].split()[0])

    header2 = lines[1].lstrip("# ")
    parts = header2.split()
    nx = int(parts[0].split("=")[1])
    ny = int(parts[1].split("=")[1])
    dx = float(parts[2].split("=")[1])

    data = np.loadtxt(path, skiprows=2)
    screen = np.fft.fftshift(data.reshape(ny, nx))

    return screen, z_pos, nx, dx


def main() -> None:
    apply_style()

    # load all backscattering screens
    screens = []
    for idx in BACK_INDICES:
        path = SCREEN_DIR / f"screen_{idx:02d}.dat"
        if not path.exists():
            continue
        screen, z_pos, nx, dx = load_screen(path)
        screens.append((screen, z_pos, nx, dx))

    if not screens:
        print("No backscattering screen data found")
        return

    n_screens = len(screens)

    # Use 2x5 layout for 10 screens
    nrows, ncols = 2, 5
    if n_screens < 10:
        ncols = (n_screens + 1) // 2
        nrows = 2

    W = column_widths_in["full"]
    fig, axes = plt.subplots(nrows, ncols, figsize=(W, W * 0.42),
                             constrained_layout=True)
    axes_flat = axes.flatten()

    # Compute global vmin/vmax from all screens
    all_data = np.concatenate([s[0].ravel() for s in screens])
    vmax = np.percentile(all_data[all_data > 0], 99.5)
    vmin = max(np.percentile(all_data[all_data > 0], 5), 1e-4 * vmax)

    norm = LogNorm(vmin=vmin, vmax=vmax)

    for i, (screen, z_pos, nx, dx) in enumerate(screens):
        if i >= nrows * ncols:
            break
        ax = axes_flat[i]

        # k-space extent
        dk = 2 * np.pi / (nx * dx)
        k_half = (nx // 2) * dk
        extent = [-k_half, k_half, -k_half, k_half]

        # Clip to positive for LogNorm
        screen_clipped = np.clip(screen, vmin, None)

        im = ax.imshow(screen_clipped, origin="lower", extent=extent,
                       cmap="magma", norm=norm, aspect="equal",
                       interpolation="bilinear")

        # Label with z position
        ax.text(0.05, 0.92, f"$z = {z_pos:.1f}$",
                transform=ax.transAxes, fontsize=6, color="white",
                ha="left", va="top",
                bbox=dict(facecolor="black", alpha=0.4, pad=1,
                          edgecolor="none", boxstyle="round,pad=0.2"))

        # Panel label
        panel_label(ax, f"({chr(97 + i)})", x=0.88, y=0.92)
        # Override color for visibility on dark background
        ax.texts[-1].set_color("white")

        # Limit k-range to meaningful region
        k_lim = min(k_half, 8.0)
        ax.set_xlim(-k_lim, k_lim)
        ax.set_ylim(-k_lim, k_lim)

        # Only bottom-row gets x-labels, only left-column gets y-labels
        if i >= ncols:
            ax.set_xlabel(r"$k_x$ (Bohr$^{-1}$)", fontsize=7)
        else:
            ax.set_xticklabels([])
        if i % ncols == 0:
            ax.set_ylabel(r"$k_y$ (Bohr$^{-1}$)", fontsize=7)
        else:
            ax.set_yticklabels([])

        ax.tick_params(labelsize=5, length=2, width=0.4)

    # Hide unused axes
    for j in range(n_screens, nrows * ncols):
        axes_flat[j].set_visible(False)

    # Shared colorbar
    cbar = fig.colorbar(im, ax=axes_flat[:n_screens].tolist(),
                        shrink=0.8, pad=0.01, aspect=30)
    cbar.set_label(r"$I(k_x, k_y)$ (arb.)", fontsize=8)
    cbar.ax.tick_params(labelsize=6)

    out = "docs/reports/report1/figures/fig_coronene_sweep.png"
    with TufteCritic.disabled():
        fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
