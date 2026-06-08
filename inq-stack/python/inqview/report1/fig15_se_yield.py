"""fig15 — Secondary electron yield from screen integrals.

Computes N_SE from integrated screen intensities on the backscattering
(entrance) and transmission (exit) sides of coronene. Single run: the
paper-replica coronene propagation.

Run:
    python -m inqview.report1.fig15_se_yield
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

from inqview.report1._shared_style import (
    apply_style,
    palette_sweep5,
    palette_regime3,
    column_widths_in,
    panel_label,
    TufteCritic,
)

SCREEN_DIR = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens"
)

# windowed screens (if available)
WINDOWED_DIR = Path(
    "ResearchProject/systems/coronene/run_propagate_paper_replica/"
    "results/screens_leed_window"
)

BACK_INDICES = list(range(10, 20))   # z > 0
FWD_INDICES = list(range(0, 10))     # z < 0


def load_screen_integral(path: Path) -> tuple[float, float]:
    """Load screen and return (z_pos, total_integrated_intensity)."""
    lines = open(path).readlines()
    z_pos = float(lines[0].split("z=")[1].split()[0])
    parts = lines[1].lstrip("# ").split()
    dx = float(parts[2].split("=")[1])
    dy = float(parts[3].split("=")[1])
    data = np.loadtxt(path, skiprows=2)
    total = data.sum() * dx * dy
    return z_pos, total


def main() -> None:
    apply_style()

    # use windowed screens if available, else main screens
    screen_dir = WINDOWED_DIR if WINDOWED_DIR.exists() else SCREEN_DIR
    label_suffix = " (windowed)" if screen_dir == WINDOWED_DIR else " (total)"

    # compute integrated intensity per screen
    back_z, back_I = [], []
    for idx in BACK_INDICES:
        path = screen_dir / f"screen_{idx:02d}.dat"
        if path.exists():
            z, I = load_screen_integral(path)
            back_z.append(z)
            back_I.append(I)

    fwd_z, fwd_I = [], []
    for idx in FWD_INDICES:
        path = screen_dir / f"screen_{idx:02d}.dat"
        if path.exists():
            z, I = load_screen_integral(path)
            fwd_z.append(abs(z))  # use |z| for comparison
            fwd_I.append(I)

    W = column_widths_in["single"]
    fig, ax = plt.subplots(figsize=(W, W * 0.75))

    c_back = palette_regime3[2]  # green
    c_fwd = "#C07020"           # warm orange

    if back_z:
        ax.plot(back_z, back_I, "o-", color=c_back, markersize=3.5,
                linewidth=0.9, label=f"Backscattering{label_suffix}")
    if fwd_z:
        ax.plot(fwd_z, fwd_I, "s-", color=c_fwd, markersize=3.5,
                linewidth=0.9, label=f"Transmission{label_suffix}")

    # total yields
    N_back = sum(back_I) if back_I else 0
    N_fwd = sum(fwd_I) if fwd_I else 0

    ax.text(0.97, 0.97,
            (f"$N_{{\\mathrm{{SE,back}}}} = {N_back:.1f}$\n"
             f"$N_{{\\mathrm{{SE,fwd}}}} = {N_fwd:.1f}$\n"
             f"$N_{{\\mathrm{{SE,total}}}} = {N_back + N_fwd:.1f}$"),
            transform=ax.transAxes, fontsize=6.5,
            va="top", ha="right",
            bbox=dict(facecolor="white", edgecolor="#b0b0b0",
                      linewidth=0.4, pad=3, alpha=0.9))

    ax.set_xlabel(r"$|z_{\mathrm{screen}}|$ (Bohr)")
    ax.set_ylabel(r"Integrated intensity (arb.)")
    ax.set_yscale("log")
    ax.legend(fontsize=6, loc="center right", frameon=True, framealpha=0.9,
              edgecolor="#b0b0b0")

    critic = TufteCritic()
    issues = critic.critique(fig)
    if issues:
        print(f"TufteCritic: {len(issues)} issue(s)")
        for iss in issues:
            print(f"  {iss}")

    out = "docs/reports/report1/figures/fig15_se_yield.png"
    fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
