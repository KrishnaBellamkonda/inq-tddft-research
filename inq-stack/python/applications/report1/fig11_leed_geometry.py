"""fig11 — LEED simulation geometry schematic.

Side-view of the simulation cell showing WP injection, target plane,
backscattering and transmission screens.
Pure matplotlib line art — no simulation data required.

Run:
    python -m applications.report1.fig11_leed_geometry
"""

from __future__ import annotations

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse, Rectangle

from applications.report1._shared_style import (
    apply_style,
    palette_sweep5,
    palette_regime3,
    column_widths_in,
    TufteCritic,
)


def _gaussian_blob(ax, cx, cz, sigma, color, label=None, alpha=0.35):
    """Draw 1σ, 2σ, 3σ concentric circles with fading opacity."""
    for n_sig, a_frac in [(3, 0.12), (2, 0.22), (1, 0.40)]:
        r = sigma * n_sig
        e = Ellipse((cz, cx), width=2*r, height=2*r,
                    facecolor=color, alpha=a_frac,
                    edgecolor=color, linewidth=0.3, linestyle="--",
                    zorder=3)
        ax.add_patch(e)
        ax.text(cz, cx - r - 0.5, rf"${n_sig}\sigma$", fontsize=4,
                ha="center", va="top", color=color, alpha=0.7)
    if label:
        ax.text(cz, cx + sigma * 3.3, label, fontsize=5.5,
                ha="center", va="bottom", color=color)


def main() -> None:
    apply_style()

    W = column_widths_in["full"]
    fig, ax = plt.subplots(figsize=(W, W * 0.32))

    # ── cell parameters (to scale) ──────────────────────────────────
    Lz = 59.90
    Lx = 34.77
    z_target = 0
    z_wp_init = 12.0
    sigma_wp = 1.0015

    # 10 backscatter + 10 forward screens (actual z-positions from config)
    z_back = [1.485, 4.914, 7.944, 11.374, 14.403, 17.833, 20.863, 24.292, 27.322, 28.952]
    z_fwd = [-1.545, -4.974, -8.004, -11.434, -14.463, -17.893, -20.923, -24.352, -27.382, -28.952]

    # cell boundary
    cell = Rectangle((-Lz/2, -Lx/2), Lz, Lx,
                     fill=False, edgecolor="#404040", linewidth=0.6,
                     linestyle="--", zorder=1)
    ax.add_patch(cell)

    # target plane — full cell height
    ax.plot([z_target, z_target], [-Lx/2, Lx/2],
            color="#404040", linewidth=1.5, zorder=2)
    ax.text(z_target, Lx/2 + 1.5, r"\textbf{target}", fontsize=6,
            ha="center", va="bottom", color="#404040")

    # WP blob
    c_wp = palette_sweep5[4]
    _gaussian_blob(ax, 0, z_wp_init, sigma_wp, c_wp,
                   label=r"WP ($\sigma$, $E$)", alpha=0.4)

    # velocity arrow
    ax.annotate("", xy=(z_wp_init - 8, 0), xytext=(z_wp_init - 1, 0),
                arrowprops=dict(arrowstyle="-|>", color=c_wp,
                                lw=1.2, mutation_scale=10),
                zorder=5)
    ax.text(z_wp_init - 4.5, -3, r"$\mathbf{k}_0$", fontsize=7,
            ha="center", va="top", color=c_wp)

    # backscattering screens — full cell height
    c_back = palette_regime3[2]
    for zs in z_back:
        ax.plot([zs, zs], [-Lx/2, Lx/2],
                color=c_back, linewidth=0.5, alpha=0.6, zorder=1)
    ax.text(z_back[4], -Lx/2 + 1.5, r"\textit{back screens}", fontsize=5,
            ha="center", va="bottom", color=c_back)

    # forward screens — full cell height
    c_fwd = "#C07020"
    for zs in z_fwd:
        ax.plot([zs, zs], [-Lx/2, Lx/2],
                color=c_fwd, linewidth=0.5, alpha=0.6, zorder=1)
    ax.text(z_fwd[4], -Lx/2 + 1.5, r"\textit{fwd screens}", fontsize=5,
            ha="center", va="bottom", color=c_fwd)

    # periodic image labels
    for side in [-Lz/2 - 1, Lz/2 + 1]:
        ax.text(side, 0, r"$\cdots$", fontsize=9, ha="center", va="center",
                color="#808080")

    # dimension annotation — L_z (horizontal, bottom)
    ax.annotate("", xy=(-Lz/2, -Lx/2 - 3), xytext=(Lz/2, -Lx/2 - 3),
                arrowprops=dict(arrowstyle="<->", color="#808080", lw=0.5))
    ax.text(0, -Lx/2 - 4.5, r"$L_z$", fontsize=6, ha="center",
            va="top", color="#808080")

    # dimension annotation — L_x (vertical, right side)
    ax.annotate("", xy=(Lz/2 + 3, -Lx/2), xytext=(Lz/2 + 3, Lx/2),
                arrowprops=dict(arrowstyle="<->", color="#808080", lw=0.5))
    ax.text(Lz/2 + 4.5, 0, r"$L_x$", fontsize=6, ha="left",
            va="center", color="#808080")

    # dimension annotation — L_y (out-of-plane, lower-right corner)
    ax.annotate("", xy=(Lz/2 + 0.5, -Lx/2 - 1.5),
                xytext=(Lz/2 + 4.5, -Lx/2 - 4.5),
                arrowprops=dict(arrowstyle="-|>", color="#808080", lw=0.5))
    ax.text(Lz/2 + 5, -Lx/2 - 5, r"$L_y$", fontsize=6, ha="left",
            va="top", color="#808080")

    ax.set_xlim(-Lz/2 - 5, Lz/2 + 8)
    ax.set_ylim(-Lx/2 - 7, Lx/2 + 4)
    ax.set_aspect("equal")
    ax.axis("off")

    out = "docs/reports/report1/figures/fig11_leed_geometry.png"
    with TufteCritic.disabled():
        fig.savefig(out, dpi=600, bbox_inches="tight", pad_inches=0.02)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
