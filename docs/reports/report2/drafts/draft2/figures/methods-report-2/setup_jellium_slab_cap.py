"""Variant of setup_jellium_slab.png with a pastel CAP overlay.

Produces setup_jellium_slab_cap.png alongside the main slab panel.
Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 setup_jellium_slab_cap.py
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]
sys.path.insert(0, str(REPO / "inq-stack/python"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from inqview.visualisation import style
from inqview import load_vti

style.apply_theme()

HERE = Path(__file__).parent

VTI = (
    "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
    "/shared_gs/slab_n102_L25x25x140_w0p5_h0p5/density_gs_system/density_gs_system.vti"
)
VMIN_ABS   = 1e-5
SIGMA_WP   = 5.0
Z_LAUNCH   = -40.0
CAP_ONSET  = 60.0

CAP_COLOR  = "#A8D4E0"   # pastel teal-cyan — complementary to pinkish-orange
CAP_ALPHA  = 0.28

LZ_SLABS = [
    (r"$L_z = 15$ Bohr",  7.5, "#5B9BD5"),
    (r"$L_z = 25$ Bohr", 12.5, "#43A047"),
    (r"$L_z = 35$ Bohr", 17.5, "#ED7D31"),
]


def _log_cbar(fig, ax, im, label):
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size="4%", pad=0.08)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return cb


def build() -> None:
    v = load_vti(VTI)
    ny = v.data.shape[1]
    n_gs_xz = v.data[:, ny // 2, :]

    XX_g, ZZ_g = np.meshgrid(v.x, v.z, indexing="ij")
    sigma_d = SIGMA_WP / np.sqrt(2.0)
    norm3d   = (2.0 * np.pi * sigma_d**2) ** (-1.5)
    n_wp_xz  = norm3d * np.exp(-(XX_g**2 + (ZZ_g - Z_LAUNCH)**2) / (2.0 * sigma_d**2))
    n_total  = np.maximum(n_gs_xz + n_wp_xz, VMIN_ABS * 0.5)

    vmax    = float(np.percentile(n_total[n_total > VMIN_ABS], 99.8))
    lognorm = LogNorm(vmin=VMIN_ABS, vmax=vmax)

    fig, ax = plt.subplots(figsize=(style.TWO_COL_W_IN, 2.2))

    im = ax.pcolormesh(v.z, v.x, n_total,
                       cmap=style.cmap_for("density"),
                       norm=lognorm, rasterized=True, shading="auto")
    _log_cbar(fig, ax, im, r"$n$ (Bohr$^{-3}$)")

    # CAP regions — pastel overlay on both sides
    zmin, zmax = v.z.min(), v.z.max()
    ax.axvspan(zmin,      -CAP_ONSET, color=CAP_COLOR, alpha=CAP_ALPHA, zorder=2, lw=0)
    ax.axvspan( CAP_ONSET, zmax,      color=CAP_COLOR, alpha=CAP_ALPHA, zorder=2, lw=0)

    # CAP boundary lines (on top of overlay)
    for sign in [+1, -1]:
        ax.axvline(sign * CAP_ONSET, color="white", lw=0.9, ls=":",
                   alpha=0.85, zorder=3, label="CAP" if sign == 1 else None)

    # slab-width indicators
    for label, half_w, col in LZ_SLABS:
        for sign in [+1, -1]:
            ax.axvline(sign * half_w, color=col, lw=1.2, ls="--", alpha=0.85,
                       zorder=3, label=label if sign == 1 else None)

    ax.set_xlabel(r"$z$ (Bohr)", fontsize=9)
    ax.set_ylabel(r"$x$ (Bohr)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.75)

    out = HERE / "setup_jellium_slab_cap.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"written: {out}")


if __name__ == "__main__":
    build()
