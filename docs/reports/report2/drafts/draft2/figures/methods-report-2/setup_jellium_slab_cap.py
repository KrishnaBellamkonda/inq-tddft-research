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
from matplotlib.patches import Patch
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

CAP_COLOR  = "#C8B8E8"   # soft lavender — neutral against inferno warm tones
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

    fig, ax = plt.subplots(figsize=(style.TWO_COL_W_IN, 2.64))
    fig.subplots_adjust(left=0.07, right=0.97, bottom=0.11, top=0.97)

    ax.pcolormesh(v.z, v.x, n_total,
                  cmap=style.cmap_for("density"),
                  norm=lognorm, rasterized=True, shading="auto")

    # CAP regions — pastel overlay on both sides
    zmin, zmax = v.z.min(), v.z.max()
    ax.axvspan(zmin,      -CAP_ONSET, color=CAP_COLOR, alpha=CAP_ALPHA, zorder=2, lw=0)
    ax.axvspan( CAP_ONSET, zmax,      color=CAP_COLOR, alpha=CAP_ALPHA, zorder=2, lw=0)

    # "CAP" text centred in each overlay region
    xmin_ax, xmax_ax = v.x.min(), v.x.max()
    y_text = xmin_ax * 0.78          # near the bottom of the x-axis range
    for z_mid in [(zmin - CAP_ONSET) / 2, (CAP_ONSET + zmax) / 2]:
        ax.text(z_mid, y_text, "CAP", color="white", fontsize=7.5,
                ha="center", va="center", fontweight="bold",
                alpha=0.90, zorder=4)

    # slab-width indicators
    for label, half_w, col in LZ_SLABS:
        for sign in [+1, -1]:
            ax.axvline(sign * half_w, color=col, lw=1.2, ls="--", alpha=0.85,
                       zorder=3, label=label if sign == 1 else None)

    ax.set_xlabel(r"$z$ (Bohr)", fontsize=9)
    ax.set_ylabel(r"$x$ (Bohr)", fontsize=9)
    ax.tick_params(labelsize=8)

    density_proxy = Patch(facecolor=plt.cm.inferno(0.97), edgecolor="none",
                          label=r"Electron density $n(x,z)$")
    handles, labels = ax.get_legend_handles_labels()
    ax.legend([density_proxy] + handles,
              [density_proxy.get_label()] + labels,
              fontsize=7, loc="upper right", framealpha=0.75)

    out = HERE / "setup_jellium_slab_cap.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"written: {out}")


if __name__ == "__main__":
    build()
