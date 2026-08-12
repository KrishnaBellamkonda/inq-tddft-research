"""Draft-2 setup panels: jellium slab + cylindrical jellium setup figures.

Produces three PNGs in this directory at the theme's default one-col size
(3.5 × 3.0 in, 600 DPI, bbox_inches=None):
  setup_jellium_slab.png
  setup_cylindrical_jellium_sweep.png
  setup_cylindrical_jellium.png

Imports all logic from the draft-1 source; only output paths and save
parameters differ.

Run:
  /local/data/public/skcb2/tddft/venv/bin/python3 build_setup_panels.py
"""
from __future__ import annotations

import os, sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[7]   # tddft root
sys.path.insert(0, str(REPO / "inq-stack/python"))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.special import erfc
from inqview.visualisation import style
from inqview import load_vti

style.apply_theme()

HERE = Path(__file__).parent


# ── helpers (identical to draft-1 source) ────────────────────────────────────

def n0_from_rs(rs: float) -> float:
    return 3.0 / (4.0 * np.pi * rs**3)


def annular_bg(X, Z, n0, r_in, r_out, w):
    d = np.abs(X)
    return n0 * 0.5 * erfc((d - r_out) / w) * 0.5 * erfc((r_in - d) / w)


def annular_bg_xy(X, Y, n0, r_in, r_out, w):
    d = np.sqrt(X**2 + Y**2)
    return n0 * 0.5 * erfc((d - r_out) / w) * 0.5 * erfc((r_in - d) / w)


def gaussian_wp(X, Z, z0, sigma_wp):
    sigma_d = sigma_wp / np.sqrt(2.0)
    norm = (2.0 * np.pi * sigma_d**2) ** (-1.5)
    return norm * np.exp(-(X**2 + (Z - z0)**2) / (2.0 * sigma_d**2))


def _log_cbar(fig, ax, im, label):
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size="4%", pad=0.08)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(label, fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return cb


def _log_norm(data, n0, floor_frac=1e-4):
    vmin = n0 * floor_frac
    vmax = float(np.percentile(data[data > vmin], 99.8))
    return LogNorm(vmin=vmin, vmax=vmax)


# ── figure 3: jellium slab setup ─────────────────────────────────────────────

def build_slab() -> None:
    VTI = (
        "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium"
        "/shared_gs/slab_n102_L25x25x140_w0p5_h0p5/density_gs_system/density_gs_system.vti"
    )
    VMIN_ABS = 1e-5
    SIGMA_WP = 5.0
    Z_LAUNCH = -40.0
    CAP_ONSET = 60.0
    LZ_SLABS = [
        (r"$L_z = 15$ Bohr",  7.5, "#5B9BD5"),
        (r"$L_z = 25$ Bohr", 12.5, "#43A047"),
        (r"$L_z = 35$ Bohr", 17.5, "#ED7D31"),
    ]

    v = load_vti(VTI)
    ny = v.data.shape[1]
    n_gs_xz = v.data[:, ny // 2, :]

    XX_g, ZZ_g = np.meshgrid(v.x, v.z, indexing="ij")
    sigma_d = SIGMA_WP / np.sqrt(2.0)
    norm3d = (2.0 * np.pi * sigma_d**2) ** (-1.5)
    n_wp_xz = norm3d * np.exp(-(XX_g**2 + (ZZ_g - Z_LAUNCH)**2) / (2.0 * sigma_d**2))
    n_total_xz = np.maximum(n_gs_xz + n_wp_xz, VMIN_ABS * 0.5)

    vmax = float(np.percentile(n_total_xz[n_total_xz > VMIN_ABS], 99.8))
    lognorm = LogNorm(vmin=VMIN_ABS, vmax=vmax)

    fig, ax = plt.subplots(figsize=(style.TWO_COL_W_IN, 2.2))
    im = ax.pcolormesh(v.z, v.x, n_total_xz,
                       cmap=style.cmap_for("density"),
                       norm=lognorm, rasterized=True, shading="auto")
    _log_cbar(fig, ax, im, r"$n$ (Bohr$^{-3}$)")

    for sign in [+1, -1]:
        ax.axvline(sign * CAP_ONSET, color="white", lw=0.9, ls=":", alpha=0.65,
                   label="CAP" if sign == 1 else None)
    for label, half_w, col in LZ_SLABS:
        for sign in [+1, -1]:
            ax.axvline(sign * half_w, color=col, lw=1.2, ls="--", alpha=0.85,
                       label=label if sign == 1 else None)

    ax.set_xlabel(r"$z$ (Bohr)", fontsize=9)
    ax.set_ylabel(r"$x$ (Bohr)", fontsize=9)
    ax.tick_params(labelsize=8)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.75)

    out = HERE / "setup_jellium_slab.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"written: {out}")


# ── figure 4: cylindrical jellium sweep ──────────────────────────────────────

def build_cylindrical_sweep() -> None:
    R_S = 3.0
    W = 0.5
    LADDER = [
        {"r_in": 10, "r_out": 14.000, "label": r"$R_{\rm in} = 10$"},
        {"r_in":  8, "r_out": 14.000, "label": r"$R_{\rm in} = 8$"},
        {"r_in":  6, "r_out": 13.986, "label": r"$R_{\rm in} = 6$"},
        {"r_in":  4, "r_out": 14.000, "label": r"$R_{\rm in} = 4$"},
        {"r_in":  0, "r_out": 13.986, "label": r"$R_{\rm in} = 0$ (filled)"},
    ]
    R_OUT_REP = 14.0
    PASTEL = ["#5BA4CF", "#63C174", "#F5A623", "#E86565", "#B79FDB"]

    n0 = n0_from_rs(R_S)
    LBOX = 35.0
    xy = np.linspace(-LBOX / 2, LBOX / 2, 600)
    XX, YY = np.meshgrid(xy, xy, indexing="ij")
    n_gs = annular_bg_xy(XX, YY, n0, r_in=10.0, r_out=R_OUT_REP, w=W)
    n_gs = np.maximum(n_gs, n0 * 1e-6)

    fig, ax = plt.subplots(figsize=(3.5, 3.5))
    ax.set_aspect("equal")
    norm = _log_norm(n_gs, n0)
    im = ax.pcolormesh(XX, YY, n_gs, cmap=style.cmap_for("density"),
                       norm=norm, rasterized=True, shading="auto")
    _log_cbar(fig, ax, im, r"$n$ (Bohr$^{-3}$)")

    theta = np.linspace(0, 2 * np.pi, 500)
    ax.plot(R_OUT_REP * np.cos(theta), R_OUT_REP * np.sin(theta),
            color="white", lw=1.0, ls=":", alpha=0.60)
    for rung, col in zip(LADDER, PASTEL):
        r_in = rung["r_in"]
        lbl = rung["label"]
        if r_in > 0:
            ax.plot(r_in * np.cos(theta), r_in * np.sin(theta),
                    color=col, lw=1.8, ls="--", label=lbl)
        else:
            ax.plot([], [], color=col, lw=1.8, ls="--", label=lbl)

    ax.plot(0, 0, "w+", markersize=9, markeredgewidth=1.6, zorder=6)
    ax.text(0.6, -1.5, "WP", color="white", fontsize=8,
            ha="left", va="top", fontweight="bold")
    ax.set_xlabel(r"$x$ (Bohr)", fontsize=9)
    ax.set_ylabel(r"$y$ (Bohr)", fontsize=9)
    ax.tick_params(labelsize=8)
    lim = LBOX / 2
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.55,
              title=r"proximity ladder", title_fontsize=7)

    out = HERE / "setup_cylindrical_jellium_sweep.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"written: {out}")


# ── figure 5: cylindrical jellium xz slice ───────────────────────────────────

def build_cylindrical_xz() -> None:
    R_S = 3.0; R_IN = 10.0; R_OUT = 14.0; LZ = 60.0; LXY = 40.0
    W_TUBE = 0.5; SIGMA_WP = 4.0; Z_LAUNCH = -28.0

    n0 = n0_from_rs(R_S)
    x = np.linspace(-LXY / 2, LXY / 2, 400)
    z = np.linspace(-LZ / 2, LZ / 2, 600)
    XX, ZZ = np.meshgrid(x, z, indexing="ij")
    n_bg = annular_bg(XX, ZZ, n0, R_IN, R_OUT, W_TUBE)
    n_wp = gaussian_wp(XX, ZZ, Z_LAUNCH, SIGMA_WP)
    n_total = np.maximum(n_bg + n_wp, n0 * 1e-6)

    fig, ax = plt.subplots(figsize=(3.8, 3.0))
    norm = _log_norm(n_total, n0)
    im = ax.pcolormesh(ZZ.T, XX.T, n_total.T,
                       cmap=style.cmap_for("density"),
                       norm=norm, rasterized=True, shading="auto")
    _log_cbar(fig, ax, im, r"$n$ (Bohr$^{-3}$)")

    for sign in [+1, -1]:
        ax.axhline(sign * R_IN,  color="white", lw=0.8, ls="--", alpha=0.75)
        ax.axhline(sign * R_OUT, color="white", lw=0.6, ls=":",  alpha=0.55)

    ax.set_xlabel(r"$z$ (Bohr)", fontsize=9)
    ax.set_ylabel(r"$x$ (Bohr)", fontsize=9)
    ax.tick_params(labelsize=8)

    out = HERE / "setup_cylindrical_jellium.png"
    fig.savefig(out, dpi=600, bbox_inches=None)
    plt.close(fig)
    print(f"written: {out}")


if __name__ == "__main__":
    build_slab()
    build_cylindrical_sweep()
    build_cylindrical_xz()
    print("All setup panels written to", HERE)
