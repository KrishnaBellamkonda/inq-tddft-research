#!/usr/bin/env python3
"""
04_leed_simulation/analysis.py
===============================
Post-processing for the coronene TDDFT WP scattering simulation.

Produces:
  results/fig1_density_snapshots.png  — density on flake plane at several times
                                        (replicates Fig. 1 of Tsubonoya et al.)
  results/fig2_leed_pattern.png       — I(x,y) = integrated density on obs. plane
                                        (replicates Fig. 2 of Tsubonoya et al.)

Usage:
    cd 04_leed_simulation
    python3 analysis.py

Requirements: numpy, matplotlib
"""

import os, glob, re
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

RESULTS = "results"

# Cell dimensions (bohr) — must match config.hpp
LX_BOHR = 34.76
LY_BOHR = 34.76


# ─────────────────────────────────────────────────────────────────────────────
# Helper: read a 2D density slice file
# ─────────────────────────────────────────────────────────────────────────────
def read_slice(path):
    """Load a 2D density slice.  Returns (array[Ny,Nx], time_au, z_bohr)."""
    with open(path) as f:
        header = f.readline()   # # t=... z=...
    m_t = re.search(r"t=([0-9.e+\-]+)", header)
    m_z = re.search(r"z=([0-9.e+\-]+)", header)
    t_au   = float(m_t.group(1)) if m_t else 0.0
    z_bohr = float(m_z.group(1)) if m_z else 0.0
    data = np.loadtxt(path, comments="#")
    return data, t_au, z_bohr


# ─────────────────────────────────────────────────────────────────────────────
# Fig. 1 replica: density snapshots on the coronene plane (z=0)
# ─────────────────────────────────────────────────────────────────────────────
def plot_fig1_snapshots():
    snap_files = sorted(glob.glob(os.path.join(RESULTS, "snapshot_t*.txt")))
    if not snap_files:
        print("  No snapshot files found — skipping Fig. 1.")
        return

    n = min(len(snap_files), 6)           # show up to 6 frames
    indices = np.linspace(0, len(snap_files)-1, n, dtype=int)

    fig, axes = plt.subplots(1, n, figsize=(3.5*n, 3.5))
    if n == 1:
        axes = [axes]

    vmax = None
    for idx in indices:
        d, _, _ = read_slice(snap_files[idx])
        vmax = max(vmax or 0, d.max())

    BOHR_TO_ANG = 0.529177
    extent = [0, LX_BOHR*BOHR_TO_ANG, 0, LY_BOHR*BOHR_TO_ANG]

    for ax, idx in zip(axes, indices):
        data, t_au, _ = read_slice(snap_files[idx])
        im = ax.imshow(data, origin="lower", vmin=0, vmax=vmax,
                       cmap="inferno", extent=extent, aspect="equal")
        ax.set_title(f"t = {t_au:.2f} a.u.", fontsize=10)
        ax.set_xlabel("x (Å)")
        ax.set_ylabel("y (Å)")

    plt.colorbar(im, ax=axes[-1], label=r"$n(x,y,0,t)$ (a.u.)")
    fig.suptitle("Coronene: electron density on flake plane (z=0)\nFig. 1 replica",
                 fontsize=12)
    plt.tight_layout()
    out = os.path.join(RESULTS, "fig1_density_snapshots.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Fig. 2 replica: LEED pattern I(x,y) on observation plane
# ─────────────────────────────────────────────────────────────────────────────
def plot_fig2_leed():
    leed_file = os.path.join(RESULTS, "leed_pattern.txt")
    if not os.path.exists(leed_file):
        print("  leed_pattern.txt not found — skipping Fig. 2.")
        return

    data = np.loadtxt(leed_file, comments="#")
    Ny, Nx = data.shape

    BOHR_TO_ANG = 0.529177
    extent = [0, LX_BOHR*BOHR_TO_ANG, 0, LY_BOHR*BOHR_TO_ANG]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))

    # Linear scale
    im1 = axes[0].imshow(data, origin="lower", cmap="hot", extent=extent)
    axes[0].set_title("LEED pattern — linear scale")
    axes[0].set_xlabel("x (Å)"); axes[0].set_ylabel("y (Å)")
    plt.colorbar(im1, ax=axes[0], label=r"$I(x,y)$ (a.u.·time)")

    # Log scale (better for diffraction spot visibility)
    data_nz = np.where(data > 0, data, data[data > 0].min() if data.any() else 1e-30)
    im2 = axes[1].imshow(data_nz, origin="lower", cmap="hot",
                          norm=LogNorm(), extent=extent)
    axes[1].set_title("LEED pattern — log scale")
    axes[1].set_xlabel("x (Å)"); axes[1].set_ylabel("y (Å)")
    plt.colorbar(im2, ax=axes[1], label=r"$I(x,y)$ (log)")

    fig.suptitle("Coronene: LEED diffraction pattern\nFig. 2 replica"
                 " (expect 6-fold D6h symmetry)", fontsize=12)
    plt.tight_layout()
    out = os.path.join(RESULTS, "fig2_leed_pattern.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
# Symmetry check: compute angular autocorrelation of LEED pattern
# (6-fold symmetry → peaks at 60°, 120°, 180°, 240°, 300°)
# ─────────────────────────────────────────────────────────────────────────────
def check_symmetry():
    leed_file = os.path.join(RESULTS, "leed_pattern.txt")
    if not os.path.exists(leed_file):
        return

    data = np.loadtxt(leed_file, comments="#")
    Ny, Nx = data.shape
    cy, cx = Ny // 2, Nx // 2

    # Build radial-averaged angular profile
    angles = np.linspace(0, 2*np.pi, 360, endpoint=False)
    r_range = np.arange(5, min(cx, cy))  # skip centre
    profile = np.zeros(360)
    for ir, r in enumerate(r_range):
        for ia, ang in enumerate(angles):
            xi = cx + r * np.cos(ang)
            yi = cy + r * np.sin(ang)
            xi, yi = int(round(xi)), int(round(yi))
            if 0 <= xi < Nx and 0 <= yi < Ny:
                profile[ia] += data[yi, xi]

    fig, ax = plt.subplots(figsize=(6, 3))
    ax.plot(np.degrees(angles), profile / profile.max())
    ax.set_xlabel("Angle (degrees)")
    ax.set_ylabel("Normalised radial intensity")
    ax.set_title("Angular profile of LEED pattern\n(6-fold symmetry → peaks at 0°,60°,120°,…)")
    ax.axhline(0.5, color="gray", linestyle="--", alpha=0.5)
    plt.tight_layout()
    out = os.path.join(RESULTS, "leed_angular_profile.png")
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out}")


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    print("=== Fig. 1: density snapshots ===")
    plot_fig1_snapshots()

    print("\n=== Fig. 2: LEED pattern ===")
    plot_fig2_leed()

    print("\n=== Symmetry check ===")
    check_symmetry()

    print("\nAll done.")
