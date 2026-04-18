#!/usr/bin/env python3
"""
make_video.py
─────────────
Produces a publication-quality MP4 video of the Gaussian wave packet
spreading in a periodic box.

Each frame shows:
  Left panel  — 2D colour map of ρ(x,y) at z = L/2  (midplane)
  Right panel — 1D density profile ρ(x) through the midline, with the
                analytical Gaussian envelope overlaid
  Inset       — σ(t) accumulated progress vs. analytical curve

References
----------
[1] Angelo, "Quantum Wave-Packet Preparation and Electron Dynamics in
    Jellium", Cavendish Lab Report, Candidate 3221L.
    Eq. (6): σ(t) = σ₀ √(1 + t²/(4σ₀⁴))
    Fig. 4 : representative snapshots of free propagation.

Usage
-----
    cd ResearchProject/systems/jellium
    python3 analysis/make_video.py

Requirements
------------
    pip install matplotlib numpy
    sudo apt install ffmpeg    # or conda install -c conda-forge ffmpeg

Output
------
    results/gaussian_spreading.mp4
"""

import os
import sys
import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.animation import FFMpegWriter
from matplotlib.colors import LogNorm
import warnings

# ─── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR  = os.path.join(SCRIPT_DIR, "..", "results")
SNAP_DIR     = os.path.join(RESULTS_DIR, "snapshots")
SIGMA_FILE   = os.path.join(RESULTS_DIR, "sigma_t.txt")
GRID_FILE    = os.path.join(RESULTS_DIR, "grid_info.txt")
OUTPUT_VIDEO = os.path.join(RESULTS_DIR, "gaussian_spreading.mp4")

FPS = 10   # frames per second (1 frame = 1 atu; 50 frames → 5 s video)

# ─── Helpers ──────────────────────────────────────────────────────────────────
def read_grid_info(path):
    info = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                try:
                    info[parts[0]] = float(parts[1])
                except (IndexError, ValueError):
                    pass
    except FileNotFoundError:
        pass
    return info


def load_snapshots(snap_dir):
    """Return sorted list of (t_index, slice_path, line_path)."""
    slices = sorted(glob.glob(os.path.join(snap_dir, "slice_*.npy")))
    snaps  = []
    for s in slices:
        idx = int(re.search(r"slice_(\d+)", s).group(1))
        ln  = os.path.join(snap_dir, f"line_{idx:04d}.npy")
        snaps.append((idx, s, ln))
    return snaps


# ─── Main ─────────────────────────────────────────────────────────────────────
def main():
    # Check prerequisites
    for req in [SIGMA_FILE, SNAP_DIR]:
        if not os.path.exists(req):
            print(f"ERROR: {req} not found. Run the simulation first.")
            sys.exit(1)

    # Load parameters
    info      = read_grid_info(GRID_FILE)
    L         = info.get("L",       40.0)
    sigma_wp  = info.get("sigma_wp", 1.0)
    k0        = info.get("k0",       1.5)
    dt        = info.get("dt",       0.04)
    snap_ev   = int(info.get("snap_every", 25))
    T_rev     = L**2 / np.pi

    # Load σ(t) data for inset
    sig_data  = np.loadtxt(SIGMA_FILE, comments="#")
    if sig_data.ndim == 1:
        sig_data = sig_data[np.newaxis, :]
    t_all     = sig_data[:, 0]
    snum_all  = sig_data[:, 1]
    sana_all  = sig_data[:, 2]

    # Load snapshot file list
    snaps = load_snapshots(SNAP_DIR)
    if not snaps:
        print("ERROR: no snapshot files found in", SNAP_DIR)
        sys.exit(1)
    n_frames = len(snaps)
    print(f"Found {n_frames} frames.")

    # Probe first frame for grid size
    idx0, s0, l0 = snaps[0]
    slice0 = np.load(s0)   # shape (ny, nx)
    line0  = np.load(l0)   # shape (nx,)
    ny, nx = slice0.shape

    # Spatial axes [bohr]
    x_axis = np.linspace(0, L, nx, endpoint=False)
    y_axis = np.linspace(0, L, ny, endpoint=False)
    dt_snap = snap_ev * dt   # time per frame [atu]

    # Colour scale: use 1st and last frame to set vmax/vmin
    slice_last = np.load(snaps[-1][1])
    rho_max    = float(slice0.max())      # initial peak
    rho_floor  = max(slice0.min(), 1e-8 * rho_max)   # floor for log scale

    # Analytical Gaussian envelope at t=0 for reference
    def gaussian_envelope(x, t):
        """Analytical |ψ(x)|² at cell midline y=L/2, z=L/2."""
        sig = sigma_wp * np.sqrt(1.0 + t**2 / (4.0 * sigma_wp**4))
        # centre moves at k0 (mod L)
        cx  = (L/2.0 + k0 * t) % L
        # minimum-image x distance
        dx  = x - cx
        dx -= L * np.round(dx / L)
        norm = 1.0 / (np.sqrt(2.0 * np.pi) * sig)
        return norm * np.exp(-dx**2 / (2.0 * sig**2))

    # ── Set up figure layout ──────────────────────────────────────────────────
    fig = plt.figure(figsize=(14, 6))
    gs  = gridspec.GridSpec(1, 2, figure=fig, wspace=0.35)

    # Left: 2D midplane colour map
    ax_2d = fig.add_subplot(gs[0, 0])
    # Right: 1D line profile
    ax_1d = fig.add_subplot(gs[0, 1])
    # Inset: σ(t) progress (top-right corner of 1D panel)
    ax_ins = ax_1d.inset_axes([0.55, 0.55, 0.43, 0.42])

    # Initial frame
    t0 = 0.0
    im = ax_2d.imshow(
        slice0,
        origin="lower",
        extent=[0, L, 0, L],
        cmap="inferno",
        norm=LogNorm(vmin=rho_floor, vmax=rho_max),
        aspect="equal",
    )
    cbar = fig.colorbar(im, ax=ax_2d, label=r"$\rho$ [a$_0^{-3}$]", fraction=0.046, pad=0.04)

    ax_2d.set_xlabel("$x$ [bohr]", fontsize=11)
    ax_2d.set_ylabel("$y$ [bohr]", fontsize=11)
    title_2d = ax_2d.set_title(
        r"$\rho(x,y,\,z{=}L/2)$" + f",  $t = {t0:.2f}$ atu", fontsize=11
    )

    # 1D line plot
    (line_plot,)   = ax_1d.plot(x_axis, line0,             "b-",  lw=1.5,
                                 label=r"$\rho(x,\,y{=}z{=}L/2)$ [INQ]")
    (env_plot,)    = ax_1d.plot(x_axis, gaussian_envelope(x_axis, t0), "r--",
                                 lw=1.5, label="Analytical envelope")
    ax_1d.set_xlabel("$x$ [bohr]", fontsize=11)
    ax_1d.set_ylabel(r"$\rho$ [a$_0^{-3}$]", fontsize=11)
    ax_1d.set_title(r"$\rho(x)$ through cell centre", fontsize=11)
    ax_1d.set_xlim(0, L)
    ax_1d.legend(fontsize=8, loc="upper left")
    ax_1d.grid(True, alpha=0.3)

    # σ(t) inset
    t_ref = np.linspace(0, t_all[-1] * 1.05, 300)
    s_ref = sigma_wp * np.sqrt(1.0 + t_ref**2 / (4.0 * sigma_wp**4))
    ax_ins.plot(t_ref, s_ref, "r-", lw=1, label="Analytical")
    (scat_ins,)   = ax_ins.plot([], [], "b.", ms=4, label="INQ")
    (vline_ins,)  = ax_ins.plot([], [], "k-", lw=1)
    ax_ins.set_xlabel("$t$ [atu]", fontsize=7)
    ax_ins.set_ylabel(r"$\sigma$ [a$_0$]", fontsize=7)
    ax_ins.tick_params(labelsize=7)
    ax_ins.legend(fontsize=6)
    ax_ins.set_xlim(0, t_all[-1] * 1.05)
    ax_ins.set_ylim(0, sana_all[-1] * 1.1)
    ax_ins.grid(True, alpha=0.25)

    plt.tight_layout()

    # ── Animation update function ─────────────────────────────────────────────
    def update(frame_idx):
        idx, s_path, l_path = snaps[frame_idx]
        t_now = idx * dt_snap     # time for this snapshot index

        # Map snapshot index → sigma data row
        row_idx = min(frame_idx, len(t_all) - 1)
        t_now   = t_all[row_idx]

        slice_arr = np.load(s_path)
        line_arr  = np.load(l_path)

        # Update 2D colour map
        im.set_data(slice_arr)
        title_2d.set_text(r"$\rho(x,y,\,z{=}L/2)$" + f",  $t = {t_now:.2f}$ atu")

        # Update colour scale to current frame's range (keep floor, adjust top)
        cur_max = max(float(slice_arr.max()), rho_floor * 2)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            im.set_norm(LogNorm(vmin=rho_floor, vmax=cur_max))

        # Update 1D profile
        line_plot.set_ydata(line_arr)
        env_plot.set_ydata(gaussian_envelope(x_axis, t_now))
        y_top = max(float(line_arr.max()) * 1.2, gaussian_envelope(x_axis, t_now).max() * 1.2)
        ax_1d.set_ylim(0, y_top)

        # Update σ(t) inset
        scat_ins.set_data(t_all[:row_idx+1], snum_all[:row_idx+1])
        vline_ins.set_data([t_now, t_now],
                           [0, sana_all[row_idx] * 1.1])

        return [im, title_2d, line_plot, env_plot, scat_ins, vline_ins]

    # ── Build animation and encode ─────────────────────────────────────────────
    from matplotlib.animation import FuncAnimation

    anim = FuncAnimation(
        fig, update,
        frames=n_frames,
        interval=1000 // FPS,
        blit=True,
    )

    writer = FFMpegWriter(
        fps=FPS,
        bitrate=3000,
        extra_args=[
            "-vcodec", "libx264",
            "-pix_fmt", "yuv420p",
            "-crf", "18",
        ],
    )

    print(f"Encoding {n_frames} frames at {FPS} fps → {OUTPUT_VIDEO}")
    anim.save(OUTPUT_VIDEO, writer=writer, dpi=120)
    plt.close(fig)
    print(f"Done. Video: {OUTPUT_VIDEO}")
    print(f"Duration: {n_frames/FPS:.1f} s  ({n_frames} frames)")


if __name__ == "__main__":
    main()
