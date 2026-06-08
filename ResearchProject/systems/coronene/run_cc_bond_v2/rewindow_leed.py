#!/usr/bin/env python3
"""Recompute LEED screen accumulations with 5σ interference-free windows.

Backscattering screens (z > 0):
  t_start = (z_WP_init - z_screen + 5σ) / v₀  (initial WP 5σ past screen)
  t_end   = t_mol + (z_screen - 5σ) / v₀       (reflected WP 5σ before screen)

Transmission screens (z < 0):
  t_start = 0
  t_end   = (z_WP_init + Lz/2 - 5σ) / v₀       (WP 5σ from cell boundary)

Loads instantaneous screen snapshots, selects those within the window,
averages, and plots the re-windowed LEED patterns.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

RUN_DIR = Path(__file__).resolve().parent
RESULTS = RUN_DIR / "results"
SNAP_DIR = RESULTS / "raw" / "screens" / "instantaneous"
OUT_DIR = RESULTS / "analysis" / "screens_5sigma"

# --- Run parameters ---
Z_WP_INIT = 12.0        # WP launch z (Bohr)
V0 = 3.834              # WP speed (Bohr/a.u.), |k₀|
SIGMA = 1.0016           # WP width (Bohr)
LX = 34.771              # cell x-dimension (Bohr)
LY = 34.771              # cell y-dimension (Bohr)
LZ = 59.904              # cell z-dimension (Bohr)
DT = 0.02                # timestep (a.u.)
SCREEN_SNAP_EVERY = 30   # steps between snapshots
N_SIGMA = 5              # clearance in σ units

T_MOL = Z_WP_INIT / V0   # WP reaches molecule at z=0
CLEARANCE = N_SIGMA * SIGMA


def load_screen_config() -> list[dict]:
    cfg_path = RESULTS / "raw" / "screens" / "screen_config.csv"
    screens = []
    with open(cfg_path) as f:
        header = f.readline().strip().split(",")
        for line in f:
            parts = line.strip().split(",")
            screens.append({
                "index": int(parts[0]),
                "z": float(parts[1]),
                "label": parts[2],
                "kind": parts[3],
            })
    return screens


def compute_5sigma_window(screen: dict) -> tuple[float, float]:
    z_s = screen["z"]
    kind = screen["kind"]

    if kind == "forward":
        t_start = 0.0
        t_end = (Z_WP_INIT + LZ / 2 - CLEARANCE) / V0
    else:
        if z_s >= Z_WP_INIT:
            t_start = 0.0
        else:
            t_start = (Z_WP_INIT - z_s + CLEARANCE) / V0

        t_end = T_MOL + (z_s - CLEARANCE) / V0
        if t_end < 0:
            t_end = 0.0

    return (max(t_start, 0.0), t_end)


def load_snapshot(screen_idx: int, step: int) -> np.ndarray | None:
    fname = SNAP_DIR / f"screen_{screen_idx:02d}_t{step:06d}.dat"
    if not fname.exists():
        return None
    lines = fname.read_text().strip().split("\n")
    data = []
    for line in lines:
        if line.startswith("#"):
            continue
        vals = [float(x) for x in line.split()]
        data.append(vals)
    arr = np.array(data)
    return np.fft.fftshift(arr)


def get_snapshot_times() -> list[tuple[int, float]]:
    files = sorted(SNAP_DIR.glob("screen_00_t*.dat"))
    result = []
    for f in files:
        m = re.search(r"t(\d+)", f.stem)
        if m:
            step = int(m.group(1))
            result.append((step, step * DT))
    return result


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    screens = load_screen_config()
    snap_times = get_snapshot_times()

    print(f"=== LEED Re-windowing with {N_SIGMA}σ clearance ===")
    print(f"  σ = {SIGMA:.4f} Bohr, clearance = {CLEARANCE:.3f} Bohr")
    print(f"  v₀ = {V0:.4f}, T_mol = {T_MOL:.3f} a.u.")
    print(f"  {len(snap_times)} time snapshots available")
    print()

    all_patterns = {}
    window_info = []

    for sc in screens:
        idx = sc["index"]
        t_start, t_end = compute_5sigma_window(sc)
        valid = t_end > t_start

        selected = [(step, t) for step, t in snap_times
                    if t >= t_start and t <= t_end]

        window_info.append({
            **sc, "t_start": t_start, "t_end": t_end,
            "n_snaps": len(selected), "valid": valid,
        })

        status = f"VALID ({len(selected)} snaps)" if valid and selected else "EMPTY"
        print(f"  Screen {idx:2d} (z={sc['z']:+8.3f}, {sc['kind']:>7s}): "
              f"t ∈ [{t_start:.2f}, {t_end:.2f}] a.u. — {status}")

        if not selected:
            all_patterns[idx] = None
            continue

        acc = None
        for step, t in selected:
            snap = load_snapshot(idx, step)
            if snap is None:
                continue
            if acc is None:
                acc = snap.copy()
            else:
                acc += snap
        if acc is not None:
            acc /= len(selected)
        all_patterns[idx] = acc

    # --- Write window summary ---
    with open(OUT_DIR / "window_summary.csv", "w") as f:
        f.write("screen,z_bohr,kind,t_start_au,t_end_au,n_snapshots,valid\n")
        for w in window_info:
            f.write(f"{w['index']},{w['z']:.4f},{w['kind']},"
                    f"{w['t_start']:.4f},{w['t_end']:.4f},"
                    f"{w['n_snaps']},{w['valid']}\n")
    print(f"\n  Saved: {OUT_DIR / 'window_summary.csv'}")

    # --- Plot individual screens ---
    for idx, pattern in all_patterns.items():
        if pattern is None:
            continue
        sc = screens[idx]
        w = window_info[idx]

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Physical extent: after fftshift, array (0,0) = (-L/2, -L/2)
        extent = [-LX/2, LX/2, -LY/2, LY/2]

        im1 = ax1.imshow(pattern, origin="lower", cmap="RdBu_r",
                         norm=mcolors.CenteredNorm(), extent=extent)
        fig.colorbar(im1, ax=ax1, label="Density")
        ax1.set_title(f"Screen {idx} (z={sc['z']:+.1f} Bohr, {sc['kind']})\n"
                      f"Linear scale, {w['n_snaps']} snapshots")

        pos = pattern.copy()
        pos[pos <= 0] = np.nan
        im2 = ax2.imshow(pos, origin="lower", cmap="inferno",
                         norm=mcolors.LogNorm(vmin=np.nanpercentile(pos, 5),
                                              vmax=np.nanpercentile(pos, 99)),
                         extent=extent)
        fig.colorbar(im2, ax=ax2, label="Density (log)")
        ax2.set_title(f"Log scale\nt ∈ [{w['t_start']:.2f}, {w['t_end']:.2f}] a.u.")

        for ax in (ax1, ax2):
            ax.set_xlabel("x / Bohr")
            ax.set_ylabel("y / Bohr")

        fig.suptitle(f"run_cc_bond: {N_SIGMA}σ IFW LEED — {sc['label']}", fontsize=12)
        fig.tight_layout()
        fig.savefig(OUT_DIR / f"screen_{idx:02d}_5sigma.png", dpi=150)
        plt.close(fig)

    # --- Grid plot: all backscattering screens ---
    back_screens = [i for i, sc in enumerate(screens)
                    if sc["kind"] == "back" and all_patterns.get(i) is not None]
    if back_screens:
        n_cols = min(5, len(back_screens))
        n_rows = (len(back_screens) + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
        if n_rows == 1:
            axes = [axes] if n_cols == 1 else list(axes)
        else:
            axes = [ax for row in axes for ax in row]

        for i, idx in enumerate(back_screens):
            ax = axes[i]
            p = all_patterns[idx]
            w = window_info[idx]
            pos = p.copy()
            pos[pos <= 0] = np.nan
            extent = [-LX/2, LX/2, -LY/2, LY/2]
            if np.any(np.isfinite(pos) & (pos > 0)):
                im = ax.imshow(pos, origin="lower", cmap="inferno",
                               norm=mcolors.LogNorm(
                                   vmin=np.nanpercentile(pos[pos > 0], 5),
                                   vmax=np.nanpercentile(pos[pos > 0], 99)),
                               extent=extent)
            else:
                im = ax.imshow(p, origin="lower", cmap="RdBu_r", extent=extent)
            ax.set_title(f"S{idx} z={screens[idx]['z']:+.1f}\n"
                         f"[{w['t_start']:.1f},{w['t_end']:.1f}] ({w['n_snaps']}sn)",
                         fontsize=9)
            ax.set_xlabel("x / Bohr", fontsize=7)
            ax.set_ylabel("y / Bohr", fontsize=7)
            ax.tick_params(labelsize=6)

        for i in range(len(back_screens), len(axes)):
            axes[i].set_visible(False)

        fig.suptitle(f"run_cc_bond: backscattering LEED — {N_SIGMA}σ IFW windows", fontsize=13)
        fig.tight_layout()
        fig.savefig(OUT_DIR / "backscattering_grid_5sigma.png", dpi=150)
        plt.close(fig)
        print(f"  Saved: backscattering_grid_5sigma.png")

    # --- Grid plot: all transmission screens ---
    fwd_screens = [i for i, sc in enumerate(screens)
                   if sc["kind"] == "forward" and all_patterns.get(i) is not None]
    if fwd_screens:
        n_cols = min(5, len(fwd_screens))
        n_rows = (len(fwd_screens) + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
        if n_rows == 1:
            axes = [axes] if n_cols == 1 else list(axes)
        else:
            axes = [ax for row in axes for ax in row]

        for i, idx in enumerate(fwd_screens):
            ax = axes[i]
            p = all_patterns[idx]
            w = window_info[idx]
            pos = p.copy()
            pos[pos <= 0] = np.nan
            if np.any(np.isfinite(pos) & (pos > 0)):
                im = ax.imshow(pos, origin="lower", cmap="inferno",
                               norm=mcolors.LogNorm(
                                   vmin=np.nanpercentile(pos[pos > 0], 5),
                                   vmax=np.nanpercentile(pos[pos > 0], 99)))
            else:
                im = ax.imshow(p, origin="lower", cmap="RdBu_r")
            ax.set_title(f"S{idx} z={screens[idx]['z']:+.1f}\n"
                         f"[{w['t_start']:.1f},{w['t_end']:.1f}] ({w['n_snaps']}sn)",
                         fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])

        for i in range(len(fwd_screens), len(axes)):
            axes[i].set_visible(False)

        fig.suptitle(f"run_cc_bond: transmission LEED — {N_SIGMA}σ IFW windows", fontsize=13)
        fig.tight_layout()
        fig.savefig(OUT_DIR / "transmission_grid_5sigma.png", dpi=150)
        plt.close(fig)
        print(f"  Saved: transmission_grid_5sigma.png")

    print("\nDone.")


if __name__ == "__main__":
    main()
