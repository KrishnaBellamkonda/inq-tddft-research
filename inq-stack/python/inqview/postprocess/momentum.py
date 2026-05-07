"""Phase: ``momentum`` — visualise the on-the-fly momentum distribution.

Reads ``results/raw/observables/momentum_distribution.csv`` produced by
``inqkit::observables::MomentumDistribution`` (long-format CSV, one row
per (step, |k|-bin) with columns ``step,time_au,k_bohr_inv,n_total,n_wp``).
The file's first line is a ``# l_bohr=... n_bins=... wp_idx=...`` header.

Outputs:

* ``analysis/observables/momentum_distribution.gif``
   — n(|k|, t) curve animated in time, with the launch k_0 marked.
* ``analysis/observables/momentum_heatmap.png``
   — 2D heatmap (time vs |k|) of n_total(|k|, t), fixed colour scale.

The ``run_summary.txt`` is parsed for ``wp_k0_bohr_inv`` to draw the
analytic launch-k marker; if absent, the marker is suppressed.
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import _common


def _read_summary_k0(results_dir: Path) -> float | None:
    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return None
    text = rs.read_text()
    m = re.search(r"wp_k0_bohr_inv\s*=\s*\S+\s+\S+\s+(\S+)", text)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            return None
    return None


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    csv_path = results_dir / "raw" / "observables" / "momentum_distribution.csv"
    if not csv_path.exists():
        return {"skipped": f"missing: {csv_path}"}

    df = pd.read_csv(csv_path, comment="#")
    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    artefacts: list[str] = []

    k0 = _read_summary_k0(results_dir)

    times = sorted(df["time_au"].unique())
    k_vals = sorted(df["k_bohr_inv"].unique())
    nT, nK = len(times), len(k_vals)

    # Pivot: rows = time, cols = |k|. n_total only for the heatmap.
    grid_total = (df.pivot_table(index="time_au", columns="k_bohr_inv",
                                 values="n_total", aggfunc="sum")
                    .reindex(index=times, columns=k_vals).to_numpy())
    grid_wp = (df.pivot_table(index="time_au", columns="k_bohr_inv",
                              values="n_wp", aggfunc="sum")
                 .reindex(index=times, columns=k_vals).to_numpy())

    # Heatmap (n_total)
    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(grid_total.T, origin="lower",
                   extent=[times[0], times[-1], k_vals[0], k_vals[-1]],
                   aspect="auto", cmap="viridis")
    if k0 is not None:
        ax.axhline(k0, color="red", lw=0.8, ls="--", label=f"k0 = {k0:.3f}")
        ax.legend(loc="upper right")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("|k| (Bohr$^{-1}$)")
    ax.set_title(f"{run_name}: one-body momentum distribution n(|k|, t)")
    fig.colorbar(im, ax=ax, label="n_total")
    fig.tight_layout()
    out = out_dir / "momentum_heatmap.png"
    if rebuild or not out.exists():
        fig.savefig(out, dpi=130)
    plt.close(fig)
    artefacts.append(str(out))

    # Animated GIF: n(|k|, t) curve, fixed y-limits.
    vmax = float(np.nanmax([grid_total.max(), grid_wp.max()]))
    fig, ax = plt.subplots(figsize=(7, 4))
    line_total, = ax.plot([], [], lw=1.6, label="n_total")
    line_wp,    = ax.plot([], [], lw=1.2, color="C3", label="n_wp")
    ax.set_xlim(k_vals[0], k_vals[-1])
    ax.set_ylim(0, vmax * 1.05 if vmax > 0 else 1)
    ax.set_xlabel("|k| (Bohr$^{-1}$)")
    ax.set_ylabel("n(|k|)")
    if k0 is not None:
        ax.axvline(k0, color="red", lw=0.8, ls="--", label=f"k0 = {k0:.3f}")
    ax.legend()
    title = ax.set_title("")

    def _frame(i):
        line_total.set_data(k_vals, grid_total[i])
        line_wp.set_data(k_vals, grid_wp[i])
        title.set_text(f"{run_name}: t = {times[i]:.3f} a.u.")
        return line_total, line_wp, title

    anim = animation.FuncAnimation(fig, _frame, frames=nT,
                                    interval=80, blit=False)
    out = out_dir / "momentum_distribution.gif"
    if rebuild or not out.exists():
        anim.save(out, writer="pillow", dpi=120)
    plt.close(fig)
    artefacts.append(str(out))

    # ---- bath-only (no_wp) variants ------------------------------------
    # n_bath = n_total - n_wp; same k-grid and time grid.
    grid_bath = grid_total - grid_wp

    fig, ax = plt.subplots(figsize=(8, 4))
    im = ax.imshow(grid_bath.T, origin="lower",
                   extent=[times[0], times[-1], k_vals[0], k_vals[-1]],
                   aspect="auto", cmap="viridis")
    if k0 is not None:
        ax.axhline(k0, color="red", lw=0.8, ls="--", label=f"k0 = {k0:.3f}")
        ax.legend(loc="upper right")
    ax.set_xlabel("time (a.u.)")
    ax.set_ylabel("|k| (Bohr$^{-1}$)")
    ax.set_title(f"{run_name}: bath-only momentum n(|k|, t) (no WP)")
    fig.colorbar(im, ax=ax, label="n_total - n_wp")
    fig.tight_layout()
    out = out_dir / "momentum_heatmap_no_wp.png"
    if rebuild or not out.exists():
        fig.savefig(out, dpi=130)
    plt.close(fig)
    artefacts.append(str(out))

    vmax_bath = float(np.nanmax(grid_bath))
    fig, ax = plt.subplots(figsize=(7, 4))
    line_bath, = ax.plot([], [], lw=1.6, color="C0", label="n_total - n_wp")
    ax.set_xlim(k_vals[0], k_vals[-1])
    ax.set_ylim(0, vmax_bath * 1.05 if vmax_bath > 0 else 1)
    ax.set_xlabel("|k| (Bohr$^{-1}$)")
    ax.set_ylabel("n_bath(|k|)")
    if k0 is not None:
        ax.axvline(k0, color="red", lw=0.8, ls="--", label=f"k0 = {k0:.3f}")
    ax.legend()
    title_bath = ax.set_title("")

    def _frame_bath(i):
        line_bath.set_data(k_vals, grid_bath[i])
        title_bath.set_text(f"{run_name}: bath-only t = {times[i]:.3f} a.u.")
        return line_bath, title_bath

    anim_bath = animation.FuncAnimation(fig, _frame_bath, frames=nT,
                                         interval=80, blit=False)
    out = out_dir / "momentum_distribution_no_wp.gif"
    if rebuild or not out.exists():
        anim_bath.save(out, writer="pillow", dpi=120)
    plt.close(fig)
    artefacts.append(str(out))

    return {"artefacts": artefacts}
