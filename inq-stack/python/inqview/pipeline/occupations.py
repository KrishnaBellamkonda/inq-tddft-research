"""Phase: ``occupations`` — RT bar GIFs of f_i(t) and Δf_i(t).

Reads ``raw/observables/occupations_vs_time.csv`` (the per-step dump
written by the new ``inqkit::observables::OccupationsWriter``) and
produces two animated bar plots, paralleling the ``state_energies``
phase but for occupations:

* ``analysis/observables/occupations_absolute.gif`` — f_i(t) ordered by
  state index, with a vertical dashed HOMO line.
* ``analysis/observables/occupations_delta.gif`` — Δf_i(t) =
  f_i(t) − f_i(0); diverging colours (red = gain, blue = loss); HOMO
  line.

Note on physics: in INQ's TDDFT propagation the KS occupations f_i are
held FROZEN (the basis evolves; the coefficients on the basis don't).
These plots therefore primarily serve as an *audit* — the bars should
not visibly change. A non-flat trace would be a numerics red flag.
The corresponding "physically dynamic" quantity is the GS-projected
occupation n_i^GS(t) = Σ_j f_j |⟨ψ_i^GS|ψ_j(t)⟩|², which lives in a
separate (planned) postprocess phase that requires the full overlap
matrix.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import _common
from .state_energies import _read_homo_index


def run(results_dir: Path, *, run_name: str, rebuild: bool, **opts) -> dict:
    # Accept either the legacy 'occupations_vs_time.csv' or the new
    # 'occupations.csv' produced by inqkit::observables::OccupationsWriter.
    csv_path = results_dir / "raw" / "observables" / "occupations_vs_time.csv"
    if not csv_path.exists():
        alt_path = results_dir / "raw" / "observables" / "occupations.csv"
        if alt_path.exists():
            csv_path = alt_path
        else:
            return {"skipped": f"missing: {csv_path} (and {alt_path})"}
    df = pd.read_csv(csv_path)
    if df.empty:
        return {"skipped": "occupations_vs_time.csv is empty"}

    # Restrict to gamma point (kpoint_index=0) for the bar plots.
    if "kpoint_index" in df.columns:
        df = df[df["kpoint_index"] == 0]
    if df.empty:
        return {"skipped": "no kpoint=0 rows"}

    homo_state_index = opts.get("homo_state_index",
                                _read_homo_index(results_dir))

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    artefacts: list[str] = []

    times  = sorted(df["time_au"].unique())
    states = sorted(df["state_index"].unique())
    nT, nS = len(times), len(states)

    F = (df.pivot_table(index="time_au", columns="state_index",
                        values="occupation")
           .reindex(index=times, columns=states).to_numpy())
    if nT < 1 or nS < 1:
        return {"skipped": "no data"}

    dF = F - F[0:1, :]

    f_min = float(np.nanmin(F))
    f_max = float(np.nanmax(F))
    df_max = float(np.nanmax(np.abs(dF))) or 1e-12  # avoid zero ylim

    # HOMO marker: rank in the (ordered-by-state-index) bar series.
    states_arr = np.array(states)
    homo_rank = None
    if homo_state_index is not None:
        try:
            homo_rank = int(np.where(states_arr == homo_state_index)[0][0])
        except IndexError:
            homo_rank = None

    def _add_homo_line(ax_):
        if homo_rank is not None:
            ax_.axvline(homo_rank + 0.5, color="black", lw=1.0, ls="--",
                        alpha=0.7,
                        label=f"HOMO (state {homo_state_index})")
            ax_.legend(loc="upper right", fontsize=8)

    # ---- Absolute occupations GIF ----
    fig, ax = plt.subplots(figsize=(9, 4))
    bars = ax.bar(np.arange(nS), F[0], color="steelblue")
    ax.set_xlabel("state index")
    ax.set_ylabel("occupation $f_i(t)$")
    pad = 0.05 * max(1.0, abs(f_max - f_min))
    ax.set_ylim(min(0.0, f_min) - pad, f_max + pad)
    _add_homo_line(ax)
    title = ax.set_title("")
    ax.grid(alpha=0.3, axis="y")

    def _frame_abs(i):
        for b, h in zip(bars, F[i]):
            b.set_height(h)
        title.set_text(f"{run_name}: f_i  t = {times[i]:.3f} a.u.")
        return list(bars) + [title]

    anim = animation.FuncAnimation(fig, _frame_abs, frames=nT,
                                    interval=80, blit=False)
    out = out_dir / "occupations_absolute.gif"
    if rebuild or not out.exists():
        anim.save(out, writer="pillow", dpi=120)
    plt.close(fig)
    artefacts.append(str(out))

    # ---- Delta occupations GIF (diverging) ----
    fig, ax = plt.subplots(figsize=(9, 4))
    bar_colors = ["red" if v > 0 else "steelblue" for v in dF[0]]
    bars = ax.bar(np.arange(nS), dF[0], color=bar_colors)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xlabel("state index")
    ax.set_ylabel(r"$\Delta f_i(t) = f_i(t) - f_i(0)$")
    ax.set_ylim(-1.05 * df_max, 1.05 * df_max)
    _add_homo_line(ax)
    title = ax.set_title("")
    ax.grid(alpha=0.3, axis="y")

    # No-offset y-axis per observables_reference §13.1 styling rule 1.
    from matplotlib.ticker import ScalarFormatter
    fmt = ScalarFormatter(useOffset=False, useMathText=True)
    fmt.set_powerlimits((-3, 3))
    ax.yaxis.set_major_formatter(fmt)

    def _frame_delta(i):
        for b, h in zip(bars, dF[i]):
            b.set_height(h)
            b.set_color("red" if h > 0 else "steelblue")
        title.set_text(f"{run_name}: Δf_i  t = {times[i]:.3f} a.u.")
        return list(bars) + [title]

    anim = animation.FuncAnimation(fig, _frame_delta, frames=nT,
                                    interval=80, blit=False)
    out = out_dir / "occupations_delta.gif"
    if rebuild or not out.exists():
        anim.save(out, writer="pillow", dpi=120)
    plt.close(fig)
    artefacts.append(str(out))

    return {"artefacts": artefacts}
