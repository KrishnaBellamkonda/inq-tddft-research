"""Phase: ``state_energies`` — KS state energy bar-plot GIFs.

Reads ``results/raw/observables/state_energies.csv`` (long-format CSV
written by ``inqkit::observables::StateEnergyWriter`` with columns
``step,time_au,kpoint_index,state_index,weight,occupation,
E_expect_ha,E_variance_ha2``).

Outputs:

* ``analysis/observables/ks_energies_absolute.gif``
   — bar plot of E_i(t) for all states ordered by initial energy,
   fixed y-limits across frames.
* ``analysis/observables/ks_energies_delta.gif``
   — diverging bar plot of dE_i(t) = E_i(t) - E_i(t=0); red = gain,
   blue = loss; symmetric fixed y-limits.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# TODO: The same import statement to be corrected?
from . import _common



def _read_wp_state_index(results_dir: Path) -> int | None:
    """Parse run_summary.txt for ``wp_state_index = <int>``."""
    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return None
    import re
    m = re.search(r"wp_state_index\s*=\s*(\d+)", rs.read_text())
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def _read_homo_index(results_dir: Path,
                     occ_threshold: float = 0.5) -> int | None:
    """Detect the HOMO state index from raw/observables/eigenvalues/occupations.csv.

    Returns the highest state_index whose occupation column is >=
    ``occ_threshold``. With INQ's spin-paired storage occupation max = 2,
    so 0.5 corresponds to "more than 1/4 filled" — generous enough to put
    the marker line at the smeared band edge."""
    occ_csv = results_dir / "raw" / "observables" / "eigenvalues" / "occupations.csv"
    if not occ_csv.exists():
        return None
    try:
        occ = pd.read_csv(occ_csv)
        filled = occ[occ["occupation"] >= occ_threshold]["state_index"]
        if filled.empty:
            return None
        return int(filled.max())
    except Exception:
        return None


def run(results_dir: Path, *, run_name: str, rebuild: bool, **opts) -> dict:
    csv_path = results_dir / "raw" / "observables" / "state_energies.csv"
    if not csv_path.exists():
        return {"skipped": f"missing: {csv_path}"}
    df = pd.read_csv(csv_path)
    wp_state_index = _read_wp_state_index(results_dir)
    homo_state_index = opts.get("homo_state_index",
                                _read_homo_index(results_dir))

    # Restrict to gamma point (kpoint_index = 0) for the bar plots.
    df = df[df["kpoint_index"] == 0]
    if df.empty:
        return {"skipped": "no kpoint=0 rows"}

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    artefacts: list[str] = []

    times = sorted(df["time_au"].unique())
    states = sorted(df["state_index"].unique())
    nT, nS = len(times), len(states)

    # E_i(t): rows = time, cols = state, ordered by initial energy.
    E = (df.pivot_table(index="time_au", columns="state_index",
                        values="E_expect_ha")
           .reindex(index=times, columns=states).to_numpy())
    if nT < 1 or nS < 1:
        return {"skipped": "no data"}

    # Order states by their initial energy.
    order = np.argsort(E[0])
    states_ord = np.array(states)[order]
    E_ord = E[:, order]
    dE = E_ord - E_ord[0:1, :]

    e_min = float(np.nanmin(E_ord))
    e_max = float(np.nanmax(E_ord))
    de_max = float(np.nanmax(np.abs(dE))) or 1.0

    # HOMO marker: rank in the (ordered-by-initial-energy) bar series.
    homo_rank = None
    if homo_state_index is not None:
        try:
            homo_rank = int(np.where(states_ord == homo_state_index)[0][0])
        except IndexError:
            homo_rank = None

    def _add_homo_line(ax_):
        if homo_rank is not None:
            ax_.axvline(homo_rank + 0.5, color="black", lw=1.0, ls="--",
                        alpha=0.7,
                        label=f"HOMO (state {homo_state_index})")
            ax_.legend(loc="upper left", fontsize=8)

    # ---- Absolute energies GIF ----
    fig, ax = plt.subplots(figsize=(9, 4))
    bars = ax.bar(np.arange(nS), E_ord[0], color="steelblue")
    ax.set_xlabel("state index (ordered by initial energy)")
    ax.set_ylabel("E_i(t) (Ha)")
    ax.set_ylim(e_min - 0.05 * abs(e_max - e_min),
                e_max + 0.05 * abs(e_max - e_min))
    _add_homo_line(ax)
    title = ax.set_title("")
    ax.grid(alpha=0.3, axis="y")

    def _frame_abs(i):
        for b, h in zip(bars, E_ord[i]):
            b.set_height(h)
        title.set_text(f"{run_name}: E_i  t = {times[i]:.3f} a.u.")
        return list(bars) + [title]

    anim = animation.FuncAnimation(fig, _frame_abs, frames=nT,
                                    interval=80, blit=False)
    out = out_dir / "ks_energies_absolute.gif"
    if rebuild or not out.exists():
        anim.save(out, writer="pillow", dpi=120)
    plt.close(fig)
    artefacts.append(str(out))

    # ---- Delta-E GIF (diverging) ----
    fig, ax = plt.subplots(figsize=(9, 4))
    bar_colors = ["red" if v > 0 else "steelblue" for v in dE[0]]
    bars = ax.bar(np.arange(nS), dE[0], color=bar_colors)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xlabel("state index (ordered by initial energy)")
    ax.set_ylabel(r"$\Delta E_i(t) = E_i(t) - E_i(0)$ (Ha)")
    ax.set_ylim(-1.05 * de_max, 1.05 * de_max)
    _add_homo_line(ax)
    title = ax.set_title("")
    ax.grid(alpha=0.3, axis="y")

    def _frame_delta(i):
        for b, h in zip(bars, dE[i]):
            b.set_height(h)
            b.set_color("red" if h > 0 else "steelblue")
        title.set_text(f"{run_name}: dE_i  t = {times[i]:.3f} a.u.")
        return list(bars) + [title]

    anim = animation.FuncAnimation(fig, _frame_delta, frames=nT,
                                    interval=80, blit=False)
    out = out_dir / "ks_energies_delta.gif"
    if rebuild or not out.exists():
        anim.save(out, writer="pillow", dpi=120)
    plt.close(fig)
    artefacts.append(str(out))

    # ---- bath-only ("no_wp") absolute and delta GIFs --------------------
    # WP is identified by state_index == wp_state_index (from run_summary).
    # If no WP marker found, no_wp variants are skipped.
    if wp_state_index is None:
        return {"artefacts": artefacts,
                "no_wp_skipped": "wp_state_index not found in run_summary.txt"}
    df_bath = df[df["state_index"] != wp_state_index]
    if not df_bath.empty:
        states_bath = sorted(df_bath["state_index"].unique())
        nS_bath = len(states_bath)
        E_bath = (df_bath.pivot_table(index="time_au",
                                       columns="state_index",
                                       values="E_expect_ha")
                         .reindex(index=times, columns=states_bath)
                         .to_numpy())
        order_b = np.argsort(E_bath[0])
        E_bath_ord = E_bath[:, order_b]
        dE_bath = E_bath_ord - E_bath_ord[0:1, :]
        eb_min = float(np.nanmin(E_bath_ord))
        eb_max = float(np.nanmax(E_bath_ord))
        deb_max = float(np.nanmax(np.abs(dE_bath))) or 1.0

        # Bath-only HOMO marker — rank in this filtered sub-series.
        bath_homo_rank = None
        if homo_state_index is not None:
            bath_states_ord = np.array(states_bath)[order_b]
            try:
                bath_homo_rank = int(
                    np.where(bath_states_ord == homo_state_index)[0][0])
            except IndexError:
                bath_homo_rank = None

        def _add_bath_homo_line(ax_):
            if bath_homo_rank is not None:
                ax_.axvline(bath_homo_rank + 0.5, color="black",
                            lw=1.0, ls="--", alpha=0.7,
                            label=f"HOMO (state {homo_state_index})")
                ax_.legend(loc="upper left", fontsize=8)

        fig, ax = plt.subplots(figsize=(9, 4))
        bars = ax.bar(np.arange(nS_bath), E_bath_ord[0], color="steelblue")
        ax.set_xlabel("bath state index (ordered by initial energy)")
        ax.set_ylabel("E_i(t) (Ha)  — bath only (no WP)")
        ax.set_ylim(eb_min - 0.05 * abs(eb_max - eb_min),
                    eb_max + 0.05 * abs(eb_max - eb_min))
        _add_bath_homo_line(ax)
        title = ax.set_title("")
        ax.grid(alpha=0.3, axis="y")

        def _frame_abs_bath(i):
            for b, h in zip(bars, E_bath_ord[i]):
                b.set_height(h)
            title.set_text(f"{run_name}: bath E_i  t = {times[i]:.3f} a.u.")
            return list(bars) + [title]

        anim = animation.FuncAnimation(fig, _frame_abs_bath, frames=nT,
                                        interval=80, blit=False)
        out = out_dir / "ks_energies_absolute_no_wp.gif"
        if rebuild or not out.exists():
            anim.save(out, writer="pillow", dpi=120)
        plt.close(fig)
        artefacts.append(str(out))

        fig, ax = plt.subplots(figsize=(9, 4))
        bar_colors = ["red" if v > 0 else "steelblue" for v in dE_bath[0]]
        bars = ax.bar(np.arange(nS_bath), dE_bath[0], color=bar_colors)
        ax.axhline(0, color="black", lw=0.6)
        ax.set_xlabel("bath state index (ordered by initial energy)")
        ax.set_ylabel(r"$\Delta E_i(t)$ (Ha)  — bath only (no WP)")
        ax.set_ylim(-1.05 * deb_max, 1.05 * deb_max)
        _add_bath_homo_line(ax)
        title = ax.set_title("")
        ax.grid(alpha=0.3, axis="y")

        def _frame_delta_bath(i):
            for b, h in zip(bars, dE_bath[i]):
                b.set_height(h)
                b.set_color("red" if h > 0 else "steelblue")
            title.set_text(f"{run_name}: bath dE_i  t = {times[i]:.3f} a.u.")
            return list(bars) + [title]

        anim = animation.FuncAnimation(fig, _frame_delta_bath, frames=nT,
                                        interval=80, blit=False)
        out = out_dir / "ks_energies_delta_no_wp.gif"
        if rebuild or not out.exists():
            anim.save(out, writer="pillow", dpi=120)
        plt.close(fig)
        artefacts.append(str(out))

    return {"artefacts": artefacts}
