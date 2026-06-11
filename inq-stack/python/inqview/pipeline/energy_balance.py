"""Phase: ``energy_balance`` — projectile/bath/unaccounted energy ledger.

Implements observables_reference §13.3: at every recorded step, compute

    ΔE_WP(t)          = E_WP(t) − E_WP(0)            single-state (occ × ε_WP)
    ΔE_bath(t)        = Σ_{i ≠ WP} f_i [ε_i(t) − ε_i(0)]   occ-weighted bath sum
    ΔE_total_obs(t)   = E_total(t) − E_total(0)      drift sanity from observables.csv
    Unaccounted(t)    = ΔE_total_obs − (ΔE_WP + ΔE_bath)   suggests excitation
                                                            into initially-empty
                                                            states (numerical
                                                            sink + WP→empty
                                                            slots transitions)

All four traces are written to ``analysis/observables/energy_balance.csv``
and rendered on a single time-axis plot ``energy_balance.png``. The plot
uses ScalarFormatter(useOffset=False) per §13.1.1 and applies the
campaign IFW highlight.

Inputs:
  - ``raw/observables/state_energies.csv`` — per-step ε_i(t) (long format
    ``step, time_au, kpoint_index, state_index, eigenvalue_ha``).
  - ``raw/observables/eigenvalues/occupations.csv`` — initial GS
    occupations f_i (cols ``state_index, occupation``).
  - ``raw/observables/observables.csv`` — ``energy_total(t)`` for drift.
  - ``run_summary.txt`` — ``wp_state_index`` to identify the WP slot.

The phase silently skips when any input is absent (e.g. free-space runs
without per-state energies).
"""

# TODO: The output of this, tells us how the jellium bath and wp orbital
# exchange energy is to be a part of the minimum set of observables to be calculated. 

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np
import pandas as pd

from . import _common


HA_TO_EV = 27.21138625


def _read_wp_index(results_dir: Path) -> int | None:
    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return None
    m = re.search(r"^\s*wp_state_index\s*=\s*(\d+)",
                  rs.read_text(), flags=re.MULTILINE)
    return int(m.group(1)) if m else None


def run(results_dir: Path, *, run_name: str, rebuild: bool, **_) -> dict:
    raw = results_dir / "raw" / "observables"
    se_csv  = raw / "state_energies.csv"
    occ_csv = raw / "eigenvalues" / "occupations.csv"
    obs_csv = raw / "observables.csv"

    for path in (se_csv, occ_csv, obs_csv):
        if not path.exists():
            return {"skipped": f"missing input: {path}"}

    wp_idx = _read_wp_index(results_dir)
    if wp_idx is None:
        return {"skipped": "no wp_state_index in run_summary.txt"}

    se  = pd.read_csv(se_csv, comment="#")
    occ = pd.read_csv(occ_csv, comment="#")
    obs = pd.read_csv(obs_csv, comment="#")

    # State energies use E_expect_ha (= <ψ_i(t) | H(t) | ψ_i(t)>); the
    # legacy column name is eigenvalue_ha. Accept either.
    e_col = next((c for c in ("E_expect_ha", "eigenvalue_ha") if c in se.columns), None)
    if e_col is None:
        return {"skipped": f"state_energies.csv has no E_expect_ha or eigenvalue_ha (cols: {set(se.columns)})"}
    if {"step", "time_au", "state_index"}.difference(se.columns):
        return {"skipped": f"state_energies.csv missing step/time_au/state_index: {set(se.columns)}"}
    if {"state_index", "occupation"}.difference(occ.columns):
        return {"skipped": f"occupations.csv missing required columns: {set(occ.columns)}"}

    # Pure ledger compute moved to the analysis layer (ADR 0003 split).
    from ..analysis.energy_balance import compute_ledger
    df = compute_ledger(se, occ, obs, wp_idx, e_col=e_col)

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    csv_out = out_dir / "energy_balance.csv"
    if _common.need_rebuild(csv_out, rebuild):
        df.to_csv(csv_out, index=False)

    # IFW highlight (single-run, per the campaign rule)
    ifw = _common.post_ifw_window_from_summary(results_dir)

    out_png = out_dir / "energy_balance.png"
    if _common.need_rebuild(out_png, rebuild):
        fig, ax = plt.subplots(figsize=(9, 5))
        if ifw is not None:
            _common.ifw_highlight(ax, ifw[0])
        ax.plot(df["time_au"], df["dE_wp_ev"],       "C3-", lw=1.6,
                label=r"$\Delta E_{\rm WP}(t)$")
        ax.plot(df["time_au"], df["dE_bath_ev"],     "C0-", lw=1.6,
                label=r"$\Delta E_{\rm bath}(t) = \sum_{i\neq{\rm WP}} f_i\,\Delta\varepsilon_i$")
        ax.plot(df["time_au"], df["dE_total_ev"],    "k-",  lw=1.2,
                label=r"$\Delta E_{\rm total}^{\rm obs}(t)$ (drift sanity)")
        ax.plot(df["time_au"], df["unaccounted_ev"], "C5--", lw=1.3,
                label=r"Unaccounted $= \Delta E_{\rm total} - (\Delta E_{\rm WP} + \Delta E_{\rm bath})$")
        ax.axhline(0.0, color="0.6", lw=0.7)
        ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False,
                                                     useMathText=True))
        ax.set_xlabel("time (a.u.)")
        ax.set_ylabel("Δ energy (eV)")
        ax.set_title(f"{run_name}: projectile / bath / unaccounted energy ledger")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best", fontsize=9)
        fig.tight_layout()
        fig.savefig(out_png, dpi=150)
        plt.close(fig)

    return {"csv": str(csv_out), "png": str(out_png),
            "wp_state_index": wp_idx,
            "dE_wp_final_ev": float(df["dE_wp_ev"].iloc[-1]),
            "dE_bath_final_ev": float(df["dE_bath_ev"].iloc[-1]),
            "unaccounted_final_ev": float(df["unaccounted_ev"].iloc[-1])}
