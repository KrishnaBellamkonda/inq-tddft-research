"""Phase: ``bath_energy`` — total KS-orbital energy excluding the WP slot.

Reads ``results/raw/observables/state_energies.csv`` (the per-orbital
``<phi_i|H|phi_i>`` time series written by
``inqkit::observables::StateEnergyWriter``) and produces:

* ``analysis/observables/bath_energy_vs_time.csv``
  — columns ``step,time_au,bath_energy_ha,bath_energy_ev,
  delta_bath_energy_ev`` (dE relative to the first snapshot).
* ``analysis/observables/bath_energy_vs_time.png``
  — line plot of bath energy vs time, with a second axis showing the
  delta from t=0 in eV.

For WP runs (``wp_state_index = N`` in ``run_summary.txt``) the WP slot
is excluded from the sum, so ``bath_energy(t)`` answers the question
"how much energy is in the bath orbitals at time t?". For classical
runs (no WP slot) the sum runs over all states — ``bath_energy(t)``
then equals the band-structure sum sum_i f_i*epsilon_i(t).

Pure Python on existing CSVs — no new C++ observable needed because
``StateEnergyWriter`` already records per-orbital ``<phi_i|H|phi_i>``.
"""

from __future__ import annotations

from pathlib import Path
import re

import matplotlib.pyplot as plt
import pandas as pd

# TODO: Is this importing convention right?
from . import _common

# TODO: Need to ensure that the 
# TODO: Build jupyter notebooks for future analysis, so that the analysis
# is together and easier to combine the answers. 

# TODO: This should also be a part of the minimum set of outcomes. 


HA_TO_EV = 27.21138625


def _read_wp_state_index(results_dir: Path) -> int | None:
    rs = results_dir / "run_summary.txt"
    if not rs.exists():
        return None
    m = re.search(r"wp_state_index\s*=\s*(\d+)", rs.read_text())
    if not m:
        return None
    try:
        return int(m.group(1))
    except ValueError:
        return None


def run(results_dir: Path, *, run_name: str, rebuild: bool, **opts) -> dict:
    csv_path = results_dir / "raw" / "observables" / "state_energies.csv"
    if not csv_path.exists():
        return {"skipped": f"missing: {csv_path}"}

    df = pd.read_csv(csv_path)
    if df.empty:
        return {"skipped": "state_energies.csv is empty"}

    wp_idx = _read_wp_state_index(results_dir)
    df_kpt0 = df[df["kpoint_index"] == 0]

    if wp_idx is not None:
        bath = df_kpt0[df_kpt0["state_index"] != wp_idx]
        wp_excluded = True
    else:
        bath = df_kpt0
        wp_excluded = False

    # Occupation-weighted band-energy sum: sum_i f_i * <phi_i|H|phi_i>
    # (the band-structure sum). For the bath alone this is the natural
    # "total energy in the bath orbitals" diagnostic.
    bath = bath.assign(
        e_weighted=bath["E_expect_ha"] * bath["occupation"])
    grouped = (bath.groupby(["step", "time_au"], as_index=False)
                   ["e_weighted"].sum()
                   .rename(columns={"e_weighted": "bath_energy_ha"}))
    grouped["bath_energy_ev"] = grouped["bath_energy_ha"] * HA_TO_EV
    e0 = grouped["bath_energy_ev"].iloc[0]
    grouped["delta_bath_energy_ev"] = grouped["bath_energy_ev"] - e0

    out_dir = _common.ensure_dir(results_dir / "analysis" / "observables")
    out_csv = out_dir / "bath_energy_vs_time.csv"
    grouped.to_csv(out_csv, index=False)

    out_png = out_dir / "bath_energy_vs_time.png"
    if _common.need_rebuild(out_png, rebuild):
        fig, ax1 = plt.subplots(figsize=(8, 4.5))
        ax1.plot(grouped["time_au"], grouped["bath_energy_ev"], "k-", lw=1.5,
                 label="bath energy (eV)")
        ax1.set_xlabel("time (a.u.)")
        ax1.set_ylabel("bath energy (eV)")
        ax2 = ax1.twinx()
        ax2.plot(grouped["time_au"], grouped["delta_bath_energy_ev"],
                 "C3-", lw=1.0, alpha=0.7,
                 label="dE_bath (eV)")
        ax2.set_ylabel("dE_bath (eV)", color="C3")
        ax2.tick_params(axis="y", labelcolor="C3")

        what = ("KS bath energy (excluding WP slot)" if wp_excluded
                else "KS band-structure sum (no WP)")
        ax1.set_title(_common.title(run_name, what))
        fig.tight_layout()
        fig.savefig(out_png, dpi=150)
        plt.close(fig)

    return {
        "wp_excluded":     wp_excluded,
        "wp_state_index":  wp_idx,
        "n_snapshots":     len(grouped),
        "bath_energy_ev":  out_csv,
        "png":             out_png,
    }
