"""Band-sum energy ledger (analysis layer; the pure compute half of the
``energy_balance`` phase, ADR 0003 split).

Complementary to ``analysis.energy_components`` (the functional KS decomposition,
IV-M07): this is the **band-sum** ledger that attributes the total-energy drift
to the WP slot vs the occupied bath, leaving an "unaccounted" residual that
flags excitation into initially-empty states.

    dE_wp(t)        = ΔE of the WP slot (occupation 1 by construction)
    dE_bath(t)      = Σ_{i≠wp} f_i · Δε_i(t)          (occ-weighted bath)
    dE_total(t)     = E_total_obs(t) − E_total_obs(0) (interp from observables)
    unaccounted(t)  = dE_total − (dE_wp + dE_bath)

Pure numpy + pandas (no matplotlib). The plotting + IFW highlight stay in
``inqview.pipeline.energy_balance``.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

HA_TO_EV = 27.211386245988

_E_COLS = ("E_expect_ha", "eigenvalue_ha")  # accept either state-energy column


def compute_ledger(se: pd.DataFrame, occ: pd.DataFrame, obs: pd.DataFrame,
                   wp_idx: int, *, e_col: str | None = None) -> pd.DataFrame:
    """Per-step band-sum energy ledger as a DataFrame.

    ``se`` = state_energies (long: step,time_au,state_index,E_expect_ha);
    ``occ`` = occupations (state_index,occupation); ``obs`` = observables
    (time_au,energy_total); ``wp_idx`` = the WP slot's state index. Columns:
    time_au, dE_{wp,bath,total,unaccounted}_{ha,ev}.
    """
    if e_col is None:
        e_col = next((c for c in _E_COLS if c in se.columns), None)
    if e_col is None:
        raise ValueError(f"state_energies has no {_E_COLS}; cols: {set(se.columns)}")

    occ_map = dict(zip(occ["state_index"].astype(int), occ["occupation"].astype(float)))
    se = se.sort_values(["state_index", "time_au"])
    first_per_state = se.groupby("state_index")[e_col].first()

    rows = []
    for t in sorted(se["time_au"].unique()):
        snap = se[se["time_au"] == t][["state_index", e_col]]
        d_state = snap[e_col].to_numpy() - first_per_state.loc[snap["state_index"]].to_numpy()
        wp_mask = snap["state_index"].to_numpy() == wp_idx
        dE_wp = float((d_state[wp_mask] * 1.0).sum())          # WP occ = 1 by construction
        bath_states = snap["state_index"].to_numpy()[~wp_mask]
        f_bath = np.array([occ_map.get(int(s), 0.0) for s in bath_states])
        dE_bath = float((d_state[~wp_mask] * f_bath).sum())
        rows.append((t, dE_wp, dE_bath))

    df = pd.DataFrame(rows, columns=["time_au", "dE_wp_ha", "dE_bath_ha"])
    if "energy_total" in obs.columns:
        # np.interp requires monotonically increasing xp — sort obs by time so a
        # restart / concatenated observables.csv can't yield silently-wrong dE_total.
        obs_s = obs.sort_values("time_au")
        ot = obs_s["time_au"].to_numpy()
        oe = obs_s["energy_total"].to_numpy()
        df["dE_total_ha"] = np.interp(df["time_au"], ot, oe) - float(oe[0])
    else:
        df["dE_total_ha"] = np.nan
    df["unaccounted_ha"] = df["dE_total_ha"] - (df["dE_wp_ha"] + df["dE_bath_ha"])
    for c in ("dE_wp", "dE_bath", "dE_total", "unaccounted"):
        df[c + "_ev"] = df[c + "_ha"] * HA_TO_EV
    return df
