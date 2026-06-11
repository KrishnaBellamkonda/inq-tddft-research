"""Analytic test for the band-sum energy ledger (B2 split).

Two states over two steps with hand-chosen energies so the ledger is known up
front: WP slot (idx 1) Δε=3, bath state (idx 0, occ 2) Δε=1, total drift 4.
  dE_wp = 3·1            = 3
  dE_bath = 2·1          = 2
  dE_total              = 4   (interp from observables)
  unaccounted = 4−(3+2) = −1
Pure numpy/pandas (deps-clean).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from inqview.analysis.energy_balance import HA_TO_EV, compute_ledger

pytestmark = pytest.mark.analysis

_SE = pd.DataFrame({
    "step":        [0, 0, 2, 2],
    "time_au":     [0.0, 0.0, 1.0, 1.0],
    "state_index": [0, 1, 0, 1],
    "E_expect_ha": [10.0, 5.0, 11.0, 8.0],     # Δε_0=1, Δε_wp=3
})
_OCC = pd.DataFrame({"state_index": [0, 1], "occupation": [2.0, 0.0]})  # WP GS occ=0
_OBS = pd.DataFrame({"time_au": [0.0, 1.0], "energy_total": [100.0, 104.0]})


def test_ledger_known_values():
    df = compute_ledger(_SE, _OCC, _OBS, wp_idx=1)
    # t=0 row is all-zero (Δ from itself)
    assert np.allclose(df["dE_wp_ha"].to_numpy(),       [0.0, 3.0])
    assert np.allclose(df["dE_bath_ha"].to_numpy(),     [0.0, 2.0])   # 2.0·1.0
    assert np.allclose(df["dE_total_ha"].to_numpy(),    [0.0, 4.0])
    assert np.allclose(df["unaccounted_ha"].to_numpy(), [0.0, -1.0])


def test_ledger_ev_columns_scale():
    df = compute_ledger(_SE, _OCC, _OBS, wp_idx=1)
    assert np.allclose(df["dE_wp_ev"].to_numpy(), df["dE_wp_ha"].to_numpy() * HA_TO_EV)
    assert np.allclose(df["unaccounted_ev"].to_numpy(),
                       df["unaccounted_ha"].to_numpy() * HA_TO_EV)


def test_ledger_missing_energy_column_raises():
    bad = _SE.rename(columns={"E_expect_ha": "something_else"})
    with pytest.raises(ValueError):
        compute_ledger(bad, _OCC, _OBS, wp_idx=1)
