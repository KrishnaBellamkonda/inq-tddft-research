"""Tests for the functional energy-flow kernel (``inqview.analysis.energy_components``).

Pure (numpy/pandas). The decisive invariant is exact by construction:
E_kin + E_H + E_xc + E_ext == E_total at every step, and Σ ΔE_component == ΔE_total.
Built on a hand-made observables table with known components (IV-M07).
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from inqview.analysis import energy_components as ec

pytestmark = pytest.mark.analysis


def _frame():
    # known components; E_total = kin + H + xc + ext(known external)
    t = np.array([0.0, 1.0, 2.0])
    E_kin = np.array([1.00, 1.20, 1.50])
    E_H = np.array([0.50, 0.45, 0.40])
    E_xc = np.array([-0.30, -0.32, -0.28])
    E_ext = np.array([-2.00, -2.05, -2.10])           # the "true" external term
    E_total = E_kin + E_H + E_xc + E_ext
    return pd.DataFrame({"time_au": t, "energy_total": E_total,
                         "energy_kinetic": E_kin, "energy_hartree": E_H,
                         "energy_xc": E_xc}), E_ext


def test_external_is_recovered_as_residual():
    df, E_ext_true = _frame()
    r = ec.compute(df)
    np.testing.assert_allclose(r.E_ext, E_ext_true, atol=1e-12)


def test_components_sum_to_total_exactly():
    df, _ = _frame()
    r = ec.compute(df)
    np.testing.assert_allclose(r.component_sum(), r.E_total, atol=1e-12)


def test_delta_components_sum_to_delta_total():
    df, _ = _frame()
    r = ec.compute(df)
    dsum = r.dE_kin + r.dE_hartree + r.dE_xc + r.dE_ext
    np.testing.assert_allclose(dsum, r.dE_total, atol=1e-12)


def test_conserved_total_gives_zero_drift():
    """If E_total is constant, ΔE_total ≈ 0 (energy-conservation sanity)."""
    df, _ = _frame()
    df["energy_total"] = df["energy_total"].iloc[0]      # force constant total
    r = ec.compute(df)
    assert abs(r.dE_total[-1]) < 1e-12


def test_redistribution_in_ev():
    df, _ = _frame()
    r = ec.compute(df)
    red = r.redistribution_ev()
    # ΔE_kin = (1.50-1.00) Ha * 27.2114 ≈ 13.6 eV
    assert red["kinetic"] == pytest.approx(0.5 * ec.HA_TO_EV, rel=1e-6)
    # the component ΔE's still sum to the total ΔE in eV
    assert (red["kinetic"] + red["hartree"] + red["xc"] + red["external"]
            == pytest.approx(red["total"], abs=1e-9))


def test_missing_column_raises():
    df, _ = _frame()
    with pytest.raises(ValueError):
        ec.compute(df.drop(columns=["energy_xc"]))


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
