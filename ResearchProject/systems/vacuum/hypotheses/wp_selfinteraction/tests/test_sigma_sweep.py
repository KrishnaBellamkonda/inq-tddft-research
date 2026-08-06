"""Known-case tests for sigma_sweep.py.

The anchor is the closed-form free Gaussian, built ON the sweep's own scaled
protocol. A synthetic sigma-sweep constructed from the analytic solution must
come back with zero excess at every sigma, and `E_PP(0) * sigma` must be the same
constant at every sigma — that constant is the gate the whole protocol rests on.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

N = 200
TAU_END = 1.875


def _sigma_free(t, sigma):
    return np.sqrt(sigma**2 / 2.0 + np.asarray(t, float) ** 2 / (2.0 * sigma**2))


def _write(root: Path, sigma: float, theory: str, *, width_factor=1.0,
           var_p_growth=0.0):
    """A synthetic sweep run built FROM the analytic solution on the scaled grid."""
    import sigma_sweep as S
    name = f"sweep_s{S.tag_of(sigma)}_{theory}"
    obs = root / name / "raw" / "observables"
    obs.mkdir(parents=True, exist_ok=True)
    dt = 0.00125 * sigma**2
    step = np.arange(N + 1)
    t = step * dt * (1500 / N)          # same t_end = 1.875 sigma^2, fewer rows
    sig = _sigma_free(t, sigma) * width_factor
    vp = 1.0 / (2.0 * sigma**2) + var_p_growth * t
    epp = S.epp0_predicted(sigma) * (_sigma_free(0.0, sigma) / _sigma_free(t, sigma))

    pd.DataFrame({"step": step, "time_au": t,
                  "px_mean": 0.0, "py_mean": 0.0, "pz_mean": 0.0,
                  "sigma_px2": vp, "sigma_py2": vp, "sigma_pz2": vp,
                  "e_kin_ha": 1.5 * vp, "norm_check": 1.0}
                 ).to_csv(obs / "wp_momentum_stats.csv", index=False)
    pd.DataFrame({"step": step, "time_au": t,
                  "sigma_x2": sig**2, "sigma_y2": sig**2, "sigma_z2": sig**2,
                  "sigma_x_circ": sig, "sigma_y_circ": sig, "sigma_z_circ": sig,
                  "norm_check": 1.0}
                 ).to_csv(obs / "wp_real_space_stats.csv", index=False)
    pd.DataFrame({"step": step, "time_au": t, "e_ss": 0.0, "e_pp": epp,
                  "e_ps": 0.0, "e_sb": 0.0, "e_pb": 0.0, "e_bb": 0.0}
                 ).to_csv(obs / "interactions.csv", index=False)
    (root / name / "run_summary.txt").write_text("run_completed = true\n")


@pytest.fixture()
def synth(tmp_path, monkeypatch):
    root = tmp_path / "results"
    monkeypatch.setenv("WPSI_RESULTS", str(root))
    sys.modules.pop("sigma_sweep", None)
    import sigma_sweep                                    # noqa: F401
    for s in (1.0, 2.0, 4.0):
        _write(root, s, "noninteracting")
        # planted excess GROWS with sigma, as the coupling ~ sigma argument says
        _write(root, s, "hartree", width_factor=1.0 + 0.10 * s,
               var_p_growth=1e-5)
        _write(root, s, "lda", width_factor=1.0 + 0.02 * s, var_p_growth=2e-6)
        _write(root, s, "sic_pzrun")                      # exact: back on the reference
    sys.modules.pop("sigma_sweep", None)
    import sigma_sweep as S
    yield S
    sys.modules.pop("sigma_sweep", None)


# ---------------------------------------------------------------------------
# the protocol gate — the constant everything rests on
# ---------------------------------------------------------------------------

def test_epp0_times_sigma_is_the_same_constant_at_every_sigma(synth):
    """The sweep's defining invariant.

    Because the BOX scales with sigma (L = 18 sigma), the Madelung offset is also
    ~1/sigma, so E_PP(0)*sigma is a pure constant. If someone rescales the box
    independently of sigma this test fails, which is the intended alarm.
    """
    S = synth
    vals = [S.epp0_predicted(s) * s for s in (1.0, 2.0, 4.0, 8.0)]
    assert np.allclose(vals, vals[0], rtol=1e-12)
    expected = 1.0 / math.sqrt(2.0 * math.pi) - S.XI_MADELUNG / 36.0
    assert vals[0] == pytest.approx(expected, rel=1e-12)
    assert vals[0] == pytest.approx(0.32013, abs=1e-5)


def test_protocol_gate_passes_on_synthetic_free_runs(synth):
    S = synth
    g = S.protocol_gate(S.load_all())
    assert len(g) == 3
    assert g.max_rel_sigma_err.max() < 1e-12          # built from the formula
    assert g.max_var_p_drift.max() < 1e-12            # var(p) exactly conserved
    assert np.allclose(g.epp0_vs_predicted, 1.0, rtol=1e-12)
    # the scaling invariant, measured rather than predicted
    assert g.epp0_x_sigma.std() < 1e-12


def test_tag_encoding_matches_the_run_directory_names(synth):
    S = synth
    assert S.tag_of(4.0) == "4p0"
    assert S.tag_of(1.0) == "1p0"


# ---------------------------------------------------------------------------
# the trend
# ---------------------------------------------------------------------------

def test_planted_excess_is_recovered_and_grows_with_sigma(synth):
    S = synth
    tab = S.sigma_table(S.load_all())
    for _, r in tab.iterrows():
        assert r["excess_hartree"] == pytest.approx(1.0 + 0.10 * r["sigma"], rel=1e-9)
        assert r["excess_lda"] == pytest.approx(1.0 + 0.02 * r["sigma"], rel=1e-9)
    assert tab.excess_hartree.is_monotonic_increasing


def test_xc_cancellation_is_computed_as_a_fraction_of_the_hartree_excess(synth):
    """1 - (lda-1)/(hartree-1). Planted 0.02s over 0.10s => 0.8 at every sigma."""
    S = synth
    tab = S.sigma_table(S.load_all())
    assert np.allclose(tab.xc_cancellation, 0.8, rtol=1e-9)


def test_sic_run_lands_back_on_the_reference(synth):
    S = synth
    tab = S.sigma_table(S.load_all())
    assert np.allclose(tab.excess_sic_pzrun, 1.0, atol=1e-12)
    assert tab.sic_residual.max() < 1e-12


# ---------------------------------------------------------------------------
# the fixed-physical-time view, and its honesty about coverage
# ---------------------------------------------------------------------------

def test_fixed_physical_time_drops_sigma_it_cannot_cover_and_says_so(synth):
    """t = 30 a.u. is beyond t_end = 1.875 sigma^2 for sigma < 4.

    Silently truncating here would read as "every sigma covered"; the dropped
    list is the deliverable.
    """
    S = synth
    tab, dropped = S.at_fixed_physical_time(S.load_all(), t_au=30.0)
    assert set(dropped) == {1.0, 2.0}
    assert list(tab.sigma) == [4.0]
    assert tab.tau.iloc[0] == pytest.approx(30.0 / 16.0)


def test_fixed_physical_time_covers_everything_at_a_short_enough_time(synth):
    S = synth
    tab, dropped = S.at_fixed_physical_time(S.load_all(), t_au=1.5)
    assert dropped == []
    assert len(tab) == 3


# ---------------------------------------------------------------------------
# the partial-run trap
# ---------------------------------------------------------------------------

def test_a_run_that_stops_short_of_tau_gives_NaN_not_a_nearest_match(synth,
                                                                    tmp_path):
    """Regression for a real mis-analysis (2026-08-02).

    A sweep analysed mid-flight had its `lda` runs still writing. `argmin` over
    tau happily returned their LAST row — at a smaller tau — and divided it by a
    COMPLETED reference at full tau. The excess came out ~40 % low, smooth and
    monotonic in sigma, and was misread as a cross-binary artefact. Truncating a
    comparison must therefore refuse, not approximate.
    """
    S = synth
    runs = S.load_all()
    ref = runs[(4.0, "noninteracting")]
    partial = runs[(4.0, "hartree")]
    half = partial.tau.size // 2
    partial = S.SweepRun(sigma=partial.sigma, theory=partial.theory,
                         t=partial.t[:half], tau=partial.tau[:half],
                         sigma_iso=partial.sigma_iso[:half],
                         var_p3d=partial.var_p3d[:half],
                         e_pp=partial.e_pp[:half], complete=False)

    # the honest answer at a tau the partial run never reached
    assert math.isnan(S._excess_at(partial, ref, TAU_END))
    # ... while a tau it DID reach still works
    assert not math.isnan(S._excess_at(partial, ref, float(partial.tau[-1])))

    # and the un-truncated run is unaffected
    assert not math.isnan(S._excess_at(runs[(4.0, "hartree")], ref, TAU_END))
