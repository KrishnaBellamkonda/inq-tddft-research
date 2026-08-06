"""Known-case tests for selfinteraction.py.

The closed-form free-particle solution is the anchor for everything here, so it
is tested first and directly: a synthetic run built FROM the analytic formula
must come back out of the loader with zero excess spreading and zero var(p)
drift. Anything that breaks that identity breaks the whole measurement, because
the self-interaction is extracted as a difference against exactly this baseline.

The second group tests the DIFFERENCE machinery: a synthetic "interacting" run
with a planted excess width must be reported with that excess and no other.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

HERE = Path(__file__).resolve().parent
PKG = HERE.parent
if str(PKG) not in sys.path:
    sys.path.insert(0, str(PKG))

SIGMA = 4.0
DT = 0.02
N = 300
HA_TO_EV = 27.211386245988


def _sigma_free(t, sigma=SIGMA):
    return np.sqrt(sigma**2 / 2.0 + np.asarray(t, float) ** 2 / (2.0 * sigma**2))


def _write_run(root: Path, name: str, *, width_factor=1.0, var_p_growth=0.0,
               e_pp0=0.0712, e_pp_decay=0.0, hartree_matches_epp=True,
               sic=False):
    """A synthetic run built FROM the analytic solution.

    `width_factor` multiplies the free width uniformly (a planted excess);
    `var_p_growth` adds a linear drift to var(p) (a planted violation of the
    free-evolution invariant). `sic=True` additionally writes a populated
    sic.csv (an active correction), as the extended binary does.
    """
    obs = root / name / "raw" / "observables"
    obs.mkdir(parents=True, exist_ok=True)
    step = np.arange(N + 1)
    t = step * DT
    sig = _sigma_free(t) * width_factor
    vp = 1.0 / (2.0 * SIGMA**2) + var_p_growth * t
    e_pp = e_pp0 - e_pp_decay * t

    pd.DataFrame({
        "step": step, "time_au": t,
        "px_mean": 0.0, "py_mean": 0.0, "pz_mean": 0.0,
        "px2_mean": vp, "py2_mean": vp, "pz2_mean": vp,
        "sigma_px2": vp, "sigma_py2": vp, "sigma_pz2": vp,
        "e_kin_ha": 1.5 * vp, "norm_check": 1.0,
    }).to_csv(obs / "wp_momentum_stats.csv", index=False)

    pd.DataFrame({
        "step": step, "time_au": t,
        "x_mean": 0.0, "y_mean": 0.0, "z_mean": 0.0,
        "x2_mean": sig**2, "y2_mean": sig**2, "z2_mean": sig**2,
        "sigma_x2": sig**2, "sigma_y2": sig**2, "sigma_z2": sig**2,
        "norm_check": 1.0,
        "x_mean_circ": 0.0, "y_mean_circ": 0.0, "z_mean_circ": 0.0,
        "R_x": 0.9, "R_y": 0.9, "R_z": 0.9,
        "sigma_x_circ": sig, "sigma_y_circ": sig, "sigma_z_circ": sig,
    }).to_csv(obs / "wp_real_space_stats.csv", index=False)

    hartree = e_pp if hartree_matches_epp else np.zeros_like(e_pp)
    pd.DataFrame({
        "step": step, "time_au": t, "total": 1.5 * vp + hartree,
        "kinetic": 1.5 * vp, "hartree": hartree, "external": 0.0,
        "non_local": 0.0, "xc": 0.0, "exact_exchange": 0.0,
        "ion": 0.0, "ion_kinetic": 0.0,
    }).to_csv(obs / "energies.csv", index=False)

    pd.DataFrame({
        "step": step, "time_au": t, "e_ss": 0.0, "e_pp": e_pp, "e_ps": 0.0,
        "e_sb": 0.0, "e_pb": 0.0, "e_bb": 0.0,
        "e_hartree_inq": hartree, "closure_pp_minus_hartree": e_pp - hartree,
        "norm_wp": 1.0,
    }).to_csv(obs / "interactions.csv", index=False)

    if sic:
        pd.DataFrame({
            "step": step, "time_au": t,
            "u_self_ha": e_pp, "exc_self_ha": -0.4 * e_pp,
            "e_corrected_ha": 1.5 * vp,          # constant when free
            "max_overlap_pre": 0.0, "norm_removed": 0.0,
            "cum_norm_removed": 1e-14,
        }).to_csv(obs / "sic.csv", index=False)
    else:
        # the extended binary writes a header-only sic.csv when no correction
        # is active; the loader must treat it as "no SIC data"
        (obs / "sic.csv").write_text(
            "step,time_au,u_self_ha,exc_self_ha,e_corrected_ha,"
            "max_overlap_pre,norm_removed,cum_norm_removed\n")

    (root / name / "run_summary.txt").write_text(
        f"run_completed = true\nwp = gaussian sigma {SIGMA} k0 0 mass 1\n")


@pytest.fixture()
def synth(tmp_path, monkeypatch):
    root = tmp_path / "results"
    # the reference: exactly free, E_PP present but not felt (hartree = 0)
    _write_run(root, "noninteracting", hartree_matches_epp=False)
    # hartree: 8 % excess width, var(p) drifts, E_PP decays as the packet spreads
    _write_run(root, "hartree", width_factor=1.08, var_p_growth=2.0e-5,
               e_pp_decay=4.0e-4)
    # lda: 12 % excess -- more than hartree, so the xc part is 4 %
    _write_run(root, "lda", width_factor=1.12, var_p_growth=3.0e-5,
               e_pp_decay=5.0e-4)
    monkeypatch.setenv("WPSI_RESULTS", str(root))
    sys.modules.pop("selfinteraction", None)
    import selfinteraction
    yield selfinteraction
    sys.modules.pop("selfinteraction", None)


# ---------------------------------------------------------------------------
# The closed-form anchor
# ---------------------------------------------------------------------------

def test_analytic_reference_formulas(synth):
    S = synth
    assert S.sigma_dens_free(0.0, 4.0) == pytest.approx(4.0 / np.sqrt(2.0))
    assert S.var_p_free(4.0) == pytest.approx(0.03125)
    # 3/(4 sigma^2) in eV -- the t=0 value of T2-T1 in the channeling analysis
    assert S.localisation_ev(4.0) == pytest.approx(1.2755, rel=1e-3)
    # the spreading law must be the DENSITY width, not the wavefunction width:
    # sigma_dens(0) = sigma/sqrt2, NOT sigma. Getting this wrong inflates every
    # excess-spreading ratio by sqrt2.
    assert S.sigma_dens_free(0.0, 4.0) < 4.0


def test_reference_run_has_zero_excess_and_passes_the_numerics_gate(synth):
    S = synth
    ref = S.load("noninteracting")
    assert np.allclose(ref.excess_spreading(), 1.0, atol=1e-12)
    g = S.numerics_gate(ref)
    assert g["passed"], g
    assert g["max_rel_sigma_error"] < 1e-12
    assert g["max_var_p_drift"] < 1e-12


def test_var_p_is_the_sharpest_gate(synth):
    """A planted var(p) drift must be caught even though the width also moved."""
    S = synth
    g = S.numerics_gate(S.load("lda"))
    assert not g["passed"]
    assert g["max_var_p_drift"] > 1e-3


# ---------------------------------------------------------------------------
# The difference machinery
# ---------------------------------------------------------------------------

def test_excess_width_is_measured_against_the_reference_run(synth):
    """Planted 8 % and 12 % excesses must come back as 8 % and 12 %."""
    S = synth
    runs = S.load_all()
    ref = runs["noninteracting"]
    assert np.allclose(S.effect(runs["hartree"], ref).sigma_ratio, 1.08, atol=1e-12)
    assert np.allclose(S.effect(runs["lda"], ref).sigma_ratio, 1.12, atol=1e-12)


def test_summary_table_splits_hartree_from_xc(synth):
    """lda - hartree must isolate the xc part: 12 % - 8 % = 4 %."""
    S = synth
    tab = S.summary_table(S.load_all())
    row = tab[tab.theory.str.startswith("xc part")].iloc[0]
    assert row["excess_width_pct"] == pytest.approx(4.0, abs=1e-9)
    h = tab[tab.theory == "hartree"].iloc[0]
    assert h["excess_width_pct"] == pytest.approx(8.0, abs=1e-9)


def test_e_pp_release_is_reported_with_the_right_sign(synth):
    """E_PP falls as the packet spreads, so 'released' must be POSITIVE."""
    S = synth
    tab = S.summary_table(S.load_all())
    h = tab[tab.theory == "hartree"].iloc[0]
    assert h["E_PP_start_eV"] > h["E_PP_end_eV"]
    assert h["E_PP_released_eV"] > 0.0


def test_closure_gate_is_applied_only_to_interacting_theories(synth):
    """INQ reports hartree = 0 for a non-interacting run BY CONSTRUCTION.

    Gating the closure on that run would fail for a reason that is not an error.
    The offline E_PP stays non-zero there and is a pure size diagnostic.
    """
    S = synth
    runs = S.load_all()
    c_ref = S.closure(runs["noninteracting"])
    assert c_ref["gated"] is False and c_ref["passed"]
    assert np.abs(runs["noninteracting"].e_pp).max() > 0    # still measured
    c_h = S.closure(runs["hartree"])
    assert c_h["gated"] is True and c_h["passed"]


def test_closure_gate_fails_on_a_broken_interacting_run(synth, tmp_path):
    S = synth
    root = tmp_path / "results"
    _write_run(root, "hartree_broken", width_factor=1.08, hartree_matches_epp=False)
    run = S.load("hartree", name="hartree_broken")
    assert not S.closure(run)["passed"]


def test_channeling_comparison_is_a_ratio_not_a_decomposition(synth):
    S = synth
    out = S.channeling_comparison(S.load_all())
    exc = S.CHANNELING_EXCESS - 1.0
    # planted 12 % vacuum excess against the channeling excess
    assert out["lda_fraction_of_channeling"] == pytest.approx(0.12 / exc, rel=1e-6)
    assert 0.0 < out["lda_fraction_of_channeling"] < 1.0


def test_channeling_constants_share_one_width_definition(synth):
    """The vacuum ratio is 3-D isotropic, so the channeling constants must be too.

    Regression guard for a real defect: `CHANNELING_EXCESS` was 1.467, the
    TRANSVERSE <r_perp>/free ratio, while everything it was divided into was a
    3-D geometric-mean ratio. Mixing the two understated the self-interaction
    fraction by ~4 points and produced no error of any kind.
    """
    S = synth
    # the 3-D and transverse numbers are genuinely different -- if someone
    # "tidies" them into one value, that is the bug coming back
    assert S.CHANNELING_EXCESS != S.CHANNELING_EXCESS_TRANSVERSE
    assert S.CHANNELING_EXCESS_SIC < S.CHANNELING_EXCESS   # SIC removes spreading

    # the directly-measured SIC fraction must be a fraction, and must agree with
    # the vacuum PREDICTION -- their agreement is the study's headline result
    out = S.channeling_comparison(S.load_all())
    assert 0.0 < out["sic_measured_fraction"] < 1.0
    assert out["sic_measured_fraction"] == pytest.approx(
        (S.CHANNELING_EXCESS - S.CHANNELING_EXCESS_SIC)
        / (S.CHANNELING_EXCESS - 1.0), rel=1e-9)


def test_wrap_indicator_is_zero_for_a_contained_packet(synth):
    S = synth
    ref = S.load("noninteracting")
    assert float(ref.wrap_indicator.max()) == pytest.approx(0.0, abs=1e-12)


# ---------------------------------------------------------------------------
# The SIC runs (intervention arm)
# ---------------------------------------------------------------------------

def test_load_all_skips_absent_theories_but_requires_the_reference(synth):
    """The fixture writes only the three original runs; the two sic_* dirs are
    absent (as they are for the `_smoke` suffix). load_all must skip them
    silently — and fail loudly if the REFERENCE itself is missing."""
    S = synth
    runs = S.load_all()
    assert set(runs) == {"noninteracting", "hartree", "lda"}
    with pytest.raises(FileNotFoundError):
        S.load_all("_no_such_suffix")


def test_headeronly_sic_csv_means_no_sic_data(synth):
    """Non-SIC runs of the extended binary write a header-only sic.csv."""
    S = synth
    run = S.load("lda")
    assert run.u_self is None and run.e_corrected is None


def test_sic_run_loads_diagnostics_and_a_perfect_correction_passes_the_gate(synth):
    """A corrected run planted ON the free solution must (a) surface its SIC
    diagnostics, (b) pass the reference's own numerics gate — the Tier V
    criterion — and (c) report ~0 excess in the summary table."""
    S = synth
    _write_run(S.RESULTS, "sic_pzrun", width_factor=1.0, e_pp_decay=4.0e-4,
               sic=True)
    run = S.load("sic_pzrun")
    assert run.u_self is not None and run.u_self[0] == pytest.approx(0.0712)
    assert run.e_corrected is not None
    # the conserved quantity: planted constant, so zero drift
    assert np.ptp(run.e_corrected) == pytest.approx(0.0, abs=1e-15)

    g = S.numerics_gate(run)
    assert g["passed"], g

    tab = S.summary_table(S.load_all())
    row = tab[tab.theory == "sic_pzrun"].iloc[0]
    assert row["excess_width_pct"] == pytest.approx(0.0, abs=1e-9)
    # the hartree/lda difference row must survive the extra rows
    assert tab[tab.theory.str.startswith("xc part")].iloc[0][
        "excess_width_pct"] == pytest.approx(4.0, abs=1e-9)
