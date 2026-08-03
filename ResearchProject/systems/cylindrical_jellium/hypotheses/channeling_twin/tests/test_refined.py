"""Known-case tests for refined.py — the refined-analysis data layer.

Two kinds of test here, and the distinction matters:

  SYNTHETIC   a run whose answer is known in closed form is written to a temp
              tree and loaded through the real loaders. These pin the ALGEBRA.
  PRODUCTION  the real channeling run is loaded and the identities are asserted
              on it. These pin that the algebra still holds on data that has
              been through the GPU, and they are skipped (not failed) if the
              results tree is absent, so the suite stays portable.

The identities under test are the ones the notebook's conclusions rest on:

  1. T2 - T1 == var(p)/2m EXACTLY. If this drifts, the drift/spread split that
     the whole study turns on is not a split of anything.
  2. T2 == INQ's own e_kin_ha. Independent confirmation that our reconstruction
     from moments matches what the engine reported.
  3. The classical energy budget closes: d(E_bath) + d(KE_proj) ~ 0.
  4. s_pintegral recovers a known constant-velocity track exactly (the
     cumulative trapezoid is exact for a linear integrand).
  5. Interaction deltas are identically zero at t = 0 by construction.
  6. fit_in_window recovers a planted slope to machine precision.
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


# ---------------------------------------------------------------------------
# Synthetic-run fixtures
# ---------------------------------------------------------------------------

N_STEPS = 200
DT = 0.02
V0 = 1.917011
Z0 = -28.0
HA_TO_EV = 27.211386
S_PLANT = 0.500           # eV/Bohr, the stopping power planted in the fake run

# CONSTANT DECELERATION is what makes the planted S exact, and the reason is
# worth stating because the obvious construction is wrong. If instead you plant
# T1 linear in t and then integrate the resulting p_z, the path is NOT V0*t --
# the projectile is slowing -- so dT1/ds is not the constant you planted and the
# fit misses by ~4e-4 relative (this test caught exactly that on first write).
#
# With dp/dt = a constant instead:
#     dT1/dt = p (dp/dt) = p a        ds/dt = p        =>  dT1/ds = a  EXACTLY
# so T1 is exactly linear in s, the OLS fit is exact, and the trapezoid path is
# exact too because p(t) is linear (the trapezoid rule is exact for linear
# integrands). Nothing is approximated anywhere in this fixture.
ACCEL = -S_PLANT / HA_TO_EV       # a.u., dp/dt


def _write_wp_run(root: Path, name: str = "wp") -> Path:
    """A wavepacket run with analytically known moments and an exact planted S.

    var(p) is given a deliberate linear growth so that T2 - T1 is NOT constant
    and the identity test has something to bite on.
    """
    obs = root / name / "raw" / "observables"
    obs.mkdir(parents=True, exist_ok=True)
    step = np.arange(N_STEPS + 1)
    t = step * DT

    pz = V0 + ACCEL * t

    var_pz = 0.03125 + 0.0004 * t          # deliberate growth
    var_px = np.full_like(t, 0.03125)
    var_py = np.full_like(t, 0.03125)
    var3 = var_px + var_py + var_pz
    e_kin = 0.5 * (pz**2 + var3)           # == T2 by construction

    pd.DataFrame({
        "step": step, "time_au": t,
        "px_mean": 0.0, "py_mean": 0.0, "pz_mean": pz,
        "px2_mean": var_px, "py2_mean": var_py, "pz2_mean": pz**2 + var_pz,
        "sigma_px2": var_px, "sigma_py2": var_py, "sigma_pz2": var_pz,
        "e_kin_ha": e_kin, "norm_check": 1.0,
    }).to_csv(obs / "wp_momentum_stats.csv", index=False)

    z = Z0 + np.concatenate([[0.0], np.cumsum(0.5 * (pz[1:] + pz[:-1]) * np.diff(t))])
    pd.DataFrame({
        "step": step, "time_au": t,
        "x_mean": 0.0, "y_mean": 0.0, "z_mean": z,
        "x2_mean": 8.0, "y2_mean": 8.0, "z2_mean": 8.0,
        "sigma_x2": 8.0, "sigma_y2": 8.0, "sigma_z2": 8.0,
        "norm_check": 1.0,
        "x_mean_circ": 0.0, "y_mean_circ": 0.0, "z_mean_circ": z,
        "R_x": 0.9, "R_y": 0.9, "R_z": 0.9,
        "sigma_x_circ": 2.8284, "sigma_y_circ": 2.8284, "sigma_z_circ": 2.8284,
    }).to_csv(obs / "wp_real_space_stats.csv", index=False)

    pd.DataFrame({
        "step": step, "time_au": t,
        "f_bore": 0.99, "f_wall": 0.01, "f_outside": 0.0,
        "r_mean": 3.5, "r2_mean": 16.0, "sigma_r": 1.85, "norm_total": 1.0,
    }).to_csv(obs / "wp_radial_occupancy.csv", index=False)

    _write_observables(obs, step, t)
    _write_interactions(obs, step, t)
    _write_momentum_distribution(obs, step, t)
    (root / name / "run_summary.txt").write_text(
        f"run_completed = true\nn_steps = {N_STEPS}\n")
    return root / name


N_BINS = 16


def _write_momentum_distribution(obs: Path, step, t):
    """LONG-FORMAT: N_BINS rows per step. The shape is the point of the test.

    Written every 10th step so the cadence differs from the scalars, matching
    production, and so nearest_slices has something to snap to.
    """
    k = np.linspace(0.05, 4.0, N_BINS)
    rows = []
    for s, tt in zip(step, t):
        if s % 10:
            continue
        n = np.exp(-((k - (V0 + ACCEL * tt)) ** 2) / (2 * 0.18**2))
        n /= n.sum()
        rows.append(pd.DataFrame({"step": s, "time_au": tt, "k_bohr_inv": k,
                                  "n_total": n * 160.0, "n_wp": n}))
    pd.concat(rows, ignore_index=True).to_csv(
        obs / "momentum_distribution.csv", index=False)


def _write_cl_run(root: Path, name: str = "classical") -> Path:
    """A classical run whose energy budget closes EXACTLY by construction.

    E_total(bath) is defined as the negative of the projectile KE change, so the
    closure residual is identically zero up to float rounding. That makes the
    closure test a test of the LOADER (does it pair the right columns, subtract
    the right baseline, convert units once) rather than of the physics — which is
    what a unit test can honestly check. The physics closure is checked on the
    production run in the second test class.
    """
    obs = root / name / "raw" / "observables"
    obs.mkdir(parents=True, exist_ok=True)
    step = np.arange(N_STEPS + 1)
    t = step * DT

    # Same constant-deceleration construction as the WP fixture, so z, vz and ke
    # are mutually exact rather than three independently-planted curves.
    vz = V0 + ACCEL * t
    z = Z0 + V0 * t + 0.5 * ACCEL * t**2
    ke_ev = 0.5 * vz**2 * HA_TO_EV
    e_bath_ha = -(ke_ev - ke_ev[0]) / HA_TO_EV       # exact anti-correlation

    pd.DataFrame({
        "step": step, "time_au": t,
        "proj_z": z, "proj_vz": vz,
        "energy_proj_ke": ke_ev / HA_TO_EV,
        "energy_proj_bg_ideal": -0.265,
        "proj_z_unwrapped": z, "n_wraps": 0,
        "proj_x": 0.0, "proj_y": 0.0,
        "force_x": 0.0, "force_y": 0.0, "force_z": -0.001,
    }).to_csv(obs / "projectile.csv", index=False)

    _write_observables(obs, step, t, energy_total=e_bath_ha)
    _write_interactions(obs, step, t)
    (root / name / "run_summary.txt").write_text(
        f"run_completed = true\nn_steps = {N_STEPS}\n")
    return root / name


def _write_observables(obs: Path, step, t, energy_total=None):
    n = len(step)
    tot = np.zeros(n) if energy_total is None else np.asarray(energy_total)
    pd.DataFrame({
        "step": step, "time_au": t, "energy_total": tot,
        "energy_kinetic": 0.0, "energy_hartree": 0.0, "energy_xc": 0.0,
        "energy_external": 0.0, "energy_nonlocal": 0.0, "energy_ion": 0.0,
        "energy_ion_kinetic": 0.0, "energy_exact_exchange": 0.0,
        "energy_nvxc": 0.0, "energy_eigenvalues": 0.0,
        "current_x": 0.0, "current_y": 0.0, "current_z": 0.0, "density_l2": 0.0,
    }).to_csv(obs / "observables.csv", index=False)


def _write_interactions(obs: Path, step, t):
    pd.DataFrame({
        "step": step, "time_au": t,
        "e_ss": 1.0 + 0.01 * t, "e_pp": 0.5, "e_ps": -0.2 - 0.005 * t,
        "e_sb": -2.0, "e_pb": -0.3, "e_bb": 3.0,
        "e_hartree_check": 1.0, "e_external_check": 1.0,
        "norm_proj": 1.0, "norm_electrons": 160.0,
    }).to_csv(obs / "interactions.csv", index=False)


def _fresh_refined():
    """Import ``refined`` with its module-level path constants re-evaluated.

    ``CHAN_WP_RESULTS`` / ``CHAN_CL_RESULTS`` are read ONCE, at
    ``channeling_stopping`` import time, into module constants. So a cached
    module keeps whichever tree it first saw, and monkeypatch restoring the
    environment at teardown does NOT undo that. Dropping both modules from
    ``sys.modules`` is what actually re-reads the environment.

    This is not hypothetical: on first write, the production tests below ran
    against the SYNTHETIC tree left cached by an earlier test and "measured" a
    0.077 eV energy loss on a run that really loses 5.13 eV. The assertion
    failed, which is the only reason it was noticed — a looser bound would have
    passed and the suite would have been silently testing nothing.
    """
    for mod in ("refined", "channeling_stopping"):
        sys.modules.pop(mod, None)
    import refined
    return refined


@pytest.fixture()
def synthetic(tmp_path, monkeypatch):
    """Point refined.py at a synthetic results tree and re-import it."""
    wp_root = tmp_path / "wp_results"
    cl_root = tmp_path / "cl_results"
    _write_wp_run(wp_root)
    _write_cl_run(cl_root)
    monkeypatch.setenv("CHAN_WP_RESULTS", str(wp_root))
    monkeypatch.setenv("CHAN_CL_RESULTS", str(cl_root))
    yield _fresh_refined()
    # Teardown: evict the synthetic-bound modules so the next test re-reads the
    # (restored) environment rather than inheriting tmp_path.
    for mod in ("refined", "channeling_stopping"):
        sys.modules.pop(mod, None)


@pytest.fixture()
def production_refined(monkeypatch):
    """``refined`` bound to the REAL results tree, whatever ran before."""
    monkeypatch.delenv("CHAN_WP_RESULTS", raising=False)
    monkeypatch.delenv("CHAN_CL_RESULTS", raising=False)
    yield _fresh_refined()
    for mod in ("refined", "channeling_stopping"):
        sys.modules.pop(mod, None)


# ---------------------------------------------------------------------------
# Algebraic identities on synthetic data
# ---------------------------------------------------------------------------

def test_t2_minus_t1_is_exactly_the_variance_term(synthetic):
    """T2 - T1 == var(p)/2m. The split the whole study rests on."""
    R = synthetic
    df = R.wp_frame()
    expect_ev = 0.5 * df["var_p3d"].to_numpy() * R.HA_TO_EV
    got = (df["T2_total_ev"] - df["T1_drift_ev"]).to_numpy()
    assert np.allclose(got, expect_ev, rtol=0, atol=1e-12)
    assert np.allclose(df["var_term_ev"].to_numpy(), expect_ev, rtol=0, atol=1e-12)


def test_t2_reconstruction_matches_inq_kinetic_energy(synthetic):
    """Our T2 from moments == the engine's own e_kin_ha, to machine precision."""
    df = synthetic.wp_frame()
    assert np.allclose(df["T2_total_ev"], df["e_kin_inq_ev"], rtol=0, atol=1e-10)


def test_label_swap_is_the_users_convention_not_the_engines(synthetic):
    """T1 must be the SMALLER (drift-only) branch; T2 the total.

    Guards the exact confusion refined.py's docstring warns about: if someone
    'fixes' the naming back to ks_stopping.py's convention, T1 becomes the total
    and this fails.
    """
    df = synthetic.wp_frame()
    assert (df["T2_total_ev"] > df["T1_drift_ev"]).all()
    assert df["var_term_ev"].min() > 0.0


def test_p_integral_path_recovers_a_known_track(synthetic):
    """Cumulative trapezoid of <p_z> is EXACT for the track that generated it.

    The synthetic z_mean_circ was written as the same trapezoid of the same pz,
    so the two position definitions must agree to rounding. On real data they do
    not have to — and where they part company is the diagnostic.
    """
    df = synthetic.wp_frame()
    assert np.allclose(df["s_centroid"], df["s_pintegral"], atol=1e-10)
    assert abs(df["s_pintegral"].iloc[0] - Z0) < 1e-12


def test_classical_energy_budget_closes(synthetic):
    """d(E_bath) + d(KE_proj) == 0 for a run built to conserve energy."""
    df = synthetic.cl_frame()
    assert np.abs(df["closure_ev"]).max() < 1e-9
    assert df["d_ke_ev"].iloc[0] == 0.0
    assert df["d_e_total_ev"].iloc[0] == 0.0
    assert df["d_ke_ev"].iloc[-1] < 0.0        # the projectile LOSES energy


def test_interaction_deltas_start_at_zero(synthetic):
    for half in ("wp", "classical"):
        df = synthetic.interactions(half)
        for term in synthetic.INTERACTION_TERMS:
            col = f"d_{term}_ev"
            assert col in df.columns, f"{col} missing for {half}"
            assert abs(df[col].iloc[0]) < 1e-12


def test_fit_recovers_a_planted_stopping_power(synthetic):
    """fit_in_window on the synthetic run returns S_PLANT."""
    R = synthetic
    df = R.wp_frame()
    fit = R.fit_in_window(df["s_pintegral"].to_numpy(),
                          df["T1_drift_ev"].to_numpy(),
                          df["t"].to_numpy(), 0.0, N_STEPS * DT)
    assert fit["S"] == pytest.approx(S_PLANT, rel=1e-6)
    assert fit["r2"] > 0.999999
    assert fit["n"] == N_STEPS + 1


def test_fit_returns_nan_not_garbage_on_an_empty_window(synthetic):
    R = synthetic
    df = R.wp_frame()
    fit = R.fit_in_window(df["s_pintegral"].to_numpy(),
                          df["T1_drift_ev"].to_numpy(),
                          df["t"].to_numpy(), 100.0, 200.0)
    assert np.isnan(fit["S"]) and fit["n"] == 0


def test_momentum_slices_keeps_the_whole_k_axis(synthetic):
    """REGRESSION: the k axis must survive loading.

    ``ks_stopping._concat_segments`` ends with
    ``drop_duplicates(subset="step", keep="last")`` — correct for a scalar
    observable, catastrophic for this long-format file, where it keeps ONE bin
    per step and discards the rest. Using it here produced distributions that
    were a single point at the Nyquist bin integrating to zero, with no error
    raised (observed 2026-08-02). This test fails if that loader is ever
    reinstated.
    """
    md = synthetic.momentum_slices()
    per_step = md.groupby("step").size()
    assert per_step.nunique() == 1
    assert per_step.iloc[0] == N_BINS, "the k axis was collapsed on load"
    # And the distribution must be a distribution: peaked near the drift
    # momentum, not pinned to an axis endpoint.
    sl = md[md.step == md.step.min()].sort_values("k")
    kpk = sl.k.to_numpy()[np.argmax(sl.n_wp.to_numpy())]
    assert abs(kpk - V0) < 0.3
    assert sl.n_wp.sum() == pytest.approx(1.0, rel=1e-9)


def test_momentum_slices_rejects_a_collapsed_k_axis(synthetic, tmp_path):
    """The guard must FAIL LOUDLY, not return a one-bin frame."""
    obs = (tmp_path / "wp_results" / "wp" / "raw" / "observables")
    md = pd.read_csv(obs / "momentum_distribution.csv")
    md.drop_duplicates(subset="step", keep="last").to_csv(
        obs / "momentum_distribution.csv", index=False)
    with pytest.raises(ValueError, match="collapsed"):
        synthetic.momentum_slices()


def test_nearest_slices_snaps_and_does_not_interpolate(synthetic):
    md = pd.DataFrame({"step": [0, 15, 30] * 3,
                       "t": [0.0, 0.3, 0.6] * 3,
                       "k": [0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3],
                       "n_wp": range(9), "n_total": range(9)})
    got = synthetic.nearest_slices(md, [0.0, 0.29, 5.0])
    assert [t for t, _ in got] == [0.0, 0.3, 0.6]
    for _, sl in got:
        assert sl["k"].is_monotonic_increasing


# ---------------------------------------------------------------------------
# The same identities, on the real production run
# ---------------------------------------------------------------------------

def _production_available() -> bool:
    """Is the real results tree present? Resolved from the PATH, not an import.

    Deliberately does not import ``refined`` — importing it here would cache a
    module bound to whatever environment happened to be live during collection,
    which is the exact failure ``_fresh_refined`` exists to prevent.
    """
    root = (PKG.parents[1] / "scripts" / "channeling_twin"
            / "wp" / "results" / "wp" / "raw" / "observables")
    return root.is_dir()


production = pytest.mark.skipif(
    not _production_available(),
    reason="production channeling results not on this filesystem")


@production
def test_production_wp_identity_and_free_reference(production_refined):
    R = production_refined
    df = R.wp_frame()
    expect = 0.5 * df["var_p3d"].to_numpy() * R.HA_TO_EV
    assert np.allclose(df["var_term_ev"], expect, atol=1e-10)
    assert np.allclose(df["T2_total_ev"], df["e_kin_inq_ev"], atol=1e-8)
    # At t = 0 the packet is the launched Gaussian: var(p) is its free value in
    # every direction, so T2 - T1 is the localisation energy 3/(4 sigma^2).
    assert df["var_term_ev"].iloc[0] == pytest.approx(R.T2_MINUS_T1_FREE_EV, rel=2e-3)
    assert df["T1_drift_ev"].iloc[0] == pytest.approx(50.0, rel=1e-4)


@production
def test_production_classical_budget_closes(production_refined):
    R = production_refined
    df = R.cl_frame()
    # The real integrator, not a constructed identity: this is the physics gate.
    # Measured 2026-08-02: max |closure| = 2.2e-5 eV over 1501 steps.
    assert np.abs(df["closure_ev"]).max() < 1e-3
    assert df["d_ke_ev"].iloc[-1] == pytest.approx(-5.1256, abs=1e-3)
    assert df["d_e_total_ev"].iloc[-1] == pytest.approx(+5.1256, abs=1e-3)


@production
def test_production_momentum_distribution_is_a_distribution(production_refined):
    """The same regression, on the real file: 128 bins per step, peaked at k0."""
    R = production_refined
    md = R.momentum_slices()
    per_step = md.groupby("step").size()
    assert per_step.nunique() == 1 and per_step.iloc[0] == 128
    sl = md[md.step == 0].sort_values("k")
    k, n = sl.k.to_numpy(), sl.n_wp.to_numpy()
    # The launched packet is a Gaussian centred on k0 = v0 = 1.917 a.u.
    assert abs(k[np.argmax(n)] - R.CS.V0) < 0.15
    assert n.sum() == pytest.approx(1.0, rel=1e-3)


@production
def test_production_unwrapped_path_starts_at_the_launch_point(production_refined):
    """Pins the proj_z_unwrapped one-step lag documented in refined.cl_frame.

    The run writes proj_z_unwrapped[i] == proj_z[i+1], so its first entry is
    already one step downrange. Our derived path must start at the true launch
    point instead, and must stay a fixed 0.038 Bohr (= v dt) behind the raw
    column — fixed, which is why no published slope was affected.
    """
    R = production_refined
    df = R.cl_frame()
    assert df["z_unwrapped"].iloc[0] == pytest.approx(-28.0, abs=1e-9)
    assert df["z_unwrapped_raw"].iloc[0] != pytest.approx(-28.0, abs=1e-6)
    lag = (df["z_unwrapped_raw"] - df["z_unwrapped"]).to_numpy()
    assert lag.min() > 0.03 and lag.max() < 0.04


# ---------------------------------------------------------------------------
# 2-D momentum map helpers and the twin comparison (added 2026-08-02)
# ---------------------------------------------------------------------------

def test_kz_asymmetry_is_exactly_neutral_on_a_symmetric_distribution(synthetic):
    """A symmetric input must give skew 0 and frac_above 0.5.

    This is the test the NAIVE implementation fails. Summing grid points with
    ``kz > mean`` returns 0.454 on the real t=0 packet -- which is EXACTLY
    symmetric -- because only ~8 k_z points carry it and the mean falls between
    them. The CDF interpolation is what makes the metric usable.
    """
    R = synthetic
    kz = np.linspace(1.0, 3.0, 21)
    p = np.exp(-((kz - 2.0) ** 2) / (2 * 0.2**2)); p /= p.sum()
    a = R.kz_asymmetry(kz, p)
    assert a["mean_kz"] == pytest.approx(2.0, abs=1e-9)
    assert a["skewness"] == pytest.approx(0.0, abs=1e-9)
    assert a["frac_above_mean"] == pytest.approx(0.5, abs=1e-3)
    assert a["median_minus_mean"] == pytest.approx(0.0, abs=1e-3)


def test_kz_asymmetry_detects_a_planted_skew(synthetic):
    """A deliberately one-sided distribution must report the right SIGN."""
    R = synthetic
    kz = np.linspace(0.0, 6.0, 121)
    # exponential tail toward HIGH k -> positive skew, mass below the mean
    p = np.where(kz >= 2.0, np.exp(-(kz - 2.0) / 0.4), 0.0); p /= p.sum()
    a = R.kz_asymmetry(kz, p)
    assert a["skewness"] > 1.0
    assert a["frac_above_mean"] < 0.5, "a right-skewed distribution has most mass BELOW the mean"
    assert a["median_minus_mean"] < 0.0


def test_kz_asymmetry_accepts_a_2d_map(synthetic):
    """Passing the full (k_z, k_perp) map must marginalise, not error."""
    R = synthetic
    kz = np.linspace(1.0, 3.0, 21)
    p1 = np.exp(-((kz - 2.0) ** 2) / (2 * 0.2**2)); p1 /= p1.sum()
    P = np.outer(p1, np.array([0.5, 0.3, 0.2]))
    assert R.kz_asymmetry(kz, P)["mean_kz"] == pytest.approx(
        R.kz_asymmetry(kz, p1)["mean_kz"], abs=1e-12)


def test_impulse_ratio_is_one_for_identical_twins(synthetic):
    """If both halves decelerate identically the ratio must be exactly 1."""
    R = synthetic
    wp, cl = R.wp_frame(), R.cl_frame()
    imp = R.impulse_comparison(wp, cl)
    # both synthetic halves use the SAME constant deceleration by construction
    good = np.isfinite(imp["impulse_ratio"]) & (np.abs(imp["dp_cl"]) > 1e-9)
    assert good.sum() > 10
    assert np.allclose(imp["impulse_ratio"][good], 1.0, atol=1e-6)


def test_impulse_ratio_scales_with_a_weakened_drag(synthetic, monkeypatch):
    """Halving the WP's impulse must show up as a ratio of 0.5, not 0.25.

    Guards against the ratio being taken on the ENERGY (which goes as p^2) by
    mistake -- that would report ~0.25 here.
    """
    R = synthetic
    wp, cl = R.wp_frame(), R.cl_frame()
    wp = wp.copy()
    wp["pz"] = wp["pz"].iloc[0] + 0.5 * (wp["pz"] - wp["pz"].iloc[0])
    imp = R.impulse_comparison(wp, cl)
    good = np.isfinite(imp["impulse_ratio"]) & (np.abs(imp["dp_cl"]) > 1e-9)
    assert np.allclose(imp["impulse_ratio"][good], 0.5, atol=1e-6)


def test_combined_projectile_coupling_is_the_sum_of_the_two_terms(synthetic):
    R = synthetic
    for half in ("wp", "classical"):
        d = R.interactions(half)
        got = R.combined_projectile_coupling(d)
        assert np.allclose(got, d["d_e_ps_ev"] + d["d_e_pb_ev"])
        assert got[0] == pytest.approx(0.0, abs=1e-12)


@production
def test_production_momentum_map_round_trips_to_recorded_moments(production_refined):
    """The 2-D map from the orbital dump must reproduce wp_momentum_stats.

    This is the check that licenses everything read off the map: if the FFT
    ordering (ifftshift) were wrong, <k_z> would be shifted by a phase ramp and
    the whole picture would be plausible and wrong.
    """
    R = production_refined
    steps = R.available_wf_steps()
    if not steps:
        pytest.skip("no wavefunction dumps")
    kz, kperp, P = R.momentum_map(step=0)
    a = R.kz_asymmetry(kz, P)
    wp = R.wp_frame()
    assert a["mean_kz"] == pytest.approx(float(wp["pz"].iloc[0]), rel=1e-6)
    assert a["sigma_kz"] ** 2 == pytest.approx(float(wp["var_pz"].iloc[0]), rel=1e-4)
    # the launched packet is EXACTLY symmetric
    assert a["skewness"] == pytest.approx(0.0, abs=1e-3)
    assert a["frac_above_mean"] == pytest.approx(0.5, abs=5e-3)


@production
def test_production_combined_coupling_agrees_far_better_than_its_parts(production_refined):
    """delta(E_PS+E_PB) matches between the halves though the terms do not.

    The physical claim: E_PS + E_PB is the projectile's coupling to the NET
    (neutral) charge density, which is gauge-clean; the split into bath and
    background is not. Pinned with real numbers so a regression in either term
    is caught.
    """
    R = production_refined
    iwp, icl = R.interactions("wp"), R.interactions("classical")
    swp = R.combined_projectile_coupling(iwp)
    scl = R.combined_projectile_coupling(icl)
    n = min(len(swp), len(scl))
    worst_sum = float(np.abs(swp[:n] - scl[:n]).max())
    worst_ps = float(np.abs(iwp["d_e_ps_ev"].to_numpy()[:n]
                            - icl["d_e_ps_ev"].to_numpy()[:n]).max())
    assert worst_sum < 0.35, f"combined coupling diverged: {worst_sum:.3f} eV"
    assert worst_ps > 2.0, "E_PS alone should differ a lot; if not, the premise changed"
    assert worst_sum < 0.2 * worst_ps
