"""Known-case tests for ks_stopping.py (bulk-jellium KS stopping engine).

Non-circular by construction: every expected value is an analytic constant
chosen BEFORE the data is synthesised, not read back from the implementation.

Run:
    venv/bin/python -m pytest \
      ResearchProject/systems/jellium/hypotheses/bulk_ks_stopping/tests/ -v
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ks_stopping import (  # noqa: E402
    HA_TO_EV, ClassicalRun, Interactions, PairPhase, StoppingFit, WPRun,
    fit_stopping, local_stopping, unwrap_periodic, _ols,
)


# ---------------------------------------------------------------------------
# Helpers — synthesise the dataclasses directly, no file I/O
# ---------------------------------------------------------------------------

def _interactions(half: str, norm: np.ndarray, t: np.ndarray | None = None,
                  e_pp: np.ndarray | None = None) -> Interactions:
    n = len(norm)
    t = np.arange(n, dtype=float) * 0.1 if t is None else t
    z = np.zeros(n)
    return Interactions(
        run_dir=Path("/nonexistent"), half=half, t=t,
        step=np.arange(n), e_ss=z.copy(), e_ps=z.copy(),
        e_pp=z.copy() if e_pp is None else e_pp,
        e_sb=z.copy(), e_pb=z.copy(), e_bb=z.copy(),
        norm=norm, proj_z=None, closure=z.copy(),
    )


# ---------------------------------------------------------------------------
# unwrap_periodic
# ---------------------------------------------------------------------------

def test_unwrap_no_crossing_is_identity():
    """A trajectory that never touches a face must come back unchanged."""
    L = 80.0
    z = np.linspace(-32.0, 30.0, 200)
    assert np.allclose(unwrap_periodic(z, L), z, atol=1e-9)


def test_unwrap_recovers_a_face_crossing():
    """A projectile leaving at +L/2 and re-entering at -L/2 is one straight line.

    Build the continuous truth first, then wrap it into (-L/2, L/2] the way the
    circular estimator reports it, and require the unwrap to invert that.
    """
    L = 80.0
    truth = np.linspace(-32.0, 55.0, 400)          # crosses +40 partway through
    wrapped = (truth + L / 2) % L - L / 2          # what the run would report
    got = unwrap_periodic(wrapped, L)
    got = got - got[0] + wrapped[0]                # pin the branch, as the loader does
    assert np.allclose(got, truth, atol=1e-9)


def test_unwrap_survives_two_crossings():
    L = 80.0
    truth = np.linspace(-32.0, 130.0, 800)         # wraps twice
    wrapped = (truth + L / 2) % L - L / 2
    got = unwrap_periodic(wrapped, L)
    got = got - got[0] + wrapped[0]
    assert np.allclose(got, truth, atol=1e-9)


# ---------------------------------------------------------------------------
# OLS
# ---------------------------------------------------------------------------

def test_ols_recovers_an_exact_line():
    x = np.linspace(0.0, 10.0, 50)
    y = 3.0 - 0.25 * x
    slope, stderr, r2, model = _ols(x, y)
    assert slope == pytest.approx(-0.25, abs=1e-12)
    assert r2 == pytest.approx(1.0, abs=1e-12)
    assert stderr == pytest.approx(0.0, abs=1e-10)
    assert np.allclose(model, y)


# ---------------------------------------------------------------------------
# fit_stopping — the headline arithmetic
# ---------------------------------------------------------------------------

def test_fit_stopping_recovers_a_planted_S():
    """Plant S = 0.5 eV/Bohr exactly and require it back.

    The sign convention is the point of this test: S = -dT/ds, so a projectile
    LOSING energy as it advances (dT/ds < 0) must give a POSITIVE S.
    """
    S_true_ev = 0.5                                  # eV/Bohr
    slope_ha = -S_true_ev / HA_TO_EV                 # dT/ds in Ha/Bohr
    t = np.linspace(0.0, 30.0, 601)
    v = 2.7111
    s = -32.0 + v * t
    T = 3.6749 + slope_ha * (s - s[0])

    fit = fit_stopping(s, T, t, 4.0, 19.0, "planted", v=np.full_like(t, v))
    assert fit.S_ev_per_bohr == pytest.approx(S_true_ev, abs=1e-9)
    assert fit.r2 == pytest.approx(1.0, abs=1e-9)
    assert fit.stderr < 1e-9
    assert fit.window_syst < 1e-9        # a pure line is window-independent
    assert fit.n_points > 100
    assert fit.t_window == (4.0, 19.0)
    assert fit.mean_v == pytest.approx(v, abs=1e-9)


def test_fit_stopping_sign_for_an_accelerating_projectile():
    """Energy GAIN along the path must report a NEGATIVE S, not |S|."""
    t = np.linspace(0.0, 30.0, 601)
    s = -32.0 + 2.7111 * t
    T = 3.6749 + (0.3 / HA_TO_EV) * (s - s[0])       # gaining 0.3 eV/Bohr
    fit = fit_stopping(s, T, t, 4.0, 19.0, "gain")
    assert fit.S_ev_per_bohr == pytest.approx(-0.3, abs=1e-9)


def test_fit_stopping_window_selects_the_right_span():
    """Only points inside [t0, t1] may enter the fit.

    Planted as a kink: S = 1.0 before t = 19 and a wildly different slope after.
    If the window leaked, the recovered S would be pulled away from 1.0.
    """
    t = np.linspace(0.0, 30.0, 601)
    v = 2.7111
    s = -32.0 + v * t
    T = np.where(
        t <= 19.0,
        3.6749 - (1.0 / HA_TO_EV) * (s - s[0]),
        3.6749 - (1.0 / HA_TO_EV) * (s - s[0]) - (5.0 / HA_TO_EV) * (s - (-32.0 + v * 19.0)),
    )
    fit = fit_stopping(s, T, t, 4.0, 19.0, "kinked", window_scan=0.0)
    assert fit.S_ev_per_bohr == pytest.approx(1.0, abs=1e-6)


def test_window_systematic_grows_with_curvature():
    """A curved T(s) must report a LARGER window systematic than a straight one.

    This is the guard that stops the systematic from being decorative: if the
    scan were not actually re-fitting, both cases would return the same number.
    """
    t = np.linspace(0.0, 30.0, 601)
    s = -32.0 + 2.7111 * t
    ds = s - s[0]

    straight = fit_stopping(s, 3.6749 - (0.5 / HA_TO_EV) * ds, t, 8.0, 19.0, "straight")
    curved = fit_stopping(
        s, 3.6749 - (0.5 / HA_TO_EV) * ds - (0.01 / HA_TO_EV) * ds**2,
        t, 8.0, 19.0, "curved")

    assert straight.window_syst < 1e-9
    assert curved.window_syst > 0.1
    assert curved.r2 < straight.r2


def test_uncertainty_combines_in_quadrature():
    f = StoppingFit(label="x", S_ev_per_bohr=1.0, stderr=0.3, window_syst=0.4)
    assert f.uncertainty == pytest.approx(0.5, abs=1e-12)


def test_fit_is_robust_to_noise():
    """With symmetric noise the recovered S must sit inside its own error bar."""
    rng = np.random.default_rng(20260730)
    S_true = 0.8
    t = np.linspace(0.0, 30.0, 601)
    s = -32.0 + 2.7111 * t
    T = 3.6749 - (S_true / HA_TO_EV) * (s - s[0])
    T = T + rng.normal(0.0, 1e-4, size=T.shape)

    fit = fit_stopping(s, T, t, 4.0, 19.0, "noisy")
    assert abs(fit.S_ev_per_bohr - S_true) < 3.0 * fit.uncertainty
    assert fit.r2 > 0.99


# ---------------------------------------------------------------------------
# The T1 - T2 localisation energy
# ---------------------------------------------------------------------------

def test_localisation_energy_at_t0_matches_the_analytic_value():
    """T1 - T2 = 3/(8 sigma^2) for a Gaussian of wavepacket width sigma.

    sigma_WP = 2 Bohr is this study's value, so the offset between the two KE
    definitions must start at 2.551 eV. Computed here from the momentum moments
    the run actually writes, exactly as the loader does.
    """
    sigma = 2.0
    k0 = 2.7111
    sigma_p2 = 1.0 / (4.0 * sigma**2)          # = 0.0625

    T1 = 0.5 * (k0**2 + 3.0 * sigma_p2)        # <p^2>/2
    T2 = 0.5 * k0**2                           # <p>^2/2

    assert (T1 - T2) == pytest.approx(3.0 / (8.0 * sigma**2), abs=1e-12)
    assert (T1 - T2) * HA_TO_EV == pytest.approx(2.551, abs=0.002)


# ---------------------------------------------------------------------------
# local_stopping — the rolling-OLS S(z)
# ---------------------------------------------------------------------------

def test_local_stopping_recovers_a_planted_constant_slope():
    """T(s) with a planted constant dT/ds must give that S everywhere inside."""
    S_true = 0.35                                  # eV/Bohr, chosen up front
    s = np.linspace(-30.0, 20.0, 400)
    T = 4.0 - (S_true / HA_TO_EV) * (s - s[0])     # Ha, exactly linear

    S = local_stopping(s, T, half_width=12)
    interior = S[12:-12]
    assert np.all(np.isfinite(interior))
    assert np.allclose(interior, S_true, atol=1e-9)


def test_local_stopping_tracks_a_known_linear_ramp_in_S():
    """If dT/ds varies linearly, the rolling slope must follow it.

    T = -(a s + b s^2/2)/HA_TO_EV gives S(s) = a + b s exactly. Checked at the
    window centre, where the centred estimator is unbiased for a linear S.
    """
    a, b = 0.20, 0.010                             # eV/Bohr, eV/Bohr^2
    s = np.linspace(-30.0, 20.0, 600)
    T = -(a * s + 0.5 * b * s**2) / HA_TO_EV

    S = local_stopping(s, T, half_width=15)
    mid = slice(60, -60)
    assert np.allclose(S[mid], a + b * s[mid], atol=1e-6)


def test_local_stopping_edges_are_filled_not_extrapolated():
    """The first/last half_width points repeat the nearest interior value.

    Guards the documented contract: an extrapolated edge would manufacture a
    spike at the start/end of every S(z) plot that looks like impact physics.
    """
    s = np.linspace(0.0, 10.0, 100)
    T = -np.exp(s / 5.0) / HA_TO_EV                # strongly curved
    hw = 10
    S = local_stopping(s, T, half_width=hw)

    assert np.all(S[:hw] == S[hw])
    assert np.all(S[-hw:] == S[-hw - 1])


def test_local_stopping_returns_all_nan_when_too_short():
    """Fewer than 2*half_width+1 samples cannot support a centred slope."""
    S = local_stopping(np.arange(5.0), np.arange(5.0), half_width=12)
    assert len(S) == 5 and np.all(np.isnan(S))


# ---------------------------------------------------------------------------
# Interactions — the clipping detector
# ---------------------------------------------------------------------------

def test_clip_index_finds_the_contiguous_tail():
    """Only the trailing run of clipped rows counts as clipping."""
    norm = np.ones(100)
    norm[80:] = 0.99                               # the projectile leaving the box
    ix = _interactions("classical", norm)
    assert ix.clip_index == 80
    assert ix.clip_time == pytest.approx(8.0)


def test_clip_index_ignores_early_discretisation_dips():
    """Rows a few 1e-9 under 1.0 at LAUNCH are grid discretisation, not clipping.

    This is the exact bug that mislabelled t=0 as the clipping onset on
    2026-08-01: a global min() over norm_proj picks the launch rows, which sit at
    0.999999996 because the Gaussian is discretised on the grid. Only the
    contiguous tail is the projectile meeting the face.
    """
    norm = np.ones(100)
    norm[:3] = 1.0 - 4.2e-9                        # measured launch deficit
    norm[90:] = 0.994                              # the real clipping
    ix = _interactions("classical", norm)
    assert ix.clip_index == 90


def test_clip_time_is_inf_when_never_clipped():
    ix = _interactions("classical", np.ones(50))
    assert ix.clip_index is None
    assert ix.clip_time == float("inf")


def test_wp_half_never_reports_clipping():
    """A wavepacket has no rigid cloud to clip, whatever its norm does."""
    ix = _interactions("wp", np.full(50, 0.8))
    assert ix.clip_index is None
    assert ix.clip_time == float("inf")


def test_in_window_truncates_at_the_clipping_onset():
    """A window extending past the onset must be cut back to it."""
    norm = np.ones(100)
    norm[60:] = 0.99                               # onset at t = 6.0
    ix = _interactions("classical", norm)

    m = ix.in_window(1.0, 9.0)                     # asks for more than is clean
    assert ix.t[m].max() <= 6.0
    assert ix.t[m].min() >= 1.0


# ---------------------------------------------------------------------------
# PairPhase — the divergence marker
# ---------------------------------------------------------------------------

def _pair(v_wp: np.ndarray, v_cl: np.ndarray, t: np.ndarray) -> PairPhase:
    n = len(t)
    z = np.zeros(n)
    wp = WPRun(run_dir=Path("/nonexistent"), box_length_z=80.0, t=t,
               step=np.arange(n), T1=z.copy(), T2=z.copy(), pz=v_wp,
               s3=z.copy(), s3_naive=z.copy(), s4=z.copy(),
               norm=np.ones(n), sigma_z=z.copy(), parseval=z.copy())
    cl = ClassicalRun(run_dir=Path("/nonexistent"), t=t, step=np.arange(n),
                      z=z.copy(), vz=v_cl, T=z.copy())
    return PairPhase(family="synthetic", wp=wp, cl=cl)


def test_divergence_fires_at_the_planted_threshold():
    """Two velocity histories separating linearly cross 5% of v0 at a known t."""
    t = np.linspace(0.0, 10.0, 1001)
    v0 = 2.0
    v_cl = np.full_like(t, v0)
    v_wp = v0 + 0.01 * t                # gap = 0.01 t; 5% of v0 = 0.1 at t = 10
    d = _pair(v_wp, v_cl, t).divergence(frac=0.05)
    assert d["t"] == pytest.approx(10.0, abs=0.02)


def test_divergence_is_nan_when_they_never_part():
    """Identical histories must report NaN, not a spurious index 0."""
    t = np.linspace(0.0, 10.0, 201)
    v = np.full_like(t, 2.0)
    d = _pair(v, v, t).divergence(frac=0.05)
    assert np.isnan(d["t"])


def test_divergence_uses_absolute_gap_in_both_directions():
    """A wavepacket that is SLOWER must trigger just as one that is faster."""
    t = np.linspace(0.0, 10.0, 1001)
    v0 = 2.0
    fast = _pair(v0 + 0.01 * t, np.full_like(t, v0), t).divergence(frac=0.05)
    slow = _pair(v0 - 0.01 * t, np.full_like(t, v0), t).divergence(frac=0.05)
    assert fast["t"] == pytest.approx(slow["t"], abs=1e-9)
