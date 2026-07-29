"""Tests for WP-integrity metrics (``inqview.analysis.wp_integrity``).

Pure numpy with analytic expectations (IV-M05): KL of known distributions, IPR
ordering (localised > delocalised) + a closed-form value, and density-weighted
variance recovering a known Gaussian width.
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.analysis import wp_integrity as wi

pytestmark = pytest.mark.analysis


# --- momentum KL -----------------------------------------------------------
def test_kl_self_zero_and_known_value():
    p = np.array([0.25, 0.25, 0.5])
    assert wi.momentum_kl(p, p) == pytest.approx(0.0, abs=1e-12)
    a = np.array([0.5, 0.5])
    b = np.array([0.25, 0.75])
    expected = 0.5 * np.log(0.5 / 0.25) + 0.5 * np.log(0.5 / 0.75)
    assert wi.momentum_kl(a, b) == pytest.approx(expected, rel=1e-10)


def test_kl_nonnegative_unnormalised_inputs():
    """Accepts raw histograms (normalises internally); KL ≥ 0."""
    rng = np.random.default_rng(1)
    for _ in range(10):
        assert wi.momentum_kl(rng.random(6), rng.random(6)) >= -1e-12


# --- IPR (localisation) ----------------------------------------------------
def test_ipr_localised_greater_than_delocalised():
    spike = np.zeros(100); spike[50] = 1.0
    uniform = np.ones(100)
    assert wi.ipr(spike) > wi.ipr(uniform)


def test_ipr_uniform_closed_form():
    """Uniform over N cells (ρ=1 each, dV=1): IPR = N·1 / (N)² = 1/N."""
    N = 50
    assert wi.ipr(np.ones(N), dV=1.0) == pytest.approx(1.0 / N, rel=1e-12)


# --- real-space spread -----------------------------------------------------
def test_variance_recovers_gaussian_width():
    """A sampled Gaussian density of std σ returns variance ≈ σ²."""
    x = np.linspace(-20, 20, 2001)
    sigma = 2.5
    rho = np.exp(-0.5 * (x / sigma) ** 2)
    assert wi.real_space_variance(rho, x) == pytest.approx(sigma ** 2, rel=1e-3)


def test_variance_zero_density_is_zero():
    x = np.linspace(-5, 5, 11)
    assert wi.real_space_variance(np.zeros_like(x), x) == 0.0


# --- KL time series (drift-from-launch + frame-to-frame rate, IV-M05) ------
def _drift_series():
    base = np.array([0.6, 0.3, 0.1])
    drift = np.array([0.1, 0.3, 0.6])
    return np.array([(1 - a) * base + a * drift for a in np.linspace(0, 1, 5)])


def test_kl_series_from_initial_starts_zero_and_rises():
    kl = wi.kl_series(_drift_series(), reference="initial")
    assert kl[0] == pytest.approx(0.0, abs=1e-12)
    assert np.all(np.diff(kl) > 0)              # drifts away from launch


def test_kl_series_previous_zero_for_constant():
    P = np.tile([0.2, 0.5, 0.3], (4, 1))
    kl = wi.kl_series(P, reference="previous")
    assert np.allclose(kl, 0.0, atol=1e-12)     # steady WP → no frame-to-frame drift


def test_kl_series_previous_positive_when_changing():
    kl = wi.kl_series(_drift_series(), reference="previous")
    assert kl[0] == pytest.approx(0.0, abs=1e-12)
    assert np.all(kl[1:] > 0)


def test_kl_series_invalid_reference_raises():
    with pytest.raises(ValueError):
        wi.kl_series(np.ones((2, 3)), reference="bogus")


# --- dataclass container ---------------------------------------------------
def test_wpintegrity_holds_series():
    t = np.array([0.0, 1.0])
    w = wi.WPIntegrity(time_au=t, kl_mom=np.array([0.0, 0.3]),
                       sigma_r=np.array([1.4, 1.8]), ipr=np.array([0.2, 0.1]))
    assert w.time_au.shape == (2,)
    assert w.ipr[1] < w.ipr[0]          # delocalised over time


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
