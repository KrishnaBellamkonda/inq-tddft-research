"""Known-case tests for the linear-response residual / form-factor kernel.

Validation-gates: every branch of formfactor_residual is checked against a
synthetic ground truth where the answer is known analytically.

Run: venv/bin/python3 -m pytest docs/campaigns/ml-patterns/tests/test_formfactor_residual.py -q
"""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))          # campaign dir (kernels package)
from kernels import formfactor_residual as FR       # noqa: E402


def _pair(sigma, M=40, Nq=30, dq=0.13, drift=0.0, excess_amp=0.0,
          noise_level=1e-6, seed=0):
    """Synthetic matched pair: n_WP = F(q)*n_cl*(1+drift*t) + high-q excess."""
    rng = np.random.default_rng(seed)
    q = dq * np.arange(Nq)
    t = np.linspace(0.0, 5.0, M)
    # classical shell amplitude: smooth-ish positive signal, decaying in q
    base = np.exp(-0.3 * (q * 1.0) ** 2) + 0.05
    ncl = np.abs(base[None, :] * (1.0 + 0.1 * rng.standard_normal((M, Nq))))
    F = FR.form_factor(q, sigma)
    nwp = F[None, :] * ncl * (1.0 + drift * t[:, None])
    if excess_amp:
        hi = q >= np.median(q[q > 0])
        nwp[:, hi] += excess_amp * base[hi][None, :]
    noise = np.full(Nq, noise_level)
    return FR.PairResult(sigma=sigma, q=q, ncl_t=ncl, nwp_t=nwp,
                         noise_cl=noise, tcommon=t)


def test_exact_form_factor_recovers_width():
    """n_WP = F(q)*n_cl exactly -> sigma_fit ~ sigma_WP, flat, no excess."""
    r = FR.residual_test(_pair(sigma=0.5), snr=3.0)
    assert abs(r["sigma_fit"] - 0.5) < 0.05, r["sigma_fit"]
    assert r["matches_sigma_wp"] is True
    assert r["t_flatness"] < 0.05, r["t_flatness"]
    assert abs(r["highq_excess_over_noise"]) < 5.0 or not np.isfinite(
        r["highq_excess_over_noise"])
    assert r["fit_r2"] > 0.98, r["fit_r2"]


def test_high_q_excess_detected():
    """Injecting a high-q excess (nonlinear/quantum fingerprint) is flagged."""
    clean = FR.residual_test(_pair(sigma=0.5, excess_amp=0.0), snr=3.0)
    dirty = FR.residual_test(_pair(sigma=0.5, excess_amp=0.2), snr=3.0)
    assert dirty["highq_excess_over_noise"] > clean["highq_excess_over_noise"]
    assert dirty["highq_excess_over_noise"] > 10.0, dirty["highq_excess_over_noise"]


def test_time_drift_breaks_flatness():
    """A t-drift in the ratio (deceleration mismatch) inflates t-flatness."""
    flat = FR.residual_test(_pair(sigma=0.5, drift=0.0), snr=3.0)
    drift = FR.residual_test(_pair(sigma=0.5, drift=0.15), snr=3.0)
    assert drift["t_flatness"] > 5.0 * flat["t_flatness"]
    assert drift["t_flatness"] > 0.1, drift["t_flatness"]


def test_fork_a_selects_sigma_wp():
    """a(sigma) with a=0.5 sigma^2 -> slope 0.5 -> selects sigma_WP."""
    per = [FR.residual_test(_pair(sigma=s), snr=3.0) for s in (0.5, 1.0, 2.0)]
    col = FR.collapse_fork_a(per)
    assert col["ok"]
    assert abs(col["slope"] - 0.5) < 0.05, col["slope"]
    assert col["selects"] == "sigma_WP"


def test_fork_a_selects_sigma_pot():
    """If the true filter used sigma_pot=sigma/sqrt2, a=0.25 sigma^2 -> slope 0.25."""
    # emulate: label sigma but build F with sigma_pot=sigma/sqrt2
    per = []
    for s in (0.5, 1.0, 2.0):
        p = _pair(sigma=s)
        # rebuild nwp using the narrower potential width
        q = p.q
        Fpot = FR.form_factor(q, s / np.sqrt(2.0))
        p.nwp_t = Fpot[None, :] * p.ncl_t
        per.append(FR.residual_test(p, snr=3.0))
    col = FR.collapse_fork_a(per)
    assert col["ok"]
    assert abs(col["slope"] - 0.25) < 0.05, col["slope"]
    assert col["selects"] == "sigma_pot"


def test_gaussian_exponent_fit_exact():
    """fit_gaussian_exponent recovers a from a clean exp(-a q^2)."""
    q = 0.13 * np.arange(1, 30)
    a_true = 0.42
    ratio = np.exp(-a_true * q**2)
    fit = FR.fit_gaussian_exponent(q, ratio)
    assert abs(fit["a"] - a_true) < 1e-6, fit["a"]
    assert fit["r2"] > 0.999


def test_radial_spectrum_localizes_plane_wave():
    """A single cosine at wavevector k puts its power in the |q|=|k| shell."""
    n = 32; dx = 0.4
    x = np.arange(n) * dx
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    kz = 5 * (2 * np.pi / (n * dx))               # 5th z-mode
    field = np.cos(kz * Z).astype(np.float32)
    q, amp, noise, count = FR.radial_spectrum(field, dx)
    peak = q[np.argmax(amp)]
    assert abs(peak - kz) < (2 * np.pi / (n * dx)), (peak, kz)
