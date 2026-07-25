"""Known-case tests for inqview.analysis.lindhard_elf (corrected Lindhard ELF).

The companion module fixes the missing-plasmon f-sum failure of
inqview.pipeline.lindhard by using the full complex-argument Lindhard function.

Gates:
  (a) Thomas-Fermi static limit: eps(q->0, 0) -> 1 + k_TF^2/q^2.
  (b) f-sum rule: int_0^inf omega Im[-1/eps] domega = (pi/2) omega_p^2 at all q
      (the gate the old module fails for q < 0.5).
  (c) Im eps non-zero only inside the e-h continuum (plus the small-eta plasmon).
  (d) S_LR(v; sigma) q-grid convergence < 0.5%.

Run:
    cd /local/data/public/skcb2/tddft && venv/bin/python3 -m pytest \
      inq-stack/tests/python/inqview/analysis/test_lindhard_elf.py -v
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.analysis import lindhard_elf as E

RS = 5.69
KF = E.kF_from_rs(RS)
WP = E.omega_p(KF)
TARGET = 0.5 * np.pi * WP ** 2


def test_thomas_fermi_static_limit():
    """eps(q->0, omega=0) -> 1 + k_TF^2/q^2."""
    kTF2 = E.k_TF(KF) ** 2
    # q -> 0 at omega=0; use eta << q vF so the regulator u = i eta/(q vF) -> 0.
    for q in (1e-2, 2e-2, 3e-2):
        eps = float(E.epsilon_rpa(q, 0.0, KF, eta=1e-9).real)
        expected = 1.0 + kTF2 / q ** 2
        assert abs(eps - expected) / expected < 5e-3, (q, eps, expected)


@pytest.mark.parametrize("q", [0.1, 0.2, 0.3, 0.5, 0.8, 1.2])
def test_f_sum_rule(q):
    """int omega Im[-1/eps] domega = (pi/2) wp^2 at all q (was broken < 0.5)."""
    eta = max(1e-3, 2e-3 * q)
    w = np.linspace(1e-5, 3.0, 600000)
    Lw = E.loss_function(np.full_like(w, q), w, KF, eta=eta)
    integ = np.trapezoid(w * Lw, w)
    assert abs(integ - TARGET) / TARGET < 1e-2, (q, integ / TARGET)


def test_imag_eps_only_in_continuum():
    """Im eps = 0 outside the e-h continuum and the plasmon (eta->0)."""
    # well above the continuum upper edge (q vF + q^2/2) and away from plasmon:
    q = 0.3
    w_far = 2.5  # >> wp and >> continuum edge
    im_eps = float(E.epsilon_rpa(q, w_far, KF, eta=1e-4).imag)
    assert abs(im_eps) < 1e-3, im_eps
    # inside continuum (low omega, single-particle) Im eps > 0:
    im_in = float(E.epsilon_rpa(q, 0.5 * q * KF, KF, eta=1e-4).imag)
    assert im_in > 0.0


def test_stopping_qgrid_convergence():
    """S_LR(v; sigma) stable < 0.5% under q-grid refinement."""
    v, sigma = 1.0, 0.5
    # adaptive omega resolution (tied to eta) is on by default; refine q only.
    s_coarse = E.stopping_power_sigma(v, KF, sigma, n_q=600)
    s_fine = E.stopping_power_sigma(v, KF, sigma, n_q=1200)
    assert s_fine > 0.0
    assert abs(s_fine - s_coarse) / s_fine < 5e-3, (s_coarse, s_fine)


def test_stopping_positive_and_monotone_lowv():
    """Low-v friction: S(v) ~ Q v, so S increasing and positive for small v."""
    sigma = 0.5
    s1 = E.stopping_power_sigma(0.2, KF, sigma)
    s2 = E.stopping_power_sigma(0.4, KF, sigma)
    assert 0.0 < s1 < s2


@pytest.mark.parametrize("v", [0.2, 0.62, 1.0, 1.94, 2.98])
def test_stopping_point_converged(v):
    """THE single point-charge reference is converged: stable < 1% under both
    qmax-margin and q-grid refinement (the natural kinematic cutoff means no
    1/sigma blow-up). Guards the one analytical curve overlaid on every plot."""
    s_a = E.stopping_power_point(v, KF, margin=2.0, n_q=2000)
    s_b = E.stopping_power_point(v, KF, margin=6.0, n_q=8000)
    assert s_b > 0.0
    assert abs(s_a - s_b) / s_b < 1e-2, (v, s_a, s_b)


def test_stopping_point_above_finite_sigma():
    """A bare point charge couples to all q, so it must stop at least as hard as
    any finite-width (form-factor-suppressed) projectile at the same v."""
    v = 1.0
    s_point = E.stopping_power_point(v, KF)
    for sigma in (0.2, 0.35, 0.5):
        assert s_point > E.stopping_power_sigma(v, KF, sigma), sigma
