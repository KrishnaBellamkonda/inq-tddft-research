"""Tests for the plasmon spectral locator (``inqview.analysis.plasmon_spectrum``).

Pure numpy with the analytic reduced system from the formula-validation dossier:
an undamped-plasmon phasor n_q(t)=A·e^{−iω_p t} → a δ-peak at ω_p, with exact
1/q² scaling. Also checks the axial-mode extraction picks the right q-mode.
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.analysis import plasmon_spectrum as ps

pytestmark = pytest.mark.analysis


def _on_bin_omega(k, nt, dt):
    return 2.0 * np.pi * k / (nt * dt)          # angular freq exactly on bin k


def test_phasor_peaks_at_omega_p_for_all_q():
    """Undamped plasmon e^{−iω_p t} → spectral peak at ω_p for every q-mode."""
    nt, dt = 200, 0.5
    omega_p = _on_bin_omega(10, nt, dt)
    t = np.arange(nt) * dt
    q_vals = np.array([0.2, 0.4, 0.6])
    n_q_t = (np.exp(-1j * omega_p * t)[:, None] * np.ones(q_vals.size))
    spec = ps.spectrum_from_nq(n_q_t, dt, q_vals, window=False)
    for j in range(q_vals.size):
        assert spec.peak_omega[j] == pytest.approx(omega_p, rel=1e-6)


def test_one_over_q_squared_weighting():
    """Equal-amplitude phasors → power equal across q, so loss·q² is constant."""
    nt, dt = 200, 0.5
    omega_p = _on_bin_omega(10, nt, dt)
    t = np.arange(nt) * dt
    q_vals = np.array([0.2, 0.4, 0.6])
    n_q_t = np.exp(-1j * omega_p * t)[:, None] * np.ones(q_vals.size)
    spec = ps.spectrum_from_nq(n_q_t, dt, q_vals, window=False)
    peak_bin = int(np.argmax(spec.power[:, 0]))
    loss_q2 = spec.loss[peak_bin, :] * q_vals ** 2     # = power, q-independent
    assert np.allclose(loss_q2, loss_q2[0], rtol=1e-9)


def test_complex_fft_separates_plus_and_minus_omega():
    """e^{−iω_p t} has weight at ONE sign of ω (IV-E01: not folded)."""
    nt, dt = 200, 0.5
    omega_p = _on_bin_omega(10, nt, dt)
    t = np.arange(nt) * dt
    n_q_t = np.exp(-1j * omega_p * t)[:, None]
    spec = ps.spectrum_from_nq(n_q_t, dt, np.array([0.3]), window=False)
    neg = spec.power[spec.omega < 0, 0].max()
    pos = spec.power[spec.omega > 0, 0].max()
    assert min(neg, pos) < 1e-6 * max(neg, pos)        # power on ONE side only


def test_axial_extraction_picks_the_right_mode():
    """δn = cos(2π·m·z/nz) along z (last axis, inqkit z-fastest layout) →
    only axial mode m is populated."""
    nx = ny = nz = 16
    Lz = 8.0
    m_true = 3
    z = np.arange(nz)
    # variation along the LAST axis (z); inqkit layout is (nx, ny, nz)
    dn = np.cos(2.0 * np.pi * m_true * z / nz)[None, None, :] * np.ones((nx, ny, nz))
    n_q_t, q_vals = ps.extract_axial_nq(dn[None, ...], Lz, m_max=6)
    amps = np.abs(n_q_t[0])
    assert int(np.argmax(amps)) == m_true - 1          # mode m_true dominates
    assert amps[m_true - 1] > 50.0 * np.delete(amps, m_true - 1).max()
    assert q_vals[m_true - 1] == pytest.approx(2.0 * np.pi * m_true / Lz)


def test_3d_binned_peak_in_correct_q_shell():
    """A z-only density modulation cos(q0·z)·cos(ω_p·t) → the |q| shell at q0
    peaks at ω_p and dominates a no-modulation shell. Analytic, on-bin."""
    N, dx = 8, 1.0
    L = N * dx
    m = 2
    q0 = 2.0 * np.pi * m / L                       # on a grid mode (axis z)
    nt, dt = 64, 0.5
    omega_p = _on_bin_omega(8, nt, dt)
    t = np.arange(nt) * dt
    z = np.arange(N) * dx
    spatial = np.broadcast_to(np.cos(q0 * z)[None, None, :], (N, N, N))
    dn = np.cos(omega_p * t)[:, None, None, None] * spatial[None, ...]

    spec = ps.spectrum_3d_binned(dn, dt, dx, n_bins=8, window=False)
    assert spec.mode == "3d_binned"
    b = int(np.argmin(np.abs(spec.q - q0)))        # shell nearest q0
    assert spec.peak_omega[b] == pytest.approx(omega_p, rel=1e-6)
    far = int(np.argmin(np.abs(spec.q - 0.3 * q0)))  # low-q shell, no modulation
    assert spec.power[:, b].max() > 100.0 * (spec.power[:, far].max() + 1e-30)


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
