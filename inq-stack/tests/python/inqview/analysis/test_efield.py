"""Tests for the FFT-Poisson electric-field kernel (``inqview.analysis.efield``).

Pure numpy on synthetic densities. Three locked known-cases:
  1. uniform density  → E ≡ 0 (G=0 removed = neutralizing background);
  2. single cosine    → E_z = −(4πA/G₀) sin(G₀z) analytic (machine precision);
  3. Gaussian charge  → isolated-Gaussian erf field (physical, loose tol).
All in native atomic units; a fourth test checks the SI scaling is one constant.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.special import erf

import inqview.analysis.efield as efield

pytestmark = pytest.mark.analysis


def test_uniform_density_gives_zero_field():
    # A uniform charge sources no field under periodic BC (G=0 removed).
    n = np.full((16, 16, 16), 0.7)
    E = efield.electric_field(n, 0.5)
    assert E.magnitude.max() < 1e-12


def test_cosine_density_matches_analytic_field():
    # n(z) = n0 + A cos(G0 z), varying along axis 2 (z), commensurate G0 = 2π/Lz.
    # δρ = −A cos(G0 z) ⇒ φ = −(4πA/G0²) cos(G0 z) ⇒ E_z = −(4πA/G0) sin(G0 z).
    nx = ny = 8
    nz = 32
    dz = 0.5
    Lz = nz * dz
    n0, A = 1.0, 0.1
    G0 = 2.0 * np.pi / Lz  # m = 1, exactly representable on the grid

    z = np.arange(nz) * dz
    n1d = n0 + A * np.cos(G0 * z)
    n = np.broadcast_to(n1d, (nx, ny, nz)).copy()  # vary along z only

    E = efield.electric_field(n, (0.5, 0.5, dz))

    ez_analytic = -(4.0 * np.pi * A / G0) * np.sin(G0 * z)
    ez_expected = np.broadcast_to(ez_analytic, (nx, ny, nz))

    assert np.allclose(E.ez, ez_expected, atol=1e-10)
    assert np.allclose(E.ex, 0.0, atol=1e-10)
    assert np.allclose(E.ey, 0.0, atol=1e-10)


def test_gaussian_charge_matches_isolated_erf_field():
    # A unit-integral Gaussian electron density; away from the centre and well
    # inside the box, the radial field magnitude matches the isolated Gaussian
    # |E|(r) = f(r)/r²,  f(r) = erf(r/(σ√2)) − √(2/π)(r/σ) exp(−r²/2σ²).
    N = 48
    dx = 0.5
    sigma = 1.5
    L = N * dx
    c = L / 2.0  # centre on a node (N even ⇒ centre at index N/2)

    ax = np.arange(N) * dx
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    r2 = (X - c) ** 2 + (Y - c) ** 2 + (Z - c) ** 2
    n = np.exp(-r2 / (2.0 * sigma ** 2))
    n /= n.sum() * dx ** 3  # normalise ∫n dV = 1 electron

    E = efield.electric_field(n, dx)
    mag = E.magnitude

    ic = N // 2
    for r in (3.0, 4.0):  # 2σ–~2.7σ, and r < L/4 = 6 (small periodic-image error)
        k = ic + int(round(r / dx))
        # sample along +x from the centre
        num = mag[k, ic, ic]
        f = erf(r / (sigma * np.sqrt(2.0))) - np.sqrt(2.0 / np.pi) * (r / sigma) * np.exp(
            -r ** 2 / (2.0 * sigma ** 2)
        )
        analytic = f / r ** 2
        assert num == pytest.approx(analytic, rel=0.05)


def test_si_units_are_a_constant_rescale():
    rng = np.random.default_rng(0)
    n = 1.0 + 0.05 * rng.standard_normal((12, 12, 12))
    Ea = efield.electric_field(n, 0.4, units="atomic")
    Es = efield.electric_field(n, 0.4, units="SI")
    assert Es.units == "SI"
    assert np.allclose(Es.ez, Ea.ez * 5.14220674763e11, rtol=1e-12)
