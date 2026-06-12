"""Corrected Lindhard RPA dielectric / loss function and Gaussian-projectile
linear-response stopping power.

Supersedes ``inqview.pipeline.lindhard`` for any f-sum-critical use. That older
module computes ``u = (omega + i*eta)/(q vF)`` but passes only ``u.real`` to the
imaginary part, so Im chi0 is exactly zero outside the e-h continuum and the
**undamped plasmon pole is missing** — its f-sum rule fails at small q
(ratio 0.01-0.13 for q < 0.5, exact only where the plasmon is Landau-damped).

This module evaluates the Lindhard function with the **full complex argument**
(complex log), so the plasmon appears as a finite-eta Lorentzian and the f-sum
rule int_0^inf omega*Im[-1/eps] domega = (pi/2) omega_p^2 holds at all q
(verified to <1e-3, 2026-06-12).

Conventions (Hartree a.u., m_e = 1):
  z = q/(2 kF),  u = (omega + i eta)/(q vF),  vF = kF
  chi0(q,omega) = -(kF/pi^2) * [ 1/2 + (1/8z)( (1-(z-u)^2) ln((z-u+1)/(z-u-1))
                                             + (1-(z+u)^2) ln((z+u+1)/(z+u-1)) ) ]
  eps = 1 - (4 pi / q^2) chi0      (Giuliani-Vignale chi0 < 0 convention)
  ELF L(q,omega) = Im[-1/eps]

Stopping power for a Gaussian-charge (width sigma) unit projectile (Option B):
  S_LR(v; sigma) = (2/(pi v^2)) int_0^inf (dq/q) exp(-q^2 sigma^2)
                                  int_0^{q v} domega  omega  L(q, omega)
The form factor exp(-q^2 sigma^2) = |V_sigma(q)/V_point(q)|^2 with
V_sigma(q) ∝ exp(-q^2 sigma^2/2) (the erf-smoothed Gaussian charge).

References (to be cited in the report; notes under docs/sources/):
  Lindhard 1954 (RPA dielectric); Lindhard & Winther 1964 (stopping integral);
  Echenique-Nieminen-Ritchie 1981 / Echenique et al. 1986 (nonlinear low-v
  benchmark, Z=-1 Barkas); Correa 2018 (rt-TDDFT stopping practice).
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "kF_from_rs", "density_from_kF", "omega_p", "k_TF",
    "chi0", "epsilon_rpa", "loss_function", "plasmon_dispersion",
    "stopping_power_sigma",
]

_PI = np.pi


def kF_from_rs(rs: float) -> float:
    return (9.0 * _PI / 4.0) ** (1.0 / 3.0) / rs


def density_from_kF(kF: float) -> float:
    return kF ** 3 / (3.0 * _PI ** 2)


def omega_p(kF: float) -> float:
    return float(np.sqrt(4.0 * _PI * density_from_kF(kF)))


def k_TF(kF: float) -> float:
    """Thomas-Fermi screening wavevector: k_TF^2 = 4 kF / pi (a.u.)."""
    return float(np.sqrt(4.0 * kF / _PI))


def chi0(q, omega, kF: float, *, eta: float = 1e-3) -> np.ndarray:
    """Complex Lindhard susceptibility with the FULL complex argument."""
    q = np.asarray(q, dtype=float)
    omega = np.asarray(omega, dtype=float)
    q, omega = np.broadcast_arrays(q, omega)
    vF = kF
    z = q / (2.0 * kF)
    z = np.where(z == 0.0, 1e-12, z)
    u = (omega + 1j * eta) / (q * vF)
    a = z - u
    b = z + u
    clog = lambda x: np.log((x + 1.0) / (x - 1.0))
    F = 0.5 + (1.0 / (8.0 * z)) * (
        (1.0 - a * a) * clog(a) + (1.0 - b * b) * clog(b)
    )
    return -(kF / _PI ** 2) * F


def epsilon_rpa(q, omega, kF: float, *, eta: float = 1e-3) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    return 1.0 - (4.0 * _PI / np.where(q == 0.0, np.inf, q * q)) * chi0(
        q, omega, kF, eta=eta
    )


def loss_function(q, omega, kF: float, *, eta: float = 1e-3) -> np.ndarray:
    """Energy-loss function L(q, omega) = Im[-1/eps_RPA]."""
    eps = epsilon_rpa(q, omega, kF, eta=eta)
    return -(1.0 / eps).imag


def plasmon_dispersion(q, kF: float) -> np.ndarray:
    """Bohm-Gross plasmon dispersion omega_pl(q) = sqrt(wp^2 + 3/5 vF^2 q^2)."""
    q = np.asarray(q, dtype=float)
    wp = omega_p(kF)
    return np.sqrt(wp ** 2 + 0.6 * kF ** 2 * q ** 2)


def stopping_power_sigma(
    v: float, kF: float, sigma: float, *,
    qmin: float = 1e-3, qmax: float | None = None,
    n_q: int = 400, n_omega: int | None = None, eta: float = 1e-2,
) -> float:
    """Gaussian-projectile linear-response stopping power S_LR(v; sigma).

    S = (2/(pi v^2)) int (dq/q) e^{-q^2 sigma^2} int_0^{qv} omega L domega.

    The inner omega-integral spans the sharp plasmon Lorentzian (width ~eta), so
    the omega resolution is tied to eta (dω ≈ eta/8) unless n_omega is forced.
    A moderate eta (1e-2) broadens the plasmon enough to integrate stably while
    leaving the f-sum and the Bragg-peak position intact.
    """
    if qmax is None:
        qmax = 2.0 * kF + 2.0 * v + 6.0 / max(sigma, 1e-3)  # form-factor support
    q_grid = np.linspace(qmin, qmax, n_q)
    inner = np.empty_like(q_grid)
    for i, q in enumerate(q_grid):
        wmax = q * v
        n_w = n_omega if n_omega is not None else max(400, int(wmax / (eta / 8.0)))
        w = np.linspace(1e-6, wmax, n_w)
        Lw = loss_function(np.full_like(w, q), w, kF, eta=eta)
        inner[i] = np.trapezoid(w * Lw, w)
    integrand = np.exp(-(q_grid ** 2) * sigma ** 2) * inner / q_grid
    return float((2.0 / (_PI * v ** 2)) * np.trapezoid(integrand, q_grid))
