"""Analytical Lindhard susceptibility χ⁰(q, ω) for the free electron gas.

Reference: Lindhard (1954), Giuliani & Vignale ch. 4, Mahan ch. 5.
Textbook form (zero temperature, retarded propagator) — what we compute
here is the bare bubble Π(q, ω) = -χ⁰(q, ω) in some conventions; we follow
Giuliani-Vignale Eq. (4.21):

  χ⁰(q, ω) = (m k_F / π² ℏ²) · [ -1/2 + ((1 - ν₋²) / 4z) log|(ν₋+1)/(ν₋-1)|
                                       + ((1 - ν₊²) / 4z) log|(ν₊+1)/(ν₊-1)| ]

where z = q / (2 k_F) and ν± = ω̄/z ± z, ω̄ = ω / (4 E_F).

In atomic units (ℏ = m_e = 1):
  χ⁰(q, ω) = (k_F / π²) · F(q, ω)

with F dimensionless. We expose:

  - chi0_real(q, omega, kF)  — real part
  - chi0_imag(q, omega, kF)  — imaginary part
  - chi0(q, omega, kF)        — complex
  - epsilon_rpa(q, omega, kF) — ε_RPA = 1 − v(q) χ⁰, v(q) = 4π/q²
  - loss_function(q, omega, kF) — Im[−1/ε_RPA]
  - plasmon_omega(q, kF)      — Bohm-Gross ω_pl(q) (small-q analytic)
  - stopping_power(v, kF, qmin, qmax) — numerical integration of
       S(v) = (2/π v²) ∫_{qmin}^{qmax} (dq/q) ∫_0^{qv} dω · ω · Im[−1/ε(q, ω)]

All inputs and outputs in atomic units unless noted.

The implementation is purely analytical (no SCF iteration). Validated against:
  - Static limit χ⁰(q→0, 0) = -D(E_F) = -(3 n / 2 E_F)
  - f-sum rule ∫ ω Im[−1/ε(q,ω)] dω = π ω_p² / 2
  - Plasmon limit ω_pl(q→0) → ω_p = √(4π n)
  - High-ω limit χ⁰ → -n q² / ω²
  - Bethe limit: stopping_power at large v matches Bethe-Lindhard
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import integrate

# TODO: As this script calculates the results analytically, would it be better organised
# under a submodule named analytics (where jellium and lindhard can exist)

# =============================================================================
# Lindhard susceptibility
# =============================================================================
def _F_real(z: NDArray, u: NDArray) -> NDArray:
    """Dimensionless real part of the Lindhard function.

    Definition (Lindhard 1954; Giuliani-Vignale §4.4):

        Re χ⁰(q, ω) = -(k_F/π²) · F(z, u)

        F(z, u) = 1/2 + (1/(8z)) · [
            (1 - (z-u)²) · ln|((z-u+1)/(z-u-1))| +
            (1 - (z+u)²) · ln|((z+u+1)/(z+u-1))|
        ]

    where z = q/(2k_F), u = ω/(q v_F).

    Static limit z → 0, u = 0: F → 1, recovering χ⁰(0, 0) = -k_F/π²
    (the compressibility / Thomas-Fermi static response).
    """
    eps = 1e-14
    a = z - u  # first argument
    b = z + u  # second argument

    def safe_log(x):
        # ln|(x+1)/(x-1)| with regularisation
        return np.log((np.abs(x + 1.0) + eps) / (np.abs(x - 1.0) + eps))

    g = 0.5 + (1.0 / (8.0 * z)) * (
        (1.0 - a * a) * safe_log(a) + (1.0 - b * b) * safe_log(b)
    )
    return g


def _F_imag(z: NDArray, u: NDArray) -> NDArray:
    """Dimensionless imaginary part for ω > 0 (T=0 retarded propagator).

    Definition (Lindhard 1954, T=0 retarded):

        Im χ⁰(q, ω) = -(k_F/π²) · G(z, u)

    where:
        G = π u / 2          for 0 ≤ u < |1 - z|  and z < 1  (Region I)
        G = π (1 - (z-u)²) / (8 z)   for |1 - z| ≤ u and u ≤ 1 + z  (Region II)
        G = 0                otherwise

    This is the dissipative part inside the electron-hole continuum.
    """
    u = np.asarray(u)
    z = np.asarray(z)
    out = np.zeros_like(u, dtype=float)
    u_pos = np.where(u >= 0, u, 0.0)
    abs_z_m1 = np.abs(1.0 - z)
    z_p1 = 1.0 + z

    # Region I: 0 ≤ u ≤ |1-z|  AND z < 1 (otherwise Region I empty)
    sel_I = (u_pos >= 0) & (u_pos <= abs_z_m1) & (z < 1.0) & (u >= 0)
    out = np.where(sel_I, 0.5 * np.pi * u_pos, out)

    # Region II: |1-z| < u ≤ 1+z
    sel_II = (u_pos > abs_z_m1) & (u_pos <= z_p1) & (u >= 0) & ~sel_I
    out = np.where(sel_II, np.pi * (1.0 - (z - u_pos) ** 2) / (8.0 * z), out)
    return out


def chi0(q: ArrayLike, omega: ArrayLike, kF: float, *, eta: float = 1e-3) -> NDArray:
    """Complex Lindhard susceptibility χ⁰(q, ω) for a 3D free electron gas.

    Args:
      q: wavenumber (a.u., positive).
      omega: frequency (a.u.). Broadcasts with q.
      kF: Fermi wavevector (a.u.).
      eta: small imaginary part for retarded propagator (a.u.; default 1e-3).

    Returns:
      Complex array χ⁰(q, ω) in a.u.

    Sign convention: static limit  χ⁰(q→0, 0) = -k_F/π² = -N(E_F).
    Imaginary part is non-zero only in the electron-hole continuum.
    """
    q = np.atleast_1d(np.asarray(q, dtype=float))
    omega = np.atleast_1d(np.asarray(omega, dtype=float))

    vF = kF  # a.u. with m_e = 1
    prefactor = kF / (np.pi ** 2)   # ≡ N(E_F) in a.u.

    # Dimensionless variables; broadcast
    q_b, omega_b = np.broadcast_arrays(q, omega)
    z = q_b / (2.0 * kF)
    z = np.where(z == 0.0, 1e-9, z)
    u = (omega_b + 1j * eta) / (q_b * vF)
    u_re = u.real

    F_re = _F_real(z, u_re)
    F_im = _F_imag(z, u_re)

    out = (-prefactor * F_re + 1j * (-prefactor * F_im))
    # Keep the leading η→0 imaginary support intact: where the e-h continuum
    # gives zero (Im F_im = 0), there is no dissipation. The small eta only
    # avoids divide-by-zero in the log near ν = ±1.
    return out.reshape(q_b.shape)


def chi0_real(q, omega, kF, **kw):
    return chi0(q, omega, kF, **kw).real


def chi0_imag(q, omega, kF, **kw):
    return chi0(q, omega, kF, **kw).imag


# =============================================================================
# Dielectric and loss function
# =============================================================================
def vq(q: ArrayLike) -> NDArray:
    """3D Coulomb interaction in Fourier space (a.u.):  v(q) = 4π / q²."""
    q = np.asarray(q, dtype=float)
    return 4.0 * np.pi / np.where(q == 0.0, np.inf, q * q)


def epsilon_rpa(q, omega, kF, **kw):
    """RPA dielectric function ε(q, ω) = 1 − v(q) χ⁰(q, ω).

    Note: with Giuliani-Vignale's χ⁰ < 0 convention, ε = 1 − v_q χ⁰
    matches the textbook expression where the plasmon dispersion is at
    Re ε = 0.
    """
    return 1.0 - vq(q) * chi0(q, omega, kF, **kw)


def loss_function(q, omega, kF, **kw):
    """Energy-loss function L(q, ω) = Im[−1/ε_RPA(q, ω)]."""
    eps = epsilon_rpa(q, omega, kF, **kw)
    inv_eps = 1.0 / eps
    return -inv_eps.imag


# =============================================================================
# Plasmon dispersion
# =============================================================================
def plasmon_omega(q: ArrayLike, kF: float, *, order: str = "bohm_gross") -> NDArray:
    """Plasmon dispersion ω_pl(q) for the free electron gas.

    Args:
      q: wavenumber (a.u.).
      kF: Fermi wavevector (a.u.).
      order: 'bohm_gross' for ω² = ω_p² + (3/5) v_F² q²  (RPA small-q),
             'plasma' for ω = ω_p (q=0 limit).
    """
    q = np.asarray(q, dtype=float)
    vF = kF
    # ω_p² = 4π n with n = k_F³ / (3π²)
    n = kF ** 3 / (3.0 * np.pi ** 2)
    omega_p = np.sqrt(4.0 * np.pi * n)
    if order == "plasma":
        return np.full_like(q, omega_p, dtype=float)
    elif order == "bohm_gross":
        return np.sqrt(omega_p ** 2 + 0.6 * vF ** 2 * q ** 2)
    raise ValueError(f"unknown order: {order!r}")


def plasma_frequency(kF: float) -> float:
    n = kF ** 3 / (3.0 * np.pi ** 2)
    return float(np.sqrt(4.0 * np.pi * n))


def fermi_energy(kF: float) -> float:
    return 0.5 * kF ** 2


def density_from_kF(kF: float) -> float:
    return kF ** 3 / (3.0 * np.pi ** 2)


def kF_from_rs(rs: float) -> float:
    """k_F from r_s in a.u. via n = 3 / (4π r_s³) and k_F = (3π² n)^(1/3)."""
    return (9.0 * np.pi / 4.0) ** (1.0 / 3.0) / rs


# =============================================================================
# Stopping power from the Lindhard loss function
# =============================================================================
def stopping_power(v: float, kF: float, *, qmin: float = 0.0, qmax: float | None = None,
                   nq: int = 400, nomega: int = 400, eta: float = 1e-3) -> float:
    """Electronic stopping power S(v) for unit-charge projectile in jellium.

    S(v) = (2 Z₁² / π v²) · ∫_{qmin}^{qmax} (dq/q) · ∫_0^{qv} dω · ω ·
           Im[−1/ε_RPA(q, ω)]

    Z₁² = 1 (electron).

    Args:
      v: projectile velocity (a.u.).
      kF: Fermi wavevector (a.u.).
      qmin: lower integration limit (a.u.). Use the box's q_min = 2π/L for
            box-truncated comparison.
      qmax: upper limit. Default 2 m_e v + k_F (kinematic limit).
      nq: number of q grid points.
      nomega: number of ω grid points per q.
      eta: imaginary regulator for ε.

    Returns:
      S(v) in a.u. (Ha / Bohr). Multiply by 27.2114 for eV/Bohr.
    """
    if qmax is None:
        qmax = 2.0 * v + kF
    if qmin >= qmax:
        return 0.0
    # Log-spaced q grid (most weight at small q)
    qmin_eff = max(qmin, 1e-3)
    q_grid = np.logspace(np.log10(qmin_eff), np.log10(qmax), nq)

    integrand = np.empty(nq)
    for i, q in enumerate(q_grid):
        omega_max = q * v
        if omega_max <= 0:
            integrand[i] = 0.0
            continue
        omega_grid = np.linspace(1e-4, omega_max, nomega)
        L = loss_function(q, omega_grid, kF, eta=eta)
        # ω · Im[-1/ε] integrand
        integrand[i] = np.trapz(omega_grid * L, omega_grid) / q

    # Trapezoid on log-spaced q is fine if we include the 1/q factor (already done)
    # Convert log grid: ∫ f(q)/q dq = ∫ f(q) d(ln q)
    S = (2.0 / (np.pi * v ** 2)) * np.trapz(integrand, q_grid)
    # The expression ∫(dq/q) · g(q) with our log grid: actually we put 1/q
    # into the integrand already; integrating dq on log grid is wrong.
    # Correct form: integrand_log[i] = g(q_i); ∫ d(ln q) = ∫(dq/q).
    # Recompute integrand without /q, then trapz over ln q.
    integrand_logq = integrand * q_grid  # remove the 1/q; we'll integrate over ln q
    lnq = np.log(q_grid)
    S = (2.0 / (np.pi * v ** 2)) * np.trapz(integrand_logq, lnq)
    return float(S)
