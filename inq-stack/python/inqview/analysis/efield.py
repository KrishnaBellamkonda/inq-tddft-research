"""Electric field from a density grid via an FFT Poisson solve (inqview.analysis.efield).

Computes E(r) = −∇φ, with φ the electrostatic potential of the charge density
ρ = −n(r), under PERIODIC boundary conditions, by an FFT Poisson solve. NATIVE
ATOMIC UNITS (Hartree a.u.: e = mₑ = ℏ = 1, Coulomb constant 1/(4πε₀) = 1, so
Poisson reads  ∇²φ = −4π ρ).

Setting the G=0 Fourier component of φ to zero removes the (otherwise divergent)
mean — i.e. it imposes a uniform **neutralizing background**. For a jellium bath
that background *is* the physical positive jellium, so the returned field is the
field of the density FLUCTUATION δn = n − ⟨n⟩ (a uniform charge sources no field
under periodic BC anyway). In atomic units this E is directly related to INQ's
Hartree potential v_H (built from +n): E = +∇v_H — a free cross-check.

Method (periodic, uniform grid):
  1. number → charge density   ρ(r) = −n(r)            (electron charge = −1 a.u.)
  2. forward FFT               ρ(r) → ρ̃(G)
  3. Poisson in k-space        φ̃(G) = 4π ρ̃(G) / |G|² ,   φ̃(0) := 0
  4. field in k-space          Ẽ(G) = −i G φ̃(G)
  5. inverse FFT               Ẽ(G) → E(r)   (real part; imag is FFT round-off)

Density layout follows the inqkit x-slowest / z-fastest convention: shape
(nx, ny, nz). ``spacing`` is a scalar dx (Bohr) or a 3-tuple (dx, dy, dz).

Units: ``units="atomic"`` (default) returns E in Hartree atomic units
(E_h / (e a₀)); ``units="SI"`` multiplies by the a.u. field constant to give V/m.

Pure numpy (deps-clean: numpy only).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# CODATA atomic unit of electric field, in V/m (E_h / (e a₀)).
_AU_FIELD_TO_SI_VM = 5.14220674763e11


@dataclass(frozen=True)
class EField:
    """Cartesian electric-field components on the density grid.

    ``ex/ey/ez`` each have shape (nx, ny, nz). ``units`` is "atomic"
    (E_h/(e a₀)) or "SI" (V/m).
    """

    ex: np.ndarray
    ey: np.ndarray
    ez: np.ndarray
    units: str

    @property
    def magnitude(self) -> np.ndarray:
        """|E|(r), same shape and units as the components."""
        return np.sqrt(self.ex ** 2 + self.ey ** 2 + self.ez ** 2)


def _spacing3(spacing) -> tuple[float, float, float]:
    if np.isscalar(spacing):
        return (float(spacing),) * 3
    s = tuple(float(v) for v in spacing)
    if len(s) != 3:
        raise ValueError("spacing must be a scalar or length-3 (dx, dy, dz)")
    return s


def electric_field(n_grid, spacing, *, units: str = "atomic") -> EField:
    """E(r) of the charge density ρ = −n via a periodic FFT Poisson solve.

    Parameters
    ----------
    n_grid : array_like, shape (nx, ny, nz)
        Electron NUMBER density in e / Bohr³ (≥ 0). x-slowest, z-fastest.
    spacing : float or (dx, dy, dz)
        Grid spacing in Bohr.
    units : {"atomic", "SI"}
        Output units. "atomic" → E_h/(e a₀); "SI" → V/m.

    Returns
    -------
    EField
        Frozen dataclass with ex, ey, ez (each (nx, ny, nz)) and a ``magnitude``.

    Notes
    -----
    The G=0 component of φ is set to zero, so the result is the field of
    δn = n − ⟨n⟩ (the neutralizing-background / jellium convention).
    """
    if units not in ("atomic", "SI"):
        raise ValueError('units must be "atomic" or "SI"')
    n = np.asarray(n_grid, dtype=float)
    if n.ndim != 3:
        raise ValueError("n_grid must be a 3D array (nx, ny, nz)")
    dx, dy, dz = _spacing3(spacing)
    nx, ny, nz = n.shape

    # 1. number density → charge density (electron charge = −1 a.u.)
    rho = -n

    # 2. forward FFT of the charge density
    rho_k = np.fft.fftn(rho)

    # reciprocal-grid wavevectors G = 2π · fftfreq  (radians / Bohr)
    gx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    gy = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
    gz = 2.0 * np.pi * np.fft.fftfreq(nz, d=dz)
    GX, GY, GZ = np.meshgrid(gx, gy, gz, indexing="ij")
    G2 = GX ** 2 + GY ** 2 + GZ ** 2
    G2[0, 0, 0] = 1.0  # avoid div-by-zero; φ̃(0) zeroed on the next line

    # 3. Poisson in Fourier space (a.u.):  ∇²φ = −4π ρ  →  φ̃ = 4π ρ̃ / |G|²
    phi_k = 4.0 * np.pi * rho_k / G2
    phi_k[0, 0, 0] = 0.0  # drop the mean → uniform neutralizing background

    # 4. field in Fourier space:  Ẽ(G) = −i G φ̃(G)
    ex_k = -1j * GX * phi_k
    ey_k = -1j * GY * phi_k
    ez_k = -1j * GZ * phi_k

    # 5. inverse FFT; ρ real ⇒ E real (imag part is FFT round-off, discarded)
    ex = np.fft.ifftn(ex_k).real
    ey = np.fft.ifftn(ey_k).real
    ez = np.fft.ifftn(ez_k).real

    if units == "SI":
        ex = ex * _AU_FIELD_TO_SI_VM
        ey = ey * _AU_FIELD_TO_SI_VM
        ez = ez * _AU_FIELD_TO_SI_VM

    return EField(ex=ex, ey=ey, ez=ez, units=units)
