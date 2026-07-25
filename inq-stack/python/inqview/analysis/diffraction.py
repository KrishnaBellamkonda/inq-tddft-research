"""diffraction.py — kinematic LEED diffraction from a real-space screen density.

A ``PlaneScreen`` (inqkit) records the *time-integrated real-space* electron
density crossing a plane at fixed ``z``:  ρ̄(x, y) = ∑_t ρ(x, y, z, t)·dt.
Loaded by :func:`inqview.io.leed.load_leed_pattern`.

The **LEED / diffraction pattern** is the kinematic (single-scattering) far-field
intensity, i.e. the squared magnitude of the 2D Fourier transform of that
real-space density:

    I(k_x, k_y) = | FFT2[ ρ̄(x, y) · w(x, y) ] |²

with an optional Hann window ``w`` to suppress the periodic-cell edge ringing
that would otherwise smear the diffraction orders.  For a periodic crystal the
peaks sit on the 2D reciprocal lattice (graphene → hexagonal); the forward
(+z) screens give the *transmission* pattern, the backward (−z) screens the
*reflection* pattern.

Deps-clean: numpy only (lives in the ``analysis`` layer; no VTK/matplotlib).

Reference for kinematic LEED intensity ∝ |FT(scatterer density)|²: standard
kinematic (single-scattering) diffraction theory, e.g. Van Hove, Weinberg &
Chan, *Low-Energy Electron Diffraction* (Springer, 1986), Ch. 2.  This is a
qualitative replica diagnostic (single-scattering, time-integrated density),
NOT a full dynamical LEED calculation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["Diffraction", "diffraction_pattern", "hann2d"]


def hann2d(ny: int, nx: int) -> np.ndarray:
    """Separable 2D Hann window of shape (ny, nx). Reduces FFT edge ringing."""
    wy = np.hanning(ny) if ny > 1 else np.ones(1)
    wx = np.hanning(nx) if nx > 1 else np.ones(1)
    return np.outer(wy, wx)


@dataclass(frozen=True)
class Diffraction:
    """Kinematic diffraction intensity in momentum space.

    intensity : (ny, nx) float — |FFT2(density·window)|², fftshifted so the
                zero-order (specular) peak is at the array centre.
    kx, ky    : 1-D momentum axes [rad / bohr], fftshifted to match.
    """

    intensity: np.ndarray
    kx: np.ndarray
    ky: np.ndarray

    @property
    def extent(self) -> tuple[float, float, float, float]:
        """(kx_min, kx_max, ky_min, ky_max) for imshow."""
        return (float(self.kx[0]), float(self.kx[-1]),
                float(self.ky[0]), float(self.ky[-1]))


def diffraction_pattern(density: np.ndarray, dx: float, dy: float,
                        *, hann: bool = True,
                        subtract_mean: bool = True) -> Diffraction:
    """Kinematic LEED pattern of a real-space screen density.

    Parameters
    ----------
    density : (ny, nx) real array — time-integrated density ρ̄(x, y).
    dx, dy  : grid spacing [bohr] along x (columns) and y (rows).
    hann    : apply a 2D Hann window before the FFT (default True) to suppress
              periodic-edge ringing that smears the diffraction orders.
    subtract_mean : remove the DC (mean) so the huge zero-order peak does not
              dominate the colour scale (default True).  The specular peak then
              reads ~0; set False to keep the absolute zero order.

    Returns
    -------
    Diffraction with fftshifted ``intensity`` and momentum axes (rad/bohr).
    """
    rho = np.asarray(density, dtype=np.float64)
    if rho.ndim != 2:
        raise ValueError(f"density must be 2-D (ny, nx); got shape {rho.shape}")
    ny, nx = rho.shape
    if subtract_mean:
        rho = rho - rho.mean()
    if hann:
        rho = rho * hann2d(ny, nx)
    amp = np.fft.fft2(rho)
    inten = np.abs(np.fft.fftshift(amp)) ** 2
    # momentum axes: fftfreq gives cycles/bohr; ×2π → rad/bohr
    kx = np.fft.fftshift(np.fft.fftfreq(nx, d=dx)) * 2.0 * np.pi
    ky = np.fft.fftshift(np.fft.fftfreq(ny, d=dy)) * 2.0 * np.pi
    return Diffraction(intensity=inten, kx=kx, ky=ky)
