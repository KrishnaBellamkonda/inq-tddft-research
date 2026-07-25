"""Known-case tests for inqview.analysis.diffraction (kinematic LEED)."""
from __future__ import annotations

import numpy as np

from inqview.analysis.diffraction import diffraction_pattern, hann2d


def test_cosine_grating_peaks_at_expected_k():
    """A real-space cosine grating of period λ must diffract to ±k=±2π/λ.

    Build cos(2π x / λ) over an INTEGER number of periods so the FFT is exact
    (no leakage). Disable the Hann window so the peak stays a single bin.
    """
    nx, ny = 64, 48
    dx = dy = 0.5                      # bohr
    n_periods = 8                      # integer => exact FFT bins
    lam = nx * dx / n_periods          # period in bohr
    x = np.arange(nx) * dx
    row = np.cos(2.0 * np.pi * x / lam)
    rho = np.tile(row, (ny, 1))        # constant along y => peaks on the kx axis

    d = diffraction_pattern(rho, dx, dy, hann=False, subtract_mean=True)

    # zero out the central (ky) column region is unnecessary: peaks are on kx.
    iy0 = np.argmin(np.abs(d.ky))      # ky ~ 0 row (cosine is y-uniform)
    line = d.intensity[iy0]
    kpk = abs(d.kx[int(np.argmax(line))])
    assert np.isclose(kpk, 2.0 * np.pi / lam, rtol=1e-6), (kpk, 2 * np.pi / lam)


def test_real_input_gives_centrosymmetric_spectrum():
    """I(k) == I(-k) for a real density. Exact under [::-1,::-1] only for ODD
    sizes (even N leaves the Nyquist bin unpaired after fftshift)."""
    rng = np.random.default_rng(1)
    rho = rng.random((31, 27))                 # odd, odd
    d = diffraction_pattern(rho, 0.5, 0.5, hann=False, subtract_mean=False)
    assert np.allclose(d.intensity, d.intensity[::-1, ::-1], atol=1e-8)


def test_dc_removed_specular_is_small():
    """subtract_mean=True should make the zero-order (specular) bin ~0."""
    rng = np.random.default_rng(0)
    rho = rng.random((32, 40)) + 5.0          # large DC offset
    d = diffraction_pattern(rho, 0.3, 0.3, hann=False, subtract_mean=True)
    iy0 = np.argmin(np.abs(d.ky)); ix0 = np.argmin(np.abs(d.kx))
    peak = d.intensity.max()
    assert d.intensity[iy0, ix0] < 1e-6 * peak


def test_axes_units_rad_per_bohr():
    """kx spacing = 2π/(nx·dx); Nyquist magnitude = π/dx."""
    nx, dx = 50, 0.4
    d = diffraction_pattern(np.zeros((10, nx)), dx, dx, hann=False)
    dk = d.kx[1] - d.kx[0]
    assert np.isclose(dk, 2.0 * np.pi / (nx * dx), rtol=1e-12)
    assert np.max(np.abs(d.kx)) <= np.pi / dx + 1e-9


def test_hann2d_shape_and_endpoints():
    w = hann2d(8, 6)
    assert w.shape == (8, 6)
    assert np.isclose(w[0].max(), 0.0, atol=1e-12)   # Hann endpoints vanish
