"""Inverse-FFT helpers for LEED screen patterns.

LEED `.dat` files store the time-integrated electron density `ρ(x,y; z_screen)`
sampled on a regular xy grid. Phase-2 fixed the FFT-shift convention on read
(`inqview.load_leed_pattern` applies `np.fft.fftshift`), so by the time the
data reaches this module it is centred:

    pattern.data.shape  = (ny, nx)
    pattern.extent_bohr = (-Lx/2, +Lx/2, -Ly/2, +Ly/2)
    pattern.dx_bohr, pattern.dy_bohr  — physical grid spacing

Two reconstruction methods are exposed:

* **Patterson autocorrelation** (`method="patterson"`):
  LEED is intensity-only, so we treat `pattern.data` as `|F(k)|²` where
  `F(k)` is the spatial Fourier transform of the projected density at
  `z = z_screen`. By the Wiener–Khinchin theorem,
  ``IFFT(|F|²) = ρ ⋆ ρ`` — the autocorrelation of the projected density.
  Pair-correlation peaks land at every interatomic separation in the
  molecule (in coronene: 1.4 / 2.4 / 2.8 Bohr for nearest- /
  second- / third-neighbour C–C). This is the physically defensible
  reconstruction when phases are unknown.

* **Phase-less amplitude** (`method="amp_only"`):
  Take ``A(k) = sqrt(|F(k)|²)`` and assign zero phase, then return
  `|IFFT(A)|²`. A coarse heuristic — phase information is lost in
  intensity-only LEED, so this is not a reliable reconstruction of the
  true projected density, but it is sometimes informative for
  visualising rough envelope.

Both methods optionally apply a 2D Hann window before the inverse FFT
to suppress ringing from the periodic boundary.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    from .. import LeedPattern


def _hann_window_2d(ny: int, nx: int) -> np.ndarray:
    wy = np.hanning(ny)
    wx = np.hanning(nx)
    return wy[:, None] * wx[None, :]


def reconstruct_real_space(
    pattern: "LeedPattern",
    *,
    method: Literal["patterson", "amp_only"] = "patterson",
    hann: bool = True,
) -> np.ndarray:
    """Inverse-FFT a LEED screen back to a real-space density estimate.

    Parameters
    ----------
    pattern : LeedPattern
        Already centred (load_leed_pattern applies np.fft.fftshift).
    method : "patterson" (default) or "amp_only".
    hann : if True, multiply the input by a 2D Hann window before
        IFFT to suppress edge ringing.

    Returns
    -------
    np.ndarray of shape (ny, nx), real, non-negative
        The reconstructed real-space density estimate. The physical extent
        matches `pattern.extent_bohr` (the IFFT preserves the array shape).
    """
    data = pattern.data.astype(np.float64, copy=True)
    if hann:
        data *= _hann_window_2d(*data.shape)

    if method == "patterson":
        # IFFT of intensity ⇒ autocorrelation. Inverse-shift so the centred
        # input matches numpy's FFT-natural convention; then ifft2; then
        # shift back so the autocorrelation is centred about the array.
        spec = np.fft.ifftshift(data)
        auto = np.fft.fftshift(np.real(np.fft.ifft2(spec)))
        # Real autocorrelation peaks at zero-shift (the array centre); return
        # the absolute value so the cosmetic floor doesn't dip below zero
        # due to numerical residuals.
        return np.abs(auto)
    elif method == "amp_only":
        amp = np.sqrt(np.maximum(data, 0.0))
        spec = np.fft.ifftshift(amp)
        rec = np.fft.fftshift(np.fft.ifft2(spec))
        return np.abs(rec) ** 2
    else:
        raise ValueError(f"unknown method {method!r}; use 'patterson' or 'amp_only'")
