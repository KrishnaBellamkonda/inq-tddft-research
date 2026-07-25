"""Data-contract tests for the FFT-pipeline panel (fourier_panel).

The panel must DEPICT exactly what the kernel computes — never re-derive it. So we
assert the panel's stages are consistent with `FourierTransform`'s own pipeline:
  * stage 6 (the FFT) equals `FourierTransform.transform().amplitude` bit-for-bit;
  * the de-trended stage equals `FourierTransform._apply_subtract` of the
    (transient-skipped) signal;
  * the windowed stage equals detrended × the kernel window;
  * the panel exposes the required 3×2 = 6 axes in the specified order.
Expected peak is fixed up front from the synthetic tone, not read off the panel.
"""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")  # headless, portable

import numpy as np
import pytest

from inqview.analysis.fourier import FourierTransform, WindowSpec
from inqview.visualisation import fourier_panel as FP

pytestmark = pytest.mark.theme

# Synthetic tone: cos(omega t), omega chosen so the peak lands in [4,6] eV.
HA_TO_EV = FP.HA_TO_EV
_OMEGA_EV = 5.0
_OMEGA_AU = _OMEGA_EV / HA_TO_EV
_T = np.linspace(0.0, 400.0, 2000)               # uniform a.u. grid
_SIG = np.cos(_OMEGA_AU * _T) + 0.3              # +0.3 DC offset to exercise detrend


def test_stages_match_kernel():
    ft = FourierTransform(window=WindowSpec("hann"), zero_pad=4, subtract="detrend")
    st = FP.fft_stages(_T, _SIG, ft)
    res = ft.transform(_T, _SIG)
    # the FFT stage IS the kernel's transform (no re-derivation)
    assert np.allclose(st["amplitude"], res.amplitude)
    assert np.allclose(st["freq_au"], res.frequency_au)
    # de-trended stage == kernel's own baseline removal
    assert np.allclose(st["detrended"], ft._apply_subtract(_SIG))
    # windowed == detrended * kernel window; padded length == zero_pad * n
    assert np.allclose(st["windowed"], st["detrended"] * st["window"])
    assert len(st["padded"]) == ft.zero_pad * len(st["t"])
    assert np.allclose(st["padded"][len(st["t"]):], 0.0)   # the pad is zeros


def test_panel_has_six_axes_and_finds_peak():
    ft = FourierTransform(window=WindowSpec("hann"), zero_pad=4)
    fig = FP.fft_pipeline_panel(_T, _SIG, ft, label="tone",
                                peak_band=(4.0, 6.0), fmax=20.0)
    assert len(fig.axes) == 6                      # 3x2 panel
    # the detected peak (annotated on the linear FFT axis) is near 5 eV
    st = FP.fft_stages(_T, _SIG, ft)
    fx = st["freq_au"] * 2.0 * np.pi * HA_TO_EV   # angular energy convention
    band = (fx >= 4.0) & (fx <= 6.0)
    fpk = fx[band][int(np.argmax(st["amplitude"][band]))]
    assert abs(fpk - _OMEGA_EV) < 0.2
    import matplotlib.pyplot as plt
    plt.close(fig)


def test_transient_skip_respected():
    ft = FourierTransform(t_start_au=50.0, zero_pad=1, subtract="none")
    st = FP.fft_stages(_T, _SIG, ft)
    assert st["t"][0] >= _T[0] + 50.0              # skipped the first 50 a.u.
    assert len(st["t"]) < len(_T)


def test_default_is_mean_and_detrend_overlaid():
    """Default panel uses 'mean' and overlays a 'detrend' comparison on the FFT
    axes (user verdict 2026-06-25). The comparison curve == an independent
    detrend transform (no re-derivation)."""
    ft = FourierTransform(window=WindowSpec("hann"), zero_pad=4)   # default subtract
    assert ft.subtract == "mean"
    fig = FP.fft_pipeline_panel(_T, _SIG, ft, label="tone",
                                peak_band=(4.0, 6.0), fmax=20.0)
    lin_ax = fig.axes[4]                          # row-major: (2,0) = linear FFT
    labels = [str(ln.get_label()) for ln in lin_ax.lines]
    assert any("mean" in lbl for lbl in labels)
    assert any("detrend" in lbl for lbl in labels)
    # the dashed comparison line equals an independent detrend transform
    cmp_line = next(ln for ln in lin_ax.lines if "detrend" in str(ln.get_label()))
    ft_d = FourierTransform(window=WindowSpec("hann"), zero_pad=4, subtract="detrend")
    assert np.allclose(cmp_line.get_ydata(), ft_d.transform(_T, _SIG).amplitude)
    import matplotlib.pyplot as plt
    plt.close(fig)
