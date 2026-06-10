"""Analytic tests for the windowed FFT kernel (``inqview.fourier``).

Pure-numpy, machine-independent (ADR 0005). Expected values are the closed-form
spectrum of a known tone, derived up front — NOT captured from code output.

Convention recovered from the signal-validation dossier
(``docs/validation/fft-normalization-validation.md``): for a real on-bin tone
``A cos(2π f0 t)`` the one-sided ``amplitude = |rfft|/n`` with interior bins ×2
returns the physical amplitude ``A`` for a **boxcar** window. A windowed case
is reduced by the window coherent gain ``mean(win)`` until the IV-E03 fix
(``/win.sum()``) lands — that single case is marked ``xfail`` below so the
baseline is green and the fix flips it red→green.

The kernel now lives at ``inqview.analysis.fourier`` (ADR 0003 restructure);
the assertions are migration-invariant — only the import line moved.
"""
from __future__ import annotations

import numpy as np
import pytest

from _signals import make_tone
from inqview.analysis.fourier import FourierTransform, WindowSpec

pytestmark = pytest.mark.analysis


def _peak(res):
    """(freq, amplitude) at the spectrum's maximum-amplitude bin."""
    k = int(np.argmax(res.amplitude))
    return float(res.frequency_au[k]), float(res.amplitude[k])


# ---------------------------------------------------------------------------
# Peak LOCATION — robust to window/normalization (the always-reliable assertion)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("window", ["boxcar", "hann", "hamming", "blackman"])
def test_peak_is_at_f0(tone, window):
    """The spectral peak sits at the tone frequency, for every window."""
    ft = FourierTransform(window=WindowSpec(window), detrend=False, zero_pad=1)
    res = ft.transform(tone.t, tone.y)
    f_peak, _ = _peak(res)
    # within half a bin (bin width = 1/(n*dt))
    bin_w = 1.0 / (tone.n * tone.dt)
    assert f_peak == pytest.approx(tone.f0, abs=0.5 * bin_w)


# ---------------------------------------------------------------------------
# Amplitude calibration — boxcar correct now; windowed pending IV-E03
# ---------------------------------------------------------------------------
def test_boxcar_amplitude_is_A(tone):
    """Boxcar one-sided amplitude recovers the physical amplitude A exactly."""
    ft = FourierTransform(window=WindowSpec("boxcar"), detrend=False, zero_pad=1)
    res = ft.transform(tone.t, tone.y)
    _, a_peak = _peak(res)
    assert a_peak == pytest.approx(tone.A, rel=1e-6)


@pytest.mark.parametrize("window", ["hann", "hamming", "blackman"])
def test_windowed_amplitude_is_A_after_coherent_gain_fix(tone, window):
    """IV-E03 FIXED: with /Σwin normalisation every window recovers A."""
    ft = FourierTransform(window=WindowSpec(window), detrend=False, zero_pad=1)
    res = ft.transform(tone.t, tone.y)
    _, a_peak = _peak(res)
    assert a_peak == pytest.approx(tone.A, rel=2e-2)


# ---------------------------------------------------------------------------
# Zero-pad only interpolates the axis — peak unmoved, amplitude preserved
# ---------------------------------------------------------------------------
def test_zero_pad_preserves_peak(tone):
    """zero_pad densifies the freq axis without changing peak f0 or amplitude."""
    base = FourierTransform(window=WindowSpec("boxcar"), detrend=False, zero_pad=1)
    padded = FourierTransform(window=WindowSpec("boxcar"), detrend=False, zero_pad=4)
    f1, a1 = _peak(base.transform(tone.t, tone.y))
    f4, a4 = _peak(padded.transform(tone.t, tone.y))
    bin_w = 1.0 / (tone.n * tone.dt)
    assert f4 == pytest.approx(f1, abs=0.5 * bin_w)
    assert a4 == pytest.approx(a1, rel=2e-2)


# ---------------------------------------------------------------------------
# DC / drift removal (todo.txt #2; pre-cursor to the IV-M12 subtract= API)
# ---------------------------------------------------------------------------
def test_dc_offset_hijacks_peak_without_detrend():
    """A constant offset with detrend=False drives the max to the ω≈0 bin.

    This is exactly the todo.txt #2 failure mode: un-subtracted spectra report a
    bogus DC 'peak'. The IV-M12 subtract= API (default canonicals) will fix it;
    here we pin the current behaviour so the fix is a visible red→green.
    """
    tn = make_tone(offset=50.0)
    ft = FourierTransform(window=WindowSpec("boxcar"), detrend=False, zero_pad=1)
    res = ft.transform(tn.t, tn.y)
    f_peak, _ = _peak(res)
    assert f_peak == pytest.approx(0.0, abs=1e-9)   # hijacked to DC


def test_detrend_removes_dc_and_recovers_f0():
    """detrend=True strips a constant+linear baseline; the true peak returns."""
    tn = make_tone(offset=50.0, drift=0.3)
    ft = FourierTransform(window=WindowSpec("boxcar"), detrend=True, zero_pad=1)
    res = ft.transform(tn.t, tn.y)
    f_peak, _ = _peak(res)
    bin_w = 1.0 / (tn.n * tn.dt)
    assert f_peak == pytest.approx(tn.f0, abs=0.5 * bin_w)


# ---------------------------------------------------------------------------
# subtract= API (IV-M12): per-baseline DC removal, peak always preserved
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("mode", ["initial", "mean", "detrend"])
def test_subtract_modes_remove_dc_and_keep_peak(mode):
    """initial/mean/detrend all strip a constant offset; peak stays at f0.

    Uses a SINE tone (zero at t=0) so initial-value subtraction is well-posed —
    it is the induced-response baseline (Δs(0)=0), and a cos tone (cos 0 = 1)
    would leave a spurious −A DC under 'initial'.
    """
    tn = make_tone(offset=50.0)
    y = 50.0 + tn.A * np.sin(2.0 * np.pi * tn.f0 * tn.t)     # induced-style signal
    ft = FourierTransform(window=WindowSpec("boxcar"), zero_pad=1, subtract=mode)
    res = ft.transform(tn.t, y)
    assert res.subtract == mode
    f_peak, a_peak = _peak(res)
    bin_w = 1.0 / (tn.n * tn.dt)
    assert f_peak == pytest.approx(tn.f0, abs=0.5 * bin_w)   # not hijacked to DC
    assert a_peak == pytest.approx(tn.A, rel=3e-2)           # amplitude intact


def test_only_detrend_suppresses_a_genuine_linear_drift():
    """Against a real drift, only detrend kills the low-frequency LEAKAGE.

    mean/detrend both zero the DC *bin* by construction; the distinguishing
    quantity is leakage into the low non-zero bins (1..10, well below the tone
    at bin 20), where the un-removed ramp shows up.
    """
    tn = make_tone(offset=0.0, drift=0.5)
    low = {}
    for mode in ("initial", "mean", "detrend"):
        ft = FourierTransform(window=WindowSpec("boxcar"), zero_pad=1, subtract=mode)
        res = ft.transform(tn.t, tn.y)
        low[mode] = float(res.amplitude[1:11].sum())        # low-band leakage
    assert low["detrend"] < low["initial"]
    assert low["detrend"] < low["mean"]


def test_subtract_none_lets_offset_hijack_peak():
    tn = make_tone(offset=50.0)
    ft = FourierTransform(window=WindowSpec("boxcar"), zero_pad=1, subtract="none")
    f_peak, _ = _peak(ft.transform(tn.t, tn.y))
    assert f_peak == pytest.approx(0.0, abs=1e-9)


def test_invalid_subtract_raises():
    with pytest.raises(ValueError):
        FourierTransform(subtract="bogus")


if __name__ == "__main__":
    import subprocess
    import sys

    sys.exit(subprocess.call(["pytest", "-v", __file__]))
