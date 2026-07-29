"""Known-case tests for ``inqview.pipeline.density_fourier.loss_locator``.

Pure-numpy, machine-independent. These pin the BUG-A/B fix the audit confirmed
(``docs/validation/loss-function-formula-validation.md``, user verdict
2026-06-25): the loss-function peak-locator uses the COMPLEX phasor ``n_q(t)``
(not ``.real``) and returns ``|n_q(ω)|²/q²`` (not the bare ``|n_q|``).

Expected values are analytic truth fixed UP FRONT, independent of the code:

* A forward density wave ``n_q(t) = e^{+i ω0 t}`` puts all its spectral weight at
  ``+ω0`` (numpy's positive-frequency half) — this is the half ``loss_locator``
  inspects, and the half where the real E15 ``n_q(t)`` peaks (Stage C → 3.53 eV).
  The peak sits AT ω0 with the full amplitude.
* Taking ``.real`` makes ``cos(ω0 t)`` which splits weight onto ±ω0; the
  positive-half amplitude is HALVED ⇒ ``|·|²`` is quartered (ratio 0.25). This is
  the BUG-A folding the audit measured as a 0.500 amplitude ratio.
* For equal-amplitude modes, ``|n_q|²/q²`` scales as ``1/q²`` (BUG-B intent):
  doubling q quarters the locator value.

Run:
    cd /local/data/public/skcb2/tddft
    venv/bin/python3 -m pytest \
        inq-stack/tests/python/inqview/pipeline/test_density_fourier_loss.py -v
"""
from __future__ import annotations

import numpy as np
import pytest

from inqview.pipeline.density_fourier import loss_locator

# Pure numpy kernel (no matplotlib/VTK pulled), so the 'analysis'-tier mark fits.
pytestmark = pytest.mark.analysis


def _grid(n=400, dt=0.5):
    t = np.arange(n) * dt
    return t, dt, n


def test_complex_phasor_peaks_at_omega0_not_folded():
    """A forward phasor e^{+iω0 t} peaks at ω0 in the positive-frequency half."""
    t, dt, n = _grid()
    omega0 = 0.30                       # a.u. (well inside the band)
    q = 0.5
    sig = np.exp(+1j * omega0 * t)      # unit-amplitude forward density wave
    n_pad = n
    spec = loss_locator(sig, np.ones(n), n_pad, q)
    freq = np.fft.rfftfreq(n_pad, d=dt)
    omega = 2.0 * np.pi * freq
    k = int(np.argmax(spec))
    bin_w = omega[1] - omega[0]
    assert omega[k] == pytest.approx(omega0, abs=0.5 * bin_w)


def test_real_part_folding_quarters_the_peak():
    """BUG-A: feeding .real (cos) instead of the complex phasor halves the
    amplitude ⇒ quarters the |·|² locator at the peak (ratio 0.25)."""
    t, dt, n = _grid()
    omega0 = 0.30
    q = 0.5
    n_pad = n
    full = loss_locator(np.exp(+1j * omega0 * t), np.ones(n), n_pad, q)
    folded = loss_locator(np.cos(omega0 * t).astype(complex), np.ones(n), n_pad, q)
    pk_full = float(full.max())
    pk_folded = float(folded.max())
    assert pk_folded / pk_full == pytest.approx(0.25, rel=0.05)


def test_loss_scales_as_one_over_q_squared():
    """BUG-B: |n_q|²/q² scales as 1/q² for equal-amplitude modes."""
    t, dt, n = _grid()
    omega0 = 0.30
    sig = np.exp(+1j * omega0 * t)
    n_pad = n
    s1 = loss_locator(sig, np.ones(n), n_pad, q=0.5).max()
    s2 = loss_locator(sig, np.ones(n), n_pad, q=1.0).max()   # q doubled
    assert s2 / s1 == pytest.approx(0.25, rel=1e-6)           # (q1/q2)^2 = 1/4


if __name__ == "__main__":
    import subprocess
    import sys
    sys.exit(subprocess.call(["pytest", "-v", __file__]))
