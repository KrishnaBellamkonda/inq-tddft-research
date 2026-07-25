"""Fourier-analysis kernels (skill-local, shippable, numpy-only).

The LOCKED project FFT pipeline for any TDDFT time-series observable
(user verdict 2026-06-25). Every step is a deliberate, audited choice — a
spectrum is never a black box. All functions take plain numpy arrays (loading is
the caller's job) so the module is portable; no scipy, no matplotlib.

THE PIPELINE (one uniform path for EVERY observable — dipole, current, energy):

    transient cut (t < t_start dropped)
      -> baseline removal: 'mean'  (canonical default; verdict 2026-06-25)
      -> window:           Hann
      -> zero-pad:         x4   (axis interpolation only — NO new information)
      -> rfft + COHERENT-GAIN normalisation ( / Sum(win), interior bins x2 )
      -> energy axis:      ANGULAR  hbar*omega = 2*pi*f * HA_TO_EV  (eV)

Baseline = 'mean' for ALL observables (the user's uniform override of the
per-observable initial/detrend split in fft-drift-removal-validation.md). The
figure standard overlays the 'detrend' result as a dashed COMPARISON so the
choice is auditable. 'mean' removes a constant offset but NOT a linear ramp; an
energy spectrum therefore carries a residual low-omega feature — peak POSITIONS
are unaffected (verified 3e-7), only the near-DC region moves.

Grounding: Harris 1978 (window coherent gain); Yabana-Bertsch (induced-response
IC Dn(0)=0, the 'initial' alternative); fft-normalization / fft-drift-removal /
loss-function-formula validation dossiers (docs/validation/).

LOSS FUNCTION: loss_locator computes |n_q(omega)|^2 / q^2 from the COMPLEX phasor
n_q(t). It is a plasmon-peak LOCATOR only — right pole positions and 1/q^2
q-trend — NOT a quantitatively faithful -Im[1/eps] (wrong line shape, area, and
absolute 4*pi normalisation off the undamped limit).
"""
from __future__ import annotations

import numpy as np

HA_TO_EV = 27.211386245988
_VALID_SUBTRACT = ("mean", "initial", "detrend", "none")


# ----------------------------------------------------------------------------
# Windows (numpy-only; the production kernel uses scipy.get_window 'periodic',
# the ~1e-2 difference is immaterial to peak position / coherent gain).
# ----------------------------------------------------------------------------
def make_window(name: str, n: int) -> np.ndarray:
    if name == "boxcar":
        return np.ones(n)
    if name == "hann":
        return np.hanning(n)
    if name == "hamming":
        return np.hamming(n)
    if name == "blackman":
        return np.blackman(n)
    raise ValueError(f"unknown window {name!r}; use boxcar|hann|hamming|blackman")


def subtract_baseline(values: np.ndarray, mode: str = "mean") -> np.ndarray:
    """Remove the chosen baseline BEFORE windowing.

    'mean'    - remove the DC level (canonical default, all observables);
    'initial' - s - s(0)   (Yabana-Bertsch induced-response IC; needs s(0) clean);
    'detrend' - remove a least-squares linear trend (the dashed comparison;
                the only mode that kills a genuine slope / conservation ramp);
    'none'    - leave it (a constant offset then hijacks the omega~0 bin).
    """
    if mode not in _VALID_SUBTRACT:
        raise ValueError(f"unknown subtract={mode!r}; valid: {_VALID_SUBTRACT}")
    v = np.asarray(values, float)
    if mode == "none":
        return v.copy()
    if mode == "mean":
        return v - float(np.mean(v))
    if mode == "initial":
        return v - v[0]
    # 'detrend' — least-squares linear fit removal (== scipy detrend type='linear')
    t = np.arange(len(v), dtype=float)
    A = np.vstack([t, np.ones_like(t)]).T
    slope, intercept = np.linalg.lstsq(A, v, rcond=None)[0]
    return v - (slope * t + intercept)


def fft_spectrum(time_au, values, *, window="hann", subtract="mean",
                 zero_pad=4, t_start_au=0.0) -> dict:
    """The LOCKED pipeline. Returns a dict with the energy axis + amplitude.

    Keys: energy_eV, omega_au, freq_au, amplitude (coherent-gain, one-sided),
    n_used, dt_au, subtract, window.
    """
    t = np.asarray(time_au, float)
    v = np.asarray(values, float)
    if t.ndim != 1 or v.ndim != 1 or len(t) != len(v):
        raise ValueError("time_au and values must be 1-D and equal length")

    if t_start_au > 0.0:
        mask = t >= (float(t[0]) + t_start_au)
        if mask.sum() < 4:
            raise ValueError(f"t_start_au={t_start_au} leaves <4 samples")
        t, v = t[mask], v[mask]

    n = len(t)
    dt = float(t[1] - t[0])
    sig = subtract_baseline(v, subtract)
    win = make_window(window, n)
    sig_win = sig * win

    zp = max(1, int(zero_pad))
    n_pad = n * zp
    if zp > 1:
        padded = np.zeros(n_pad, float)
        padded[:n] = sig_win
    else:
        padded = sig_win

    raw = np.fft.rfft(padded)
    freq = np.fft.rfftfreq(n_pad, d=dt)
    # COHERENT GAIN (Harris 1978): /Sum(win) so a unit tone returns ~1 for ANY
    # window; interior bins x2 for the one-sided spectrum. Zero-pad only
    # interpolates the axis, so Sum uses the length-n window (calibration is
    # independent of zero_pad).
    amp = np.abs(raw) / win.sum()
    amp[1:-1] *= 2.0

    omega_au = 2.0 * np.pi * freq                 # ANGULAR frequency
    return dict(energy_eV=omega_au * HA_TO_EV, omega_au=omega_au, freq_au=freq,
                amplitude=amp, n_used=n, dt_au=dt, subtract=subtract, window=window)


def peak_in_band(energy_eV, amplitude, lo, hi):
    """(E_peak, amp) of the largest bin inside [lo, hi] eV.

    ALWAYS search a PHYSICAL band — a global argmax lands on the DC-dominated
    omega~0 bin (the recurring trap; the audit's Stage B finding).
    """
    e = np.asarray(energy_eV, float)
    a = np.asarray(amplitude, float)
    sel = (e >= lo) & (e <= hi)
    if not sel.any():
        return None, None
    idx = np.where(sel)[0]
    k = idx[int(np.argmax(a[idx]))]
    return float(e[k]), float(a[k])


def resolution_eV(t_total_au: float) -> float:
    """Spectral resolution Delta(hbar*omega) = 2*pi / T  in eV (bold-flag it
    when comparing peak shifts smaller than this)."""
    return 2.0 * np.pi / float(t_total_au) * HA_TO_EV


def loss_locator(n_q_t, win, n_pad, q):
    """Loss-function PEAK-LOCATOR |n_q(omega)|^2 / q^2 (positive-freq half).

    Uses the COMPLEX phasor n_q(t) (NOT .real — taking the real part folds the
    +/- omega lobes and halves the amplitude, ratio 0.500), and returns
    |n_q|^2/q^2 (NOT the bare |n_q|). LOCATOR ONLY — see module docstring caveat.
    """
    sig = np.asarray(n_q_t) * win
    full = np.fft.fft(sig, n=n_pad)               # complex FFT (BUG-A fix)
    half = full[: n_pad // 2 + 1]
    return (np.abs(half) ** 2) / (q ** 2)         # |n_q|^2/q^2 (BUG-B fix)


def pipeline_stages(time_au, values, *, window="hann", subtract="mean",
                    zero_pad=4, t_start_au=0.0) -> dict:
    """The six panel stages: raw | baseline-removed | windowed | padded |
    spectrum(energy_eV, amplitude). Mirrors fft_spectrum so the figure depicts
    exactly the locked pipeline (see the figure-standard in SKILL.md)."""
    t = np.asarray(time_au, float)
    v = np.asarray(values, float)
    if t_start_au > 0.0:
        mask = t >= (float(t[0]) + t_start_au)
        t, v = t[mask], v[mask]
    n = len(t)
    dt = float(t[1] - t[0])
    base = subtract_baseline(v, subtract)
    win = make_window(window, n)
    windowed = base * win
    zp = max(1, int(zero_pad))
    padded = np.zeros(n * zp, float)
    padded[:n] = windowed
    t_pad = t[0] + dt * np.arange(len(padded))
    spec = fft_spectrum(time_au, values, window=window, subtract=subtract,
                        zero_pad=zero_pad, t_start_au=t_start_au)
    return dict(t=t, raw=v, baseline_removed=base, window=win, windowed=windowed,
                t_pad=t_pad, padded=padded, energy_eV=spec["energy_eV"],
                amplitude=spec["amplitude"])


# ----------------------------------------------------------------------------
# Self-test — known analytic answers, fixed up front (run: python fourier_kernel.py)
# ----------------------------------------------------------------------------
def _selftest():
    rng_n, dt = 2000, 0.5
    t = np.arange(rng_n) * dt
    k0 = 40                                        # integer => on-bin tone
    f0 = k0 / (rng_n * dt)
    A = 1.7
    tone = A * np.cos(2.0 * np.pi * f0 * t)

    # 1. coherent gain: every window recovers A at the peak (boxcar exact).
    for w, tol in (("boxcar", 1e-6), ("hann", 2e-2), ("hamming", 3e-2),
                   ("blackman", 3e-2)):
        s = fft_spectrum(t, tone, window=w, subtract="none", zero_pad=1)
        k = int(np.argmax(s["amplitude"]))
        assert abs(s["amplitude"][k] - A) <= tol * max(1.0, A) + 1e-6, (w, s["amplitude"][k])

    # 2. mean removal strips a constant offset; the true peak returns (not DC).
    s = fft_spectrum(t, tone + 50.0, window="boxcar", subtract="mean", zero_pad=1)
    kpk = int(np.argmax(s["amplitude"]))
    assert s["freq_au"][kpk] == max(s["freq_au"][kpk], 0.0) and kpk != 0

    # 3. zero-pad keeps the peak energy invariant.
    s1 = fft_spectrum(t, tone, window="hann", subtract="none", zero_pad=1)
    s4 = fft_spectrum(t, tone, window="hann", subtract="none", zero_pad=4)
    e1, _ = peak_in_band(s1["energy_eV"], s1["amplitude"], 0.01, 1e9)
    e4, _ = peak_in_band(s4["energy_eV"], s4["amplitude"], 0.01, 1e9)
    assert abs(e1 - e4) < 0.05 * e1

    # 4. ANGULAR energy convention: a tone at hbar*omega = 5 eV peaks at 5 eV.
    omega0_au = 5.0 / HA_TO_EV
    f5 = omega0_au / (2.0 * np.pi)
    sig5 = np.cos(2.0 * np.pi * f5 * t)
    s5 = fft_spectrum(t, sig5, window="hann", subtract="mean", zero_pad=4)
    epk, _ = peak_in_band(s5["energy_eV"], s5["amplitude"], 3.0, 7.0)
    assert abs(epk - 5.0) < 0.1, epk

    # 5. loss_locator: forward phasor peaks in +freq half; .real folding quarters
    #    the |.|^2 peak (ratio 0.25); |n_q|^2/q^2 scales as 1/q^2.
    om = 0.30
    win = np.ones(rng_n)
    full = loss_locator(np.exp(+1j * om * t), win, rng_n, q=0.5)
    fold = loss_locator(np.cos(om * t).astype(complex), win, rng_n, q=0.5)
    assert abs(fold.max() / full.max() - 0.25) < 0.05
    s_q1 = loss_locator(np.exp(+1j * om * t), win, rng_n, q=0.5).max()
    s_q2 = loss_locator(np.exp(+1j * om * t), win, rng_n, q=1.0).max()
    assert abs(s_q2 / s_q1 - 0.25) < 1e-6

    # 6. resolution helper
    assert abs(resolution_eV(t[-1]) - 2.0 * np.pi / t[-1] * HA_TO_EV) < 1e-9

    print("fourier_kernel._selftest: all assertions passed")


if __name__ == "__main__":
    _selftest()
