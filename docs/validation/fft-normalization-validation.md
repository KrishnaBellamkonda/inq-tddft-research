# FFT normalization validation (independent)

Independent signal-processing review of the amplitude/power normalization
convention in `inq-stack/python/inqview/fourier.py`
(`FourierTransform.transform()`). Fresh-context derivation, cross-checked
against a pure-numpy numerical experiment and authoritative DSP references.
No accept/reject verdict is given — verdict lines at the end are left blank
for the user.

All numerical claims below were reproduced with
`/local/data/public/skcb2/tddft/venv/bin/python3` (numpy + scipy), not
asserted from memory.

---

## Code convention as-implemented

For a real time series `values` sampled uniformly at `dt` (length `n`), the
transform does (lines 145–196 of `fourier.py`):

1. §13.6 transient cut: drop samples with `t < t_start_au` (optional).
2. Linear detrend: `scipy.signal.detrend(values, type="linear")` (default on).
3. Window: `sig_win = sig * win`, where `win = get_window(name, n)`
   (default Hann), length `n`.
4. Zero-pad: `n_pad = n * zero_pad` (default `zero_pad = 4`); the windowed
   signal occupies `[:n]`, the rest is zeros.
5. `raw = np.fft.rfft(sig_padded)`; `freq = np.fft.rfftfreq(n_pad, d=dt)`.
6. Normalization:
   - `amplitude = np.abs(raw) / n`   ← divide by the **original** length `n`,
     not `n_pad`.
   - `amplitude[1:-1] *= 2.0`        ← one-sided interior-bin doubling
     (DC and Nyquist not doubled).
   - `power = amplitude ** 2`.
7. Optional Gaussian smoothing of `amplitude` (frequency-bin domain).

The single missing ingredient relative to a *calibrated* amplitude spectrum
is a **window coherent-gain correction**. Everything else (interior
doubling, division by original `n`, DC/Nyquist exclusion from doubling) is
the standard one-sided convention and is correct.

---

## Q1 one-sided amplitude

**Setup.** `x(t) = A cos(2π f₀ t)`, sampled at `n` points, **no window**
(boxcar), `zero_pad = 1`, `f₀` exactly on a bin (integer number of cycles in
the record).

**Derivation.** The DFT of a real cosine that lands exactly on bin `k₀` is
two delta-like spikes, at `+k₀` and `−k₀` (the latter aliased to `n−k₀`),
each of magnitude `|X[k₀]| = A·n/2`. `np.fft.rfft` returns only the
non-negative-frequency half (bins `0 … n/2`), so the `+f₀` spike alone
appears with `|raw[k₀]| = A·n/2`.

The code then computes:

```
amplitude[k₀] = |raw[k₀]| / n        = (A n / 2) / n = A/2
amplitude[k₀] *= 2  (interior bin)   = A
```

So with a **boxcar** window the code returns exactly **A** — the physical
amplitude. **Verified numerically:** `A = 3.0`, `n = 1000`, boxcar →
peak amplitude `3.00000` at the correct bin.

**Exactly-correct convention to recover A** (real one-sided amplitude
spectrum, tone on a bin, rectangular window):
`A = (2/n)·|rfft| ` for interior bins, `(1/n)·|rfft|` for the DC and Nyquist
bins. That is precisely what the code does. The interior×2 + ÷n combination
is the right convention for recovering the peak amplitude of a real tone.

Caveat (not a bug, but a known DFT property): if `f₀` does **not** fall on a
bin, scalloping loss reduces the recovered peak (up to ~3.9 dB ≈ ×0.64 for a
boxcar worst case at the half-bin). Windowing widens the main lobe and
greatly reduces this scalloping loss — one reason a Hann/Kaiser window is
preferred for off-bin tones — but windowing introduces the coherent-gain
deficit of Q2.

---

## Q2 window coherent gain (amplitude) vs sum(win²) (PSD)

**Coherent gain.** Harris (1978) defines the *coherent gain* of a window as
the sum of its coefficients, `Σ win`, and tabulates it normalized by its
maximum `N`, i.e. `mean(win) = Σ win / n`. For a unit-amplitude on-bin tone,
windowing scales the recovered DFT spike by exactly this factor: a coherent
sinusoid is multiplied pointwise by the window and the DFT of the product at
the tone bin is `(A·n/2)·mean(win)` (interior). Values:
`mean(win) ≈ 1.0` (boxcar), `≈ 0.5` (Hann).

**Does the current code under-report by mean(win)?** Yes. The code divides
only by `n`, never by `Σ win`, so for a windowed tone it returns
`A·mean(win)`, not `A`. **Verified numerically** (`A = 3.0`, on-bin):

| window  | returned peak | mean(win) | peak / mean(win) |
|---------|---------------|-----------|------------------|
| boxcar  | 3.00000       | 1.00000   | 3.00000 (= A)    |
| hann    | 1.50000       | 0.50000   | 3.00000 (= A)    |

So a Hann-windowed unit tone is returned as **≈ 0.5**, a boxcar-windowed
unit tone as **≈ 1.0**. The two are not on the same amplitude scale — a
calibration error of `mean(win)` whenever a non-rectangular window is used.

**Correct amplitude fix:** divide by the window's coherent gain instead of
(or in addition to) `n`:

```
amplitude = np.abs(raw) / win.sum()   # = (|raw|/n) / win.mean()
amplitude[1:-1] *= 2.0
```

Equivalently, keep the existing `/n` and add `amplitude /= win.mean()`. This
restores `A` for **any** window (boxcar `mean=1` is unchanged, so the fix is
backward-compatible for the rectangular case).

**Amplitude vs PSD — the distinction must be explicit.**

- **Amplitude spectrum** (recovering the height `A` of coherent tones): use
  the **coherent gain** `Σ win` (Harris 1978). This is what `fourier.py`'s
  `amplitude`/`power = amplitude²` is implicitly trying to be — a *peak
  amplitude* and its square, appropriate for line spectra (discrete tones,
  oscillation amplitudes).

- **Power spectral density** (variance per Hz of broadband / noise-like
  signals): use the **equivalent noise bandwidth**, driven by `Σ win²`.
  Harris: `ENBW = N · Σ(win²) / (Σ win)²`. The one-sided PSD is

  ```
  PSD = (2 / (fs · Σ(win²))) · |raw|²        # interior bins
  ```

  with `fs = 1/dt`, which makes `Σ PSD · Δf ≈ var(signal)` independent of the
  window. The `Σ win²` normalization is *wrong* for recovering a tone's
  amplitude, and `Σ win` (coherent gain) is *wrong* for a calibrated noise
  PSD. They differ by exactly the normalized ENBW factor.

**Which does this codebase want?** For TDDFT oscillation amplitudes / dipole
or current response peaks (line-like features), the **amplitude / coherent-
gain** convention is the right target, so the recommended fix is the
coherent-gain correction. If the user ever wants a calibrated *power
spectral density* of a broadband signal, that is a separate `sum(win²)`/ENBW
normalization and should be a distinct method, not a patch to `amplitude`.

---

## Q3 zero-pad

Dividing by the **original** `n` (not `n_pad`) is **correct**, and is the
key reason zero-padding here only interpolates the frequency axis without
changing the peak height.

**Why.** Zero-padding appends zeros; it adds no signal energy and does not
change the value of the underlying continuous-frequency transform — it only
samples that transform on a denser grid (`rfftfreq(n_pad, dt)` has
`n·zero_pad/2 + 1` bins). The amplitude of an on-bin tone is set by the
`n` non-zero samples; the correct denominator is therefore `n`, the number
of *real* (non-zero) data points. Dividing by `n_pad` would shrink every
amplitude by `1/zero_pad` — a spurious, pad-dependent attenuation.

**Verified numerically:** with `zero_pad ∈ {1,4}` the recovered on-bin peak
amplitude is unchanged (the extra bins merely interpolate the lobe shape).
Conclusion: `÷ n` is right; `÷ n_pad` would be a bug.

Minor note: zero-padding plus interior×2 makes the discrete *sum* of `power`
across the denser grid no longer a clean Parseval estimate (the lobe is now
oversampled, so a naive bin sum over-counts). Parseval-type checks (Q4)
should be done at `zero_pad = 1`, or with an explicit `Δf` integration that
accounts for the oversampling. Peak-amplitude reads are unaffected.

---

## Q4 Parseval

**Identity that should hold (zero_pad = 1).** Let `s = win * detrend(values)`
be the windowed signal and `raw = rfft(s)`. Discrete Parseval for the DFT:

```
Σ_{m=0}^{n-1} s[m]²  =  (1/n) Σ_{k=0}^{n-1} |X[k]|²
```

In one-sided `rfft` form with the same interior-doubling used for power:

```
mean(s²)  =  Σ_k P_onesided[k],
   where P_onesided[k] = (|raw[k]| / n)²,  and interior bins ×2
   (DC and Nyquist not doubled).
```

So **the sum of the one-sided power equals the mean-square (variance, since
detrend removes the mean) of the windowed signal**. **Verified numerically**
(boxcar, on-bin tone `A=3`): `Σ P_onesided = 4.5` and `mean(s²) = 4.5`
(ratio `1.0` to machine precision; note `A²/2 = 4.5`).

**With the coherent-gain amplitude fix in place**, the `amplitude` array is
scaled up by `1/mean(win)`, so `power = amplitude²` is scaled by
`1/mean(win)²`. The Parseval check must then carry the window factor
explicitly. The clean, window-independent statement is to phrase Parseval on
the **windowed** signal `s` before any coherent-gain rescaling:

```
Σ_k (|raw[k]|/n)²  (interior ×2)  ==  mean(s²)  ==  mean( (win·detrend(values))² )
```

This is the relation a unit test should assert. To relate it back to the
*unwindowed* variance you would divide the RHS by `mean(win²)` (the ENBW /
power-correction factor), which is exactly why amplitude calibration uses
`Σ win` but power/variance calibration uses `Σ win²`.

---

## Q5 cross-run amplitude comparability

The code's TODO (lines 8–10) flags exactly this: comparing an INQ spectrum
against a QBall / `QuantumKickExtension` spectrum. For the **amplitudes**
(not just peak *positions*) to be comparable, every per-window / per-length
scale factor must be equalized:

1. **Window coherent gain (`mean(win)`).** This is the dominant pitfall. If
   one run uses Hann and the other boxcar, their amplitudes differ by ~2×
   purely from `mean(win)` (Q2). Either use the *same* window for both runs,
   or apply the coherent-gain correction so each is on the absolute `A`
   scale. Without the fix, "INQ peak is half the QBall peak" can be a pure
   windowing artefact.

2. **Detrending.** Both must use the same detrend setting; linear detrend on
   one and not the other shifts low-frequency content and can change a
   near-DC peak's height.

3. **Transient cut (`t_start_au`) and record length `n`.** Different `n` (or
   different cut) changes frequency resolution `Δf = 1/(n·dt)` and the
   scalloping/leakage of any off-bin tone. For amplitude comparison the
   safest route is to (a) match the *windowed coherent gain* via the fix,
   and (b) report whether each tone is on- or off-bin (resolution-limited).
   Identical `n` and identical `dt` remove this entirely.

4. **Sampling interval `dt` / units.** The amplitude convention here is
   independent of `dt` (it divides by `n`, a count, not by `fs`), so a
   *peak amplitude* read is `dt`-independent — good. But a *PSD* read would
   be `dt`-dependent (`1/fs` factor), so do not mix the two between runs.

5. **Same observable scaling / units.** Ensure the INQ and QBall columns are
   in the same physical units before comparing magnitudes (atomic units vs
   SI, dipole per electron vs total, etc.) — outside `fourier.py` but
   essential for amplitude (not just peak) comparison.

**Bottom line:** with the coherent-gain fix applied and identical units,
amplitudes become window-independent and directly comparable even if the two
runs used different windows; without it, the two runs must use the *same*
window *and* same `n`/`dt`/detrend for amplitudes to mean the same thing.

---

## Recommended fix + portable analytic test (numpy pseudocode, expected values, tolerances)

### Minimal corrected normalization

Replace lines 179–180 of `fourier.py`:

```python
# current
amplitude = np.abs(raw) / n
amplitude[1:-1] *= 2.0
```

with the coherent-gain-corrected amplitude:

```python
# corrected: calibrated one-sided amplitude spectrum (Harris 1978 coherent gain)
cg = win.sum()                    # coherent gain = Σ win  (= n * mean(win))
amplitude = np.abs(raw) / cg      # divide by window coherent gain, not by n
amplitude[1:-1] *= 2.0            # one-sided interior doubling (DC/Nyquist excluded)
power = amplitude ** 2
```

Notes:
- `win.sum()` reduces to `n` for boxcar, so the rectangular case is
  unchanged — backward compatible.
- This recovers the true tone amplitude `A` for *any* window.
- If a calibrated *power spectral density* is ever needed, add a **separate**
  method using `2 / (fs * np.sum(win**2)) * np.abs(raw)**2` (ENBW / `Σ win²`
  normalization) — do not overload `amplitude`.

### Portable pure-numpy unit test

```python
import numpy as np

def one_sided_amplitude(values, dt, win, zero_pad=1):
    """Calibrated one-sided amplitude spectrum (coherent-gain corrected)."""
    n = len(values)
    sig = values - values.mean()              # (test uses explicit drift removal below)
    sig_win = sig * win
    n_pad = n * zero_pad
    padded = np.zeros(n_pad); padded[:n] = sig_win
    raw = np.fft.rfft(padded)
    freq = np.fft.rfftfreq(n_pad, d=dt)
    amp = np.abs(raw) / win.sum()             # <-- coherent-gain normalization
    amp[1:-1] *= 2.0
    return freq, amp

def test_amplitude_calibration():
    n, dt, A = 1000, 0.1, 3.0
    cycles = 50                               # integer -> f0 lands exactly on a bin
    f0 = cycles / (n * dt)
    t = np.arange(n) * dt
    drift = 0.7 + 0.002 * t                   # linear drift, must be removed by detrend
    x = A * np.cos(2 * np.pi * f0 * t) + drift

    # detrend exactly cancels the linear drift before windowing
    p = np.polyfit(t, x, 1)
    x_dt = x - np.polyval(p, t)

    for wname, win in [("boxcar", np.ones(n)), ("hann", np.hanning(n))]:
        freq, amp = one_sided_amplitude(x_dt, dt, win, zero_pad=1)
        i = np.argmax(amp)
        # (a) peak at f0
        assert np.isclose(freq[i], f0, rtol=0, atol=1.0 / (n * dt) / 2), wname
        # (b) recovered amplitude == A for BOTH windows after coherent-gain fix
        assert np.isclose(amp[i], A, rtol=2e-2, atol=1e-3), (wname, amp[i])

    # (c) Parseval on the WINDOWED signal (boxcar, zero_pad=1):
    win = np.ones(n)
    s = x_dt * win
    raw = np.fft.rfft(s)
    P = (np.abs(raw) / n) ** 2                # NOTE: /n here (raw Parseval form)
    P[1:-1] *= 2.0
    assert np.isclose(P.sum(), np.mean(s ** 2), rtol=1e-6, atol=1e-9)

    print("FFT amplitude-calibration + Parseval tests passed.")

test_amplitude_calibration()
```

### Expected values and tolerances

| Check | Expected | Tolerance (machine-independent) |
|---|---|---|
| (a) peak frequency | `f0 = 0.5` (here) | `atol = Δf/2 = 1/(2 n dt)` — within half a bin |
| (b) boxcar peak amp | `A` (e.g. 3.0) | `rtol = 2e-2, atol = 1e-3` |
| (b) hann peak amp **after fix** | `A` (e.g. 3.0) | `rtol = 2e-2, atol = 1e-3` |
| (b) hann peak amp **before fix** | `A·mean(win) ≈ A/2` | (documents the bug; not asserted as correct) |
| (c) Parseval (windowed, raw `/n` form) | `Σ P_onesided == mean(s²)` | `rtol = 1e-6, atol = 1e-9` |

Rationale for tolerances: on-bin tones are recovered to machine precision in
principle, but the `2e-2` relative tolerance on amplitude leaves head-room
for tiny leakage from the finite record / detrend interaction and keeps the
test portable across numpy/scipy versions and BLAS backends (no bit-exact
assumption). The Parseval identity is an algebraic DFT property and holds to
~`1e-12`; `1e-6` is a safe, portable bound. The frequency tolerance of half
a bin is the correct resolution-limited criterion. Use `zero_pad = 1` for
the Parseval check (see Q3).

**Reproduced result for the as-implemented code** (before fix), confirming
the table in Q2: boxcar peak `3.00000`, hann peak `1.50000`,
`peak/mean(win) = 3.00000` for both; Parseval ratio `1.0`.

---

## Sources cited

- **F. J. Harris**, "On the use of windows for harmonic analysis with the
  discrete Fourier transform," *Proc. IEEE* **66**(1), 51–83 (1978). The
  authoritative window/figure-of-merit paper. Defines *coherent gain* =
  `Σ win` (Table 1 lists it normalized by `N`, i.e. `mean(win)`), *processing
  gain* = reciprocal of normalized ENBW, and equivalent noise bandwidth
  `ENBW = N · Σ(win²) / (Σ win)²`. PDF mirror:
  <https://www.cs.cmu.edu/afs/cs/user/bhiksha/WWW/courses/dsp/spring2013/WWW/schedule/readings/windows_comparison2_harris.pdf>
- **A. V. Oppenheim & R. W. Schafer**, *Discrete-Time Signal Processing*
  (Prentice Hall). Standard reference for the DFT, Parseval's theorem, and
  the one-sided real-spectrum conventions (rfft, interior-bin doubling,
  DC/Nyquist handling).
- **numpy.fft documentation** — definitions of `np.fft.rfft` /
  `np.fft.rfftfreq` and the unnormalized forward-transform convention
  (`raw[k] = Σ_m s[m] e^{-2πi k m / n}`), which is why the `1/n` (or
  `1/Σ win`) division is applied by hand:
  <https://numpy.org/doc/stable/reference/routines.fft.html>
- **MATLAB `enbw` reference** (corroborating the ENBW = `N·Σ(win²)/(Σ win)²`
  definition used above):
  <https://www.mathworks.com/help/signal/ref/enbw.html>
- GaussianWaves, "Window function – figure of merits" (secondary, corroborating
  coherent-gain vs ENBW summary; low-trust, used only for cross-checking the
  primary Harris definitions):
  <https://www.gaussianwaves.com/2020/09/window-function-figure-of-merits/>

All numerical results above were produced with
`/local/data/public/skcb2/tddft/venv/bin/python3` (numpy + scipy), not taken
from any source.

---

## Verdict (LEAVE BLANK — user fills)

- Current convention correct?  Y/N ___
- Apply coherent-gain fix?     Y/N ___
