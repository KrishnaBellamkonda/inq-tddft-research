# Report: spectra-preprocessing comparison — QBall `td_kicks` vs. inqview

**Author:** project session, 2026-05-05
**Status:** comparative report, no new computations yet — actionable
diff against the existing `inqview/fourier.py`.
**Sources used:**
- `QuantumKickExtension/qball-codebase/Li/td_kicks/analyse.py`
- `QuantumKickExtension/qball-codebase/Li/td_kicks/analyse_pub.py`
- `QuantumKickExtension/qball-codebase/Li/td_kicks/fft_peak_attribution.py`
- `QuantumKickExtension/qball-codebase/Li/td_kicks/fft_windowing_comparison.py`
- `inq-stack/python/inqview/fourier.py`
- `inq-stack/python/inqview/postprocess/observables.py`
- `docs/observables_reference.md §11` (extended preprocessed spectra)
- `docs/sources/correa-2018-electronic-stopping-power.md` (Correa 2018
  §6, transient discussion)

## 1. Summary verdict

The QBall reference scripts and the inqview Fourier pipeline agree on
**windowing** (Hann), **zero-padding** (8× in QBall, 4× default in
inqview, both configurable), and **normalisation** (one-sided
amplitude with interior-bin doubling). They **disagree** on:

1. **Transient exclusion before the FFT** — QBall removes the first
   half of the time series implicitly (via `dE − dE[N//2:].mean()`)
   and treats the second half as the steady-state plateau; inqview
   currently uses the *full* trace with optional `scipy.signal.detrend`
   linear-detrend.
2. **Detrending strategy** — QBall: subtract plateau mean (DC) and
   keep the time series otherwise as-is; inqview: subtract a linear
   least-squares fit of the *whole* trace.
3. **Energy-component spectra coverage** — QBall produces FFTs of
   `<etotal>`, energy components (`<ekin>`, `<ecoul>`, `<exc>`,
   `<enl>`), and electronic momentum; inqview currently only FFTs
   `energy_total`, `current_*`, `dipole_*`. No per-component energy
   spectrum.
4. **Multi-window comparison harness** — QBall has
   `fft_windowing_comparison.py` (Rectangular vs. Hann vs. Blackman)
   for distinguishing real peaks from leakage; inqview has no
   equivalent.
5. **Peak attribution** — QBall's `fft_peak_attribution.py` finds
   peaks via `scipy.signal.find_peaks` and labels them against a
   hard-coded physics dictionary (plasmon, second harmonic,
   Drude tail, …). inqview just plots the spectrum.

The most consequential of these is **#1 (transient exclusion)** —
this is the one Correa 2018 §6 explicitly recommends, and the one the
user has flagged. The others are stylistic improvements that can
follow.

## 2. Detailed comparison table

| Step | QBall (`analyse.py` lines 246–271; `analyse_pub.py` `_plot_fft_domain`; `fft_peak_attribution.py::compute_fft`) | inqview (`fourier.py::FourierTransform.transform`) |
|---|---|---|
| Input series | `dE_uc(t)` parsed from `Li.54.out` (long format) | `pandas.DataFrame` column from `observables.csv` |
| **Transient handling** | `dE_osc = dE − dE[N//2:].mean()` (subtracts plateau mean = second half) | `scipy.signal.detrend(values, type='linear')` if `detrend=True` (default), full trace |
| Window | `np.hanning(N)` (Hann) | `WindowSpec("hann")` default; configurable to `boxcar/hamming/blackman/tukey/kaiser/flattop` |
| Zero-pad multiplier | 8× (`n=N*8`) | 4× default (`zero_pad=4`); configurable |
| FFT | `np.fft.rfft` | `np.fft.rfft` |
| Frequency axis | `rfftfreq(N*8, d=dt) * 2π * Ha2eV` (eV) | `rfftfreq(n_pad, d=dt_au)` (a.u.); eV conversion at plot time |
| Output | `|FFT|²` (power) normalised to peak in 0–20 eV mask | one-sided amplitude `|raw|/N`, doubled at interior bins; `power = amplitude²` |
| Optional smoothing | none | `gaussian_filter1d(σ_bins)` post-FFT, off by default |
| Multi-window comparison | `fft_windowing_comparison.py` runs Rectangular/Hann/Blackman side-by-side and prints peak/FWHM/sidelobe table | not implemented |
| Peak finding | `scipy.signal.find_peaks(threshold=0.05, min_dist=1.5 eV)` | not implemented |
| Peak attribution | `attribute_peak()` hard-codes Li bulk plasmon (6.0–7.0 eV), 2nd harmonic, Drude, etc. | not implemented |
| Components covered | `etotal`, `ekin`, `ecoul`, `exc`, `enl`, `Pe_x` (electronic momentum), and ratios | `energy_total`, `current_{x,y,z}`, `dipole_{x,y,z}` |

## 3. Why QBall's "subtract plateau-mean" works (and what we are missing)

QBall's preprocessing line `dE_osc = dE − dE[N//2:].mean()` does two
things at once:

1. **Removes DC.** `dE` plateaus at some non-zero value $E_\infty$
   reflecting the steady-state energy absorbed by the bath; the FFT
   would otherwise have a giant 0-Hz spike.
2. **Suppresses the transient.** During the first half of the run,
   `dE − dE[N//2:].mean()` is negative (the system is still building
   up to its plateau); during the second half, the residual is the
   steady-state oscillation around the plateau. The Hann window then
   downweights both endpoints, but the **effective transient
   suppression comes from the choice of subtraction baseline**, not
   from the window.

This is the same idea as Correa 2018 §6 "discard the transient" — but
implemented as a simple plateau-mean detrend rather than as an
explicit `t_start_au` cutoff. In practice the two are nearly
equivalent if the transient is much shorter than $N/2$ time steps.

**Inqview's `scipy.signal.detrend(type='linear')` does NOT achieve
this:** it fits a linear trend $a + b t$ over the *whole* series and
subtracts it. That removes a slope but **leaves the transient in
place** — and worse, it lets the transient bend the linear fit
itself, distorting the steady-state region.

## 4. Recommended changes to inqview

### 4.1 Add `t_start_au` to `FourierTransform` (highest priority — fixes the user's complaint about transient contamination)

```python
@dataclass
class FourierTransform:
    window: WindowSpec | None = None
    detrend: bool = True
    zero_pad: int = 4
    smooth_sigma_bins: float = 0.0
    t_start_au: float = 0.0     # NEW — drop samples with t < t_start_au
    detrend_strategy: str = "linear"  # NEW — "linear" | "plateau_mean" | "none"
```

Update `transform()`:

```python
n = len(time_au)
if self.t_start_au > 0:
    mask = time_au >= self.t_start_au
    time_au = time_au[mask]
    values  = values[mask]
    n = len(time_au)

if self.detrend_strategy == "linear":
    sig = scipy_detrend(values, type="linear")
elif self.detrend_strategy == "plateau_mean":
    sig = values - values[n // 2:].mean()        # QBall-style
elif self.detrend_strategy == "none":
    sig = values.copy()
else:
    raise ValueError(...)
```

### 4.2 Default to QBall-style preprocessing for jellium WP runs

Inqview's default is `detrend_strategy="linear", t_start_au=0`. For
jellium WP runs at WP_EKIN ≤ 5 eV, the recommended default (per the
new `observables_reference.md §13.6`) is

```python
FourierTransform(detrend_strategy="plateau_mean", t_start_au=5.0)
```

— the 5 a.u. cutoff matches the WP-injection shake-up timescale at
$\sigma_r = 5$ Bohr, $v = 0.33$ a.u.

### 4.3 Add per-component energy spectra

Mirror `transform_energy()` with `transform_energy_component(component)`:

```python
def transform_energy_component(self, observables, component: str) -> FourierResult:
    return self.transform_column(observables, f"energy_{component}")
```

Then `observables.py::_plot_spectra` emits one spectrum per component
in `analysis/observables/spectra/energy_components/`.

### 4.4 Add a windowing comparison harness

Port QBall's `fft_windowing_comparison.py` to a function
`inqview.fourier.compare_windows(time_au, values, windows=("boxcar",
"hann", "blackman"))` returning a `dict[name, FourierResult]`. Useful
for distinguishing real peaks from leakage.

### 4.5 Add peak finding + optional attribution

Port QBall's `find_peaks` step into `FourierResult.peaks(threshold,
min_dist_au)`. Attribution is project-specific (Li-plasmon hard-codes
won't apply to jellium); leave attribution to per-project caller code.

## 5. Effect on the existing journal entries

The two L=50 / E=1.5 eV entries currently flag their FFT spectra as
"transient included" and defer the proper plot to a follow-up. Once
4.1 and 4.2 land, the spectra should be regenerated with
`FourierTransform(detrend_strategy="plateau_mean", t_start_au=5.0)`
and the journal entries' §6 ("Spectra — caveat") sections updated to
point at the new plots. The L=30 entry has the same caveat.

## 6. Test plan (per development-feedback-loop rule)

For each implementation step:

1. **`t_start_au` cutoff** — known-case test: pass a sinusoid at
   $\omega_0$ with $t \in [0, T]$, with the first 20 % corrupted by a
   spike. With `t_start_au = 0.2 T`, the resulting spectrum should
   peak at $\omega_0$ with no spike-induced sidelobes; without the
   cutoff, the spike dominates the low-frequency end. Compare both
   numerically.
2. **`plateau_mean` detrend** — known-case: a sinusoid offset by a
   ramp. `linear` removes the ramp but leaves the offset; `plateau_mean`
   removes the offset (= mean of second half) and leaves the ramp.
   Confirm both behave as documented.
3. **Per-component spectra** — round-trip: feed a synthetic
   `observables.csv` with `energy_kinetic = sin(ω_kt)` and
   `energy_hartree = cos(ω_Ht)`; verify the two component spectra peak
   at the right frequencies independently.
4. **Windowing harness** — sanity: feed a pure sinusoid; Rectangular
   should have wide sidelobes, Blackman narrow ones, Hann in between.
   Match QBall's `fft_windowing_comparison.py` output qualitatively.

## 7. Conclusion

The current inqview FFT pipeline is structurally sound (Hann window,
zero-padding, one-sided normalisation match QBall) but **defaults to
the wrong detrending strategy for our use case**. The two minimal
changes — `t_start_au` cutoff and `detrend_strategy="plateau_mean"`
default for jellium runs — should be implemented before any
spectra are quoted in reports, and the existing journal entries'
spectra should be regenerated. Per-component energy spectra and a
windowing-comparison harness are useful follow-ups but not blockers
for the L=50 papers.
