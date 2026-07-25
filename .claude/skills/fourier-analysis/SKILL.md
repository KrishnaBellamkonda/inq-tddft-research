---
name: fourier-analysis
description: Use when Fourier-transforming any TDDFT time-series observable (dipole, current, energy, density n_q) to find spectral peaks — plasmon resonances, response frequencies, loss-function peaks. Applies the one locked pipeline (mean baseline → Hann → 4× zero-pad → coherent-gain rfft, angular ħω axis), searches peaks in a physical band (never global argmax), emits the 6-stage diagnostic panel, and uses |n_q|²/q² only as a peak-LOCATOR with the −Im[1/ε] caveat.
---

# Fourier-transform analysis (TDDFT observables)

Turn a time-series `s(t)` into a spectrum and read OFF a physical frequency
(plasmon ħω_p, a response line, a loss-function peak) — with every processing
step a deliberate, audited choice so the spectrum is never a black box.

All kernels live skill-locally in **`fourier_kernel.py`** (numpy-only, portable;
run `python fourier_kernel.py` for the self-test). The production code is
`inqview.analysis.fourier.FourierTransform` (same pipeline) and
`inqview.visualisation.fourier_panel` (the figure standard) — reach those by this
method, not as black boxes.

## The locked pipeline (one uniform path for EVERY observable)

```
transient cut (drop t < t_start)
  → baseline:  'mean'          ← canonical default, ALL observables (verdict 2026-06-25)
  → window:    Hann
  → zero-pad:  ×4              ← axis interpolation only, NO new information
  → rfft + COHERENT-GAIN norm  ← / Σwin, interior bins ×2  (unit tone → ~1, any window)
  → energy axis: ANGULAR  ħω = 2π·f · HA_TO_EV   (eV)
```

`fft_spectrum(t, values, window='hann', subtract='mean', zero_pad=4, t_start_au=…)`
runs all of it and returns `energy_eV` + `amplitude`.

## Step-by-step — what each step is and WHY

1. **Transient cut.** Drop the injection/kick shake-up (`t < t_start`) BEFORE
   anything else, so the spectrum reflects the steady-state response (per
   observables_reference §13.6). Apply the cut first; baseline removal then sees
   only the steady part.

2. **Baseline removal = `mean` (uniform default).** Subtract the DC level from
   every observable — dipole, current, AND energy. This is the user's deliberate
   **override** (2026-06-25) of the fft-drift-removal dossier's per-observable
   split (`initial` for dipole/current, `detrend` for energy). One deterministic
   rule beats per-case branching.
   - **Honest caveat:** `mean` removes a constant but NOT a linear ramp. An energy
     signal's conservation drift is a *slope*, so an energy spectrum keeps a
     residual low-ω feature. **Peak positions are unaffected** (verified 3e-7) —
     only the near-DC region moves. If a genuine slope dominates, `detrend` is the
     remedy, which is exactly why it rides along as the comparison ↓.
   - **Always show `detrend` as a dashed COMPARISON** on the spectrum (the panel
     does this automatically) so the baseline choice is auditable.
   - Without subtraction a constant offset hijacks the ω≈0 bin → a bogus DC "peak".

3. **Window = Hann.** Tapers the endpoints to suppress spectral leakage from the
   finite record. Hann is the default; flattop only for amplitude-accurate work,
   Kaiser/Tukey for special trade-offs (production kernel exposes them).

4. **Zero-pad ×4.** Densifies the frequency axis (smoother interpolation between
   true bins) — it adds **no information** and does **not** improve resolution.
   Real resolution is set by the record length T (next item).

5. **Coherent-gain normalisation.** Divide by `Σwin` (Harris 1978), not by `n`,
   and ×2 the interior one-sided bins, so a unit-amplitude tone returns ~1.0 for
   *any* window. (Boxcar has `Σwin=n`, so this is backward-compatible.)

6. **Angular energy axis.** `fft_spectrum` returns `freq_au` (ordinary
   `rfftfreq`, cycles/a.u.) AND `energy_eV = 2π·freq_au·HA_TO_EV`. **Report in
   `energy_eV`.** Forgetting the `2π` is the recurring bug (a 5 eV tone misreads as
   ~0.8 eV·… — off by 2π); the kernel + panel bake it in.

## Peak attribution — the one hard rule

**ALWAYS locate the peak inside a PHYSICAL band, never by global argmax.** A
global max lands on the DC-dominated ω≈0 bin (the audit's Stage B trap). Use
`peak_in_band(energy_eV, amplitude, lo, hi)` with a band motivated by physics
(e.g. the plasmon `ħω_p = √(3/r_s³)·HA_TO_EV` ± a window). Cross-checks before
trusting a peak:
- **Resolution (informational, NOT a gate).** `Δ(ħω) = 2π/T` (`resolution_eV(T)`)
  is the bin width — annotate it for context if useful, but it NEVER blocks or
  withholds a loss function / spectrum. When one is requested, produce it
  regardless of whether `T` resolves every feature; a sub-bin shift is not a
  reason to refuse or to distrust by default. (User decision, 2026-06-26,
  `feedback_fourier_loss_function_gate`.)
- **Baseline-invariance.** The peak energy must not move between `mean` and the
  `detrend` comparison (it shouldn't — verified 3e-7). If it does, suspect a
  drift artefact and investigate.
- **Window-invariance.** Peak position robust across Hann/Hamming/Blackman.

## Loss function `|n_q(ω)|²/q²` — LOCATOR ONLY

For density modes `n_q(t)` (jellium plasmon), `loss_locator(n_q_t, win, n_pad, q)`
returns `|n_q(ω)|²/q²` from the **complex** phasor on the positive-frequency half.
Two fixes are baked in (audit-confirmed on E15, applied to
`pipeline/density_fourier.py` via `code-test`):
- **Use the complex `n_q(t)`, not `.real`** — taking the real part folds ±ω and
  halves the amplitude (ratio 0.500).
- **`|n_q|²/q²`, not the bare `|n_q|`.**

**CAVEAT (keep in every caption, bold):** this is a plasmon-peak **locator** —
right pole positions and the right `1/q²` q-trend — but **not** a quantitatively
faithful `−Im[1/ε]`: off the undamped limit the line shape (Lorentzian² vs
Lorentzian), the spectral area, and the absolute `4π` Coulomb normalisation are
all wrong. (`docs/validation/loss-function-formula-validation.md`, verdict
2026-06-25.) Do not read amplitudes/areas as a true loss function.

## Figure standard — the 6-stage panel (mandatory)

Every FFT of a signal ships a 3×2 diagnostic panel so the spectrum is auditable:

```
row 1:  raw signal            |  baseline-removed (mean)
row 2:  windowed signal       |  zero-padded signal
row 3:  |FFT| linear scale     |  |FFT| log scale   (mean solid + detrend dashed)
```

Render with `inqview.visualisation.fourier_panel.fft_pipeline_panel(t, values,
peak_band=…, fmax=…)` (it overlays the detrend comparison and annotates the
band peak + `Δω=2π/T`). `fourier_kernel.pipeline_stages(...)` returns the same
stage arrays if you build the figure by hand. This standard is also wired into
the `run-notebook` and `notebook-making` skills.

## Known-answer anchors (validate before trusting a new pipeline)

| Stage | Input | Known answer |
|---|---|---|
| Synthetic | unit on-bin tone | coherent gain → ~1.0 every window |
| Synthetic | complex phasor vs `.real` | `.real` folding ratio 0.500 (→ 0.25 in power) |
| Synthetic | equal-amplitude `n_q` modes | `|n_q|²/q²` ∝ 1/q² |
| QKE | Li v=0.0626 multi-k energy series | plasmon **≈6.48 eV** (paper 6.5), baseline-invariant |
| E15 jellium | `n_q` long run | `ω_p ≈ 3.47 eV`, `Δω ≈ 0.09 eV` |

Worked example notebook (three-stage audit, executed):
`docs/validation/stopping-power-extraction/fourier_analysis.ipynb`.

## Validation status

- `fourier_kernel._selftest()` — 6 known-answer assertions PASS (coherent gain,
  mean-removal, zero-pad invariance, angular 5 eV convention, loss_locator
  folding 0.25 + 1/q², resolution helper).
- Production parity tested: `inq-stack/tests/python/inqview/analysis/test_fourier.py`
  (19/19, default=mean), `…/pipeline/test_density_fourier_loss.py` (3/3, BUG-A/B),
  `…/visualisation/test_fourier_panel.py` (4/4, panel + detrend overlay).
  Catalogue: `docs/validation/test-catalogue.md`.

## Sources & dossiers
- `docs/validation/fft-normalization-validation.md` — coherent gain (verdict: correct).
- `docs/validation/fft-drift-removal-validation.md` — baseline (verdict: uniform `mean` override).
- `docs/validation/loss-function-formula-validation.md` — `|n_q|²/q²` locator (verdict: accepted as locator).
- Harris 1978 (windows); Yabana–Bertsch (induced-response IC); Giuliani–Vignale §4–5 (longitudinal loss function).
- `feedback_fourier_loss_function_gate` — loss-function production is UNCONDITIONAL (no resolution gate); `Δω` is an optional informational annotation, never a block.
