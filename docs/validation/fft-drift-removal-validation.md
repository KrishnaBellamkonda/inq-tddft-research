# FFT drift-removal / baseline-subtraction validation (independent)

Independent signal-processing + TDDFT review of the **baseline-subtraction
(DC-removal) convention** applied before FFT-ing real-time TDDFT observables
in `inq-stack/python/inqview/postprocess/observables.py`
(`_extended_spectra`, `_build_variants`) and `inq-stack/python/inqview/fourier.py`.

This is a fresh-context derivation, cross-checked against a pure-numpy
numerical experiment and authoritative DSP / TDDFT references. **No
accept/reject verdict is given** — the verdict lines at the end are left
blank for the user.

The companion dossier
`docs/validation/fft-normalization-validation.md` already settled the
**window / amplitude-normalization** convention (coherent-gain, interior-bin
doubling, `÷n` vs `÷n_pad`, Parseval). That convention is **assumed fixed
here and not re-derived** — this dossier is solely about *what baseline to
subtract before windowing*. The two are orthogonal: normalization sets the
vertical scale of a peak; baseline subtraction sets what goes into the ω≈0
region.

All numerical claims below were reproduced with
`/local/data/public/skcb2/tddft/venv/bin/python3` (numpy), not asserted from
memory. Script: `/tmp/drift_test2.py` (regenerable; see Portable test).

---

## Candidates and the physics question

Real-time TDDFT writes observables sampled uniformly at `dt`:
total energy `E(t)`, current `J_{x,y,z}(t)`, dipole `d_{x,y,z}(t)`, in two
settings — (a) weak δ-kick linear response (optical spectra), (b) Gaussian
wave-packet (WP) injected at `t=0`, scattering off a jellium electron bath.
Before FFT a baseline `b(t)` must be removed: `s'(t) = s(t) − b(t)`. Three
candidates (the project's existing `raw_subtracted / mean_subtracted /
detrended` triple, `docs/observables_reference.md` §11):

| Name | `b(t)` | Removes |
|---|---|---|
| **initial-value** | `s(0)` (a constant) | a constant offset *only* |
| **mean** | `⟨s⟩_t` (a constant) | a constant offset *only* (the DC bin exactly) |
| **linear detrend** | least-squares `a + b·t` | a constant offset **and** a linear slope |

The existing kernel defaults to **linear detrend** (`fourier.py` line ~22,
`scipy.signal.detrend(values, type="linear")`).

**The physics question.** What is the correct *baseline* of each observable —
the value the FFT should treat as "zero response"? This is not a free
numerical knob: the baseline encodes a physical claim about the initial
condition and about which slow components are physics vs artefact.

**Two facts that organize the whole analysis** (both verified numerically,
Q1/Q2 below):

1. A nonzero baseline under a finite window does **not** stay at ω=0. It
   convolves with the window transform and **leaks into the lowest non-zero
   bins**, where the most interesting low-energy TDDFT physics
   (single-particle electron–hole gap, low plasmon) also lives. Removing it
   is mandatory, not cosmetic.
2. `initial` and `mean` both subtract only a **constant**; they are powerless
   against a genuine **linear slope**. Only `detrend` removes a slope. The
   distinction between candidates therefore lives entirely in (i) which
   constant, and (ii) whether you also remove a slope.

---

## Q1 induced-response baseline (current/dipole)

### The standard optical-spectrum recipe

The canonical real-time TDDFT absorption recipe (Yabana & Bertsch, PRB 54,
4484 (1996); Marques & Gross, *Time-Dependent Density Functional Theory*,
Lect. Notes Phys. 706): apply a weak δ-kick `δv(r,t) = −κ·x·δ(t)` to the
ground state, propagate, record the **induced dipole**
`d(t) = ⟨Ψ(t)|x̂|Ψ(t)⟩ − ⟨Ψ(0)|x̂|Ψ(0)⟩`, and Fourier-transform it. The
dynamic polarizability is `α(ω) = d̃(ω)/κ` (up to the one-sided / damping
convention) and the dipole strength function / absorption cross section is
`S(ω) ∝ ω·Im α(ω)`.

The object that is transformed is **already the induced quantity** `d(t) −
d(0)`: by construction the system is in the ground state for `t < 0`, the
kick acts at `t=0`, and the *induced* response satisfies `Δd(0⁻) = 0`. The
linear-response function `χ(t−t')` is **causal** and **starts from zero** —
this is the defining initial condition of a response function (Kubo;
Fetter & Walecka, *Quantum Theory of Many-Particle Systems*, ch. on linear
response).

### Why initial-value subtraction is the *physically correct* baseline here

`s(t) − s(0)` is, definitionally, the induced response. It enforces the
correct initial condition `Δs(0) = 0` exactly and introduces **no**
time-averaged quantity into the signal.

- **Mean subtraction** `s(t) − ⟨s⟩_t` mixes the response with a
  *window-dependent* time average. `⟨s⟩_t` is not a physical quantity: it
  depends on the propagation length `T`, the transient cut, and on how much
  of the oscillation happens to fit in the record. Using it as "zero" shifts
  the response by a non-physical constant. For a *symmetric, fully-resolved*
  oscillation `⟨s⟩ → s(0)` and the two agree, but for a truncated or
  asymmetric record they differ.
- **Linear detrend** removes a least-squares slope. If the genuine physical
  response contains slow, quasi-monotonic content (e.g. a real slow induced
  polarization building up over the window, or a low-frequency mode whose
  period exceeds the record), detrend will **attribute that physics to the
  baseline and delete it**. This is the central hazard: detrend is a
  *bandstop at ω≈0 plus a slope removal*, and a genuine response near ω≈0 is
  collateral damage.

### Numerical confirmation (induced sine, `s(0)=0`)

Synthetic induced response `s(t) = A sin(2π f₀ t)` (sine → `s(0)=0` exactly,
matching the response-function IC), `A=3`, Hann window, `zero_pad=1`:

| signal | `none` | `initial` | `mean` | `detrend` |
|---|---|---|---|---|
| **pure induced** | peak 3.000 @ f₀, DC≈2e-7 | 3.000, 2e-7 | 3.000, 2e-7 | 3.000, 2e-7 |

For a *correctly-prepared induced response, all four methods are identical*
— and all four give the right answer. This is the key result: when the
signal already obeys `s(0)=0`, initial-value subtraction is a **no-op**, so
it is the *minimally invasive* choice — it can never delete physics because
it only ever removes a single constant `s(0)` that is *already* the right
zero. Peak frequency is unshifted (spread `0.0` across methods) and peak
amplitude matches `A` to `3e-7`.

**Inference:** for δ-kick optical-response dipole/current, initial-value
subtraction is the physically principled default — it *is* the textbook
recipe `d(t)−d(0)`. In exact arithmetic with a perfectly GS-prepared system
it coincides with mean and detrend on the pure tone, so the choice only
*matters* once a numerical baseline or drift contaminates the signal — which
is Q2.

---

## Q2 DC leakage: quasi-stationary vs genuine drift

### The leakage mechanism

A nonzero baseline `b` over a finite record of length `T` is a *windowed
constant*. Its transform is the window's spectral kernel `W(ω)` centred at
ω=0, with sidelobes of width ~`1/T` and Hann sidelobe roll-off. So a residual
constant does **not** live only in bin 0 — it spills into bins `±1, ±2, …`,
i.e. exactly the lowest physical frequencies. A residual **slope** is worse:
its transform is `~1/(iω)`-like, with a `1/ω` tail that decays slowly across
many low bins.

### (i) Quasi-stationary signal (induced tone + constant offset)

`s(t) = A sin(2π f₀ t) + c`, `c=5`:

| method | peak | DC bin |
|---|---|---|
| `none` | **0.010** (wrong!) amp 5.01 | 5.0 |
| `initial` | 0.500 (f₀) amp 3.000 | 1.9e-7 |
| `mean` | 0.500 amp 3.000 | 1.2e-6 |
| `detrend` | 0.500 amp 3.000 | 1.2e-6 |

With **no** subtraction the constant dominates and the reported "peak" is the
lowest bin — a pure artefact. All three subtractions remove the constant
perfectly and recover f₀ at the right amplitude. For a constant offset,
**initial = mean = detrend** are equally good (all remove a constant);
initial-value is the cleanest because here `s(0)=c` exactly so it cancels to
~1e-7.

### (ii) Genuine linear drift (induced tone + slope)

`s(t) = A sin(2π f₀ t) + g·t`, `g=0.01`:

| method | peak | DC bin (residual leakage) |
|---|---|---|
| `none` | 0.500 amp 3.000 | **1.20** |
| `initial` | 0.500 amp 3.000 | **0.4995** |
| `mean` | 0.500 amp 3.000 | **1.9e-7** |
| `detrend` | 0.500 amp 3.000 | **1.9e-7** |

Here the methods **diverge sharply**:

- `initial` subtracts only `s(0)=0` (the slope contributes 0 at t=0) and so
  leaves essentially the *entire* slope in place → large DC leakage (0.50).
  Initial-value subtraction **cannot** suppress a drift.
- `mean` removes the *mean* of the ramp, halving its excursion and — for a
  symmetric ramp over the window — happens to zero the DC bin, but it leaves
  the *slope* (a sawtooth-like residual) which still leaks into low bins via
  the window edges (here the DC bin is clean because the ramp is centred, but
  the low non-zero bins are not — mean is a partial fix only).
- `detrend` removes the slope entirely → DC leakage ~1e-7. **Only detrend
  fully suppresses a genuine linear drift.**

(The "both" case — offset *and* drift — gives the same ordering: only
`detrend` cleans it; `none` peaks at the wrong bin, `initial` leaves the
slope.)

### When linear detrend is right vs when it destroys physics

- **Right tool when** the slow trend is *known to be an artefact* — numerical
  energy non-conservation, slow norm leakage, a CAP-induced monotonic drain —
  AND the physics of interest is genuinely oscillatory at `f > 1/T`. Then the
  slope is noise and removing it is correct.
- **Destroys physics when** the signal genuinely contains low-frequency
  response with period `≳ T` (a slow collective mode, a near-DC induced
  polarization, a build-up that *is* the physics). Detrend will fit that slow
  physics as a "trend" and delete it. A response function legitimately has
  weight approaching ω=0; aggressively flattening ω≈0 is not free.

**Inference:** detrend is the most *robust against drift* but the most
*aggressive against low-ω physics*. Initial-value is the most *faithful to
the response IC* but *defenceless against drift*. Mean sits between and has no
clean physical interpretation in either limit.

---

## Q3 energy spectrum baseline

`E(t)` is categorically different from dipole/current. The **total-energy
drift is a numerical-conservation artefact**, not a physical response: an
exactly-conserving propagator gives `E(t)=const` (the system is autonomous
after the kick/injection — no external work in the linear-response window). Any
slope in `E(t)` is the *integrator's* error (finite `dt`, incomplete SCF,
basis incompleteness), and the *physical* content of an energy spectrum is
the **oscillation about the conserved value** — e.g. coherent breathing of
the energy components, or the energy signature of an excited mode.

Therefore:

- **Initial-value `E(t)−E(0)`** is a *poor* energy baseline: `E(0)` is a
  single noisy sample, and (per Q2-ii) it cannot remove the numerical drift,
  which is precisely the artefact one must kill for the energy spectrum to be
  meaningful. It leaves the drift's `1/ω` leakage all over the low bins.
- **Mean `E(t)−⟨E⟩`** removes the DC level robustly (a many-sample average,
  not one noisy point) and, for a small symmetric drift, suppresses the DC
  bin. It is a defensible default *if* the drift is negligible after the
  transient cut.
- **Linear detrend** is the most defensible when there is a *measurable*
  monotonic energy drift, because that drift is by construction an artefact —
  there is no physical slow energy ramp in an autonomous post-kick system, so
  removing the slope cannot delete physics. This is the one observable where
  the "detrend deletes physics" worry of Q2 essentially does not apply.

**Inference:** for the **energy spectrum**, `detrend` is the most defensible
default (its only risk — deleting genuine slow physics — is absent because
there is no genuine slow energy ramp), with `mean` an acceptable alternative
when the post-transient drift is already `< ~1 mHa` (the project's energy-
conservation bar, `.claude/rules/jellium-base-run-spec.md`). `initial` is the
weakest choice for energy. This *differs* from the current/dipole conclusion
and is the reason the recommendation is **per-observable-type**.

---

## Q4 WP-scattering specifics + finite-amplitude caveat

In a WP run the **injection at t=0 is the stimulus**; everything after is
response. Read through a linear-response lens this makes the WP run formally
analogous to the δ-kick: the perturbation is localized at `t=0`, and `s(t) −
s(0)` is the induced dynamical response with the correct IC `Δs(0)=0`. On that
reading:

- **Initial-value subtraction** is the most consistent baseline — it matches
  the response-function IC and does not inject a time-averaged quantity.
- **Mean subtraction** *biases the baseline by mixing in the stimulus*: the
  mean `⟨s⟩_t` is taken over a window that *includes* the strong injection
  transient and the collision, so `⟨s⟩` is dominated by the stimulus, not by
  a meaningful "rest" level. Subtracting it shifts the response by a stimulus-
  contaminated constant.

### The finite-amplitude caveat (important)

A WP run is **not** a weak perturbation. The injected WP is a finite-
amplitude object (a real added orbital / density), not an infinitesimal kick.
Linear-response intuition is therefore only **approximate**:

- The "response" is not guaranteed `∝` stimulus, so `χ(t)` reasoning is
  heuristic. `α(ω)=d̃(ω)/κ` has no clean analogue; one reads *peak positions*
  (excitation energies of the bath that the WP couples to) more safely than
  *absolute amplitudes*.
- The injection transient (first ~1–2 plasmon periods of shake-up) is a
  genuine non-linear stimulus artefact and is **already excluded** by the
  §13.6 transient cut (`t_start_au`, default 5.0 a.u. for jellium WP_EKIN ≤ 5
  eV; Correa 2018 Sec. 6). The baseline question is about *what remains after*
  that cut.
- After the transient cut, the WP is mid-flight; there can be a **genuine slow
  drift** in `J_z`/`d_z` as the packet physically translates across the box
  (a real, low-frequency physical motion), *and* a numerical drift. Here Q2-ii
  bites: `initial` will leave that drift in (good if the drift is physics you
  want; bad if it swamps the inelastic-scattering oscillations you are after),
  while `detrend` removes it (good for isolating the scattering spectrum; bad
  if the slow translation *is* the signal).

**Inference:** for WP current/dipole spectra, initial-value subtraction is the
most *interpretation-consistent* baseline (`J_z − J_z(0)` is more meaningful
than `J_z − ⟨J_z⟩`, exactly as the project TODO argues), but because WP runs
carry a real packet-translation drift, the **detrended** variant remains the
most useful for *isolating the oscillatory inelastic-scattering content* — and
the project's existing practice of emitting all three with a `_compare`
overlay (peaks surviving in `detrended` are the most trustworthy) is sound.
Do **not** trust absolute amplitudes from a WP run as if it were linear
response (finite-amplitude caveat).

---

## Recommendation (canonical default per observable type; kernel API)

Grounded in Q1–Q4 and the project's prior reasoning
(`docs/todo_later.md` "Choice of drift-removal method", §11):

### Per-observable-type canonical default

| Observable | Canonical default | Rationale | Keep available |
|---|---|---|---|
| **dipole `d`, current `J`** (δ-kick **and** WP) | **`initial`** (`s−s(0)`) | Matches the response-function IC `Δs(0)=0`; is the literal Yabana–Bertsch recipe; injects no window-dependent average; minimally invasive (no-op on a clean induced signal). | `detrend` as the drift-robust comparison variant; `mean` for diagnostics. |
| **energy `E`** | **`detrend`** | The only slow component is the numerical-conservation artefact; there is no genuine slow energy ramp, so detrend's "deletes-physics" risk is absent. | `mean` when post-transient drift `< ~1 mHa`. |

This split is deliberate and is the main finding: **`initial` for the
response-like observables (dipole/current), `detrend` for the artefact-like
observable (energy).** Mean-subtraction has no regime where it is uniquely
best and is recommended only as a diagnostic.

Caveats to attach wherever a spectrum is reported:
- Always apply the §13.6 transient cut *before* baseline subtraction.
- For WP runs, flag absolute amplitudes as not-linear-response (Q4).
- For dipole/current, when a clearly-artefactual drift survives the transient
  cut, fall back to `detrend` and *say so in the caption* — initial-value
  cannot remove a slope (Q2-ii).
- Baseline choice does **not** move peak frequency or (for an on-bin tone)
  peak amplitude — verified to `3e-7` (Q1). It only changes the ω≈0 region.
  So peak *positions* are robust to the choice; only near-DC features depend
  on it.

### Kernel API

**Yes — expose `subtract={'initial','mean','detrend','none'}`.** The kernel
should not hard-code one policy because the correct choice is
observable-dependent (above) and run-dependent (drift present or not).
Recommended signature behaviour:

```python
def fft(values, dt, *, subtract='detrend', window='hann', zero_pad=4,
        t_start_au=None): ...
```

- Keep `'none'` for callers that pre-process themselves and for the Parseval
  / calibration unit tests (which must see the raw windowed signal).
- **Default:** keep the current `'detrend'` as the *library* default (safest
  blanket choice — it never lets a baseline masquerade as a peak), BUT have
  the **postprocess layer** (`observables.py::_build_variants`) set the
  *per-column canonical* variant per the table above and mark it canonical
  (emit it as `spectrum_<col>.png` alongside the variant grid), i.e. `initial`
  for dipole/current columns, `detrend` for energy columns. This keeps the
  generic kernel conservative while making the *scientifically reported*
  spectrum the physically-correct one per observable.

---

## Portable test (numpy pseudocode + expected values + tolerances)

Pure-numpy, no project imports; machine-independent tolerances. Tests that
the chosen subtraction (a) removes the baseline, (b) leaves the f₀ peak
**unshifted** with **correct amplitude**, (c) **suppresses the ω≈0 bin**.

```python
import numpy as np

def one_sided_amplitude(values, dt, win, zero_pad=1):
    """Coherent-gain-normalized one-sided amplitude spectrum.
    (Normalization per docs/validation/fft-normalization-validation.md.)"""
    n = len(values)
    sig = values * win
    n_pad = n * zero_pad
    padded = np.zeros(n_pad); padded[:n] = sig
    raw = np.fft.rfft(padded)
    freq = np.fft.rfftfreq(n_pad, d=dt)
    amp = np.abs(raw) / win.sum()       # coherent gain (= n for boxcar)
    amp[1:-1] *= 2.0                    # one-sided interior doubling
    return freq, amp

def subtract(s, t, mode):
    if mode == 'none':    return s.copy()
    if mode == 'initial': return s - s[0]
    if mode == 'mean':    return s - s.mean()
    if mode == 'detrend':
        p = np.polyfit(t, s, 1); return s - np.polyval(p, t)
    raise ValueError(mode)

def test_drift_removal():
    n, dt, A = 1000, 0.1, 3.0
    f0 = 50 / (n * dt)                  # 50 integer cycles -> exactly on a bin
    t  = np.arange(n) * dt
    df = 1.0 / (n * dt)                 # bin spacing (= 0.01 here)
    win = np.hanning(n)

    # Induced response: SINE so s(0)=0 (correct response-function IC).
    tone = A * np.sin(2 * np.pi * f0 * t)
    cases = {
        'constant_offset':   tone + 5.0,        # (a) constant
        'linear_drift':      tone + 0.01 * t,   # (b) drift
        'both':              tone + 5.0 + 0.01 * t,
    }
    # which subtraction is expected to fully clean each case:
    full_clean = {
        'constant_offset':   ['initial', 'mean', 'detrend'],
        'linear_drift':      ['detrend'],            # only detrend kills a slope
        'both':              ['detrend'],
    }
    for cname, sig in cases.items():
        for mode in ['initial', 'mean', 'detrend']:
            s = subtract(sig, t, mode)
            freq, amp = one_sided_amplitude(s, dt, win, zero_pad=1)
            ipk = np.argmax(amp)
            # (b) peak unshifted: within half a bin of f0
            assert abs(freq[ipk] - f0) <= df / 2 + 1e-12, (cname, mode, freq[ipk])
            # (b) amplitude correct (Hann, coherent-gain normalized) to ~2%
            assert np.isclose(amp[ipk], A, rtol=2e-2, atol=1e-3), (cname, mode, amp[ipk])
            # (c) DC suppression for the methods expected to fully clean
            if mode in full_clean[cname]:
                assert amp[0] < 1e-4, (cname, mode, amp[0])

    # initial-value is a NO-OP on a clean induced signal (s(0)=0): minimally invasive
    f1, a1 = one_sided_amplitude(subtract(tone, t, 'initial'), dt, win)
    f0_, a0 = one_sided_amplitude(subtract(tone, t, 'none'),    dt, win)
    assert np.allclose(a1, a0, atol=1e-10)

    # initial-value CANNOT remove a slope (documents Q2-ii; not an error):
    s = subtract(cases['linear_drift'], t, 'initial')
    _, amp = one_sided_amplitude(s, dt, win)
    assert amp[0] > 1e-1, amp[0]        # large residual DC leakage, as expected

    print("drift-removal tests passed.")

test_drift_removal()
```

### Expected values and tolerances

| Check | Expected | Tolerance (machine-independent) |
|---|---|---|
| peak frequency, every (case, method) | `f0 = 0.5` | `|Δf| ≤ Δf/2 = 1/(2·n·dt)` (half a bin) |
| peak amplitude, every (case, method) | `A = 3.0` | `rtol = 2e-2, atol = 1e-3` (head-room for finite-record leakage; portable across numpy/BLAS) |
| DC bin, constant offset, {initial,mean,detrend} | `≈ 0` | `< 1e-4` (verified ≈ 2e-7–1e-6) |
| DC bin, linear drift / both, **detrend only** | `≈ 0` | `< 1e-4` (verified ≈ 2e-7) |
| DC bin, linear drift, **initial** | large residual | `> 1e-1` (verified ≈ 0.50) — documents that initial cannot kill a slope |
| `initial` on clean induced tone | identical to `none` | `atol = 1e-10` — proves minimal invasiveness |

**Reproduced numerically** (`/local/data/public/skcb2/tddft/venv/bin/python3`,
Hann, `zero_pad=1`, `A=3`, `f0=0.5`):

- pure induced sine: all of {none, initial, mean, detrend} → peak 3.00000 @
  0.5000, DC ≈ 1.9e-7; peak-freq spread `0.0`, peak-amp spread `2.8e-7`.
- induced + constant offset 5: `none` peaks at wrong bin 0.0100 (amp 5.01);
  initial/mean/detrend → 3.00000 @ 0.5, DC ≈ 2e-7–1e-6.
- induced + linear drift: initial DC = **0.4995** (slope survives), mean DC =
  1.9e-7, detrend DC = 1.9e-7; all peak 3.00000 @ 0.5.
- induced + both: only `detrend` fully cleans; `none` peaks at wrong bin.

These reproduce the tables in Q1–Q2 and back every recommendation above.

---

## Sources cited

- **K. Yabana & G. F. Bertsch**, "Time-dependent local-density approximation
  in real time," *Phys. Rev. B* **54**, 4484 (1996). The real-time TDDFT
  optical-response recipe: weak impulsive kick → propagate induced dipole
  `d(t)` → Fourier transform → dipole strength function / absorption. The
  transformed object is the *induced* dipole `d(t)−d(0)`, fixing the response
  initial condition `Δd(0)=0`.
- **M. A. L. Marques & E. K. U. Gross** (and the broader TDDFT community
  reviews, *Time-Dependent Density Functional Theory*, Lect. Notes Phys.
  **706**, Springer 2006). Linear-response TDDFT, the δ-kick spectrum recipe,
  `α(ω)=d̃(ω)/κ`, `S(ω)∝ω·Im α(ω)`.
- **A. L. Fetter & J. D. Walecka**, *Quantum Theory of Many-Particle
  Systems*; **R. Kubo**, linear-response theory. Causal response functions
  `χ(t−t')` start from zero — the IC that justifies initial-value subtraction
  for induced quantities.
- **F. J. Harris**, "On the use of windows for harmonic analysis with the
  DFT," *Proc. IEEE* **66**(1), 51–83 (1978). Window spectral kernels and
  sidelobe behaviour — the mechanism by which a residual baseline leaks from
  ω=0 into the low non-zero bins.
- **A. V. Oppenheim & R. W. Schafer**, *Discrete-Time Signal Processing*.
  DFT, leakage, and the fact that a windowed constant/ramp has a `W(ω)` /
  `1/ω`-like low-frequency footprint; detrending as a high-pass-at-DC
  operation.
- **A. A. Correa**, "Calculating electronic stopping power in materials from
  first principles," *Comput. Mater. Sci.* **150**, 291 (2018), Sec. 6 —
  definition of the post-stimulus transient region (cited in
  `docs/observables_reference.md` §13.6 for the `t_start_au` cut applied
  *before* baseline subtraction).
- Project prior reasoning (not external authority, recorded for continuity):
  `docs/todo_later.md` "Choice of drift-removal method for spectra";
  `docs/observables_reference.md` §11 (variant table) and §13.6 (transient
  cut); `docs/validation/fft-normalization-validation.md` (the assumed-fixed
  amplitude-normalization convention).

All numerical results above were produced with
`/local/data/public/skcb2/tddft/venv/bin/python3` (numpy), not taken from any
source. Web search was not required: every quantitative claim is reproduced
locally, and the TDDFT recipe and DSP leakage facts are standard textbook
results from the cited references.

---

## Verdict (LEAVE BLANK — user fills)

- canonical default for **current/dipole** spectra: ____
- canonical default for **energy** spectra: ____
- expose `subtract={'initial','mean','detrend','none'}` param?  Y/N ____
- per-observable canonical variant emitted as the headline `spectrum_<col>.png`?  Y/N ____
