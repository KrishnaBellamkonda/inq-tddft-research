# inqview findings (bugs + method decisions)

Mirror of `inqkit-errors.md` for the Python side. Bugs are recorded but
**fixed only after** the characterization suite guards behaviour (ADR 0005) —
each fix is a separate red→green step. `IV-E` = inqview error/finding.

## Bugs (deferred — characterization first)

### IV-E01 — loss function FFTs real part only  [CONFIRMED, to fix]
`postprocess/density_fourier.py:182` does `np.fft.fft(sig.real, ...)` on the
complex phasor `n_q(t)`, contradicting its own comment (lines 179–181). Folds
±ω together, halves amplitude, loses propagation direction. Confirmed by
formula-validation agent. The portable test (undamped-plasmon real-part
variant) is designed to catch it. NOTE: report1 `make_fig_loss_function.py`
reads the complex CSV columns, so the published figure was likely unaffected —
verify before reuse. Status: **confirmed; fix during migration (separate
red→green)**.

### IV-E02 — quantity is power, not "loss function"  [RELABEL, not a code bug]
`density_fourier.py` returns/plots `|n_q(ω)|` (and the memory called
`|n_q|²/q²` "the loss function"). Per validation, the squared/weighted quantity
is a plasmon **peak-locator** (|χ|²-like), NOT −Im[1/ε]. Resolution is a
**relabel** (→ `PlasmonSpectrum`, see IV-M04), not a numeric fix: compute
`|n_q(ω)|²` honestly, test the peak LOCATION, drop the "loss function" label.
Memory `reference_loss_function_method` updated with the caveat.

## Method decisions (locked via grilling)

### IV-E03 — fourier.py window coherent-gain not corrected  [CONFIRMED, fix known]
`fourier.py` divides by `n` but never by `Σ win`, so every windowed amplitude
is under-reported by `mean(win)` (Hann unit tone → ≈0.5, boxcar → ≈1.0).
Confirmed by signal-validation agent (`docs/validation/fft-normalization-validation.md`),
which numerically verified: boxcar `|rfft|/n`+interior×2 is correct (Q1), zero-pad
by original `n` is correct (Q3), Parseval holds (Q4). **Fix:** `amplitude =
|raw| / win.sum()` (Harris-1978 coherent gain) — backward-compatible (boxcar
`Σwin=n`). Use `Σwin²`/ENBW only for a calibrated PSD (separate convention).
Skews QKE-vs-INQ amplitude comparison (fourier.py:8) until fixed. Portable test
(sinusoid A,f0 → peak@f0 + amplitude==A for boxcar AND Hann + Parseval) in the
dossier. Status: **confirmed; fix + test during migration**.
fourier.py:12 "used for" → observable spectra (energy/dipole/current →
plasmon/absorption), `plot_spectrum`. Core `transform()` → `analysis`.

### IV-M10 — validation strategy: code-review at end (2026-06-10)
User chose: write all tests, then ONE user-triggered `/code-review ultra` pass
on the whole suite at the end (the assistant cannot launch it — billed/
user-triggered). Per-test formula/test-validation agents NOT run going forward
(though 2 were run this session: loss-fn, FFT). **Mitigation (mandatory
practice):** every expected value is still **derived analytically up front**,
never captured from current code output — this is what prevents a circular
green test that the code-quality review would not catch.

### IV-M11 — free-WP simulation tests: separate per library (2026-06-10)
- **inqkit** — ✅ BUILT + PASSED (ctest 3.5 s, first run, 0 failures, 2026-06-10)
  `test_free_wp_engine.cpp`: non-interacting theory (kinetic-only H), inject
  Gaussian WP (σ₀, k₀), propagate, assert vs ANALYTIC free-particle evolution:
  `σ_r(t)=σ₀√(1+(t/τ)²)` (τ=mσ₀²/ℏ), `⟨p⟩(t)=k₀` conserved, centroid
  `z(t)=z₀+(k₀/m)t` ballistic, norm=1, energy conserved. Exercises
  wavepacket-inject + wp_real_space_stats + wp_momentum_stats + center_of_density
  end-to-end. (Extends the "complete" inqkit suite with its first integration
  test.)
- **inqview** — build its free-space-WP fixture **independently** (not the same
  run output); post-process → assert the same analytic laws. User chose
  isolation over a shared single run (duplicated but decoupled).

### IV-M12 — FFT drift-removal/subtraction baseline: RESOLVED (2026-06-10)
todo.txt #1/#2 + todo_later "Choice of drift-removal method", settled by the
signal+TDDFT validation agent (`docs/validation/fft-drift-removal-validation.md`)
and user sign-off. **Locked:** `FourierTransform` exposes
`subtract={'initial','mean','detrend','none'}`, library default `'detrend'`
(conservative); the **postprocess layer sets the canonical per column —
`'initial'` (s−s(0)) for dipole/current** (Yabana–Bertsch induced-response IC),
**`'detrend'` for energy** (its only slow term is the numerical-conservation
artefact). `mean` is diagnostic-only (never uniquely best). Fixes the #2 bug:
with no subtraction a constant offset hijacks the "peak" to the lowest bin.
Verified numerically: baseline choice never moves peak freq/amplitude (only
ω≈0); against a genuine linear drift only detrend fully suppresses DC leakage.
Caveat: WP runs are finite-amplitude → absolute amplitudes untrustworthy
regardless. Portable test (tone+offset+drift → peak unshifted, baseline
removed, ω≈0 suppressed) in the dossier. Combine with the IV-E03 coherent-gain
fix in the fourier kernel.

### IV-M09 — features: band-structure deferred (multi-k), gs-projected kept (2026-06-10)
- **Band structure** (`orbitals_per_kpoint.py:30`): wanted for **multi-k-point
  runs** (QuantumKickExtension/QBall crystalline Li/diamond/Al), NOT the Γ-only
  jellium. **Deferred to a FUTURE TODO** — a true ε vs k-path band structure
  needs Brillouin-zone k-sampling the jellium runs don't have. Build when
  multi-k runs are post-processed. (Existing `bands_summary.png` ε-vs-k-index
  suffices for now.)
- **gs_projected_occupations** (`:55`): KEEP (not cut — catalogue's "remove"
  was wrong). It is "the cleanest TDDFT definition of excitation into orbital
  i" (`n_i^GS(t)=Σ_j f_j(0)|⟨ψ_i^GS|ψ_j(t)⟩|²`). Ensure the heatmaps are
  actually produced; make it a jellium-WP extra. Test: **t=0 identity**
  `n_i^GS(0)=f_i(0)` (overlap = identity at t=0).

### IV-M08 — pipeline stays strictly sequential (2026-06-10)
Dispatcher remains sequential. Document that phases are independent except
`summary` (builds `run_summary.txt`, read by others for `wp_state_index`/`dt`).
Parallelism (process-pool over the now-pure `analysis` phases; viz not
thread-safe) is noted as FUTURE work, not built. User chose simplest.

### IV-M07 — energy analysis: functional-component flow (2026-06-10)
New `analysis/energy_components.py` becomes THE primary energy analysis,
sidestepping band-sum double-counting. Uses `observables.csv` (`energy_total,
energy_kinetic, energy_hartree, energy_xc`) + implied `E_ext = E_total −
(E_kin+E_H+E_xc)` (residual = electron–ion). `compute(observables) →
EnergyComponents(t, E_kin, E_H, E_xc, E_ext, E_total, ΔE_* from t0)`. Answers
"which store gained the energy". Visualisation (3): initial-vs-final grouped
**bar** breakdown; `ΔE_component(t)` **time-series** lines; a **GIF** animating
the breakdown. Test invariants: `E_kin+E_H+E_xc+E_ext == E_total` (exact),
`Σ ΔE_component == ΔE_total`, `ΔE_total` drift within tolerance. Old band-sum
ledger (`energy_balance.py`, `bath_energy.py`) KEPT as **secondary, caveated**
(band sum `Σf_iε_i` ≠ total energy; "Unaccounted" conflates excitation with
Hartree/xc double-counting — relabel honestly), jellium-WP extra. **FUTURE
(deferred):** orbital-level WP×component cross-decomposition — needs saved
wavefunctions; kinetic per-orbital splittable, Hartree/xc NOT cleanly
attributable. User: "interesting but presume not too meaningful; future."

### IV-M05 — WP-integrity metric: momentum-KL + real-space spread/IPR (2026-06-10)
New `analysis/wp_integrity.py` → `WPIntegrity(t, kl_mom, sigma_r, ipr)`:
- `kl_mom` = existing momentum `KL(P_t‖P_0)` (translation-invariant → genuine
  momentum-scattering/preservation measure);
- `sigma_r(t)` = real-space WP spread (reuse `wp_real_space_stats`);
- `ipr` = `(∫ρ²)/(∫ρ)²` localisation.
Resolves kl_divergence.py:49. Test: **free-space WP** → `kl_mom≈0` (no
scattering), `sigma_r(t)=σ₀√(1+(t/τ)²)` analytic Gaussian spreading, known IPR
decay (ties to the free-space-WP integration fixture, ADR 0005).
Sub-resolutions (kl_divergence): line 40 `_common` import → splits at package
boundary (path helpers→pipeline, plot helpers→visualisation); line 46 KL is
ALREADY a time series (negligible cost; optional frame-to-frame `KL(P_t‖P_{t−Δ})`
drift-rate variant); line 53 contour viz → the `(k,t)` momentum carpet is the
existing `grid` array straight to `contourf` (visualisation feature).

### IV-M04 — loss function renamed to PlasmonSpectrum, peak-locator semantics (2026-06-10)
Resolves IV-E02 + the formula-validation finding. The kernel
(`analysis/plasmon_spectrum.py`) computes `|n_q(ω)|²` (axial + 3d_binned,
IV-M01) and returns `PlasmonSpectrum(q, omega, power, peak_omega(q))` — NOT
labelled a loss function. The usable output is the peak position ω_p(q).
Portable test: undamped-plasmon phasor `n_q(t)=A·e^{−iω_p t}` → δ-peak at ω_p
(exact location) + 1/q² scaling (`rtol 1e-6`), plus a real-part variant that
catches IV-E01. Reduced-system numpy pseudocode + expected values are in
`docs/validation/loss-function-formula-validation.md`. Proper −Im[1/ε] via
real-time χ from a weak kick was considered and deferred (heavier; needs the
drive spectrum). Memory `reference_loss_function_method` updated.

### IV-M03 — minimum observable set: one global core + extras (2026-06-10)
`pipeline` defines a single global `PHASES_MINIMUM = (summary, observables,
density)` that every run produces; run-type essentials (jellium wake/bath/
trajectory, coronene LEED screens) are added per call as opt-in extras. User
chose this over per-system presets — simpler, flat. Note (flagged, not
actioned): jellium wake/bath and LEED screens are NOT in the global core, so
callers must add them explicitly for those run types.

### IV-M02 — COD recomputed in Python, not reused from inqkit CSV (2026-06-10)
`inqview.analysis.center_of_density.compute(density)` computes COD from saved
VTI with the **correct node convention** (`ix·dx`, matching
`wp_real_space_stats`), returning `COD(x,y,z)` and a `CODComparison(wp, total,
bath)`. It does NOT consume inqkit's `center_of_density.csv` (which carries the
deferred E04 `(ix+0.5)·dx` half-cell offset). A test asserting
`python_cod − inqkit_csv_cod == dx/2` turns E04 into a documented cross-check
instead of an inherited bug. inqview is post-processing, so this needs no
production change. Ties: inqkit E04, memory `reference_canonical_bath_density`.

### IV-M01 — loss-function q-sampling: BOTH modes (2026-06-10)
The loss-function analysis kernel supports `mode='axial'` (q∥z,
q=(0,0,2πm/L_z) — longitudinal/plasmon, matches current data) and
`mode='3d_binned'` (all (kx,ky,kz) binned by |q|, direction-averaged +
isotropy check). Two code paths, two golden sets. Rationale: jellium is
isotropic so 3d_binned gives more |q| samples + a consistency check, while
axial is the simpler longitudinal-only estimator; the user wants both
exposed. Formula `L(q,ω)=|n_q(ω)|²/q²` under independent validation before
the test locks (formula-validation agent dispatched 2026-06-10).
