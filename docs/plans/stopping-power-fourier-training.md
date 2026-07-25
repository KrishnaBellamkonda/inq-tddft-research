# Plan: hands-on training — stopping power & Fourier-transform analysis

Source prompt: `docs/prompts/check_logic/check_stopping_power_calculation.md`.
Goal: the user regains hands-on control of two formula-bearing derived
observables they feel were "outsourced" — electronic **stopping power S(v)** from
a classical-projectile track, and **Fourier-transform analysis** of TDDFT
time-series — building each into a *replicable, deterministic* workflow (and a
skill if warranted). Notebook-driven, tight feedback loop.

## Scope — two tasks (do task 1 fully first)
1. **Stopping power extraction** (this notebook first).
2. **Fourier-transform analysis** (second notebook; `fourier.py` already carries
   TODOs doubting its windowing/detrend/convenience methods — that is the audit
   target).

## Method shape (agreed in grilling)
For each task: **(b) understand-in-context, built up from first principles so the
existing kernel is a destination not a black box → critically stress-test the
method → rebuild** a clean deterministic version. Reaching the existing code via
an independent derivation keeps the judgement honest (mirrors the CONTEXT.md
formula/test-validation independence ethos).

## Deliverables / placement
- `docs/validation/stopping-power-extraction/stopping_power_extraction.ipynb`
- `docs/validation/stopping-power-extraction/fourier_analysis.ipynb` (task 2)
- `docs/validation/stopping-power-extraction/README.md` — ties them; records the
  locked method.
- Executed against the **`inqview-venv`** kernel; outputs embedded; user views /
  re-runs in VSCode. Figures as PNG, canonical theme.

## Task-1 concrete example
- First example: **`v3p0`** run (`run_sv_sigma0p5/results/v3p0/`), v₀=3.0, nearly
  constant v (clean single slope). Stress-test runs held in reserve: `v0p8`
  (Barkas sweep into noise), `v0p6` (sub-v_F, noisiest).
- r_s=5.69 jellium, m=m_e, free Ehrenfest, σ=0.5 erf electron. Track columns:
  `step,time_au,x,y,z,vx,vy,vz` (motion along z).

## Task-1 first-principles ladder (section-by-section; revise after each)
1. Load & look — `v(t)`, `z(t)`, `KE_proj(s)`; also `observables.csv` energies.
2. Define S = −dE_proj/dx = −dKE/ds; one global linear fit → single number
   (Ha/Bohr and eV/Å).
3. Transient — motivate discarding `transient_bohr`; re-fit.
4. Cross-check — S′ = +dE_electrons/ds same window; E_total drift = integrator
   health.
5. Why bin by v(t) — v≈const here so one slope suffices; decelerating runs need a
   *local* slope (bridge to the kernel's design).
6. Arrive at (b) — call `inqview.analysis.stopping_extract.stopping_vs_v`; show it
   reproduces the hand number. No black box.
Then: critique phase (window size, transient length, finite-diff vs fit, KE-vs-
E_elec disagreement, low-v breakdown on v0p8) → rebuild.

## Status
- [x] Grilling: scope, method shape, location, example, ladder, loop granularity.
- [x] Section 1 — load & look (v3p0).
- [x] Sections 2–6 + critique (§7a–e) + clean rebuild (§8) BUILT & EXECUTED.
  - Deterministic builder: `build_stopping_notebook.py` (preserves §1, appends
    §2–8, executes via `inqview-venv` kernel). Notebook 36 cells, **0 errors**.
  - Numbers (grounded, v3p0): global S=0.0074 Ha/Bohr=0.38 eV/Å (no cut);
    transient-cut(3.0) S=0.0082=0.42 eV/Å; kernel mean 0.008186; rebuild 0.008171
    — three routes agree <1 %. Electronic cross-check 8.5 % (5 coarse pts);
    integrator drift 0.02 % of KE0.
  - **Honest stress-test findings (overturned my first-draft prose):**
    (1) §3 transient has NO plateau — S keeps rising with cut (0.42→0.48 over
    3→8 Bohr) → the ONE real systematic (~15 %); kernel default 3.0 likely
    under-estimates. (2) §7c finite-diff ≈ windowed (roughness ratio ≈1) — tracks
    are smooth, derivative method barely matters. (3) §7e NO low-v breakdown
    (zero sign-flips, roughness ~1e-3) — instead a clean **Bragg peak**: S rises
    to ~2.6 eV/Å near v~v_F (0.4–0.8 a.u.) and turns over below v_F (v0p6).
- [ ] **PENDING: user critique of the stopping notebook** (interactive loop).
- [ ] README locking the method (after critique).
- [x] Task 2 — Fourier-audit notebook BUILT & EXECUTED.
  - Deterministic builder: `build_fourier_notebook.py` → `fourier_analysis.ipynb`
    (31 cells, **0 errors**). Three-stage known-answer audit.
  - **All known answers reproduced** (verified before building):
    A1 coherent gain → 1.0000 for every window; A3 zero-pad peak-invariant;
    A4 Test-1 peak at ω_p + `L·q²` q-independent (1/q² exact); A5 Test-2 BUG-A
    real-part halving ratio **0.500**; B dipole_x plasmon **6.446 eV** (journal
    6.480, within the 0.276 eV bin) and **baseline-invariant**; C E15 n_q m=1
    complex FFT **3.533 eV** (stored 3.53 ≈ ω_p).
  - **Audit verdict-evidence:** `fourier.py` window/baseline/zero-pad machinery
    is sound and does NOT skew the QKE plasmon position — the real pitfall is
    naive global-argmax on a DC-dominated spectrum (→ always search a physical
    band). `density_fourier.py` BUG-A (line 182 `.real`) and BUG-B (line 183
    plots `|FFT|` not `|FFT|²/q²`) both **confirmed on real data**; corrected
    `loss_function` drop-in provided in C4 (apply via code-test in Task 3).
  - Dossier verdict lines left BLANK for the user (verification-user-owns-verdict).
- [ ] **PENDING: user critique of the Fourier notebook** (interactive loop).
- [ ] Task 3 — fill dossier verdicts + apply density_fourier BUG-A/B fixes (code-test).

## 2026-06-24 — METHOD CORRECTION (user redirection, locked-in-progress)

The user reviewed the whole stopping notebook and **overruled its primary method.**
The canonical electronic stopping power for a **classical projectile** is:

1. **PRIMARY — total-system-energy regression.** `dE_total(x) = E_total(t) −
   E_total(0)` plotted vs traversal distance `x`; discard the transient; **linear-fit
   the steady-state region**; the **slope = S**, with **error bars = `linregress`
   stderr**. Grounded: Correa 2018 Eq. (10), `S = ⟨dE/dt⟩/v` (≡ `dE/dx`); transient
   exclusion = Correa Fig. 8 (`docs/sources/correa-2018-electronic-stopping-power.md`).
   `E_total` = `observables.csv:energy_total` (electronic total energy; the projectile
   is external classical, so its energy gain = projectile loss).
2. **SANITY (a) — projectile deceleration**, `−dKE/dx` (`KE=½mv²` from the track).
3. **SANITY (b) — ∫F·v dt** cumulative work. CAVEAT (honest): the track records only
   `x,v` (no force column) so `F=m·dv/dt` ⇒ ∫F·v dt ≡ ΔKE analytically — it shares the
   kinetic channel, so it is a *deposition-profile / discretisation* check, NOT a third
   independent physics channel. Flagged to the user.
   Every extraction runs all three; a large primary-vs-sanity deviation is **reported
   to the user** for investigation (do not silently average).

**Transient detection (research done; method to lock).** Correa eyeballs it; no
published auto-detector. Physical scales (r_s=5.69): k_F=v_F=0.337, ω_p=0.128 a.u.,
**τ_p≈49 a.u. (longer than the runs)** ⇒ no plasmon-period floor possible; transient
set by e-h response ~ screening-length/v (λ_TF≈1.53 Bohr). Recommended detector =
**slope-convergence plateau** of `S(x_c)` vs cut, + physical floor + user gates:
transient < 40% of run (else flag run-time suspect), typical < 20%, fallback ~20%.

**Data for the S(v)-vs-Lindhard test.** Classical-Gaussian σ-sweep runs
(`run_classical_n162_L50_sv_sigma0p15/0p25/0p35/3p0` + `run_sv_sigma0p5`), velocities
v0p2..v3p0; `obs_every=3` ⇒ 84–101 energy samples/run (fit-viable, unlike the old
7-row v3p0 anchor). σ in run name = projectile **charge** std (σ_q); σ_WP=√2·σ_q, so
σ_q=0.35 ⇔ σ_WP=0.5 ("σ=0.5 → 0.350"). Lindhard reference: existing
`inqview.analysis.lindhard_elf.stopping_power_point` (point charge) and
`stopping_power_sigma` (Gaussian form factor). The σ-sweep **report script already
implements Method A + B** (`hypotheses/06_sigma_convergence/sigma_sweep_report.py`,
fixed 20% transient) — the new notebook rebuilds this cleanly with the principled
transient + 3 channels + the skill kernel.

- [x] **Built `transient_method_comparison.ipynb`** (17 cells, 0 err, 5 figs;
      builder `build_transient_comparison.py`). Compares M1 fixed-20% vs M2
      slope-plateau agent (user-supplied spec, verbatim), both free-intercept.
      **Result:** ΔE_total(x) is CONVEX (slope rises monotonically, no plateau in
      the 2% band) because wake-formation ~1/ω_p≈8 a.u. ≳ run length (6–20 a.u.).
      → **M2 reports S on only 1/30 runs** (and that one is a degenerate 7-pt
      flat-layout case); M1 reports on all 30, r²≈0.99, M1/LR≈1.05 at v3p0.
      Energy conservation 0.7%; 3 channels agree to 0.5%. Separate obs: M1
      undershoots Lindhard near the Bragg peak (v0.6–0.8: M1/LR≈0.5–0.67).
- [x] **Built `p5_classical_transient_comparison.ipynb`** (13 cells, 0 err, 4 figs;
      builder `build_p5_classical_comparison.py`) on the **localised-jellium Phase-5
      classical slab** run (user-identified: longer, cleaner). r_s≈4 Na-like slab,
      25 Bohr thick, classical σ_pot=0.35, v=2.711, CAP η=−0.5 but **N conserved**
      (234→233.78) so `E_total(t)−E_total(0)` is a clean **bath** signal = 23.3 eV
      (region ΔE_bath=26.5 eV, ΔKE_ion=24.7 eV — all the user's "20–30 eV").
      Reference `S = ΔE_bath/25 = 0.93 eV/Bohr`.
      **Result (different failure mode than the bulk runs):** signal is a localised-
      deposit **sigmoid** → the dominant choice is the **upper bound x_T (slab exit)**,
      not the 20% entry transient. M1 fixed-20% full=0.615 (underestimates 34% via
      post-exit flat); @slab-exit=0.815. M2 = **no_plateau** again (S(x₀) never flat),
      BUT its **endpoint check correctly fired `endpoint_contaminated`** (end-slope
      0.23 ≪ mid 0.87) — it detects the slab exit M1 ignores.
      NB earlier confusion: bulk `b2_classical_E100` is a CAP run that drains 97% of
      the bath → its raw `energy_total` is unusable (peaks +20 Ha); that was the wrong
      run. Phase-5 slab is the right one.
- [x] **ENCODED `stopping-power-extraction` skill** (`.claude/skills/stopping-power-extraction/`,
      self-contained: `SKILL.md` + `stopping_power.py` kernel + portable `_selftest`).
      Locked method (user, 2026-06-24): **window [x0, xT]** (both bounds first-class)
      + run-geometry split — **Method A** slope-fit over [x0,xT] (fixed-20% OR agent
      x0; endpoint detection for xT) for continuous/bulk runs; **Method B**
      `S=[E_total(t_f)−E_total(t0)]/L_z` for localised slabs, with a **mandatory
      convergence gate** on E_total(t_f). Guards: N(t)≈const (CAP-drainage kills
      E_total), ΔE_total≈ΔKE_ion. Self-test passes; reproduces Phase-5 S=0.931 eV/Bohr
      (gate flags marginal not_converged). Campaign frontmatter task-4 → done.
- [ ] Task 5 — encode the `fourier-analysis` skill (the gating deliverable).
- [ ] **PENDING USER REVIEW** of the stopping skill (training-loop critique).
- [ ] (superseded) transient-detection algorithm lock — folded into the decision above.
- [ ] Rewrite `stopping_power_extraction.ipynb`: primary = energy regression; one run
      walked in detail (dE_total(x), fit+errors, transient plateau, 3 channels).
- [ ] New S(v)-vs-Lindhard notebook driven by the (drafted) skill kernel across all
      (σ, v); error bars; deviation flags.

## 2026-06-25 — Latest classical run (qsp_phase2) + periodic-wrap truncation

User redirection: build a p5-style notebook for the **latest classical run**
(`localised_jellium/scripts/qsp_phase2/classical/results/p2_classical`) with the
time series **truncated before the projectile loops the periodic cell and re-enters
the slab**. Run: cell 50×50×70 (z∈[−35,+35] periodic), slab faces ±12.5, CAP |z|∈[25,35]
η=−0.7, launch z₀=−22, v₀=2.711, 2000 steps × dt 0.02 = 40 a.u.

- **Truncation calc (locked with user).** Cut where the projectile, having crossed
  the top wall (z=+35) and wrapped, exits the first (lower) CAP at physical z=−25 ≡
  **unwrapped z=+45** (path 67 Bohr from launch). Quick v₀-const estimate = 24.7 a.u.
  (step ~1236); **empirical track (decelerating) = 31.44 a.u. = step 1572** — the
  chosen cut. Keep steps 0…1572 (158/201 obs rows, 1573/2001 track rows).
- **Built `p2_classical_truncated_stopping.ipynb`** (builder
  `build_p2_classical_truncated.py`, 17 cells, **0 errors**, executed on inqview-venv).
  Imports the **skill kernel** `stopping_power.py` directly → the notebook is a live
  validation of the shipped skill on an unseen run.
- **Headline (user's choice 2026-06-25): skill Method B as-is, trust the flag.**
  - Full (un-cut) Method B = 0.104 Ha/Bohr (70.6 eV deposit) — catastrophic wrap-image
    re-plough. Truncated Method B = **0.0125 Ha/Bohr = 0.644 eV/Å, status NOT_CONVERGED**
    (tail 66% of deposit in last 15%). Reported as a LOWER BOUND.
  - Guards pass: N drained 0.024%; dKE_ion (+0.380 Ha) ≈ dE_electronic (+0.313 Ha).
- **Trainee flag recorded in notebook §6 (NOT overriding the headline):** the
  `not_converged` fires because the projectile KE has a large *reversible* excursion
  across the slab (KE_min 1.11 Ha at centre, recovers to ~3.5 Ha on exit), not an
  unfinished transit — so the skill's documented remedy "extend the run" is **inverted**
  (extending worsens the wrap). Equal-potential **slab-face** cross-check (z −12.5↔+12.5):
  S = 0.95 eV/Å (KE) / 0.99 (E_total), both channels agree ~4%, matches the independent
  `--measured-s 0.018632`. → candidate skill refinement: equal-potential-face window for
  a charged projectile in a localised slab.
- [ ] **PENDING: user critique of `p2_classical_truncated_stopping.ipynb`** (training loop).

### Skill simplification (user directive 2026-06-25)

User locked a simpler, deterministic two-branch rule and asked to encode it:
- **Localised slab → Method B `dE_total/L_slab` is PRIMARY**; kinetic / ∫F·v /
  equal-potential-face are **sanity checks only** (never averaged in).
- **All other (bulk/continuous) → fixed 20%-of-simulation-TIME transient cut** +
  free-intercept slope fit, as the default. Agent slope-plateau detector demoted to
  an optional steady-state *probe*, not the headline.

Encoded:
- `stopping_power.py`: added `fixed_time_fraction(t, x, E, frac=0.20)` (cuts first
  20% of time, fits ΔE(x) on the remainder); `_selftest()` extended (recovers known
  slope, `basis='time'`, `t_cut` correct) — **passes**.
- `SKILL.md`: §1 decision table now PRIMARY-vs-sanity; §2 leads with fixed-20%-time;
  §3 Method B = slab primary + the `not_converged` reversible-well reading caveat;
  **§3a periodic-wrap truncation** (the p2 lesson); §4 sanity-channels reframed as
  cross-check-only (+ equal-potential slab-face channel); §5 reporting = headline is
  geometry-fixed; validation-status lists the p2 worked example.

## 2026-06-25 — FFT-pipeline panel as a figure standard (user directive)

User: "for each FFT processing of a signal, a panel of figures" — row 1 raw |
de-trended; row 2 windowed | zero-padded; row 3 |FFT| linear | log. Add to the
skill(s); update the latest jellium_slab run-notebooks.

- **`inqview.visualisation.fourier_panel`** (new): `fft_pipeline_panel()` +
  `fft_stages()`. Derives every stage from `analysis.fourier.FourierTransform`
  itself (`_apply_subtract` / `window.build` / zero-pad) and reuses `transform()`
  for the authoritative FFT, so the panel always matches the kernel. Energy axis =
  **angular** ℏω=2πf·E_h (matches `omega_eV` in the blessed `fourier_analysis.ipynb`
  + `pipeline/spectral_weight`). Peak-in-band locator + Δω=2π/τ resolution note.
- **Test + catalogue row:** `test_fourier_panel.py` (3/3 PASS): stage-6 == `transform()`,
  detrended/windowed/pad identities, 6 axes, 5 eV tone peak, transient-skip.
- **Skills:** `run-notebook` (Collective-response row → mandatory FFT-pipeline panel)
  and `notebook-making` (new figure-standard subsection) updated.
- **Run-notebook builder** auto-emits the panel for the most dynamic dipole component
  under "Collective response" (plasmon band + resolution annotated).
- [x] Rebuilt p2cl/p2wp jellium_slab run-notebooks (60 / 79 cells, **0 errors** each);
      both embed `fft_pipeline_dipole_z.png` under "Collective response".

## 2026-06-25 — Verdicts delivered + code changes (task 3, via code-test)

User delivered the three dossier verdicts (verbally; verdict-prose edits to the
dossier files were declined — the lines remain BLANK pending the user's own
wording / a terser version to approve):

1. **loss-function:** formula accepted **as a peak-locator** with the bold caveat;
   use the **complex** `n_q` (BUG-A) and `|n_q|²/q²` (BUG-B). reduced-test accepted.
2. **fft-drift-removal (OVERRIDE):** single **uniform `mean`** baseline for ALL
   observables incl. energy (no exception — accepts residual low-ω ramp on energy;
   peaks unaffected). Pipeline = mean → Hann → 4× pad → coherent-gain rfft. Expose
   the `subtract=` param. Headline = the **mean** curve with **detrend overlaid
   dashed** as the comparison. Flip the **library default** to mean.
3. **fft-normalization:** coherent-gain convention correct; fix accepted (already
   shipped lines 210–211).

**Code changes (all tested, full suite 171 passed / 5 xfail, deps-clean intact):**
- `analysis/fourier.py`: default `subtract` flipped **detrend → mean**
  (back-compat: explicit `detrend=`/`subtract=` honoured; `detrend: bool|None`).
  Test `test_fourier.py` +2 cases (19/19).
- `pipeline/density_fourier.py`: extracted `loss_locator()` — complex FFT (BUG-A)
  + `|n_q|²/q²` (BUG-B), bold "peak-locator not −Im[1/ε]" caveat; CSV col +
  plot ylabel/title updated. New `test_density_fourier_loss.py` (3/3: phasor peaks
  in +freq half, `.real` folding ratio 0.25, 1/q² scaling).
- `visualisation/fourier_panel.py`: baseline stage relabelled "baseline-removed";
  FFT axes now **overlay a dashed detrend comparison** vs the mean primary.
  `test_fourier_panel.py` +1 case (4/4).
- `docs/validation/test-catalogue.md`: 3 rows added/updated.
- [x] **Dossier verdicts filled** (user authorised "fill sensibly", 2026-06-25):
      loss-function (accepted as locator + complex/`|n_q|²/q²`), fft-drift-removal
      (uniform `mean` OVERRIDE), fft-normalization (coherent-gain accepted).
- [x] **Task 5 — `fourier-analysis` skill ENCODED** (`.claude/skills/fourier-analysis/`):
      `SKILL.md` (locked pipeline mean→Hann→4×→coherent-gain, angular ħω, peak-in-band
      rule, 6-stage panel standard, loss-locator caveat, known-answer anchors) +
      self-contained numpy `fourier_kernel.py` with `_selftest()` (6 assertions PASS).
- [x] Rebuilt p2cl/p2wp run-notebooks (57 / 74 cells, **0 errors** each); both
      `fft_pipeline_dipole_z.png` panels regenerated (137 / 147 KB) showing the
      `mean` primary + dashed `detrend` comparison on the FFT axes.
- [ ] **PENDING: user review** of the code changes, the skill, and the rebuilt notebooks.
- [ ] **Gate decision (user):** both skills exist + 3 verdicts filled = the
      strict gate-release condition is MET → campaign `status` may flip to `done`,
      releasing `quantum-kick-extension` and `cap-jellium-loss-function`. Held for
      user sign-off (task 1 stopping-notebook rewrite still open in plan).
