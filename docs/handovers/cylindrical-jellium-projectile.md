# Handover: Cylindrical (annular) jellium tube — projectile-down-bore stopping

Rolling handover for the `cylindrical-jellium-projectile` campaign
(`docs/campaigns/cylindrical_jellium/cylindrical_jellium_projectile.md`).
Branch: `cylindrical-jellium` (session), off `overnight-gaussian-classical`.

---

## Milestone: 2026-06-30 (i) — stopping-skill reporting mandate, dual-S, FFT via fourier skill

Three further user directives (round 2):
1. **Energy-deposit method IS the stopping power; KE is ONLY a sanity check —
   state it every time.** Codified in `.claude/skills/stopping-power-extraction/
   SKILL.md` §5 (new "Mandatory deliverables": the primary-method plot WITH the S
   result annotated on it, PLUS the KE sanity metric clearly labelled as a check).
   `per_run.py` stopping plot relabeled: panel 2 = "PRIMARY — defined method
   S=dE_total/dx" (bold), panel 3 = "SANITY CHECK: KE method vs primary".
2. **r_s=2: report BOTH the KE-based and the authentic ΔE_total S.** The per-run
   header now prints `S(primary ΔE_total) … ; kinetic cross-check S … (ratio …)`
   for every classical run; r_s=2 shows e.g. rs2_v0p15 S=0.0136 (flagged r²=0.23)
   vs S_KE=0.0067. Both visible; verdict still the user's.
3. **All FFT-driven observables go through the `fourier-analysis` skill.** The
   pipeline's raw `fft_*`/`spectra/` plots are NO LONGER embedded (excluded in
   `per_run.collect_pipeline_figs`). New `per_run.fft_panels()` emits the audited
   6-stage `fft_pipeline_panel` (Hann + mean baseline + ×4 zero-pad + coherent gain
   + ANGULAR ħω, detrend overlaid) for `current_z` (hydrovoltaic signal) and
   `energy_total`, with the plasmon band ħω_p=√(3/r_s³)·27.211 eV shaded and
   Δω=2π/τ≈5.7 eV annotated (τ=30 a.u. coarse — informational, not a gate).
   **Same fix applied to the shared `run-notebook` builder**
   (`.claude/skills/run-notebook/run_notebook_builder.py`): raw pipeline FFT figs
   excluded; the audited panel now covers the most-dynamic dipole + current_z +
   energy_total; baseline switched `initial`→`mean` (skill canonical default). So
   EVERY future notebook routes FFT through the skill.

Rebuilt + verified: notebook 0 errors, 267 cells, 170/170 images resolve, 10
audited FFT sections (20 panels), 0 leftover raw pipeline FFT plots. skill
self-tests pass.

---

## Milestone: 2026-06-30 (h) — stopping made skill-compliant + ion overlay + r_s=2 flag

User feedback on (g): (1) stopping power not extracted correctly — must use the
`stopping-power-extraction` skill for classical projectiles; (2) overlay the moving
ion on the density plots, from each run's OWN data.

Reworked `hypotheses/annular_sv/per_run.py`:
- `stopping_analysis()` now imports the skill kernels
  (`.claude/skills/stopping-power-extraction/stopping_power.py`). Geometry =
  continuous traversal (periodic tube) → **Method A**. PRIMARY = free-intercept
  slope of the **electronic deposit `ΔE_total(x)`** (`energy_total` is electronic;
  it RISES +0.045 Ha as the ion loses KE − verified `ΔE_total ≈ −ΔKE_ion`). Window =
  early `v≥0.85·v0` (light-projectile rule overrides the skill's default 20%-time
  cut). Guards run: `conservation_guard` N(t) (0% drained, no CAP) + energy
  conservation. Kinetic channel `−dKE_ion/dx` = independent cross-check. FLAGS
  channel divergence >10% / r²<0.8 — never averages.
- Density GIFs + carpets now **overlay the moving projectile** (cyan marker+trail
  on GIFs, z(t) line on carpets): classical ion from `electron_track.csv` (x,z),
  WP centroid from `wp_real_space_stats.csv` (x_mean,z_mean). Per-run, so always the
  right run's data. Verified cyan marker pixels present in classical + WP GIFs.
- Fixed a real bug: per-run `launch_z` is read from each `run_summary.txt`
  (rs6=−23, rs4=−13, rs2=−4) — was hardcoded −23.

**KEY SCIENTIFIC FINDING (the skill guards surfaced it).** Channel agreement
(ΔE_total vs −dKE_ion) by wall density:
- r_s=6 (3 runs): ratio 1.01–1.05, r²≥0.90 — CLEAN.
- r_s=4: mostly clean; rs4_v0p15 marginal (12% divergence).
- r_s=2 (3 runs, short L_z=10 tube): ALL FLAGGED — rs2_v0p15 ratio 2.04 r²=0.23,
  rs2_v0p30 ratio 1.18 r²=0.69, rs2_v0p45 ratio 1.11. The dense short-cell runs give
  an unreliable ΔE_total slope (tiny traversal: rs2_v0p15 moves only ~2.5 Bohr).

OPEN DECISION (user owns the verdict): the sweep headline `Sv_results.csv`/β uses
the **ke_ion** channel (rs2_v0p15 S=0.0067) while the skill PRIMARY is ΔE_total
(0.0136, flagged). For r_s=6,4 they agree so β is robust; r_s=2 is the question —
accept ke_ion, accept ΔE_total, or rerun r_s=2 with a longer L_z. NOT silently
resolved. Notebook now shows both channels + flags per run.

Verified: notebook 0 errors, 277 cells, 180/180 images resolve; skill `_selftest`
passes.

---

## Milestone: 2026-06-30 (g) — per-run deep-dive sections added to report notebook

User wanted, per projectile run, the **matrix of density visualisations** + the
**high-value observables** (`docs/notes/plots_examples.md`) embedded as a section
each in `hypotheses/annular_sv/annular_sv_report.ipynb`, so every run is fully
inspectable. DONE + verified.

What was built (all in `hypotheses/annular_sv/`):
- `per_run.py` — TUBE-aware per-run generator (self-contained; uses canonical
  `inqview.load_vti`, physical order). Produces: density matrix GIFs
  {density, Δn=n(t)−n(0), Δn=n(t+dt)−n(t)} × {total[, wp, bath]} with **vertical
  wall-radius markers x=±5,±13** (NOT slab faces — there are none/no CAP);
  z–t carpets; initial-drag stopping plot (light-projectile method). Pipeline
  observables via `inqview.pipeline` phases `observables, momentum, kl_divergence`
  (energy decomposition, induced current+FFT, dipole, momentum incl/excl WP, KL).
- `build_per_run_figs.py` — driver: loops the 10 runs, writes
  `per_run_manifest.json` (figure paths RELATIVE to `hypotheses/annular_sv/`).
- `build_report.py` — extended: `per_run_cells()` reads the manifest and splices
  one `## Run <label>` subsection per run (path-referenced `<img>` → small .ipynb).

Why NOT the slab `run-notebook` builder / `make_density_gif_battery`: both are
slab-geometry (slab faces ±12.5, CAP, S=ΔE/L_z). This is a PERIODIC annular tube —
slab markers + slab stopping would be physically wrong. Reused only the
geometry-agnostic ideas, kept tube-correct stopping (initial drag).

Verified:
- Tube geometry sanity: wall density peaks at |x|=8.8 Bohr (centre of [5,13]),
  bore hollow → no centre↔edge swap (vti-coordinate rule).
- `per_run.py` stopping S(rs6,v0.30)=0.00395 Ha/Bohr == campaign `Sv_results.csv`
  (0.00394) → initial-drag extraction consistent.
- Notebook executes **0 errors**, 277 cells, all 10 runs sectioned; classical →
  3-GIF matrix, WP → full 9-GIF matrix (Total/WP/Bath) + Momentum + KL.
- **180 embedded images, 0 missing** (all paths resolve).

Honest data-limitation notes carried in the notebook (NOT fabricated): runs did
not store `state_energies.csv` → no KS eigen-energy bar-GIFs; `eigenvalues.csv`
not retrofitted → no GS KS-excitation decomposition; no E-field pipeline phase;
only 1D |k| momentum → no 2D (k_z,k_⊥) scattering map.

Plan: `docs/plans/annular-sv-per-run-sections.md`.

---

## Milestone: 2026-06-28 (f) — BUG FOUND + FIXED (light-electron deceleration)

First autonomous run (run1) FAILED at the pilot gate and emailed two PHASE-FAILED
notices. Diagnosis (real physics, not a code crash):

- The projectile is a LIGHT electron (m_e), so at v=0.30 its KE is only 0.045 Ha.
  The stopping force decelerates it ~85% (vz 0.300→0.045) within ~6 Bohr,
  depositing all its KE (ΔKE_ion −0.044 ≈ ΔE_system +0.045 Ha; energy conserved).
  The electron STOPS long before a 5-plasma-period wake can form.
- Three code bugs this exposed: (1) `n_steps_for` sized runs to 13k+ steps
  (5·2π/ω_p wake criterion) — the electron stops by ~1500 steps and then sits;
  (2) `extract_S` regressed over the full post-20% window → averaged S over v from
  0.25 down to 0.04, NOT S at v0; (3) the pilot gate ABORTED the campaign on the
  (by-design) 85% v-drift, and the finalizer then polled uselessly for 9.5 h.

FIXES (all in `scripts/annular_sv/orchestrate.py`, mirrored in `finalize.py` +
`hypotheses/annular_sv/build_report.py`):
- `n_steps_for` → `ceil(max(30, 100·v)/dt)` (1500–2250 steps): capture the
  initial-drag window + deceleration sweep, not the wake.
- `extract_S` → S(v0) = INITIAL drag = −d(KE_ion)/ds over the early
  near-constant-velocity window (vz ≥ 0.85·v0; widen to 0.70/0.50 if sparse).
  Uses the per-step track (hi-res). VALIDATED on the run1 pilot data:
  **S(r_s=6,v=0.30) = 0.00394 ± 0.00005 Ha/Bohr** (v_mean=0.286, 540 pts).
- pilot gate no longer aborts on v-drift; it requires a clean initial-drag S
  (finite, ≥30 pts) instead. Now PASSES.
- Killed the stuck finalizer; reused the validated pilot run as the rs6_v0p30
  production point (copied → `annular_sv/rs6_v0p30`).

RELAUNCHED 2026-06-28 22:51: orchestrate.py (PID 483897) + finalize.py (PID
484307), both GPUs. Pilot PASS; production fanning out (rs6_v0p15 GPU0, rs6_v0p45
GPU1 — both stepping, no OOM). GPU1 shared with another user (~18.5 GB; my runs
~1-2 GB fit the 5.5 GB free). ETA ~2-2.5 h → WP rung + notebook → emails.
Note: S(v) here = instantaneous stopping at launch v0 (the decelerating projectile
sweeps a v-range; initial slope = friction force at v0) — scientifically the right
quantity, consistent with the locked free-Ehrenfest design.

---

## Milestone: 2026-06-28 (e) — EXECUTION STARTED (autonomous, both GPUs)

User: "use both the GPUs and run orchestrate using python the cylindrical
localised jellium campaign." `status: running`. Acting as the executing agent.

**Preflight re-verified:** both GPUs FREE (cudaMemGetInfo probe: 23815/24062 MB
free each; nvidia-smi down = NVML mismatch, compute fine). Using venv python.

**Phase 0 — DONE.**
- Scaffold created (ADR-0007): `ResearchProject/systems/cylindrical_jellium/`
  with `shared_gs/ shared/{configs,pseudopotentials}/ scripts/annular_sv/
  {gs,classical,wp}/ annular_sv/ hypotheses/annular_sv/tests/`. Sweep = `annular_sv`.
- Literature grounding: 8 source notes in `docs/sources/`. Corrections: the real
  TDDFT-PENN paper is Matias 2025 (arXiv:2505.23396, verified) + Penn 1987
  (jellium-of-varying-r_s ⇒ different materials); arXiv:1805.01377 was the WRONG
  paper (Si, no Penn); Netz NOT a Nature-2022 author; Seguí–Arista 2007 =
  cylindrical-tube stopping; Lindhard 1976 grounds charge-even electron-as-cation.
  A few citation details flagged for re-verification (honest, not invented).

**Phase 1a — DONE (both ingredients gated before any production run).**
- (i) **Annular `background_shape`** added to
  `inq-stack/include/inqkit/jellium/localised_background.hpp`
  (enum `annulus`; `inner_radius` param; branch
  `mask = ½erfc((d−R_out)/w)·[1−½erfc((d−R_in)/w)]`, d ⟂ tube axis).
  GATES: `formula-validation` agent VERDICT **CONFIRM** (all 6 invariants +
  inner-complement erfc identity; caveats R_in=0/thin-wall both outside our
  regime); `code-test` **6/6 pass** (engine test T0.4 neutrality, T0.5 bore-carve
  additivity, T0.6 erfc-smoothed charge — `test_localised_background_engine.cpp`);
  catalogue rows added to `docs/validation/test-catalogue.md`.
- (ii) **Projectile UPF VERIFIED (no regeneration needed).** Key finding: INQ's
  pseudopod reader IGNORES `is_coulomb` (not referenced in inq/src,
  pseudopotential.hpp, or upf2.hpp — `upf2.hpp:235` reads PP_LOCAL verbatim ×0.5
  Ry→Ha). The stale `is_coulomb="T"` header is harmless; INQ uses the tabulated
  Gaussian. `electron_gaussian_wpsigma0p5.upf` verified by DATA
  (`hypotheses/annular_sv/tests/verify_projectile_upf.py`): V(0)=+2.257 Ha
  repulsive, Z=1.000, **σ_pot=0.3536** (fit rms 1.7e-5 Ha = σ_WP/√2 for σ_WP=0.5).
  Cutoff guard PASS (classical E_cut=537 eV≫3 eV; WP aliased-tail 0.00%). Copied
  to `cylindrical_jellium/shared/pseudopotentials/`. (sigma0p35 = equivalent
  alternative, σ_pot=0.350.)

**Phase 1b — IN PROGRESS.** Written: `shared/configs/annular_tube.hpp` (locked
geometry table); `scripts/annular_sv/gs/run.cpp` (pure-env annulus GS; writes GS
density + n₊ VTIs for the slice-emails). r_s=6 GS build+run LAUNCHED on GPU 0
(`scripts/annular_sv/gs/build_run_rs6.log`). Per-density: r_s=6 L_z48 N24,
r_s=4 L_z28 N48, r_s=2 L_z10 N136 (n0=N/V_annulus, exact neutrality).

**Phase 1b — DONE.** All 3 GS validated (battery 3/3 PASS): rs6 E=−3.78 Ha,
rs4 E=−19.05 Ha, rs2 E=−352.4 Ha; neutrality 0.00%, wall plateau 1–7% of n0,
bore 0.03–0.12·n0, symmetry 0.2–3.4%. Checkpoints in `shared_gs/tube_rs{6,4,2}`.
Slice-emails SENT (6 PNGs, all 3 densities). Validator:
`hypotheses/annular_sv/tests/validate_gs.py`; plots in
`hypotheses/annular_sv/gs_validation/`.

**Phase 2/3 — RUNNING AUTONOMOUSLY (detached, both GPUs).** Classical glide
binary built (`scripts/annular_sv/classical/run`; merges sv_ladder glide +
annulus bg, Ehrenfest, track + current_z + density VTI + induced delta).
200-step smoke PASSED: ~1.0 s/step, energy finite (drift +1.4e-3), projectile
DECELERATES (vz 0.300→0.297 over 4 a.u. — stopping signal present), wall current
+ induced density responding. Orchestrator launched:
`nohup orchestrate.py` (PID 3163983, log `scripts/annular_sv/orchestrate.log`).
Runs short pilot gate → 9-run sweep across GPU 0+1 (~12.4 h wall, ~89k steps)
→ S(v)/β(r_s) extraction → email. Idempotent resume (skips run_completed=true);
per-phase failure emails. N_STEPS per run = max(1.5 L_z/v, 5·2π/ω_p)/dt.

**Phase 4/5 — BUILT + QUEUED autonomously (do NOT touch the production run).**
- P4 WP run.cpp written (`scripts/annular_sv/wp/run.cpp`; WP electron σ_WP=0.5,
  k0=0.30 injected on-axis, annulus bg, periodic, momentum+real-space stats).
  Binary BUILT with `CUDA_VISIBLE_DEVICES=""` (compile-only, no GPU contention).
  Matched classical ghost = the existing rs6_v0p30 production run.
- P5 notebook builder `hypotheses/annular_sv/build_report.py` (executable cells:
  S(v)/β, induced wall current, wake xz linear+log, WP-vs-classical; canonical
  theme; load_vti no-fftshift). Syntax + all 6 cells validated; structure builds.
- **`scripts/annular_sv/finalize.py`** launched DETACHED (PID 3185966, log
  `finalize.log`): polls until all 9 production runs complete, then runs the WP
  rung (GPU 0, 6h timeout = injection-deadlock guard) + executes the notebook +
  emails. Separate process; does NOT modify the running orchestrator.

**Engine fact (verified, reusable):** INQ's pseudopod reader IGNORES `is_coulomb`
(`upf2.hpp:235` reads PP_LOCAL verbatim ×0.5 Ry→Ha; flag absent from inq/src,
pseudopotential.hpp, upf2.hpp). The `electron_gaussian_*.upf` "Coulomb" headers
are STALE but the Gaussian PP_LOCAL data is what INQ uses — verify UPFs by DATA
(V(r) shape), never by header.

**State of GPUs:** both busy (production). WP binary built without touching them.
Two detached Python processes own the campaign to completion + emails.

**Caveat flagged:** σ_pot=0.354 < dx=0.5 → the Gaussian potential is coarsely
sampled (~0.7 grid pts/σ). Cutoff guard (momentum-based) passes; the P2 pilot
energy-drift gate is the real test of this resolution choice.

---

## Milestone: 2026-06-27 (d) — PLAN LOCKED FINAL

User: "Lock the plan as final." Design is FROZEN. The campaign prompt carries a
`LOCKED FINAL` banner; execute AS-IS (re-verify <preflight>, run Phases 0→5).
No further design changes without an explicit user instruction.

Final-state verification (passed): `status: ready`; 7 phase-tasks; 10 section
tags; preflight intact; phases P0–P5 all present; no placeholders/contradictions;
frontmatter parses + INDEX regenerates clean.

Five defaults stand as locked (user did not override): grounding = P0 (non-running);
dx/L_xy convergence kept in GS battery; stationarity at end of P1; Phase 1 split
1a/1b; GS slice-emails non-blocking. NOT committed to git (offered).

---

## Milestone: 2026-06-27 (c) — projectile UPF gap fixed + Ehrenfest pinned

User probed "where do we make the Gaussian projectile?" — exposed an over-
optimistic reuse assumption. Verification found the `electron_gaussian_*` UPFs all
carry STALE Coulomb headers (copied from the antiproton base), so the asset cannot
be trusted by filename/header, and the σ-convention file was not cleanly pinned.
No generator script found in-repo.

Fixes written into the prompt:
- **Projectile = classical electron** (−1, species "H"), Gaussian radial potential.
  NOT a positive H⁺ (user confirmed "call it electron to avoid confusion").
- **Dynamics = free EHRENFEST** (user: "ensure we are using ehrenfest") — projectile
  is an Ehrenfest ion under the real KS force, velocity evolves self-consistently;
  made explicit in <resolved_decisions> (projectile + propagation) and <preflight>.
  Propagation = real-time LDA, ETRS electrons + Ehrenfest ion, dt=0.020.
- **Phase 1a is now "Ingredients"**: (i) annular shape (code-test+formula-validation
  +catalogue) AND (ii) projectile Gaussian UPF **VERIFY-then-REGENERATE** — verify
  candidate `electron_gaussian_wpsigma0p5.upf` (charge −1; finite repulsive V(0);
  σ_pot≈0.354; cutoff guard); regenerate + check in a generator if it fails.
- σ convention: label σ_WP=0.5 ⇒ verified UPF target σ_pot≈0.354.
- Status stays `ready`; still NO runs/code executed.

---

## Milestone: 2026-06-27 (b) — plan reorganised into validated phases

Restructured the campaign `<tasks>` + frontmatter into a 6-phase complexity
ladder (user request). GS is the **first running phase**; grounding/scaffold is
non-running prep (P0).

- **P0** Grounding & scaffold (non-gating) → **P1** Static system + ROBUST GS
  (1a annular shape coded+pre-gated; 1b GS battery + slice-emails) → **P2**
  Propagation validation (stationarity + classical-glide smoke) → **P3**
  Production S(v) sweep → **P4** Quantum rung → **P5** Synthesis notebook.
- **Robust GS battery (P1b):** SCF convergence; electron neutrality; radial
  profile n(d) [flat interior=n₀, Friedel at BOTH surfaces, spill-out, small bore
  density]; cylindrical symmetry; energy sanity vs bulk LDA; grid/box convergence
  spot-check (r_s=6: dx 0.5↔0.4, L_xy 40↔48); stationarity (no-projectile drift).
- **GS slice-emails added to P1b** (user request): per density, xz/yz/xy slices of
  BOTH n₊ (prescribed background) and n (converged electron density), via
  `inqview.load_vti` (physical order, no fftshift), emailed through the
  `email-notifications` skill (four-part + attached PNGs). The autonomous agent
  PROCEEDS on the numeric battery; the email is the user's visual record / optional
  manual checkpoint — it does NOT block.
- Status stays `ready`; still NO runs/code executed.

---

## Milestone: 2026-06-27 — design LOCKED, campaign promoted draft → ready

### Status
The rough draft was advanced to an **autonomy-ready** prompt via the `campaigns`
skill (Mode A grilling, all forks user-locked). The campaign is `status: ready`:
a fresh agent can execute it end-to-end. **No runs launched, no code written yet.**
The single command to start it autonomously: run the campaign prompt file above.

### Locked decisions (user, 2026-06-27)
- **Scope = 2 rungs.** Rung 1 = classical-electron r_s sweep → S(v) → β(r_s).
  Rung 2 = electron WP vs matched classical electron ghost at r_s=6 (quantum
  check). Deferred to later campaigns: off-axis projectile, multi-r_s concentric
  layers (Penn), the literal positive-cation variant.
- **Projectile = electron throughout** (charge −1, mass m_e, free Ehrenfest →
  genuinely decelerates). Gaussian erf-smoothed via the EXISTING asset
  `electron_gaussian_wpsigma0p5.upf` (σ_WP=0.5 ⇒ σ_pot≈0.354 per the √2 rule).
  On-axis (x=y=0) +z launch. Velocity env-driven (sv_ladder pattern) — one build
  per density serves all 3 velocities. **No new projectile code/asset.**
- **Geometry = annular jellium tube, axis ∥ z, PERIODIC infinite tube.**
  R_in=5, R_out=13 (8 Bohr wall), L_xy≈40, edge_width w≈1, dx=0.5.
  Per-density **L_z ≈ {48, 28, 10} Bohr** sized to 2×wake length at v=0.45
  (fixed transverse geometry keeps S(r_s) comparable; L_z is free because S is
  per-unit-length on a z-uniform tube).
- **r_s = {6, 4, 2}** → **N ≈ {24, 48, 136}** (round ∫n₊ to nearest even, then
  set n₀=N/V_annulus for exact neutrality). r_s=6 = cheap smoke + WP point;
  r_s=2 = expensive tail, run last. **Hard cap N ≤ ~260** (user).
- **Velocity = common absolute v = {0.15, 0.30, 0.45} a.u.** (user accepted that
  r_s=6's top point is mildly supersonic, v_F=0.32). 3×3 = 9 classical runs +
  1 WP pair (2 runs) at r_s=6.
- **dt = 0.020 a.u.** (4× the old Coulomb-forced 0.005; Gaussian smoothing
  removes the singularity). Smoke test validates energy conservation.
- **New code = ONLY the annular `background_shape`** (one enum value + one
  `else if` branch reusing `background_mask`): `mask = background_mask(d,R_out,w)·
  (1 − background_mask(d,R_in,w))`, `d=√(x²+y²)`, z-uniform. Pre-gated
  (code-test + formula-validation + catalogue row).
- **S extraction** = ΔE_system 20%-transient regression vs path (continuous-glide
  method, `stopping-power-extraction` skill).
- **Placement** = NEW system `ResearchProject/systems/cylindrical_jellium/`
  (full ADR-0007 layout; configs `#include` jellium shared base headers).

### Verified engine facts (source line-refs)
- `background_shape { slab, sphere, box }` — no annulus
  (`inq-stack/include/inqkit/jellium/localised_background.hpp:52`). Mask helper
  `½erfc((d−R)/w)` (`:70`).
- `background_perturbation.hpp:61-67` applies `v_bg=−poisson(n₊)` via a
  shape-agnostic `potential()` callback — annular n₊ plugs straight in.
- Existing run templates: `ResearchProject/systems/localised_jellium/scripts/
  {fullsuite_classical,fullsuite_wp,qsp_phase3/{gs,classical,wp}}/run.cpp`.
- Projectile recipe (classical electron, m_e, env-driven v, Gaussian UPF):
  `ResearchProject/systems/jellium/shared/configs/sv_ladder_L50_sigma0p5.hpp`
  and `electron_proj_E100_L50_cubic_sigma1.hpp`.
- Gaussian electron UPFs present:
  `ResearchProject/systems/jellium/shared/pseudopotentials/
  electron_gaussian_wpsigma0p5.upf` (σ_WP=0.5).

### Done / partial / not done
- DONE: full design grill; all decisions locked; campaign promoted to `ready`;
  index refreshed.
- NOT DONE (the autonomous run): annular shape code+tests; 3 GS; smoke; 9
  classical runs; WP rung; analysis notebook. All specified in the prompt.

### Open risks / caveats (carried into guard rails)
- r_s=6 supersonic wake at v=0.45 may be longer-lived than 2λ; smoke test MUST
  verify a clean F_z plateau across the z-period wrap. Headroom: N=24 at L_z=48,
  so L_z can rise to ~96 (N≈48) for the dilute runs if contaminated.
- r_s=6 is a small electron gas (~24 e) — finite-size/shell effects; results
  PROVISIONAL, cross-check the trend against r_s=4.
- Electron-as-proxy-for-cation rests on S being charge-even at leading order
  (Barkas is the odd correction) — grounding task must cite this.
