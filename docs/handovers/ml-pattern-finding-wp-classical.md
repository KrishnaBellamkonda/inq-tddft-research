# Handover — ml-patterns: pattern-finding in WP/classical runs

Campaign: `docs/campaigns/ml-patterns/pattern-finding-in-wp-classical-runs.md`
(`id: ml-pattern-finding-wp-classical`, status `draft`, 1/8).
Authored via grill-with-docs + campaigns Mode A (Stage 0). Last updated 2026-06-30.

## Goal
Use ML as an INTERPRETABLE discovery tool to find/explain spatial & dynamical
differences in the **induced density** n(r,t) between a **classical** and a
**quantum (wavepacket)** projectile at matched (σ, v) in jellium — physics the
scalar stopping power S(v) cannot see. Discovery aim, complexity ladder
(static → dynamic), SIE/vacuum control mandatory.

## DONE (locked, verified)
- All run-DATABASE design decisions locked (see campaign `<resolved_decisions>`):
  new from-scratch rich DB, wide CSV (`docs/run_database.csv`) + JSON mirror
  (`docs/run_database.json`), all 6 systems + `relevant_to_induced_density` flag,
  deepest parse (run_summary→config→run.cpp), full enriched schema, NULL token,
  twin linkage. Scripts → `docs/campaigns/ml-patterns/`. Builder=Opus,
  validator=Sonnet (100% scripted checks + stratified ~40–60 re-parse + all
  compound-line/fallback runs; ≤2 fix rounds).
- Scientific framing locked: aim = discovery; scope = static+dynamic ladder;
  4 target signatures (exchange/xc-hole, diffraction, wake phase, form-factor);
  vacuum-WP SIE control + linear-response reference.
- CONTEXT.md: added "Run inventory: catalogue vs database (2026-06-30)" glossary
  (catalogue vs database, twin_run_id, relevant_to_induced_density).
- Stage-0 rough draft + open-questions written to the campaign file. INDEX
  regenerated (25 campaigns).

## IN FLIGHT (background agents — await completion notifications)
- **DB builder** (Opus, agentId `a4c28ac0e527c8494`): writing
  `docs/campaigns/ml-patterns/build_run_database.py`, emitting the CSV/JSON +
  `run_database_columns.json` data dictionary. UNVERIFIED until it returns.
- **Deep research** (Opus, agentId `a3767209ada7f6d41`): cited synthesis →
  `docs/campaigns/ml-patterns/research/ml_induced_density_research.md` +
  `docs/sources/` notes. Answers: precedent, methods-per-rung, 4-signature
  physics grounding, adjacent questions, pitfalls. UNVERIFIED until it returns.

## NOT DONE / NEXT
1. When builder returns → dispatch the **Sonnet validator** over its output;
   run the ≤2-round fix loop; confirm DB clean.
2. Fold DB inventory into **Stage 2** (which matched twins actually exist;
   VTI cadence/resolution adequacy for Rung-2 dynamics).
3. Fold research into **Stage 1** (sharpen the single falsifiable primary
   signature + success/failure thresholds) and method selection per rung.
4. Continue grill: Stages 1→4 (matrix, grounding, validation/guard rails),
   then Stage 5 autonomy checklist before `status: ready`.

## 2026-06-30 milestone — science framing locked + DB feasibility

**Deep research DONE** (`docs/campaigns/ml-patterns/research/ml_induced_density_research.md`
+ 3 source notes). Verdict: gap real (nobody applies interpretable ML to the
induced-density FIELD; Ward 2024 = scalar surrogate to invert). Of the 4
signatures, wake (iii) + form-factor (iv) are ROBUST; exchange (i) + diffraction
(ii) are SIE-confounded → exploratory only.

**Science decisions LOCKED (see campaign `<resolved_decisions>`):**
- Primary target = wake (iii) + form-factor (iv); discovery axis = nonlinear
  residual vs Lindhard; falsification gate λ(v)=2πv/ω_p via DMD.
- Substrate = **point-classical (ONCV) vs finite-σ WP at matched ENERGY/velocity**;
  references = Lindhard + vacuum-WP (SIE). **No new runs** — analysis-only campaign.
- Observable = bath response n_bath(t)−n_bath^GS. Normalisation ladder = GS →
  rigid-motion → Lindhard → vacuum-WP. Never fftshift a VTI.
- Methods: Rung 1 POD on Δn = n_WP−n_classical + persistent homology; Rung 2
  DMD/Koopman (windowed early-v) + SINDy. jellium-vs-slab scope STILL OPEN.

**DB feasibility (reconciled):** initial twin linker under-counted. Recovering σ
from UPF + matching on velocity + point-as-σ→0: **jellium 172, localised 26**
matched WP↔classical pairs with density. Rung-1 bulk grid: σ_WP∈{0.5,1,3,5,8} ×
E∈{20…1500}; clean cuts — form-factor @E=100 (σ sweep), wake @σ=5 (E sweep).

**DB-revision bundle LOCKED (one builder pass AFTER validator, then re-validate):**
1. `classical_potential_form` {gaussian,coulombic} by PP_LOCAL V(r) inspection.
2. Rework twin linkage + `match_type` {point_vs_wp, sigma_matched_gauss, exact}.
3. Fill classical `sigma_wp_bohr` (from UPF) + always-fill `velocity_au` (v↔E).
(Dropped: separate classical_gaussian_width — redundant w/ filled σ + σ_pot.)

## IN FLIGHT
- **Sonnet validator** (agentId `a6580e0219137c1d0`) running on the v1 DB
  (581×133). Output → `docs/campaigns/ml-patterns/run_database_validation.md`.

## NEXT (in order)
1. Validator returns → consolidate its fixes + the 3-item revision bundle → send
   builder (a4c28ac0e527c8494) ONE revision pass → re-validate.
2. Re-pose the jellium-first vs slab-first **system-scope** decision with the
   corrected inventory (user deferred it pending the pair reconciliation).
3. Resume grill: Rung-1 method/normalisation specifics → notebook output contract
   → autonomy adaptation (analysis-only) → Stage-5 checklist → status: ready.

## Open questions
See campaign `<open_questions>`: system scope, Rung-1 method/normalisation
specifics, notebook contract, new-code pre-gating (POD/DMD/normaliser kernels).

## 2026-07-01 milestone — grill COMPLETE; campaign assembled; T0 validated

**DB (T0) DONE.** Revision applied + re-validated (round 2 **PASS**, 0 blocker / 0
major / 3 graphene-only cosmetic minors). `docs/run_database.{csv,json}` 581×137,
validated, idempotent. `twin_run_ids`/`match_type`/`classical_potential_form` live;
jellium 77 + localised 15 + cylindrical 6 best-twins. Graphene legitimately 0 twins
(σ √2-mismatch). KNOWN MINOR: graphene seeded `run_cl_centroid_s{1,2,3}`
`energy_ev`=design vs `velocity_au`=seeded (don't match graphene on energy) — not
fixed (graphene out of scope).

**All grill decisions LOCKED** (campaign `<resolved_decisions>`): falsification =
parameter-free overlays (exp(−q²σ_pot²/2), ω_p) at ±20% on HELD-OUT cells; code =
campaign-local; autonomy = fully-autonomous Python orchestrator + ≤4-try agentic
loop tuning agreement on a PRE-REGISTERED CALIBRATION split, verdict from HELD-OUT,
all attempts logged (ADR 0011); notebooks = per-rung + synthesis, auto-built.

**Artefacts written:** full campaign prompt (8 tasks, all template sections),
`docs/adr/0011-held-out-split-anti-phacking.md`, CONTEXT.md glossary, INDEX (1/8).

**READY.** User reviewed the prompt, requested 3 fixes + pinned split (all applied),
and approved. `status: ready` set 2026-07-01; campaign in the INDEX "Ready" group.

Applied at review: (1) split PINNED — form-factor calib σ_WP∈{1,5}/held-out
{0.5,3,8}; wake calib=even-velocity-index/held-out=odd. (2) T2 headline metric =
q-ratio R(q)=n_WP(q)/n_classical(q) (POD demoted to structural support). (3) divide
by ACTUAL F_ONCV(q) (not 1) — T1 must measure it + bound the q-range where ≈1.

**NEXT (execution, after ready):** the autonomous executor builds T1 kernels
(POD/DMD/normaliser, pre-gated + measures F_ONCV) → runs T2–T7 per the ladder. No
GPU; analysis-only. NOT yet committed — user to confirm the two scoped commits.

## 2026-07-01 milestone — autonomous execution: T1 done, kernels validated, T2–T7 orchestrated

**Executor session (autonomous, CPU-only, no INQ runs).** Built all campaign-local
code under `docs/campaigns/ml-patterns/`:
- `kernels/pod.py` (truncated + randomized/power-iteration SVD POD),
  `kernels/dmd.py` (exact windowed DMD/Koopman), `kernels/formfactor.py`
  (F_WP, F_ONCV-from-UPF, radial q-spectrum, R(q)), `kernels/normaliser.py`
  (VTI series loader via `inqview.load_vti`, GS subtraction ladder, co-gridding),
  `kernels/celldb.py` (cell resolver — matches WITHIN-cut on r_s/L/dx/velocity,
  NOT the DB best_twin which crossed densities), `kernels/pipeline.py` (T2 R(q)
  agreement + T3 wake DMD).
- `tests/test_kernels.py` — 12 known-case code-tests, ALL PASS.
- `orchestrate.py` — idempotent/resumable Python orchestrator, per-phase
  try/except + 4-part Gmail + auto-built notebooks under `notebooks/`.

**T1 DONE (verified).**
- 3 independent `formula-validation` agents returned **CONFIRM** (POD, DMD,
  form-factor) — given only formula+source.
- 12/12 code-tests pass. Catalogue rows added to `docs/validation/test-catalogue.md`.
- **F_ONCV(q)** computed from the ACTUAL `electron-ONCV-1.2.upf` local potential
  (Coulomb-tail-subtracted radial FT; the potential is REPULSIVE, V→+1/r):
  **F_ONCV ≈ 1 within 5% for q ≤ 1.9 1/Bohr** (within 2% for q ≤ 1.2). The T2
  prediction reduces to exp(−q²σ_pot²/2) only inside this window. Artifact:
  `artifacts/foncv.json`, `artifacts/T1_foncv.png`.

**DATA REALITY (important).** The form-factor cut bath-only density (`density_system`
= n_total−n_wp, 162 e) is CLEAN for the held-out σ_WP∈{0.5,3,8} (`_wf` runs,
hundreds of frames) — so the **held-out verdict rests on clean data**. Calibration
σ_WP∈{1,5} at E=100 lack `_wf` bath-only runs; the resolver falls back (logged
per-cell `wp_method`), acceptable since calibration only tunes config (ADR 0011).
Wake cut σ_WP=5: 3 calib + 3 held-out energies matched; 4 energies skipped+logged
(no matched classical). Cells resolve at r_s=5.69 / L50 / dx0.4.

**SMOKE (pre-T2) PASSED.** Full T2 pipeline on σ_WP=3 held-out → finite, physical
R(q); robust pipeline (shell-smoothing + prediction-floor window + σ_eff slope
cross-check) gives σ=0.5 71% within ±20%, σ=3 45% — calibration loop tunes this.

**IN FLIGHT.** `orchestrate.py` running T1→T7 in the background (log:
`artifacts/orchestrate.log`, per-phase `artifacts/T*_result.json`, status
`artifacts/phase_status.json`, notebooks `notebooks/*.ipynb`). Idempotent — re-run
`venv/bin/python3 docs/campaigns/ml-patterns/orchestrate.py` to resume; a phase
with an existing result JSON is skipped (`--force` to redo).

**Verdict logic.** T2: per held-out σ "agrees" if ≥50% of its valid q-window is
within ±20%; CONFIRM if ≥2/3 held-out σ agree, REFUTE if 0/3 with non-trivial
windows, else INCONCLUSIVE. T3: per held-out energy "agrees" if |ω_DMD−ω_p|/ω_p
≤0.20 AND Nyquist (dt<π/ω_p); CONFIRM if majority agree. All ≤4 calibration
attempts logged in each result JSON + notebook.

**Grounding note (T1 citation check).** The `[abstract-only]` method citations in
the research doc (DMD-on-electron-phonon, POD/DMD & persistent-homology reviews,
contrastive, β-VAE) were NOT re-fetched this session; the kernels rest on the
PRIMARY method sources (Schmid 2010 / Tu 2014 for DMD; Brunton & Kutz / Halko 2011
for POD; Jackson for the form factor) which the formula-validation agents
confirmed. The SIE≈7 eV figure remains labelled project-brief / externally
unverified (used only in the caveated T6).

**NOT committed** (per instructions — user owns commits).

**NEXT.** When the orchestrator finishes: read `artifacts/phase_status.json` +
each `T*_result.json` for verdicts; flip frontmatter `done:` per completed phase
and append the final verdicts here.

## 2026-07-03 milestone — bulk-only PDE-discovery REDO designed (grill-with-docs)

**Prior run reviewed.** T1–T7 executed 2026-07-01: T2 form-factor **CONFIRM**
(2/3 held-out), T3 wake **INCONCLUSIVE** (1/3), T4 slab **INCONCLUSIVE**, T5 SINDy
thin/unvalidated (2-mode latent ODE, no forward test, no interpretation). The
governing-equation goal the user cares about was the weakest phase → redo.

**Redo fully specified (grilled 2026-07-03, all decisions user-locked).** Extends
this campaign in place; new tasks **T8–T14** in the frontmatter + a
`# 2026-07-03 — Bulk-only PDE-discovery redo` section in the campaign file. Locks:
- **Scope:** pure bulk jellium only (drop slab + other systems).
- **Two-track:** A = form-factor + wake gates re-run clean on bulk; B (headline) =
  discover a governing **field PDE** (weak-form SINDy/PDE-FIND) + latent-ODE
  support, **broad agnostic library**, physics named post-hoc.
- **Separate-then-compare:** `PDE_classical` (coulombic point sweep) and `PDE_WP`
  (σ=5/σ=1 sweep) discovered independently, compared afterwards.
- **Three validation walls** (ADR **0012**, extends 0011): pinned calib/held-out
  cell split + temporal forward-prediction + bootstrap coefficient stability.
- **Engine:** Python orchestrator extending `orchestrate.py`, spawns discovery +
  skeptic + interpreter + judge agents, per-phase Gmail, **hard 12 h cap**.
- **Termination:** convergence gate + plateau stop + 12 h fallback; honest
  best-so-far on non-convergence.

**Cells pinned from the validated DB (verified 2026-07-03):** form-factor E=100
σ∈{0.5,1,3,5,8}; classical + WP velocity sweeps E∈{20,25,50,100,300,600} at 125³,
190–457 bath frames; σ=3@E=25 has ~10 001 frames (forward-prediction cell).
Track-B split: calib E∈{20,50,300}, held-out E∈{25,100,600}. Data is adequate for
∇²/∇·/∂ₜ (Nyquist-safe; ω_p≈3.5 eV).

**Artefacts written this session:** CONTEXT.md glossary "Bulk-jellium PDE-discovery
redo (2026-07-03)"; `docs/adr/0012-agnostic-pde-discovery-three-walls.md`; campaign
frontmatter T8–T14 + design section.

## 2026-07-04 milestone — REDO build complete + validated on real data (paused for launch)

User reinforced: **fully autonomous, hands-off**. Built + validated the discovery
stack; smoke-tested on real bulk-jellium data; **paused before the 12 h launch**
per user instruction ("build kernel + orchestrator, smoke-test, pause for my go").

**Built (all campaign-local under `docs/campaigns/ml-patterns/`):**
- `kernels/pdefind.py` — PDE-FIND (STRidge, Rudy et al. 2017): broad agnostic
  library u^p·∂^d u, scale-invariant thresholding, order-1/2 targets, forward
  integration (RK4 / (u,v)-leapfrog with sub-stepping), bootstrap stability,
  post-hoc rule-based physics interpreter, `mask_to_terms` (admitted-eq scoring).
- `kernels/discovery.py` — cell → axial-reduced induced field n(z,t) → three-wall
  discovery (POD-denoise + field nondimensionalisation + forward-predict the
  ADMITTED equation + bootstrap). `load_cell_axial`, `discover_cell`, `coeff_vector`.
- `orchestrate_pde.py` — autonomous T8→T14 driver: idempotent/resumable, **hard
  12 h wall-clock cap** (checked per phase AND per cell), per-cell + per-phase
  try/except, plateau-stop refine loop tuned on PINNED CALIBRATION only, verdicts
  from PINNED HELD-OUT, per-phase 4-part Gmail, auto notebook. `--smoke`/`--force`/
  `--hours`/`--no-email` flags.
- `tests/test_pdefind.py` — 6 known-case code-tests (advection, diffusion,
  2nd-order wave/plasma, forward-predict, bootstrap, noise) — **6/6 PASS**.

**Pre-gates cleared (T9 done):** `formula-validation` agent → **CONFIRM** (STRidge
faithful to Rudy 2017; 2nd-diff / leapfrog / RK4 correct; column↔name mapping
consistent; no sign/off-by-one). 6/6 code-tests. 7 catalogue rows added to
`docs/validation/test-catalogue.md`. Existing kernels 12/12 still pass.

**Real-data smoke (E=100, σ=5, v=2.71, ω_p=3.47 eV):**
- **WP order-2: `u_tt = 7.1·u_xx`** — a WAVE equation, c=√7.1=2.67 ≈ projectile
  velocity v=2.71 (wake propagates at the projectile speed); `u_xx` the ONLY
  bootstrap-stable term; forward-predicts held-out time rel-L2 **0.50** → validated.
- Classical order-2 fails the forward wall (noisier field) → a real classical-vs-WP
  DIFFERENCE, which the 12 h loop will characterise per cell/order.
- Full orchestrator `--smoke` chain ran T8→T14 in ~5 min: T11 classical REFUTE
  (1 cell), **T12 wavepacket CONFIRM**, T13 compare, T14 synthesis. Idempotency
  verified (re-run skips cached). Smoke artifacts cleared for a clean real launch.

**4 real bugs found + fixed via real data** (each a genuine methodological gain):
(1) `subsample_indices` linspace→**uniform stride** (non-uniform dt sawtooth in
u_t); (2) STRidge **normalise the target b too** (1e-5 field else culls all terms);
(3) **nondimensionalise the field** (nonlinear terms had 1e10 coeffs); (4)
forward-score the **admitted** (bootstrap-filtered) eq + **sub-stepped** integrator
(stiff u_xxx blew up explicit integration → false rejections).

**Artefacts written 2026-07-03/04:** CONTEXT glossary; ADR 0012; campaign T8–T14
+ design section; `kernels/pdefind.py`, `kernels/discovery.py`, `orchestrate_pde.py`,
`tests/test_pdefind.py`; 7 test-catalogue rows.

## 2026-07-04 milestone — SCIENTIFIC PANEL OVERTURNS the WP CONFIRM (artifact)

A 4-expert panel (opus, openings→rebuttal→judge) deliberated on the T12 WP
CONFIRM. **Verdict: the "WP wave equation u_tt = c²·u_zz with c≈v" is an
ARTIFACT, ~0.9 confidence, over-determined by two independent causes:**

1. **Data-provenance bug (E4, file-grounded + re-ran the kernel).** The Track-B
   velocity-sweep cells resolve to the NON-`_wf` jellium runs
   (`run_wp_n162_L50_E{20..600}`), where `density_system == density_total` is
   **WP-INCLUDED** (measured ∫n dV = 163 e vs GS 162 → dN = +1). So the axial
   "wake" field is dominated by the σ=5 wavepacket BLOB itself, translating at v.
   Any rigid f(z−vt) satisfies u_tt=v²u_zz *identically* → c≈v is blob kinematics,
   not collective physics. `discovery.load_cell_axial` subtracts only the GS, NOT
   the mandatory rigid-projectile-motion subtraction (campaign guard rail).
   **Proof:** re-running discovery at E=100 on bath-only `_wf` runs (σ3, σ8) gives
   **u_tt = 0** — remove the blob, the wave equation vanishes.
2. **Physics box-limit (E2/E3, independent of the bug).** Wake wavelength
   λ = 2πv/ω_p = 66/130/320 Bohr at the three velocities, ALL ≫ L=50 Bohr; and
   transit ≪ plasma period (49 a.u.). Under one wavelength fits and no steady wake
   forms → even a correctly-subtracted bath gives u_tt≈0. Resolvability ≠ existence
   (f-sum rule guarantees the wake exists, not that it is a findable PDE in L=50 at
   v≫v_F). c/v=0.87 at low v = projectile DECELERATION (v_blob/v_initial), NOT
   Bohm-Gross dispersion (panel withdrew the dispersion reading).

**Classical REFUTE:** real in that no linear wave PDE *should* be found, but cause
CONFOUNDED (genuine nonlinear/high-q screening vs transverse-MEAN dilution of the
sharp on-axis point charge) — NOT attributable without an on-axis line cut (only
plane-collapsed traces were loaded). The T13 "ratio 5" is a cross-velocity artifact
(WP E100 u_xx / classical E25 u_xx) — invalid.

**Co-moving-frame re-discovery (my proposed next test) is INSUFFICIENT** — a rigid
Gaussian blob collapses perfectly in ζ=z−vt "for the wrong reason"; the ∫dN
diagnostic + bath-only subtraction is the real discriminator.

**HARD DATA BOUNDS (E4, from disk):** `_wf` bath-only WP runs exist ONLY at E100
(σ0.5/3/8); NO `density_wp` saved across the velocity sweep → a corrected bath-only
velocity-sweep re-analysis is **IMPOSSIBLE on existing data**. Heavy-ion run and
Lz=150 long-box runs are EMPTY (0 VTIs). The ONE mineable run where λ<L:
`run_plasmon_n162_L50_E3p4_varyv` (v≈0.5, λ≈25 Bohr, 2589 VTIs).

**PANEL'S DECISIVE NEXT TEST:** (zero-cost, existing data) re-mine
`run_plasmon_n162_L50_E3p4_varyv` with (i) bath-only n_total−n_wp (blob removed) +
(ii) an ON-AXIS LINE CUT (not transverse mean); add a Klein-Gordon −ω_p²u term to
the library and test whether it SURVIVES. Resolves both the wake question (λ<L
here) and the classical line-cut-vs-dilution fork. **Real fix needs NEW data:**
HEAVY (m≈1836) constant-velocity projectile, v≈1.0 (λ≈49 Bohr), box L≥300 Bohr,
≥5 plasma periods (≥250 a.u.), `density_wp` written.

**CODE FIXES OWED before any re-run:** (a) `celldb`/`discovery` must resolve
bath-only fields (prefer `_wf`; else `total − wp`); reject WP-included
`density_system` for the sweep. (b) add rigid-projectile-motion subtraction to
`load_cell_axial`. (c) switch axial reduction to an on-axis line cut (or add it
alongside the transverse mean). (d) add optional Klein-Gordon −ω_p²u library term.

**Open questions for the user (from the panel):** (1) does the varyv run store
per-frame n_wp / `_wf` bath? (2) is a HEAVY constant-velocity projectile in scope,
or is the study committed to the light free-Ehrenfest electron (if light-only, a
resolvable collective-wake PDE is likely unreachable → reframe as near-field
screening, not a wake equation)? (3) on-axis line cut vs full-3D residual for the
classical discriminator? (4) make "plasmon resolved vs not" explicit via the
−ω_p²u term surviving?

The kernel (pdefind.py) itself is SOUND (formula-validated, 6/6 tests); the failure
was in DATA SELECTION + REDUCTION, not the discovery math. The panel did its job:
caught a plausible-but-wrong headline before it entered the thesis.

---

## 2026-07-04 — RE-ANALYSIS (user pushed back on the panel) + constructive POD/DMD result

**User pushback:** "I don't think this conclusion is right, try analysing again."
Re-analysed empirically, addressing the panel's own weak point (E4 re-ran with the
BLOB-TUNED config; a config tuned on contaminated data will find nothing in a weak
clean bath).

**Provenance re-verified independently (my measurement):** velocity-sweep run
`run_wp_n162_L50_E100` density_system = 163 e (WP-INCLUDED); `..._sigma3_wf`
total−wp = 162 e (clean, blob=1e). Artifact confirmed a 3rd way.

**Clean-bath discovery, re-tuned FROM SCRATCH (both transverse-mean AND on-axis
line cut):**
- E100 (v=2.71, λ=130≫L): no validated low-order PDE.
- `run_base_n162_L50_E1p5` (v=0.33≈v_F, **λ=16<L=50**, Bragg regime, full
  density_rt_wp): **still no validated PDE.** (Note: older runs use dir names
  `density_rt_total`/`density_rt_wp`/`density_rt_delta`, NOT density_total — the DB
  dir columns are stale for these; VTIs ARE on disk.)
- Conclusion: the governing-DE framing does not survive on clean data in ANY
  testable regime. BUT the induced bath is NOT null — on-axis line cut ptp ~1.2e-3
  (30× the transverse mean). Real structure exists; it is just not a low-order 1-D
  DE (runs <1 plasma period → no restoring dynamics; 1-D axial reduction collapses
  the 3-D wake; coarse grid on the low-v run dz=1.0).

**CONSTRUCTIVE RESULT (user chose: characterise the true 3-D bath via POD/DMD).**
`/tmp/bath_compare.py` → `artifacts/bath_pod_dmd_compare.png`. Blob-free 3-D bath
(n_total−n_wp−GS), classical point vs WP σ=3 at v=2.71:
- **WP: POD rank ~1 (93% energy in mode 1)**, DMD dominant ~7.4 eV (low).
- **Classical: POD rank ~4 (62% in mode 1)**, DMD dominant ~105 eV (high, e-h).
- Interpretation (matches panel E2 f-sum-rule): WP σ low-pass-filters → coherent
  low-rank low-freq collective-like response; point charge couples to all q →
  higher-rank, high-freq incoherent single-particle e-h. **A real, defensible
  classical↔quantum difference in the induced density.**
- Caveats: non-stationary runs (DMD growths >0, classical 105 eV is approximate /
  likely e-h+noise; POD-rank difference is the robust finding). One v, one σ.

**WHERE IT LEFT OFF:** POD/DMD on the true 3-D bath is the right artifact-free tool.
NEXT (if continued): sweep velocity × σ, tabulate POD-rank & DMD-spectrum
classical vs WP, using ONLY runs with density_wp (bath reconstructable) or classical
(bath = total−GS). The `orchestrate_pde.py` PDE-discovery path is DEPRECATED as a
headline (artifact); repurpose to a POD/DMD bath-structure sweep. Code fixes still
owed if any PDE work resumes (bath-only resolver, rigid-motion subtraction,
on-axis line cut option).

## 2026-07-04 — BATH-STRUCTURE SWEEP (POD/DMD, artifact-free) — the real result

`bath_structure_sweep.py` → `artifacts/bathstruct_{sigma,velocity}_sweep.png` +
`bathstruct_summary.json`. POD/DMD on the TRUE blob-free bath across a σ-sweep
(fixed v=2.71) and a classical velocity sweep. Emailed.

**σ-SWEEP @ v=2.71 (headline, all 125³, clean bath) — sharp, monotone:**
| σ_WP | POD rank(90%) | leading-mode energy | DMD dominant |
|---|---|---|---|
| 0 (classical point) | **4** | 0.62 | ~210 eV (high, e-h) |
| 0.5 | **1** | 0.94 | ~10 eV |
| 3 | **1** | 0.93 | ~7 eV |
| 8 | **1** | 0.99 | ~11 eV |

**Interpretation (matches panel E2 f-sum-rule):** the classical POINT charge couples
to all q → higher-rank (4-mode), incoherent, high-frequency single-particle e-h
response. ANY finite-σ WP low-pass-filters → a near-rank-1, highly coherent
(0.92-0.99), low-frequency (collective-scale ~7-11 eV) induced bath. **A clean,
defensible classical↔quantum difference in the induced density.** This is the
constructive replacement for the (artifact) governing-PDE headline.

**CLASSICAL velocity sweep:** POD rank 3→4 (weak rise with v); DMD dominant freq
rises strongly with v (225→613 eV) = e-h continuum energy scaling with projectile
velocity/q. Physically sensible.

**WP velocity points (mixed σ/grid, caveated):** E1p5 (v=0.33,σ5) rank 1 lead 0.92
DMD 1.6 eV (near ω_p scale); E100 (v=2.71,σ3) rank 1. `run_plasmon_..._varyv`
(v=0.5) rank 12 lead 0.25 = OUTLIER (coarse dt=4.0, atypical varyv run) — exclude.

**Caveats:** runs non-stationary → DMD frequencies approximate (POD-rank/coherence
is the robust descriptor); DMD growths >0. Single density (r_s=5.69). To make it
thesis-grade: a matched-σ WP velocity sweep with clean bath (needs new runs with
density_wp written), and repeat at ≥2 densities.

**STATUS: the campaign's real deliverable is now the bath-structure (POD/DMD)
classical↔quantum contrast, NOT a governing PDE.** The PDE-discovery kernel/
orchestrator remain validated + shippable but are DEPRECATED as the headline
(blob-artifact). NOT committed (user owns commits).

---

**(historical) PAUSED — awaiting user GO to launch the full 12 h PDE run.** —
superseded: the PDE run was launched, its result overturned as an artifact, and
the campaign pivoted to the POD/DMD bath-structure sweep above. To launch:
`venv/bin/python3 docs/campaigns/ml-patterns/orchestrate_pde.py` (add `&`/nohup or
background). Est. real runtime ~30–90 min (3 calib + 3 held-out cells × 4 configs ×
2 projectiles, 220 frames, bootstrap 15); the 12 h cap is a safety net. Emails
per phase. NOT launched, NOT committed (user owns commits — two scoped commits owed:
docs = campaign+ADR+CONTEXT+catalogue; code = kernels+orchestrator+tests).

---

## Milestone 2026-07-06 — linear-response residual / form-factor test (panel-chosen)

**User ask:** "conjure up some other technique to figure out the specific
differences between classical and WP induced density." Then: research online
(Floquet suggested) → run the scientific-panel workflow → decide a plan → run
autonomously (hands-off).

**Research + panel (opus×9, default depth).** Verdict transcript in
`.../subagents/workflows/wf_37fcc894-a12`. The panel REJECTED Floquet/Koopman,
HAVOK/Hankel-DMD (just restate POD/DMD coherent-vs-incoherent + sample-starved on
~100 frames), optimal transport, wavelet+transfer-entropy. It also killed every
FREQUENCY-domain method: E4 verified on disk the matched runs are only T≈4.8–17
a.u. (0.1–0.5 plasma periods) → Δω≈10–36 eV ≫ ω_p=3.5 eV, so χ(q,ω)/S(q,ω)/Floquet
quasi-energies are NOT extractable on existing data.

**Chosen technique — time-domain linear-response residual.** Null: n_ind=χ(q,ω)·
V_ext(q,ω); χ is a medium property → identical for both projectiles, so the WP is
a low-pass-filtered point charge, F(q)=exp(−q²σ²/2), and frame-by-frame
n_WP(q,t)=F(q)·n_cl(q,t). The ratio cancels χ in the TIME domain (no ω-bin needed).
d'Alembert-safe (magnitude cancels rigid f(z−vt)). Deeper than POD/DMD.

**Built + validated (this session):**
- `kernels/formfactor_residual.py` (pure numpy): radial_spectrum (3-D |q| shells +
  noise), axial_spectrum (1-D q_z), resample_time, form_factor,
  fit_gaussian_exponent, residual_test, collapse_fork_a. Fit uses the CONTIGUOUS
  DESCENDING ARM of |R| (auto-stops at the noise-floor turnaround).
- `tests/test_formfactor_residual.py`: 7 known-case tests, ALL PASS (exact-F
  recovery, high-q excess detect, t-drift flag, Fork-A σ_WP vs σ_pot selection,
  exponent fit, radial-spectrum localization). Rows added to test-catalogue.
- `linres_residual_test.py`: autonomous runner (streams frames, idempotent per-σ
  JSON, figures, 4-part email, per-pair try/except). CPU, no new runs.
- Plan: `docs/plans/linres-residual-classical-vs-wp.md`.

**RESULT (full run, maxf=100, artifacts/linres_residual_summary.json + 4 figures,
email sent):**
- **σ=0.5 (the ONLY SNR-adequate pair; F keeps signal to q≈4, ~22 shells):** the
  induced-density ratio |R(q)| tracks a Gaussian filter of width **σ_fit≈0.62 —
  near σ_WP=0.5, not σ_pot=0.35** (single-point leans σ_WP), high-q excess ≈0.15σ
  (no significant nonlinear fingerprint in-window). BUT |R(q,t)| is NOT flat over
  the full 4.8 a.u. window (t_flatness=0.43): **the static-linear-filter null is
  REJECTED.** Early-vs-late split: BOTH halves individually flat (0.19 each) while
  the full window is 0.43 → a monotone LEVEL drift, not noise → **DECELERATION**
  (light WP slows under free Ehrenfest; the classical ghost holds v). decel=True.
- **σ=3, σ=8: SNR-DEAD (correctly excluded).** Their form factor e-folds within
  ~1–4 shells (q_e=√2/σ = 0.47 / 0.18 a.u.); the descending-arm fit hits the
  ~5–10% total−wp blob-subtraction floor, giving meaningless σ_fit=1.46 / 1.26.
  Adequacy gate: a ≥ 0.15·σ² (a floored fit under-estimates a).
- **Fork A: INCONCLUSIVE on existing data** (only 1 SNR-adequate σ; need ≥2). The
  √2 trap cannot be resolved from these short/narrow-box runs; what data exist lean
  σ_WP.

**Physical answer to "the specific difference":** instantaneously (fixed v) the WP
induced bath density = the classical point-charge induced density through a Gaussian
low-pass filter F(q)=exp(−q²σ²/2) with σ≈σ_WP (same medium χ, only the coupling
differs — deeper than POD/DMD and d'Alembert-safe); DYNAMICALLY the two diverge
because the light WP decelerates while the classical ghost holds v (measured
light-projectile effect). Whether a genuine quantum/nonlinear high-q fingerprint
also exists, and the exact filter width, need the new run below.

**HELD for user GO (expensive-sim gate — user owns launches):** the panel's single
DECISIVE new run — matched classical+WP pair at v≈1–2, **L≥100 Bohr** (dq≈0.063,
puts q*=ω_p/v on-grid), **≥3 plasma periods ≈150 a.u.** (Δω≈1.1 eV), writing both
`_wf` bath and `density_wp`. Alone answers plasmon-vs-Doppler: is the WP's ~7–11 eV
DMD mode real plasmon coupling (flat 3.6 eV) or the kinematic co-moving Doppler line
ω=q_min·v≈9.3 eV. Spec'd, NOT launched.

**Open (panel Qs for user):** (a) accept the ~70% "linear-response dominates" prior
the panel leaned on, or frame the test to CHALLENGE it? (b) authorize the one new
L≥100 low-v long-time pair? NOT committed (user owns commits).
