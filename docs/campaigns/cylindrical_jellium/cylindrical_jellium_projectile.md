---
# LOCKED FINAL 2026-06-27 — design frozen by the user. Execute AS-IS: re-verify the
# <preflight> block, then run Phases 0→5 in order. Do NOT re-open design decisions
# (geometry, r_s/velocity grids, projectile, phasing) without an explicit user instruction.
id: cylindrical-jellium-projectile
area: cylindrical_jellium
title: "Annular jellium tube — electron projectile down the bore, S(v) vs wall r_s"
status: running
hypothesis: "A charged projectile gliding on-axis down the hollow bore of a PERIODIC annular jellium tube experiences a measurable electronic stopping power S(v) — from wake/image coupling to the wall electron gas — whose low-velocity friction coefficient beta(r_s)=dS/dv|_(v->0) varies monotonically with the wall density r_s, giving a TDDFT handle on the wall-electron coupling that underlies flow-induced current / quantum friction in carbon nanotubes."
handover: docs/handovers/cylindrical-jellium-projectile.md
tasks:
  - { name: "Phase 0 - Grounding & scaffold: pin TDDFT-PENN/hydrovoltaic/image-stopping/charge-even citations + create systems/cylindrical_jellium/ skeleton (non-running prep)", done: true }
  - { name: "Phase 1a - Ingredients pre-gated: annular background_shape (coded) + classical-electron Gaussian UPF (verify-then-regenerate, sigma_pot~0.354, repulsive)", done: true }
  - { name: "Phase 1b - ROBUST GS for r_s={6,4,2}: full validation battery + xz/yz/xy slices of n+ and n EMAILED for user confirmation", done: true }
  - { name: "Phase 2 - Propagation validation: GS stationarity (no-projectile) + single classical-glide smoke gate at r_s=6", done: false }
  - { name: "Phase 3 - Production S(v) sweep: v={0.15,0.30,0.45} x r_s={6,4,2} (9 runs) -> beta(r_s) with uncertainty", done: false }
  - { name: "Phase 4 - Quantum rung at r_s=6: electron WP (sigma_WP=0.5) vs matched classical ghost", done: false }
  - { name: "Phase 5 - Synthesis notebook (.ipynb) + handover/status update", done: false }
blocked_reason: ""
---

# Annular jellium tube — electron projectile down the bore, S(v) vs wall r_s

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository. You are executing this campaign AUTONOMOUSLY, end-to-end, with no
user in the loop. Re-verify the <preflight> block before burning any GPU.
</identity>

<north_star>
Toward a quantum (TDDFT) model of nanotube hydrovoltaics. In systems where a
material layer and water are in close contact, water flow induces an electronic
current and vice versa (the dipole nature of water + wall-electron coupling).
Existing TDDFT treatments do not fully account for quantum effects. Water in
carbon nanotubes emulates such systems; we model the wall as JELLIUM of variable
density r_s ("different materials = different FEGs", the TDDFT-PENN idea). This
campaign is the FIRST rung: a single annular jellium tube with a charged
projectile gliding down the bore, measuring its electronic stopping power S(v)
as a function of wall density r_s. Complexity (off-axis, concentric multi-r_s
Penn layers) is added in LATER campaigns.
</north_star>

<description>
**System.** An annular jellium cylinder: positive background fills the region
between two concentric cylinders (R_in < d < R_out, d=√(x²+y²)) with a HOLLOW
bore and the tube axis ∥ z. The cell is PERIODIC in all directions; along z the
tube is translationally invariant ("infinite tube"). A charged projectile
(an electron — see <resolved_decisions>) glides ON-AXIS (x=y=0) at +z velocity v.

**Mechanism.** The projectile is in vacuum (the bore) and never ploughs through
the gas — it loses energy only by polarising the wall: the induced wall charge
lags the moving projectile, forming a WAKE whose field retards it (the cylindrical
analogue of a charge moving parallel to a metal surface; Echenique–Ritchie image
stopping). Because the tube is infinite/periodic, the drag reaches a STEADY STATE
and S = −dE/dz is a clean per-unit-length quantity. At low v, S ∝ v
(electronic friction); the slope β(r_s)=dS/dv|_(v→0) is the quantity in the
hypothesis and the TDDFT handle on the wall coupling.

**Decision this informs.** Whether a TDDFT jellium-tube model can quantify the
quantum wall-electron coupling (stopping / induced current) underlying nanotube
hydrovoltaics — the first rung toward a quantum model.

**Falsifiability.** SUCCESS: S(v) is measurable and approximately linear at low v
for each r_s, and β(r_s) is monotonic and resolved beyond its uncertainty across
r_s={6,4,2}. FAILURE: S is within noise of zero (no measurable bore coupling at
this geometry), OR β(r_s) shows no resolved r_s dependence — either falsifies the
"measurable, density-dependent wall coupling" claim and is itself a reportable
result.
</description>

<observables_set>
Reuse the ADR-0006 minimal/maximal set (`minimum_observable_set.hpp`) wired by the
existing localised-jellium run templates. Per run:
- **Classical projectile track** (z, v, F_z) every step → S via −dE/dz (PRIMARY
  stopping channel). Energy ledger `energy_total`(t) logged (≤50-step cadence)
  for the ΔE_system regression.
- **Induced density** `density_delta` (raw + coarse) and **wake** structure
  in/around the wall; **density VTI at the 300-frame cadence**
  (WRITE_EVERY=max(1,round(N_STEPS/300))).
- **Integrated current in the wall** as the projectile passes — the hydrovoltaic
  "flow → induced current" signature (existing wired observable).
- WP rung (r_s=6): WP momentum distribution + integrity, in addition to the above.
- **Cylindrically-averaged radial profile** n(d) (d=√(x²+y²)) for GS validation
  and wake visualisation. This is a NOTEBOOK-LEVEL numpy reduction over the VTI
  (load via `inqview.load_vti`, NEVER fftshift — VTIs are physical-order), NOT a
  new pre-gated kernel. IF you instead promote it to a reusable `inqkit`/`inqview`
  kernel, it becomes pre-gated (code-test + formula-validation + catalogue row).
</observables_set>

<resolved_decisions>
All LOCKED 2026-06-27 (user grill). Values + one-line justifications; engine
claims carry source line-refs.

**Geometry (fixed transverse, per-density axial).**
- Annular shape, axis ∥ z, cell periodic. R_in=5, R_out=13 Bohr (8 Bohr wall —
  thick enough for a bulk-like wall response, sized so N(r_s=2)≤260).
- Transverse box L_xy≈40 (R_out=13 ⇒ ≥7 Bohr vacuum margin each side beyond the
  wall; erfc tail w≈1 decays in ~3 Bohr ≪ margin, so wall tails fit and periodic
  transverse images do not overlap). dx=0.5 Bohr. edge_width w≈1.0 Bohr.
- **PERIODIC infinite tube** along z. Per-density L_z sized to ≥ 2× the v=0.45
  wake length λ=2πv/ω_p (ω_p=√(3/r_s³)): **L_z = {48, 28, 10} Bohr for
  r_s={6,4,2}**. Justification: S is per-unit-length on a z-uniform tube, so L_z
  only needs to exceed the wake length (no self-wake periodic-image overlap) and
  does NOT break S(r_s) comparability.
- Engine: `background_shape` has only {slab,sphere,box}
  (`inq-stack/include/inqkit/jellium/localised_background.hpp:52`) — the annulus
  is new (see <new_code>). v_bg injected shape-agnostically via
  `background_perturbation.hpp:61-67`.

**Density sweep + electron counts.**
- r_s = {6, 4, 2} (dilute→dense; emulates different "materials"). n₀=3/(4π r_s³).
- **N ≈ {24, 48, 136}** (= n₀·π(R_out²−R_in²)·L_z, **rounded to the nearest EVEN
  integer**; then set n₀=N/V_annulus so ∫n₊=N EXACTLY — the exact-neutrality
  requirement, `localised_background.hpp:15`). Hard cap **N ≤ ~260** (user).
- r_s=6 is the CHEAP smoke + WP density (run first); r_s=2 is the expensive tail
  (run last).

**Projectile (classical electron, EHRENFEST).**
- Charge −1, mass = m_e (`PROJ_MASS_AMU=1.0/1822.8885`), carried as the fictitious
  "H" species (project idiom; not a positive H⁺).
- **Dynamics = free EHRENFEST** — the projectile moves under the real KS force and
  its velocity evolves self-consistently (genuinely decelerates); it is NOT a
  prescribed fixed-velocity trajectory. On-axis launch (x=y=0), +z, near the −z
  face with a small margin; periodic-z wrap is fine (z-uniform tube).
- **Gaussian erf-smoothed radial potential** (NOT a bare-Coulomb psp), σ labelled =
  σ_WP=0.5 Bohr ⇒ potential width σ_pot=σ_WP/√2≈0.354 (the √2 rule,
  `[[reference_sigma_matching_convention]]`). Candidate asset
  `…/pseudopotentials/electron_gaussian_wpsigma0p5.upf` — BUT the existing
  `electron_gaussian_*` files carry STALE Coulomb headers and CANNOT be trusted by
  filename/header. The exact UPF is therefore SELECTED + VERIFIED (and regenerated
  if it fails) in the Phase-1 projectile-UPF step (see <tasks>). Verified target:
  charge −1; finite REPULSIVE core V(0)>0 (true Gaussian, no Coulomb singularity);
  σ_pot≈0.354; passes the cutoff/aliasing guard.
- Velocity ENV-DRIVEN per run (sv_ladder pattern: `PROJ_V0`/`SV_N_STEPS`/
  `SV_WRITE_EVERY`/`SV_OUT_SUBDIR`), so ONE build per density-GS serves all 3
  velocities (`sv_ladder_L50_sigma0p5.hpp`).

**Run matrix.**
- Velocities: common absolute **v = {0.15, 0.30, 0.45} a.u.** across all r_s
  (direct S(v) overlay; NOTE r_s=6 v_F=0.32 ⇒ v=0.45 is mildly supersonic — the
  top point may leave the linear regime there; fit β from the lower points if so).
- Rung 1 (production): 3 v × 3 r_s = **9 classical runs** → S(v) and β(r_s).
- Rung 2 (quantum check): at **r_s=6**, one **electron WP (σ_WP=0.5, matched KE)**
  vs its **matched classical electron ghost** (same UPF) — 2 runs.

**Propagation.** Real-time TDDFT (LDA), ETRS integrator for the electrons +
**Ehrenfest** ion dynamics for the classical projectile (the projectile is an
Ehrenfest ion, not a clamped trajectory). dt=0.020 a.u. (4× the old Coulomb-forced
0.005; Gaussian smoothing removes the singularity —
`electron_proj_E100_L50_cubic_sigma1.hpp`, `sv_ladder_L50_sigma0p5.hpp`).
N_STEPS sized per run (see <guard_rails>).

**Stopping extraction.** ΔE_system(t)=`energy_total`(t)−`energy_total`(t₀),
discard first **20%** transient, linear-regress vs projectile path s=|z−z₀|;
gradient = S, stderr = uncertainty (`stopping-power-extraction` skill,
continuous-glide method). Report S, β at **2 s.f.** (3 s.f. only for genuine
near-equalities).

**File placement (ADR-0007).** NEW system
`ResearchProject/systems/cylindrical_jellium/`:
`shared_gs/` (3 GS: tube_rs6/4/2), `shared/configs/` (annular Cfg headers,
`#include` jellium base headers), `scripts/<sweep>/` (build-once run.cpp +
dispatcher, adapted from localised_jellium `fullsuite_classical|wp` +
`qsp_phase3/{gs,classical,wp}`), `<sweep>/<run>/` (outputs; logs gitignored),
`hypotheses/<sweep>/` (combined CSVs, build scripts, study .ipynb, tests/).
</resolved_decisions>

<new_code>
**ONLY new code: the annular `background_shape`.** Everything else (classical
electron projectile, Gaussian UPF, WP injection, env-driven velocity sweep,
S-extraction) already exists and is validated — adapt, do not rewrite.

Add to `inq-stack/include/inqkit/jellium/localised_background.hpp`:
- `enum class background_shape { slab, sphere, box, annulus };`
- params: reuse `center` (tube axis through center in x,y), `half_width` → R_out,
  a new `inner_radius` → R_in, `slab_axis` → the tube axis (=2 for z), `edge_width`
  → w. (Keep the change minimal and POD-only for the device lambda.)
- branch (with `d = √((x−c_x)²+(y−c_y)²)` in the plane ⟂ the tube axis, z free):
  `mask = background_mask(d, R_out, w) * (1.0 - background_mask(d, R_in, w));`
  This is `½erfc((d−R_out)/w)·½erfc((R_in−d)/w)` — 1 inside the annulus,
  erfc-softened at BOTH radial edges, uniform along the axis.

**Pre-gate BEFORE any GS/run** (`code-test` + `formula-validation` agent +
`docs/validation/test-catalogue.md` row):
- Neutrality: ∫n₊ = n₀·π(R_out²−R_in²)·L_z to grid tolerance.
- Interior: n₊ = n₀ for R_in+3w < d < R_out−3w.
- Both radial edges: mask = ½ at d=R_in and d=R_out; → 0 for d≪R_in and d≫R_out.
- Axial uniformity: n₊ independent of z.
- `formula-validation` checks the mask-as-implemented against the erfc annulus
  formula, given ONLY formula + source (not the test).
</new_code>

<guard_rails>
**Pilot-first (BLOCKING).** Before the 9-run sweep, run ONE smoke: r_s=6, v=0.30,
short (N_STEPS for ~1.5·L_z path). PASS criteria (all required):
- exit 0; no NaN; energy not complex.
- electron number conserved < 1% over the run.
- total-energy drift bounded (no monotonic blow-up; report the drift).
- **F_z(t) reaches a clean STEADY PLATEAU** and the plateau is stable as the
  projectile crosses the z-period boundary (confirms L_z ≥ wake; no self-wake
  periodic-image artefact). If contaminated, INCREASE L_z for the dilute runs
  (headroom: N=24 at L_z=48 → L_z can rise to ~96, N≈48, still ≤260) and re-smoke.
- projectile velocity drift < ~10% over the post-transient window (else shorten
  the window / flag — the light electron decelerates).
Only on a full PASS proceed to production.

**Pre-launch cutoff/aliasing guard (MANDATORY).** Run `cutoff_guard.py`
(tddft-simulations §2b, `[[reference_cutoff_aliasing_guard]]`) for the Gaussian
projectile at σ_pot≈0.354 (σ_p=1/(√2·σ_WP)=1.41 for σ_WP=0.5 — the √2 trap) and
dx=0.5. BLOCK if aliased tail > 2%.

**N_STEPS sizing per run.** t_total = max(1.5·L_z / v, 5·2π/ω_p) [floor lets the
wake form]; N_STEPS = ceil(t_total/dt), dt=0.020; WRITE_EVERY =
max(1, round(N_STEPS/300)) [300-frame cadence]. Keep the post-transient window
long enough for a stable plateau but short enough for <10% v-drift.

**Boundary (periodic infinite tube).** Transverse: wall tails fit L_xy (verified
in <resolved_decisions>). Axial: periodic, L_z ≥ 2×wake (per-density sizing).
The 4σ/1σ axial launch-stop rule does NOT apply (no single-traversal cap — the
tube is z-periodic); the binding axial constraint is the wake length.

**Abort conditions.** NaN / complex energy → abort that run, log, continue the
sweep. GPU occupied by ANOTHER user → wait / warn (see <preflight> mechanics);
GPU is the default (NVML "driver mismatch" does NOT block compute —
`[[reference_gpu_driver_mismatch]]`; verify via cudaMemGetInfo probe).

**PROVISIONAL caveats.** r_s=6 is a small gas (~24 e) → finite-size/shell effects;
cross-check the β trend against r_s=4. Electron-as-cation-proxy rests on
charge-even S at leading order (Barkas = odd correction) — cite in grounding.
Results PROVISIONAL until the notebook + validation are complete.
</guard_rails>

<tasks>
**Phased complexity ladder.** Execute phases IN ORDER; each phase's gate must
pass before the next. The cheapest density (r_s=6) leads every compute phase.
Flip the matching frontmatter `done` flag and update the handover as each phase
completes. Phase overview:

  P0 Grounding & scaffold (non-running) → P1 Static system + ROBUST GS (first
  running phase) → P2 Propagation validation → P3 Production S(v) sweep →
  P4 Quantum rung → P5 Synthesis.

---

**PHASE 0 — Grounding & scaffold** (light, non-gating, may overlap P1).
- `literature-review`: pin the exact TDDFT-PENN paper (proton in jellium spheres
  of varying r_s + Penn/optical-ELF averaging; candidates arXiv:2505.23396,
  arXiv:1805.01377); hydrovoltaic/quantum-friction core (Kavokine–Bocquet Nature
  602 (2022); PRX 13, 011019 & 011020 (2023); Král–Shapiro PRL 86, 131 (2001));
  cylindrical/parallel image-stopping (Echenique–Ritchie; Arista cylindrical-
  channel if found); the charge-even-S / Barkas leading-order claim → `docs/sources/`.
- Create the `ResearchProject/systems/cylindrical_jellium/` skeleton (ADR-0007).
- *Gate:* citations exist (incl. charge-even justification); folder scaffolded.

**PHASE 1 — Static system + ROBUST ground state** (the first RUNNING phase).

  *1a — Ingredients (pre-gate EVERYTHING before any run).*
   (i) **Annular `background_shape`** (new code) per <new_code>; `code-test` +
       `formula-validation` agent + `test-catalogue.md` row. *Gate:* all four
       background invariants pass (∫n₊=N; flat interior=n₀; both erfc edges;
       z-uniform); formula-validation confirms the mask.
   (ii) **Projectile Gaussian UPF — VERIFY-then-REGENERATE** (the
       `electron_gaussian_*` files carry stale Coulomb headers; do NOT trust the
       filename). Verify the candidate `electron_gaussian_wpsigma0p5.upf`: charge
       −1; plot V(r) and confirm a finite REPULSIVE core (true erf-Gaussian, no
       Coulomb singularity); potential width σ_pot≈0.354 (=σ_WP/√2 for σ_WP=0.5);
       passes the cutoff/aliasing guard (`[[reference_cutoff_aliasing_guard]]`).
       IF any check fails, REGENERATE a clean Gaussian UPF from a documented
       generator and CHECK THE GENERATOR INTO THE REPO (own the provenance).
       *Gate:* a verified −1 Gaussian projectile UPF exists with provenance and
       passes the cutoff guard. (Needed for Phases 2–4, not for the GS.)

  *1b — Robust GS + validation battery + SLICE EMAILS.* Build & converge the GS
  for each density (adapt localised_jellium `qsp_phase3/gs/run.cpp`), order
  **r_s=6 → 4 → 2**. Run the FULL battery per density:
    1. SCF convergence (monotone energy, residual <1e-6 Ha, no charge sloshing).
    2. Electron neutrality ∫n_elec = N (<1%).
    3. **Radial density profile n(d)** (cylindrically averaged, d=√(x²+y²)):
       flat plateau ≈ n₀ in the wall; **Friedel oscillations at BOTH surfaces**
       (inner bore + outer, period ≈ π/k_F); smooth exponential spill-out;
       **small on-axis bore density** (tails only, no pileup — the projectile
       flies there).
    4. Cylindrical symmetry: density angularly isotropic (cubic grid imprints no
       x/y anisotropy).
    5. Energy sanity: E/N plausible; interior energy density vs bulk-jellium LDA
       E/N(r_s); surface energy ≥ 0.
    6. **Grid/box convergence spot-check (r_s=6 only):** GS energy + profile
       stable vs dx (0.5 vs 0.4) and vs transverse vacuum L_xy (40 vs 48).
    7. **Stationarity:** propagate the GS with NO projectile briefly → density &
       energy constant (a true GS is a fixed point of the propagator; bridges P2).
  - **GS SLICE EMAILS (user-requested validation deliverable).** For each
    density, load the GS VTI via `inqview.load_vti` (PHYSICAL order — NEVER
    `np.fft.fftshift`; `[[reference_vti_coordinate_mapping]]`) and render, via the
    canonical theme, **xz (y=0), yz (x=0), and xy (z=mid) slices of BOTH** the
    prescribed background **n₊** and the converged electron density **n** (linear
    + log; shared colorbar where n₊/n are compared). The xy slice shows the
    annulus + hollow bore face-on; xz/yz confirm z-uniformity, the two wall bands,
    Friedel oscillations and spill-out. **Email each density's slices via the
    `email-notifications` skill** (mandatory four-part: hypothesis reminder → what
    was done [GS of the r_s=X annular tube converged; battery results] → what the
    plots show [bore/wall/Friedel/spill-out called out] → conclusion [is the
    intended tube geometry replicated + battery PASS/FAIL]); attach the slice PNGs.
  - *Autonomy note:* the agent PROCEEDS on the NUMERIC battery (the slice email is
    the user's visual record / optional manual checkpoint — the run does not block
    waiting on a reply). Do NOT Read/preview the PNGs yourself
    (`[[feedback_no_image_preview]]`).
  - *Gate:* battery passes for all 3 densities; GS checkpointed to
    `shared_gs/tube_rs{6,4,2}`; slice emails sent.

**PHASE 2 — Propagation validation** (dynamics smoke, r_s=6).
- Single classical glide (r_s=6, v=0.30; adapt `fullsuite_classical/run.cpp`).
  Apply the pilot gate + cutoff guard in <guard_rails> (exit0; no NaN/complex E;
  N conserved <1%; energy drift bounded; clean F_z plateau across the z-wrap;
  v-drift <10%). *Gate:* all criteria PASS (or L_z bumped + re-smoked to PASS).

**PHASE 3 — Production S(v) sweep.**
- v={0.15,0.30,0.45} × r_s={6,4,2} = 9 runs, env-driven v (one build per
  density-GS), dispatcher with per-phase Gmail (`email-notifications`). Order
  r_s=6 → 4 → 2. Extract S per run (ΔE_system 20%-transient regression); fit
  β(r_s)=dS/dv. Per-run `analyse.py` → REPORT.md. *Gate:* S(v) for each r_s with
  β(r_s) ± uncertainty; β(r_s) monotonicity assessed vs its error.

**PHASE 4 — Quantum rung** (r_s=6).
- Electron WP (σ_WP=0.5, matched KE) vs its matched classical electron ghost
  (adapt `fullsuite_wp/run.cpp` + `qsp_phase3/wp`). Compare S and wake. *Gate:*
  WP-vs-classical S compared; quantum effect quantified (charge-even caveat noted).

**PHASE 5 — Synthesis & notebook.**
- Executed `.ipynb` (`notebook-making`/`run-notebook`) in `hypotheses/<sweep>/`:
  S(v) + β(r_s); induced wall current (flow→current signature); wake (radial+axial,
  shared-colorbar, linear+log); WP-vs-classical panel; hydrovoltaic/quantum-
  friction connection. Canonical theme. *Gate:* executed notebook exists; handover
  + frontmatter updated; `status: done`.

(Deferred to LATER campaigns: off-axis radial-offset sweep; positive-cation
variant via a positive Gaussian UPF; concentric multi-r_s Penn layers.)
</tasks>

<rules>
- `inq/` is IMMUTABLE — the annular shape goes in `inqkit/localised_background.hpp`.
- New observable/kernel → code-test + formula-validation + catalogue row BEFORE
  expensive runs (only the annular shape qualifies; the radial profile stays
  notebook-level unless promoted).
- VTIs are PHYSICAL-order — load via `inqview.load_vti`, NEVER `np.fft.fftshift`
  (`[[reference_vti_coordinate_mapping]]`).
- σ is ALWAYS σ_WP=0.5 in any label/axis/caption; σ_pot=0.354 only in a methods
  footnote (`[[reference_sigma_matching_convention]]`).
- GPU is the default; verify via cudaMemGetInfo probe; warn if another user
  occupies it. Use the project venv for all Python.
- Report S, β at 2 s.f. (3 s.f. only for genuine near-equalities).
- Ground every physical claim (wall coupling, Penn averaging, charge-even S) per
  `literature-review`; label inferences "Inference:".
- Per-run `analyse.py` (full inqview pipeline → REPORT.md) for every completed run.
</rules>

<preflight>
Re-verify EVERY box from this prompt alone BEFORE burning GPU; if any fails, STOP
and surface it rather than running.
- [ ] Intent self-contained: falsifiable hypothesis + SUCCESS/FAILURE criteria
      (see <description>); every task has an unambiguous done-criterion.
- [ ] Setup reproducible, zero guessing: geometry (R_in=5, R_out=13, L_xy≈40,
      L_z={48,28,10}, dx=0.5, w≈1); r_s={6,4,2} ⇒ N≈{24,48,136} (even, exact
      neutrality); projectile (classical electron, m_e, **EHRENFEST**, Gaussian
      UPF σ_pot≈0.354 **VERIFIED in Phase 1a** — never trusted by filename);
      v={0.15,0.30,0.45}; dt=0.020 (real-time LDA, ETRS electrons + Ehrenfest ion);
      GS = Phase 1b (validated+checkpointed); file placement = new
      `systems/cylindrical_jellium/` (ADR-0007).
- [ ] New code/assets pre-gated: annular `background_shape` → code-test +
      formula-validation + catalogue row; projectile Gaussian UPF →
      verify-then-regenerate (stale-Coulomb-header trap) — BOTH before runs.
- [ ] Validation & guard rails: BLOCKING pilot (r_s=6, v=0.30) with numeric gate;
      cutoff/aliasing guard (σ_p=1.41, the √2 trap); periodic-tube boundary =
      L_z≥2×wake (NOT the 4σ/1σ rule); abort on NaN/complex E / GPU-occupied;
      PROVISIONAL caveats (r_s=6 small gas; charge-even/Barkas) named.
- [ ] Autonomous mechanics: GPU via cudaMemGetInfo probe (NVML broken; GPU default;
      warn if occupied); env-driven velocity sweep (one build/density); dispatcher
      concurrency + per-phase Gmail; 300-frame VTI cadence; per-run analyse.py;
      auto-built notebook; handover pointer present; agent flips frontmatter
      done/status.
- [ ] Grounding: every scientific/numerical choice cited or labelled "Inference:";
      engine claims carry source line-refs (localised_background.hpp:52/70,
      background_perturbation.hpp:61-67, the localised_jellium run.cpp templates).
</preflight>
</output>
