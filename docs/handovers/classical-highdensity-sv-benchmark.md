# Handover — High-density classical S(v) benchmark (localised jellium slab)

Campaign: `docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md`
(id `classical-highdensity-sv`, area `localised_jellium`, `status: draft`).
Branch: `overnight-gaussian-classical`. Started 2026-07-21 (grill-with-docs).

## Goal (one line)

Produce a clean **classical electronic-stopping benchmark curve S(v)** for a
denser localised jellium slab (r_s≈4.2), CAP-free with a z-open box so the
projectile exits and E_electronic plateaus — the like-for-like reference a later
WP run is compared against. Carry BOTH stopping definitions (E_absorbed/L headline;
full energy-component ledger data-collected for the still-TBD decomposition formula).

## Design LOCKED via grill (2026-07-21, 8 decisions)

1. **Oscillation scope = both.** (a) ΔE>0 rise = diagnosed CAP artifact
   (`energy-oscillation-diagnosis.md`) → **no CAP**. (b) non-plateau =
   stuck/wrapped projectile → **z-open exit + energy conservation**.
2. **Projectile rep = Gaussian charge perturbation** (`periodicity(2)`), ghost UPF
   dropped; **+1 ghost-UPF run kept as a documented failure contrast** only.
3. **Drive = Ehrenfest primary.**
4. **Slab = 25 Bohr thick, 35×35 transverse, N=100 ⇒ r_s≈4.2** (Na-like). Box
   Lz≈80–90 (lock at Phase 0).
5. **Mass = 1 electron (exact WP match), high-velocity-only** (transits a dense
   25-Bohr slab; a mass-1 electron at v=1.3 would STOP inside).
6. **Velocity grid = transit-floor + 5 up** (floor ~1.8–2.0 from pilot;
   ~{2.0,2.5,3.0,3.5,4.0,4.5}). Samples the Bethe tail ABOVE the Lindhard peak
   (peak unreachable for a transiting mass-1 electron). Prior r_s=5.68 points do
   NOT belong on this curve.
7. **σ_WP=0.5** (σ_pot≈0.354), matches all prior WP work.
8. **E_absorbed extraction is trivially clean:** no CAP ⇒ energy conserved ⇒
   post-exit plateau exact; `E_absorbed = E_electronic(plateau) − E_electronic(0)`,
   `S = E_absorbed/25`. No coupling subtraction (the old +68 eV was the mid-transit
   transient). Lindhard/bulk overlay is an eyeball comparison in the notebook,
   **never a gate** (user, 2026-07-21).

## Phase ladder (accepted)

- Manual-gated: **Phase 0 GS** → **Phase 1 vacuum-exit test** → **Phase 2
  single-transit pilot @ v≈2** (proves plateau + finds v_min; ghost-UPF contrast
  here). One-page dashboard per gate; user inspects between rungs.
- Autonomous: **Phase 3** 6-velocity sweep (per-run analyse.py + run-notebook,
  density GIF on top, step-by-step stopping section) → **Phase 4** synthesis
  phase-notebook (S(v) + full component ledger, WP-overlay-ready).

## Key engine facts (verified, file:line)

- `moving_gaussian_projectile_perturbation.hpp` — projectile IS a moving Gaussian
  charge perturbation, NOT a ghost UPF (`proj_dyn/run.cpp:4-8`).
- `inq/src/systems/cell.hpp` — `periodicity(2)` = slab (x,y periodic, z open);
  Gaussian leaves box via finite-grid clip + mixed-BC Poisson, no wrap.
- `projectile.hpp:50-56` velocity-Verlet Ehrenfest; `interaction_energies.hpp`
  emits the 6 pairwise Coulomb terms (one Poisson solve/step).
- `observables_writer.hpp` already supports the full energy decomposition.

## Status (2026-07-21)

DONE: full design grill; campaign file + `<preflight>`; CONTEXT.md glossary
(5 terms); this handover. Campaign is `draft` — NOT yet `ready` (no runs).

NOT DONE (the whole execution): all 6 tasks. Reuse/fork
`scripts/classical_slab_stopping/{run.cpp,orchestrate.py}` +
`scripts/localised_jellium_dynamics/proj_dyn/run.cpp`. Known code gap:
`run_notebook_builder.py` keys off `electron_track.csv`/ion-track columns but
these runs emit `projectile.csv` — adapt it (task 2) before Phase 3.

OPEN (offered, not yet created): **ADR 0013** for the benchmark definition
(mass-1/high-v-only + z-open + CAP-free + E_absorbed/L + σ_WP=0.5) — hard to
reverse, surprising, real trade-off (WP-mass-match vs reaching the Lindhard peak).
Awaiting user confirm.

## Milestone: 2026-07-21 — Phase 0 + Phase 1 LAUNCHED (parallel, one per GPU)

Both GPUs verified free (23.8 GB each; nvidia-smi down/NVML cosmetic, checked via
`shared/bin/gpu-status` + `vacuum/gpu_probe`).

- **Config header WRITTEN:** `shared/configs/slab_n100_L35x35x85.hpp`
  (35×35×85, dx=0.5, N=100, slab half-width 12.5, EDGE_WIDTH=1.0, EXTRA_STATES=24,
  n0=3.2653e-3 ⇒ **r_s=4.18**; projectile fields σ_WP=0.5/σ_pot=0.35355, mass 1).
- **Phase 0 GS (GPU 0):** run.cpp written at
  `scripts/classical_highdensity_sv/gs/run.cpp` (clone of `pbc_open_z/gs/run.cpp`,
  new config, `periodicity(2)`, default `inq/`, checkpoint →
  `shared_gs/slab_n100_L35x35x85_dx0p5_per2/`, GS density VTI dump). Background
  agent building+running via `inq-run` (TMPDIR set), will produce
  `hypotheses/classical_highdensity_sv/phase0_gs/phase0_gs_dashboard.png` +
  `phase0_summary.md` (n(z) profile, xz slice, ∫n dV, symmetry).
- **Phase 1 vacuum-exit (GPU 1):** background agent writing a minimal dump binary
  at `scripts/classical_highdensity_sv/vac_exit/run.cpp` (empty `periodicity(2)`
  box, sweeps a Gaussian charge center across the +42.5 far face, dumps n_proj
  VTIs + `exit_scan.csv`). Tests CLIP-vs-WRAP. Dashboard →
  `hypotheses/classical_highdensity_sv/phase1_vac_exit/`.
- Both are MANUAL GATES: user inspects the dashboards before advancing. Agents do
  NOT flip frontmatter.

## Milestone: 2026-07-21 — Phase 0 + Phase 1 COMPLETE + verification notebooks

- **Phase 0 GS DONE** (GPU 0): run_completed=true, r_s=4.181, E_GS=207.18 Ha,
  num_states=74, ∫n dV=100.000, interior/n0=1.005, spill-out≈1.06 Bohr, symmetry
  4.3e-5 (fitted centre = +dz/2; naive 11.5% is a half-cell sampling artefact).
  Checkpoint `shared_gs/slab_n100_L35x35x85_dx0p5_per2/`.
- **Phase 1 vacuum-exit DONE** (GPU 1): PASS. ∫n_proj=1.00 interior → 0.94(z=42) →
  0.06(z=43) → 1e-105(z=50); **wrap witness = 0 at ALL 53 positions**. Extended the
  binary to also dump φ_proj = poisson(n_proj) (the actual perturbation potential)
  over a full −44→+50 transit sweep (53 n_proj + 53 φ_proj VTIs).
- **Verification notebooks (executed, 0 errors, for the user's manual gate):**
  - `hypotheses/classical_highdensity_sv/phase0_gs/phase0_gs_verify.ipynb`
    (+ `build_gs_notebook.py`): run_summary hard checks, n(z) profile, xz slice,
    ∫n dV, symmetry.
  - `hypotheses/classical_highdensity_sv/phase1_vac_exit/phase1_vac_exit.ipynb`
    (+ `build_vac_exit_notebook.py`): clip/no-wrap + φ_peak curves, on-axis n_proj(z)
    & φ_proj(z) lineouts at several projectile positions, and an embedded GIF
    `perturbation_along_z.gif` of φ_proj(x,z) moving out of the box.
- Both remain MANUAL GATES — user inspects the notebooks and decides. Frontmatter
  NOT yet flipped (awaiting user acceptance of Phase 0 + Phase 1).

## Milestone: 2026-07-21 — periodicity(2) locked explicit; 2 new phases added

- **periodicity(2) confirmed for ALL runs** (GS + dynamics); stated explicitly in
  the campaign `<resolved_decisions>` + `<rules>`.
- **Potential physics clarified (verified from φ data):** φ_proj = poisson(n_proj)
  is a smooth screened-Coulomb well (erf(r/√2σ)/r free-space form), NOT a Gaussian
  (φ(z=3)=0.22 vs Gaussian 1e-16). In periodicity(2) it's the in-cell screened +
  gauge potential (deviates from free-space, goes slightly negative far out —
  charged-cell convention). Perturbation = correct in-cell potential; ghost UPF =
  free-space erf/r tabulation + r_cut aliasing ⇒ UPF approximates the perturbation.
- **New Phase 1b (vac_dynamic) — BUILDING** (GPU 0, `scripts/classical_highdensity_sv/
  vac_dynamic/run.cpp`): REAL propagation, const-v (v=2) Projectile launched z=−30,
  runs to proj_z=+127.5 (= far face +42.5 + Lz=85 beyond), dumps φ_proj VTIs +
  φ_peak(t). Validates the dynamic moving-perturbation tracking. Notebook + GIF to
  follow.
- **New Phase 2 (dynamics/Ehrenfest validation) — LOCKED A+B+C** (user 2026-07-21),
  design in campaign `<dynamics_validation>`: (A) analytic two-Gaussian force vs
  `projectile_force_z` (Catch2 unit test); (B) energy-conserving trajectory vs a
  Python Newton-ODE (velocity-Verlet check); (C) perturbation-vs-pseudopotential
  overlay. NOT yet built.
- Campaign tasks now 8 (was 6); INDEX regenerated (0/8).

## Milestone: 2026-07-21 — Phase 2 Test A PASS; vac_dynamic FAR-macro fix

- **Phase 2 Test A (analytic force) — PASS.** `force_test/run.cpp` (finite
  periodicity(0) box, fixed Gaussian source σ_s=0.5, projectile σ_pot=0.354 swept
  z=1..14) vs the closed-form two-Gaussian force. Results: source_norm=1.0000,
  **median E_num/E_ana=1.0001, median F_num/F_ana=1.0012**, max shape dev 1.0e-2.
  `projectile_force_z` reproduces the analytic Coulomb force to ~0.1%; atomic-unit
  Poisson convention confirmed. Notebook `phase2_force_test/phase2_force_test.ipynb`
  (executed). Tests B (energy-conserving trajectory) + C (perturbation-vs-UPF) NOT
  yet built.
- **GOTCHA (fixed): `FAR` is a reserved macro** (legacy FAR/NEAR, pulled by an INQ
  header) → `<<FAR<<` preprocessed to `<<<<`, breaking the vac_dynamic build.
  Renamed the variable to `ZFAR`. AVOID `FAR`/`NEAR` as identifiers in any run.cpp.
- **Phase 1b (vac_dynamic φ(t)) — DONE.** Real propagation, const-v=2, Projectile
  −30 → **+127.6** (≥ Lz beyond the +42.5 face), v≡2, **fitted dz/dt=2.0000**
  (perturbation tracks the Projectile), φ_peak 2.12 → **0** at exit. 282 frames.
  Notebook `phase1b_vac_dynamic/phase1b_vac_dynamic.ipynb` (executed) with embedded
  `phi_vs_time.gif` (φ_proj moving through the box in time). NOTE: both notebook
  builders needed a run_summary parse fix (multiple `key = val` on one line → use
  regex `grab()`, not split-on-first-'='); fixed in both.
- **Awaiting user gate** on Phase 1b + Phase 2 Test A before Phase 2 B/C.

## Milestone: 2026-07-22 — analytic INQ-native force (user directive: no differences)

User directive: the projectile force must use INQ's EXACT native Ehrenfest HF
formula (analytic gradient), not the finite difference, and be proven to equal the
force INQ computes for a real pseudopotential ion.

- **INQ native Ehrenfest mapped** (Explore, file:line): ions move only with
  `ion_dynamics=EHRENFEST` (`ionic/propagator.hpp:106`); velocity-Verlet
  (`velocity_verlet.hpp`); **local force = −∫(V_long+V_short)·∇n dr** density-gradient
  form (`forces_stress.hpp:182-187`), analytic ∇n via `operations::gradient`; public
  entry `observables::forces_stress{ions,electrons,ham,energy}.forces` (Ha/Bohr).
  Our old `Projectile` matched the *physics* (HF force + velocity-Verlet) but used
  FD for the force and a per-step (not intra-ETRS) ion move.
- **WRITTEN:** `inqkit/dynamics/projectile_force.hpp` now has
  `projectile_force_analytic(density, cell, center, sigma_pot)` (+ `_z`) = INQ's exact
  integrand −∫V_proj·∇n (V_proj=poisson(gaussian_density), analytic gradient, covariant
  reduce → volume_element·to_cartesian). May need compile fixes (delegated).
- **Ghost UPF caveat (from inq source):** existing `electron_gaussian_*.upf` were
  extended with a flipped-sign −1/r tail — NOT clean erf/r. Generate a fresh clean one
  via `inqview.io.gaussian_psp.generate_gaussian_psp` for the native-ion test.
- **LAUNCHED (background agent):** builds 2 tests — (1) analytic force vs closed-form
  two-Gaussian; (2) DECISIVE: `projectile_force_analytic` vs INQ **native**
  `forces_stress` on a real ghost-UPF ion of the same V_loc (same density) → proves
  perturbation-force == pseudopotential-ion-force. Output →
  `hypotheses/classical_highdensity_sv/phase2_native_force/`.
- After this passes: wire analytic force into the production run.cpp fork, then Tests
  B (integrator vs Newton-ODE) + C (perturbation vs native-Ehrenfest trajectory).

## Milestone: 2026-07-22 — Phase 2 native-force PASS, Test A/B PASS, Test C running

- **Native-force test (decisive) — PASS.** `projectile_force_analytic` vs INQ native
  `forces_stress` on a clean ghost-UPF ion: **<0.1% at dx=0.4** (1.5% at dx=0.5,
  collapses on refinement → pure discretization). Sign correct. Verified from raw CSVs
  (`phase2_native_force/force_vs_native_dx0p4.csv`). Clean ghost UPF generated
  (`ghost_sigma0p354.upf`). Force header compiled clean, no fix. Campaign
  resolved_decisions updated: production run.cpp MUST use `projectile_force_analytic`.
- **Test A (closed-form) — PASS:** analytic force median 0.11% vs two-Gaussian (beats
  FD 0.33%). **Test B (integrator) — PASS:** velocity-Verlet vs DOP853 Newton-ODE
  max|Δz|=4.5e-6 Bohr, energy conserved 4.1e-7 (host-only, `integrator_test/`).
- **Test C (LAUNCHED, background):** perturbation (analytic force) vs INQ NATIVE
  Ehrenfest (`ion_dynamics=EHRENFEST`, real ghost ion) — same finite-box He-source
  setup, dx=0.4, compare z(t)/vz(t)/energy → quantifies the remaining intra-step
  ordering difference end-to-end. Output → `hypotheses/classical_highdensity_sv/
  phase2_native_ehrenfest/`. Also checks whether native Ehrenfest even moves a
  z_valence=0 ghost ion (decisive fact if not).
- REMAINING after C: wire `projectile_force_analytic` into the production run.cpp fork
  (task 5), then Phase 3 pilot.

## Milestone: 2026-07-22 — Phase 3 pilot LAUNCHED autonomously (user away, review later)

User: "automatically move onto the next phase, make both run notebooks, I'll review
later today." Two background agents now running concurrently (one GPU each):
- **Test C** (phase2_native_ehrenfest): native-Ehrenfest ghost ion vs perturbation, small He system.
- **Phase 3 pilot** (background agent): builds `scripts/classical_highdensity_sv/pilot/run.cpp`
  (perturbation Ehrenfest, mass-1, v=2, **analytic force**, r_s=4.18 slab, periodicity 2,
  full ledger, ~300 frames, checkpoint) + `pilot_native/run.cpp` (native-Ehrenfest clean
  ghost-UPF ion, same slab). cutoff_guard gate (classical, 54 eV, dx=0.5 → 537 eV PASS).
  Produces TWO run-notebooks (density GIF at top, E-plateau + E_absorbed, stopping-power
  section broken into steps, native-vs-perturbation overlay) → `hypotheses/
  classical_highdensity_sv/{pilot,pilot_native}/`. Robust: Run A + notebook ship even if
  Run B (native) fails. GS is dx=0.5 (perturbation vs native ∇n differ ~1.5% there,
  0.07% at dx=0.4) — fine for the pilot; dx=0.4 flagged for the production sweep.
- Frontmatter task-done flags NOT yet flipped (await pilot results + user review).

## Milestone: 2026-07-22 — Test C DONE; Phase 2 validation COMPLETE (all 4 pass)

**Test C (native Ehrenfest ghost ion vs perturbation) — PASS, verified from raw CSVs.**
Identical shot (mass-1 projectile at a fixed He source, approach-reflect-return, dx=0.4):
- Native `ehrenfest()` DOES move the z_valence=0 ghost ion (z −6→reflect→−5.47).
- Agreement: **max|Δz|=4.68e-3 Bohr (0.11%), max|Δvz|=2.1e-3 (0.17%)**.
- **dt-halving: discrepancy IDENTICAL at dt=0.02 and 0.01 (4.683e-3, ratio 1.00 not 2)**
  ⇒ dt-INDEPENDENT ⇒ NOT the intra-step ordering; it's a fixed ~0.02 eV
  GS-representation offset (perturbation projectile absent from the GS, ghost is in it).
  The true ordering error sits below this floor.
- **Verdict: the perturbation projectile faithfully replicates INQ native Ehrenfest**
  (force = INQ's exact formula, velocity-Verlet, sub-percent trajectory). User's
  "no differences" requirement met to measurable precision.
- Files: `scripts/classical_highdensity_sv/phase2_native_ehrenfest/{c1_native,c2_pert}/run.cpp`;
  `hypotheses/classical_highdensity_sv/phase2_native_ehrenfest/native_ehrenfest_comparison.{md,png}`.

**Phase 2 scorecard (all PASS):** A closed-form 0.11% | native-force <0.1% | B integrator
Δz=4.5e-6 | C full-Ehrenfest 0.11%. Analytic force `projectile_force_analytic` is the
validated production force. Pilot (Phase 3) still running (~2 h); notebooks pending.

## Milestone: 2026-07-22 — Phase 3 pilot Run A DONE; central aim MET; E_absorbed GAUGE FLAW found; SWEEP HELD

Run A (perturbation Ehrenfest, v=2, analytic force, r_s=4.18 slab, 1600 steps) COMPLETE.
Verified from raw CSVs. Notebook `hypotheses/classical_highdensity_sv/pilot/
pilot_run_notebook.ipynb` (executed, density GIF embedded) + `PILOT_SUMMARY.md`.

- **CENTRAL AIM MET:** transit (v 2.0→1.40, proj_z −30→+70.9, did NOT stop) + clean
  z-open exit + **E_total plateaus FLAT** (last-15% std=0.0000 eV, no oscillation).
  Energy conserved (drift 0.87 eV). z-open + CAP-free design WORKS at a real velocity.
- **⚠️ E_absorbed = ΔE_total is CHARGED-CELL GAUGE-BROKEN** (locked decision 8 WRONG):
  ΔE_total=+445.7 eV → S=17.8 (unphysical) because U_proj_bg/E_PS/E_PB swing ±418.8 eV
  as the projectile moves, while −ΔKE_proj=27.8 eV. (Predicted by
  `reference_charged_cell_hartree_convention`.)
- **Gauge-clean classical S = 0.93 eV/Bohr** (KE loss across slab 23.1 eV / 25; in-slab
  −dKE/ds=0.98; mean v=1.80). Between Lindhard-point 0.57 and bulk σ=0.5 0.94 — sensible.
  But KE-loss is NOT usable for the WP → the classical↔quantum shared metric is broken.
- **SWEEP HELD.** Decision needed (user): define a gauge-clean E_absorbed (slab-only
  electronic energy via the pairwise `interactions.csv` ledger, OR KE-loss S with a
  documented caveat, OR neutralise the cell). Do NOT run the 6-velocity sweep until
  resolved — it would produce gauge-broken headline S.
- **Run B (native ghost) DONE — DISQUALIFIED for this geometry** (verified from
  native.csv). Native Ehrenfest DOES move the z_valence=0 ghost, BUT the ghost's
  unscreened 1/r tail + INQ's ion-force seeing slab electrons-but-not-the-background-
  perturbation gives a spurious static vacuum force: ghost stalls ~14 Bohr short of the
  slab (z −30 → oscillates −58.5..−26.2, never reaches near face −12.5), reverses to
  vz=−2.55, E_total swings 229 eV (not conserved). ⇒ **the perturbation projectile
  (Run A) is REQUIRED**; native ghost unusable in the jellium-background-as-perturbation
  setup. (Consistent with Test C, where the compatible finite-box/no-background setup
  gave native≈perturbation to 0.11% — a setup incompatibility, not a native flaw.)
- **BOTH run-notebooks DONE (executed, density GIFs embedded):**
  `hypotheses/classical_highdensity_sv/pilot/pilot_run_notebook.ipynb` (Run A: central
  aim + gauge diagnosis + clean S=0.93) and `pilot_native/pilot_native_run_notebook.ipynb`
  (Run B: disqualified + native-vs-perturbation z(t) overlay). PILOT_SUMMARY.md complete.
- **Awaiting user:** review both notebooks + decide the gauge-clean E_absorbed (option 1
  recommended: slab-only energy via the pairwise ledger). Sweep held.

## Milestone: 2026-07-22 — gauge "flaw" RESOLVED (baseline error); Definition 2 works; sweep UNBLOCKED

Triggered by the user asking to plot the raw `energy_total`. The +445 eV "gauge-broken
E_absorbed" was a **BASELINE ERROR**, not a broken metric. The charged-cell gauge only
contaminates E_total WHILE the projectile is in the box; both clean endpoints are
NEUTRAL-cell:
- **E_absorbed = E_total(plateau=208.1715 Ha) − E_GS(207.1832 Ha) = 0.988 Ha = 26.9 eV**
- **S = E_absorbed / 25 = 1.08 eV/Bohr** (gauge-clean), matches −ΔKE_proj=27.8 eV (within
  0.86 eV drift) and in-slab KE-loss S=0.93.
- The wrong baseline was E_total(0) (projectile at z=−30 in a CHARGED cell, U_proj_bg/E_PS
  ±419 eV → the 445 eV artifact).
- **Definition 2 = [E_total(plateau) − E_GS]/L_slab, read after full exit — WP-transferable**
  (WP plateau−E_GS is also neutral-cell). Classical↔quantum shared metric INTACT.
- **SWEEP UNBLOCKED.** Pilot notebook corrected (raw energy_total plot added + baseline
  resolution + Def-2 headline S=1.08). Locked-decision-8 refined: E_absorbed is clean with
  the GS baseline (not the t=0 baseline).
- Next: run the 6-velocity sweep (S via GS-baseline plateau; carry pairwise ledger);
  dx=0.4 for the fast points (non-blocking refinement). Awaiting user go.

## Milestone: 2026-07-22 — abrupt-exit diagnosis + DIRECT-potential perturbation (A/B running)

User asked: is charge introduced or just potential? Diagnosis (verified from pilot data):
the abrupt E_total change at exit is `norm_proj` (in-cell Gaussian CHARGE) dropping 1→0 as
the narrow σ=0.35 Gaussian clips at the +42.5 face (proj_z 41.3→44); per-step ΔE_total
spikes −0.77→+40.6 eV/step there, then flat. Cause: the perturbation builds a charge
density and `poisson::solve`s it — the periodic Poisson carries an implicit G=0 neutralizing
background (+ transverse images) whose offset lurches as the clipped charge crosses the box.
`energy_total` IS raw INQ `data.energy().total()` (real_time_session.hpp:57).

- **NEW CODE (add potential directly, no charge/Poisson/background):**
  `inqkit/jellium/gaussian_potential.hpp` (direct erf(|r-R|/√2σ)/|r-R| field),
  `inqkit/dynamics/moving_gaussian_projectile_potential.hpp` (perturbation),
  `projectile_force_direct[_z]` in `projectile_force.hpp` (matching force). Gradient (force)
  ≈ Poisson version; only the charge-dependent background offset is removed.
- **A/B run LAUNCHED (background):** `scripts/classical_highdensity_sv/pilot_direct/run.cpp`
  = pilot clone with 2 swaps (direct perturbation + direct force). Compares to the Poisson
  pilot: is the exit transient gone? does S change? Output →
  `hypotheses/classical_highdensity_sv/pilot_direct/` (A/B notebook + DIRECT_SUMMARY.md).
  Poisson pilot ref: S_def2=1.08, S_KEloss=0.93.
- Expected: near-field dominates → S ≈ same; far-field transverse-image differences may shift
  slightly; exit transient should vanish (no charge to clip). Decide adoption from the A/B.

## Resume recipe

Read this handover + the campaign file. The design spine is fully locked; next
concrete step is **Phase 0 GS** (manual gate) — build the GS at
35×35×Lz, 25-slab, N=100, r_s≈4.2, `periodicity(2)`, then the one-page dashboard.
Do NOT advance a gate without user inspection (manual phases). `cutoff_guard.py`
before any dynamic launch. Every INQ build sets
`TMPDIR=/local/data/public/skcb2/tddft/.build_tmp`.

## Milestone: 2026-07-22 — Phase 3 pilot EXECUTING (autonomous, user away)

Autonomous Phase-3 pilot agent. Both runs built + launched detached on **GPU 0**
(Test C occupied GPU 1 initially — used the other card per instruction).

- **STEP 1 cutoff_guard — PASS.** `--spacing 0.5 --kind classical --energy-ev 54`
  → E_cut=537 eV ≥ 1.10·E_kin=59 eV. Cleared for launch.
- **Ghost UPF verified by DATA (not header):** `force_vs_native/ghost_sigma0p354.upf`
  is a CLEAN r·V→+2.0 (=+1/r in Ha after UPF Rydberg convention), z_valence=0 — the
  validated clean erf/r ghost (header text is stale antiproton boilerplate; ignore it).
  Copied into `pilot_native/` for Run B.
- **RUN A (primary) — `scripts/classical_highdensity_sv/pilot/run.cpp`:** clone of
  `classical_slab_stopping/run.cpp`, config→slab_n100 (35×35×85, HALF=12.5, N=100,
  EDGE_W=1.0, dx=0.5, per2), GS=`shared_gs/slab_n100_L35x35x85_dx0p5_per2`. Force line
  REPLACED with the INQ-native analytic HF: `Fz = projectile_force_analytic_z(density,
  cell, center, σ_pot) − projectile_force_analytic_z(nplus, cell, center, σ_pot)`.
  mass 1, charge −1, σ_pot=0.35355, launch_z=−30, K0=2 (v=2), Ehrenfest,
  N_STEPS=1600, dt=0.04, SAVE_EVERY=5 (~320 frames), final checkpoint+resume. Compiled
  clean (only unused-var warning). **RUNNING** on GPU 0 — step ~21 at ~3.9 s/step
  (slowed by Run B compile sharing the card); energies finite, current_z building, no
  NaN. Full ledger (observables/projectile/interactions.csv) + conservation column
  E_elec+KE_proj+U_proj_bg wired.
- **RUN B (contrast) — `scripts/classical_highdensity_sv/pilot_native/run.cpp`:** real
  ghost-UPF ion (mass 1 via `.mass(1/1822.8885)`), native Ehrenfest
  (`options::real_time{}.ehrenfest()`), SAME slab (slab as background perturbation +
  ghost ion), OWN GS recomputed with ghost present (pert GS can't be reused). launch_z
  =−30, v=2, N_STEPS=1600, dt=0.04, density frames. Records z(t)/vz(t)+ledger; detects
  whether native Ehrenfest moves a z_valence=0 ion (decisive fact if not). **BUILDING**
  (its own INQ tree) on GPU 0; will start GS+RT after compile.
- **Notebook builder — `hypotheses/classical_highdensity_sv/build_pilot_notebook.py`**
  (parametrized `{pilot|pilot_native}`): density-GIF at TOP (total + induced Δn, via
  `_slice_stack`/`_save_gif`, load_vti, base64-embedded), z(t)/vz(t), E_elec plateau +
  E_absorbed, pairwise ledger + conservation flat-check, stopping-power section broken
  into Method-1 (E_abs/L_slab, L=25) + Method-2 (initial-drag −dKE/ds, vz≥0.85·v0
  windows), Lindhard bulk eyeball (NON-gating), and the native-vs-perturbation z/vz
  overlay in the pilot_native notebook. Executes via nbclient kernel `inqview-venv`.
- **Robustness:** Run A + its notebook ship independently of Run B. If native Ehrenfest
  can't move the z_valence=0 ghost, Run B records that and exits 0; pilot_native notebook
  still builds (overlay cell degrades gracefully).
- Frontmatter task-done flags NOT flipped (left for main session per instruction).

**Next (this agent):** on Run A completion → build+execute pilot notebook, verify
transit+plateau+S. On Run B completion → build+execute pilot_native notebook + overlay.
Write `hypotheses/classical_highdensity_sv/pilot/PILOT_SUMMARY.md`.

## Milestone 2026-07-22 — DIRECT-potential A/B test (RUN A')

**Goal.** A/B-test whether the Poisson exit transient is an artifact of the
charge/Poisson/neutralizing-background construction by re-running the Phase-3 pilot
with a DIRECT free-space erf/r projectile potential (new inqkit headers by the user).

**RUN A' — `scripts/classical_highdensity_sv/pilot_direct/run.cpp`.** Clone of
`pilot/run.cpp` with EXACTLY two functional changes (verified by diff of code-only lines):
  (a) perturbation → `inqkit::dynamics::moving_gaussian_projectile_potential proj_pert(proj,
      SIGMA_POT)` (adds `gaussian_potential` = erf(|r−R|/(√2σ))/|r−R| directly; no charge,
      no Poisson, no G=0 background) instead of `moving_gaussian_projectile_perturbation`.
  (b) BOTH Ehrenfest force calls → `projectile_force_direct_z(...)` (same HF integrand
      −∫V_proj·∇n but V_proj = the direct erf/r field), keeping the
      `force_direct(electrons.density) − force_direct(background)` structure.
Everything else identical: config slab_n100_L35x35x85, GS
`shared_gs/slab_n100_L35x35x85_dx0p5_per2`, per2, mass 1, v=2/K0=2, launch_z=−30,
N_STEPS=1600, dt=0.04, SAVE_EVERY=5, full ledger + checkpoint/resume.

**Build.** `inq-run` (TMPDIR=/local/.../.build_tmp), fresh CMake tree, GPU sm_80.
The new inqkit headers compiled CLEAN (only an unused-var warning for HA) — no inq/
edit needed. Launched detached (setsid nohup) on **GPU 0** (both cards free at launch).

**Poisson reference (measured from completed `pilot`, VALIDATED analysis logic):**
  - Exit transient: smooth ~−0.8 eV/step until proj_z≈41.5, then a sharp POSITIVE spike
    **peaking at +40.62 eV/step at proj_z=42.51** (= Lz/2, the far BOX face, NOT the slab
    face) — the G=0-background lurch as the clipped periodic charge crosses the box edge.
    E_total jumps ~+21 eV (186.9→208.2 Ha) over ~50 steps, then dead-flat plateau.
  - S_def2 (plateau−E_GS)/25 = **1.076 eV/Bohr** (E_GS=207.18322156141 Ha); S_KEloss
    (−ΔKE_proj across slab)/25 = **0.925 eV/Bohr**. Both match the prompt's stated Poisson
    values (1.08 / 0.93 eV/Bohr).
  - **Baseline subtlety (KEY):** Poisson E_total(0)=191.79 Ha (projectile at z=−30 already
    imposes the −background offset) ≠ clean-slab E_GS=207.18. DIRECT E_total(0)=210.33 Ha
    (energy_external +18.54 Ha=+504 eV higher — the direct erf/r carries NO neutralizing
    background). So `plateau−E_total(0)` is offset-contaminated for POISSON (gives 17.8
    eV/Bohr, spurious) but is the physically correct measure for DIRECT (no offset swing).
    Report S_def2 per-run: DIRECT uses own baseline, POISSON uses E_GS. S_KEloss is
    offset-free for both → the cleanest cross-run comparison.

**A/B notebook builder — `hypotheses/classical_highdensity_sv/pilot_direct/build_ab_notebook.py`**
(executes to `pilot_direct_ab_notebook.ipynb`, kernel `inqview-venv` = explicit venv python,
newly registered): DIRECT density GIF at top (base64), per-step ΔE_total vs proj_z overlay
(headline transient check, full + far-face zoom), raw E_total(t) overlay (+own-t0 shape),
S table + bar chart (S_def2 own/E_GS + S_KEloss for both), z(t)/vz(t) overlay with
max/RMS agreement. Analysis functions pre-validated against the Poisson run (peak +40.62
@ 42.51, S_def2 1.076, S_KEloss 0.925 all reproduced). GIF helper smoke-tested on DIRECT
frames.

**Status at this write:** RUN A' RUNNING on GPU 0, ~step 300/1600, proj_z≈−6 (inside slab),
energies finite & smooth, no NaN, frames writing (~2.3 s/step → ~62 min total). Completion
monitor armed (exits on `run_completed=true` or proc death). PENDING: run completion →
execute A/B notebook → `pilot_direct/DIRECT_SUMMARY.md` (transient gone? S_direct vs S_poisson?
trajectory agreement? adopt-for-sweep verdict). Campaign frontmatter NOT flipped.
