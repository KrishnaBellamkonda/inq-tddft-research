---
id: classical-highdensity-sv
area: localised_jellium
title: High-density classical S(v) benchmark for the localised jellium slab (z-open, CAP-free)
status: draft
hypothesis: "A mass-1 classical electron (Gaussian charge, σ_WP=0.5) launched at high velocity through a denser 25-Bohr jellium slab (r_s≈4.2, N=100) transits cleanly under z-open (periodicity 2) boundaries with NO CAP, so E_electronic plateaus after the projectile exits; the deposited energy E_absorbed=E_electronic(plateau)−E_electronic(0) yields a clean stopping-power benchmark S(v)=E_absorbed/L_slab over a 6-point high-velocity grid, with the full pairwise energy ledger (E_PP,E_PS,E_SS,E_SB,E_PB,E_BB,KE,E_xc) recorded so the (still-TBD) energy-decomposition stopping definition is computable post-hoc — establishing the classical baseline a later WP run is compared against."
handover: docs/handovers/classical-highdensity-sv-benchmark.md
tasks:
  - { name: "Phase 0 — GS 35×35×~85, 25-Bohr slab, N=100, r_s≈4.2, periodicity(2); density profile + energy validated (MANUAL GATE)", done: false }
  - { name: "Phase 1a — vacuum exit test (static): Gaussian charge clipped at z-open face, no wrap (MANUAL GATE)", done: false }
  - { name: "Phase 1b — vacuum DYNAMIC exit: real propagation, moving perturbation tracks the Projectile per step, φ_proj(t) leaves box >= Lz beyond face (MANUAL GATE)", done: false }
  - { name: "Phase 2 — dynamics/Ehrenfest validation: independent analytic force test (two-Gaussian) + energy-conservation trajectory + perturbation-vs-pseudopotential comparison (one approximates the other) (MANUAL GATE)", done: false }
  - { name: "Code — run.cpp periodicity(2)+mass-1+full pairwise ledger emit; run_notebook_builder reads projectile.csv + step-by-step stopping section (code-test)", done: false }
  - { name: "Phase 3 — single-transit pilot at v≈2 (transit-floor find): transit + clean plateau + E_absorbed extractable + full ledger + cutoff_guard PASS; ghost-UPF failure contrast (MANUAL GATE = central aim)", done: false }
  - { name: "Phase 4 — autonomous 6-velocity sweep (floor+5 up); per-run analyse.py + run-notebook (density GIF top, stopping section)", done: false }
  - { name: "Phase 5 — synthesis phase-notebook: S(v)=E_absorbed/L + full component ledger staged for Def-1; Lindhard/bulk eyeball overlay (NON-gating); WP-overlay-ready", done: false }
blocked_reason: ""
---

# High-density classical S(v) benchmark for the localised jellium slab

<identity>
You are a scientific computing researcher running first-principles rt-TDDFT in INQ.
You adhere to this repository's rules, skills, and workflows (tddft-simulations,
stopping-power-extraction, twin-run-analysis, run-notebook, notebook-making,
code-test, simulation-validation, handover-update, scientific-grounding).
σ-convention is UNIFIED: σ means the wavepacket σ_WP; the classical Gaussian
charge std = σ_pot = σ_WP/√2 (`.claude/rules/sigma-wp-convention.md`, CONTEXT.md).
</identity>

<description>
This campaign establishes a **classical electronic-stopping benchmark curve S(v)**
for a denser localised jellium slab, as the like-for-like reference a later
wavepacket (WP) run is compared against. It exists because prior classical
baselines were either bulk/point-charge (ADR 0010 geometry mismatch) or at the
lower r_s=5.68 density with periodic-wrap contamination; this campaign fixes the
projectile-exit and boundary problems and raises the density.

Two stopping-power **definitions** are carried through every run (both computed
classically here, to benchmark both quantum definitions later):

- **Definition 2 — localised-slab energy deposit (HEADLINE):**
  `S(v) = E_absorbed / L_slab`, with `E_absorbed = E_electronic(plateau after
  projectile fully exits) − E_electronic(0)` and `L_slab = 25 Bohr`
  (stopping-power-extraction skill, Method B; Correa 2018).
- **Definition 1 — energy-component decomposition (DATA-COLLECT ONLY, formula
  TBD):** the user is still deriving the closed stopping formula from the pairwise
  Coulomb + kinetic + xc components. This campaign does NOT lock that formula; it
  guarantees every run **emits the full ledger at high cadence** (E_PP, E_PS,
  E_SS, E_SB, E_PB, E_BB, KE_total, E_xc, E_hartree, E_external) so the formula is
  computable post-hoc without re-running.

**Central aim / success criterion:** at the Phase-2 pilot, a mass-1 classical
electron at high v **transits the slab, fully exits, and E_electronic plateaus**
(no CAP, energy conserved) → E_absorbed is a well-defined constant → a clean
S(v). If the pilot plateaus, the campaign's core problem (the historical
non-plateau / stuck-projectile / wrap-around) is solved and the 6-velocity sweep
proceeds. **Failure** = no plateau (projectile stops inside, or wake reaches the
z-walls and re-couples) at any transiting velocity.

**Two oscillation problems, both resolved by construction (not by density):**
1. The ΔE_total>0 "oscillation" was diagnosed as a **non-Hermitian CAP ledger
   artifact** (handover `energy-oscillation-diagnosis.md`, 2026-07-13, conf 0.90;
   CAP-off ⇒ conserved to <0.02 eV). → **No CAP is used.**
2. The non-plateau was the projectile getting **stuck at / wrapping around** the
   box edge. → **z-open `periodicity(2)`** lets the Gaussian potential leave the
   box; energy conservation (no CAP) then makes the post-exit plateau exact.
   Higher density (r_s≈4.2) is kept only for **wake localisation** (shorter
   plasmon wavelength → excitation stays near the slab), NOT as an oscillation
   fix.
</description>

<observables_set>
Per run (all velocities), high cadence, ADR-0006 set + full decomposition:
- **Total & components (every step or every ~2 steps):** `energy_total`,
  `energy_kinetic`, `energy_hartree`, `energy_xc`, `energy_external`
  (+ `energy_nonlocal`, `energy_ion` if present). Emitted by
  `observables_writer.hpp` (already supports full decomposition).
- **Pairwise Coulomb ledger (every step):** `E_PP, E_PS, E_SS, E_SB, E_PB, E_BB`
  from `interaction_energies.hpp` (one Poisson solve/step) → `interactions.csv`.
  Required for Definition-1 post-hoc formula.
- **Projectile track (every step):** `projectile.csv` cols
  `step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal`.
- **Density VTIs:** ~300 frames per run (boundary-rule cadence) for the mandatory
  density GIF (mid-y xz field, physical order — never fftshift a VTI).
- **Conservation check column:** `E_electronic + energy_proj_ke + U_proj_bg`
  (should be flat while projectile is in-box; correctness cross-check).
NO new observable/kernel is introduced — all exist. (If Definition-1's eventual
formula needs a new kernel, that is a SEPARATE code-test-gated task, not here.)
</observables_set>

<resolved_decisions>
All locked via grill-with-docs (2026-07-21). Engine facts carry file:line refs.

**Geometry / system.**
- Transverse box **35×35 Bohr**, slab **thickness 25 Bohr** (axis z), **N=100**
  ⇒ slab volume 30,625 Bohr³ ⇒ n_e=0.00327 ⇒ **r_s≈4.2** (≈Na; a real bump from
  the r_s=5.68 baseline). v_F=1.92/r_s≈0.46.
- Box **z-length ≈ 80–90 Bohr** (slab 25 + launch vacuum + exit vacuum + wake
  room; the old 25-slab runs used Lz≈90) — lock the exact value at Phase 0 so the
  wake never reaches the z-faces at the fastest v.
- **Boundary: `periodicity(2)` for EVERY run in this campaign (GS + all dynamics)**
  — x,y periodic (infinite slab), **z open/finite** (`inq/src/systems/cell.hpp`:
  periodicity 2 = slab; `.periodicity(2)`). This is the single setting that lets
  the moving Gaussian charge LEAVE the box (finite-grid clip + mixed-BC Poisson),
  NO wraparound — the fix for the old stuck/wrapped-projectile non-plateau.
  Verified in isolation by Phase 1a (vacuum static: wrap witness = 0 at all 53
  positions; norm clips 1→0 across the +42.5 face). The RT and the GS it loads
  MUST share periodicity(2) (a per-2 RT on a per-3 GS gives a spurious t=0 kick).
- **Projectile POTENTIAL (clarification, verified).** φ_proj = poisson(n_proj) is
  a smooth **screened-Coulomb well** (finite at centre, erf(r/(√2σ))/r form in free
  space), NOT a Gaussian. In the periodicity(2) cell it is the in-cell potential
  (laterally screened + G=0 gauge; deviates from free-space erf/r, goes slightly
  negative far out — the charged-cell convention). The perturbation computes THIS
  correct in-cell potential each step; a ghost UPF instead tabulates the free-space
  erf/r and suffers r_cut aliasing ⇒ **the UPF approximates the perturbation**
  (quantified in Phase 2).
- **NO CAP.** Avoids the diagnosed non-Hermitian ledger artifact; energy
  conserved ⇒ exact post-exit plateau.

**Projectile.**
- Rigid **Gaussian charge −1**, realised as the moving perturbation
  `inqkit::dynamics::moving_gaussian_projectile_perturbation` (NOT a ghost UPF;
  `moving_gaussian_projectile_perturbation.hpp`). Width **σ_WP=0.5 ⇒ σ_pot≈0.354**
  (sigma-wp-convention).
- **Mass = 1** (classical electron, exact WP mass-match).
- **Drive = Ehrenfest** (self-consistent Hellmann–Feynman force, velocity-Verlet;
  `projectile.hpp:50-56`). Mass-1 ⇒ decelerates, so the velocity grid is
  restricted to the **high-velocity transit window** (below).
- **Force = INQ-NATIVE analytic gradient** (user directive 2026-07-22, "no
  differences"): the projectile force uses `inqkit::dynamics::projectile_force_analytic`
  = INQ's exact local HF integrand −∫V_proj·∇n (`forces_stress.hpp:182-187`), NOT the
  finite difference. VALIDATED: matches INQ's native `forces_stress` on a ghost-UPF ion
  of the same V_loc to **<0.1%** at dx=0.4 (Phase 2 native-force test). The production
  run.cpp fork MUST call `projectile_force_analytic` (not the FD `projectile_force_z`).
- **+1 ghost-UPF run** at the pilot velocity as a DOCUMENTED FAILURE contrast
  (wrap/aliasing) — an evidence artifact, NOT a benchmark point. (Reintroduces the
  `electron_gaussian_*.upf` r_cut-aliasing pathology; memory
  `reference_ghost_upf_tail_aliasing`.)

**Velocity grid (6 points).**
- **Transit floor** v_min found by the Phase-2 pilot (est. ~1.8–2.0 at r_s≈4.2:
  a mass-1 electron must not stop inside 25 Bohr). Grid = floor + 5 upward, e.g.
  ~{2.0, 2.5, 3.0, 3.5, 4.0, 4.5}.
- This samples the **Bethe-like high-velocity tail** (above the Lindhard peak at
  v≈0.5–1); the peak is out of reach for a transiting mass-1 electron — a known,
  accepted limitation. Push v_min as low as transit allows.
- Prior r_s=5.68 datapoints are **NOT** on this curve (different density) — fresh
  curve.

**Propagator / duration.**
- ETRS or CN per prior localised-jellium practice; **dt=0.04** (validate energy
  conservation per velocity). Real-valued (no complex CAP) ⇒ default `inq/`.
- N_STEPS per velocity sized so the projectile center reaches ≳2 Bohr past the far
  slab/box edge (full exit ⇒ potential contribution →0) THEN a short plateau tail.
  Faster v ⇒ fewer steps; size per run, not a fixed count.
- **dx=0.5** provisional; σ_pot≈0.354 and momentum transfers ~2v at the fast end
  can exceed the k=0.5 Nyquist ⇒ **`cutoff_guard.py` is a mandatory per-velocity
  pre-launch gate**; drop dx to ~0.35–0.4 for the fastest points if it fails.

**Extraction (Definition 2, headline).**
- `E_absorbed = E_electronic(plateau) − E_electronic(0)`, plateau read AFTER the
  projectile Gaussian is fully out of the box. No coupling subtraction needed
  (energy conserved; the mid-transit +68 eV transient of the old baseline is not
  the deposit). `S(v) = E_absorbed / 25`. Report with mean in-slab velocity.
- Secondary classical cross-check (NOT the headline, valid classically only):
  initial-drag `S(v0) = −d(KE_proj)/ds` over the early v≥0.85·v0 window
  (light-projectile rule). Never used for the WP comparison (memory
  `feedback_quantum_stopping_not_from_projectile_ke`).

**File placement (ADR-0007; jellium grandfathered-flat but this sweep is foldered).**
- Scripts: `ResearchProject/systems/localised_jellium/scripts/classical_highdensity_sv/`
  (run.cpp forked from `classical_slab_stopping/run.cpp` + `scripts/localised_jellium_dynamics/proj_dyn/run.cpp`; orchestrate.py forked).
- GS: `shared_gs/slab_n100_L35x35x85_rs4p2_per2/` (name at Phase 0).
- Runs (outputs only, logs gitignored):
  `classical_highdensity_sv/{vac_exit,pilot_vXpY,ghostupf_vXpY,v0..v5}/`.
- Analysis: `hypotheses/classical_highdensity_sv/` (combined CSVs, build_*.py,
  study `.ipynb`, `tests/`).
</resolved_decisions>

<guard_rails>
- **Correctness-only hard gates (BLOCK):** NaN / complex / non-finite energy;
  missing or unvalidated GS; `cutoff_guard.py` aliased-tail BLOCK; energy NOT
  conserved with the projectile in-box (conservation column drifts) — indicates a
  numerics/dx problem, not physics.
- **Central-aim gate (Phase 2, numeric):** projectile transits (proj_z crosses
  far slab face with v>0) AND E_electronic plateaus after exit (|dE/dt| in the
  final 15% < 5% of E_absorbed) AND finite S. NOT gated on velocity drift (mass-1
  is supposed to decelerate) and NOT gated on agreement with Lindhard (comparison
  only, never a block — memory `feedback_fourier_loss_function_gate`).
- **Wake-localisation check (not a hard gate, informs box size):** excited density
  at |z| near the box faces ≈ 0 at run end; if not, extend Lz (never add a CAP).
- **Boundary/cadence:** ~300-frame VTI cadence; density GIF (mid-y xz, physical
  order, LINEAR|LOG, slab faces marked) embedded at the TOP of every run notebook
  (`.claude/rules/notebook-density-gif.md`).
- **Checkpoint-don't-block:** every run final-checkpointed + resumable
  (`*_RESUME=1`); on projected budget overrun, WARN (email) + proceed, never
  self-block (`.claude/rules/checkpoint-dont-block.md`).
- **PROVISIONAL caveats:** v_min is pilot-measured (est ~2); dx may need reduction
  at fast v; curve is the high-v Bethe tail (peak unreachable); Definition-1
  formula is NOT locked (data-collect only).
</guard_rails>

<tasks>
1. **Phase 0 — GS (MANUAL GATE).** Converge GS for 35×35×Lz, 25-Bohr slab,
   N=100, r_s≈4.2, `periodicity(2)`; lock Lz. Dashboard: n(z) profile, 2D slice,
   energy, ∫n dV, symmetry. Done: GS converged, checkpoint saved. **User inspects.**
   (tddft-simulations, simulation-validation.) — DONE 2026-07-21 (r_s=4.18).
2. **Phase 1a — vacuum exit test, static (MANUAL GATE).** Moving Gaussian charge
   swept across the z-open far face; verify clip (∫n 1→0) + no wrap (wrap witness
   ≈0). Done: clean clip confirmed. — DONE 2026-07-21 (PASS).
3. **Phase 1b — vacuum DYNAMIC exit (MANUAL GATE).** REAL propagation: the
   `Projectile` advances by velocity-Verlet each step, the moving perturbation
   re-solves φ_proj at each new position; run until proj_z ≥ +Lz/2 + Lz (≥ Lz
   beyond the face). Dashboard/notebook: φ_proj(t) animation + φ_peak(t) + z(t)
   linear. Validates the DYNAMIC perturbation-tracking through `propagate`. Done:
   φ_proj tracks the projectile and leaves smoothly. **User inspects.**
   (`scripts/classical_highdensity_sv/vac_dynamic/run.cpp`.)
4. **Phase 2 — dynamics / Ehrenfest validation (MANUAL GATE).** An INDEPENDENT
   check that the `Projectile` + Ehrenfest force + velocity-Verlet are correct,
   then a perturbation-vs-pseudopotential cross-check. See `<dynamics_validation>`
   for the locked test design. Done: analytic force matched, energy conserved,
   perturbation ≈ pseudopotential. **User inspects.** (code-test, simulation-validation.)
5. **Code — run.cpp + notebook builder (code-test).** Fork run.cpp to:
   `periodicity(2)`, mass-1 Ehrenfest Gaussian charge, full pairwise ledger +
   full energy decomposition emit, final checkpoint + `*_RESUME`. Adapt
   `run_notebook_builder.py` to read `projectile.csv` + a **step-by-step stopping
   section**. Done: unit test green, builder produces a notebook.
6. **Phase 3 — single-transit pilot at v≈2 (MANUAL GATE = central aim).** Full
   system, one velocity, Ehrenfest mass-1; `cutoff_guard.py` PASS first. Dashboard:
   z(t)/v(t) transit + exit, **E_electronic(t) plateau + E_absorbed**, full ledger,
   density GIF, step-by-step stopping calc, **ghost-UPF failure contrast**, Lindhard
   eyeball (NON-gating). Done: plateau + finite S + v_min measured. **User inspects.**
7. **Phase 4 — autonomous 6-velocity sweep.** Python orchestrator (idempotent
   resume, per-phase Gmail, cudaMemGetInfo probe, 2-GPU), grid = floor + 5 up. Each
   run → `analyse.py` + `run-notebook` (density GIF top, stopping section). Done:
   6 clean plateaus, 6 notebooks, S(v) points.
8. **Phase 5 — synthesis phase-notebook.** S(v)=E_absorbed/L across 6 points +
   full component ledger staged for Definition-1; Lindhard/bulk eyeball overlay
   (never a gate); WP-overlay-ready. Done: phase notebook + S(v) figure +
   S_summary.csv; frontmatter all done, status→done, handover + INDEX updated.
   (notebook-making.)
</tasks>

<dynamics_validation>
Phase 2 test design — LOCKED (user, 2026-07-21): A+B+C all three. Purpose: an
independent, closed-form check that the classical `Projectile` machinery
(`projectile.hpp` velocity-Verlet + `projectile_force.hpp` Hellmann-Feynman force)
is correct, then a demonstration that the perturbation and the pseudopotential
projectile give the same physics (one approximates the other).

- **Test A — analytic force (static, THE independent test).** The projectile force
  is F_z = −d/dR ∫ n_proj(·−R)·φ_drag with φ_drag = poisson(n_source). Feed a KNOWN
  fixed Gaussian source charge (width σ_s, charge Q_s) and compare INQ's
  `projectile_force_z` (finite-difference of the on-grid Poisson integral) to the
  **closed-form two-Gaussian Coulomb force** F(d) = −d/dd [Q_p Q_s erf(d/√(2(σ_p²+σ_s²)))/d]
  across a range of separations d. PASS = agreement to grid tolerance. This is a
  pure computation (no dynamics), best as a Catch2 unit test in
  `inq-stack/tests/include/inqkit/dynamics/`.
- **Test B — integrator + energy conservation (dynamic).** Launch the projectile
  toward a fixed repulsive Gaussian; verify (i) total energy ½mv² + U(d) conserved
  along the trajectory, (ii) the trajectory matches an independent Python
  Newton-ODE integration using the analytic F(d) (turning point / z(t)). Validates
  velocity-Verlet + force sign together, against a non-INQ reference.
- **Test C — perturbation vs pseudopotential.** Run the SAME scenario with the
  moving Gaussian-charge perturbation and with the ghost-UPF pseudopotential;
  overlay z(t)/force(d). Agreement (modulo the UPF's r_cut aliasing) demonstrates
  one approximates the other and justifies the perturbation as the clean primary.
</dynamics_validation>

<rules>
- ALWAYS run CAP-free; NEVER add a CAP to "clean up" the wake (extend Lz instead).
- ALWAYS `periodicity(2)` for RT (and the GS it loads) — a per-2 RT on a per-3 GS
  gives a spurious t=0 kick, and vice-versa.
- ALWAYS emit the full pairwise ledger + decomposition every run (Definition-1 is
  data-collect; do not drop columns to save space).
- NEVER gate on velocity drift or on Lindhard agreement. Gate only on
  correctness + plateau existence.
- NEVER report the projectile-KE-drag S as the headline or transfer it to the WP.
- NEVER `np.fft.fftshift` a VTI (vti-coordinate-mapping rule); load via
  `inqview.load_vti`.
- Density GIF at the TOP of every run notebook (notebook-density-gif rule).
- `TMPDIR=/local/data/public/skcb2/tddft/.build_tmp` for every INQ build (/tmp is
  9.8 G and fills).
</rules>

<preflight>
- [ ] Intent self-contained: hypothesis + success (Phase-2 plateau + finite S) /
      failure (no plateau at any transiting v) criteria; each task has a
      done-criterion.
- [ ] Setup reproducible: 35×35×Lz, 25-slab, N=100, r_s≈4.2, periodicity(2), no
      CAP; GS = Phase-0 task; Ehrenfest mass-1 Gaussian charge σ_WP=0.5; dt=0.04;
      per-velocity N_STEPS sized to full exit; observables + cadence enumerated;
      file placement per ADR-0007.
- [ ] New code pre-gated: run.cpp fork + notebook-builder adaptation → code-test
      (ledger closure + zero-force path) BEFORE the sweep. No new observable kernel.
- [ ] Validation & guard rails: cutoff_guard per velocity; correctness-only hard
      gates; Phase-2 numeric plateau gate; wake-localisation check; boundary +
      300-frame cadence; PROVISIONAL caveats (v_min pilot-measured; dx; Bethe tail;
      Def-1 formula TBD) named.
- [ ] Autonomous mechanics: cudaMemGetInfo GPU probe (NVML cosmetic; GPU default;
      warn if occupied by another user); Python orchestrator (idempotent resume,
      per-phase Gmail, one-shot retry); per-run + phase notebooks auto-built;
      checkpoint-don't-block; handover + frontmatter updated by the agent.
- [ ] Grounding: Method B = Correa 2018 (docs/sources/correa-2018-...); CAP
      diagnosis = energy-oscillation-diagnosis handover; engine claims carry
      file:line (cell.hpp periodicity, moving_gaussian_projectile_perturbation.hpp,
      projectile.hpp, interaction_energies.hpp).
</preflight>
