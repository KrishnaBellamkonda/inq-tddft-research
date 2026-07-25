---
# Campaign frontmatter — single source of truth for the INDEX (build_index.py).
# id NEVER changes. The EXECUTING agent flips task `done` flags + bumps `status`.
id: nazarov-gross-validation
area: nazarov_gross_comparison
title: "Nazarov–Gross validation — fixed-velocity projectile-mass sweeps (null + slow branches)"
status: ready
hypothesis: "Nazarov & Gross (arXiv:2510.26222): same-charge projectiles of
  different mass moving at the same velocity experience different friction in an
  electron liquid — but ONLY in the slow (sub-Fermi-velocity) regime. Testable
  two ways on the r_s=4.0 localised slab: (NULL) at v=2.711 a.u. = 5.6·v_F the
  stopping of m ∈ {0.5, 0.71, 1(p3), 1.41} wavepackets is flat within the
  spreading/SIE noise floor; (SLOW, next phase — piloted here) at v=0.25 a.u. =
  0.52·v_F the stopping splits measurably and mass-ordered per S = Q(m)·v."
handover: docs/handovers/nazarov-gross-comparison.md
tasks:
  - { name: "Phase 1a — two new slab GS (n234, h=0.35 and h=0.40), background-consistent, converged", done: true }
  - { name: "Phase 1b — smoke gate at h=0.40 (worst-case m=1.41 rung, ~60 steps): energy bounded, no NaN, measured step-cost fits budget (re-planned from h=0.35, 2026-07-12)", done: false }
  - { name: "Phase 1c — null branch: 3 WP runs m ∈ {0.5, 0.71, 1.41} at v=2.711, h=0.40, p3-identical otherwise (m=2.2 dropped: aliasing BLOCK at 0.40)", done: false }
  - { name: "Phase 1d — slow pilots m ∈ {1, 10} at v=0.25, h=0.40 (budget-gated; sized to remaining wall, deferred if none left)", done: false }
  - { name: "Phase 1e — analysis: S(m) null-branch figure via the p3 retained-energy ledger method + spreading-systematic check + pilot initial-drag verdicts, executed notebook", done: false }
  - { name: "Phase 2 — slow-branch production ladder (design gated on 1d pilot verdicts; NOT yet authored)", done: false }
blocked_reason: ""
---

# Nazarov–Gross validation — fixed-velocity projectile-mass sweeps

<identity>
You are a scientific computing researcher working on first-principles TDDFT
simulations. σ always means the wavepacket width σ_WP. NEVER edit `inq/`; the
per-state mass fork lives ONLY in `inq-study`. Report numbers at 2 s.f.
(3 s.f. for near-equalities).
</identity>

<description>
**Why.** Nazarov & Gross 2025 (source note:
`docs/sources/nazarov-gross-2025-quantum-projectile-stopping.md`) claim the
friction felt by a projectile in an electron liquid depends on its MASS at fixed
charge and velocity — a purely quantum effect (the classical/heavy limit
recovers mass-independent Lindhard stopping). Their result is derived in the
slow-projectile (friction) limit. The regime boundary is the Fermi velocity of
the target: the r_s = 3.996 slab (`shared/configs/slab_n234_L50.hpp`) has
k_F = v_F ≈ 0.48 a.u., E_F ≈ 3.1 eV.

**The two-branch experiment (locked 2026-07-11 grill).**
- **Null branch (this phase):** clone the p3 baseline
  (`scripts/fullsuite_wp/results/p3_wp`, analysed in
  `hypotheses/qsp_phase2/quantum_stopping_ledger_p3_26-6-26.ipynb`) at its exact
  velocity v = 2.711 a.u. (= 5.6·v_F, HIGH-velocity regime) for 4 projectile
  masses. NG predicts NO appreciable S difference here — this is the control /
  noise-floor half of the test. An observed splitting is an artefact channel
  (mass-dependent packet spreading ∝ 1/m, SIE, grid), not NG friction.
- **Slow branch (next phase; 2 pilots here):** v = 0.25 a.u. ≈ 0.52·v_F, where
  NG's mass-dependent friction is the predicted signal. The pilots de-risk it:
  is the m=1 zero-point-dominated packet (Δv/v ≈ 5.7) extractable at all, and
  does a heavy rung (m=10, Δv/v ≈ 0.57) give a clean initial-drag window?

**Decisions this informs.** Whether NG's mass-dependent stopping is visible in
real-time mean-field TDDFT with free-Gaussian projectiles; whether the slow
production ladder (Phase 2) is feasible; what the artefact floor is at high v.

**Falsifiable outcomes.**
- Null branch CONSISTENT with NG iff max|S(m) − S̄|/S̄ across the 4 same-grid
  rungs is within the artefact floor estimated from the spreading-systematic
  check (task 1e). A splitting exceeding that floor at high v FALSIFIES the
  flat-null expectation and must be reported as an artefact-channel finding
  (not NG friction) with the spreading channel as prime suspect.
- Pilots PASS iff each yields a finite initial-drag S over ≥ ~30 early-window
  points at v ≥ 0.85·v₀ (light-projectile rule — NEVER gate on v-drift).
  m=1 failing while m=10 passes is itself a reportable pilot verdict (redesign
  the light rungs: wider σ or eigenstate-width packet) — it does NOT block m=10.
</description>

<locked_parameters>
Everything inherited from p3 verbatim unless stated (user-locked 2026-07-11):
- **Target:** localised slab, N=234, 50³ Bohr periodic box, slab half-width 12.5
  along z, r_s = 3.996, LDA, ETRS, gamma-only, temperature 0.00862 eV, no CAP.
  Cfg: `shared/configs/slab_n234_L50.hpp`.
- **Projectile:** σ_WP = 0.5 Bohr free Gaussian, charge −1, injected into the
  last extra state and orthogonalised against occupied; mass set per-state via
  the PROVEN inq-study fork `electrons.inverse_mass()[0][wp_idx]` (validated in
  the muon-mass-fork campaign Phases 1–3; bath stays mass-1 electrons, so the
  GS is mass-independent and shared per spacing).
- **Sweep invariant: VELOCITY** (fixed-velocity mass sweep — see CONTEXT.md).
  k₀ = m·v and E = m·v²/2 per rung.
- **Null branch:** v = 2.7110633401 (p3's k₀ for m=1), masses {0.5, 0.71, 1.41};
  m=1 anchored by the EXISTING p3 run (NOT rerun; cross-grid caveat below).
  launch_z = −23, dt = 0.02, N_STEPS = 880, WRITE_EVERY = 10 (p3-identical).
- **Spacing (RE-PLANNED 2026-07-12, user-locked):** h = 0.40 for ALL runs on the
  converged GS40. The original h = 0.35 plan was EXECUTED and MEASURED unusable
  on the 24 GB GPUs: ~260 s/step effective (memory-thrash + 143=11×13 FFT radix;
  evidence `scripts/nazarov_gross/wp/smoke.log`) ⇒ ~60 h per 880-step run vs
  ~4 h predicted. Energy conservation was FINE — the failure was cost, not
  physics. Aliasing ceiling at h=0.40 is m ≈ 1.8 (2%-tail tier), so the locked
  m=2.2 rung is IMPOSSIBLE and was dropped (user chose 3 surviving rungs over
  adding a new m=1.8 rung). Heavier rungs must not be added at this spacing.
- **Cross-grid caveat (accepted):** p3 (h=0.5) anchors m=1; flatness is judged
  primarily among the 3 same-grid h=0.40 runs, with p3 as a secondary anchor.
- **Slow pilots:** v = 0.25, masses {1, 10}, h = 0.40 (its own GS; aliasing
  trivial: m=10 ⇒ k₀ = 2.5, clean 3σ). launch_z = −13.5 (1 Bohr = 2σ outside the
  slab face — the slow packet must not waste the run approaching), dt = 0.02,
  N_STEPS sized to the remaining wall budget (floor 800 steps ≈ 16 a.u.; defer
  the pilots entirely if even the floor does not fit).
- **Compute budget: 14 h wall × 2 GPUs** (clock reset at the 2026-07-12
  relaunch). Sequence: smoke (h=0.40, budget gate) → round 1 (m0.5 ∥ m0.71) →
  round 2 (m1.41 ∥ pilot m=10) → pilot m=1 (budget-sized tail). Estimated
  ~11 h at the predicted ~11 s/step; the smoke MEASURES the true cost and the
  null rounds do not launch if two sequential rounds cannot finish in budget.
- **S extraction:** null branch uses the SAME retained-energy ledger method as
  `quantum_stopping_ledger_p3_26-6-26.ipynb` (comparability with the anchor);
  pilots use the initial-drag window (S = −dKE/ds over v ≥ 0.85·v₀, per
  `.claude/rules/light-projectile-stopping.md`).
</locked_parameters>

<observables_set>
The p3 canonical suite, unchanged (fullsuite_wp writer set): density VTIs
total/system/gs_system/wp + wavefunction_wp, density_delta(+coarse);
observables.csv (energies/current/dipole/L2), state_energies, occupations,
eigenvalues; WP: momentum_distribution, wp_momentum_stats,
**wp_real_space_stats** (σ_ρ(t) — REQUIRED for the spreading-systematic check),
overlap + overlap_full, electron_number. NO new observable/kernel is introduced.
</observables_set>

<guard_rails>
- **Smoke gate before the expensive runs:** ~60 steps at h=0.40 with the
  worst-case rung (m=1.41, largest k₀): completes, no NaN, |E_total drift|
  bounded (< 1e-3 Ha over the smoke), WP norm ≈ 1 after injection, AND the
  measured step cost fits the remaining budget (HARD gate on the null rounds —
  added 2026-07-12 after the h=0.35 cost blow-up). Failure ⇒
  status blocked + email + STOP.
- **Cutoff guard** (`.claude/skills/tddft-simulations/cutoff_guard.py` tiers):
  every rung checked pre-launch with σ_p = 1/(√2·σ_WP) = 1.414 (NOT 1/(2σ)).
  BLOCK if aliased tail > 2%.
- **Light-projectile rule:** NEVER abort on velocity drift; slow pilots are
  EXPECTED to decelerate/stop. Gate only on "a clean initial-drag slope exists".
- **Boundary/cadence:** launch positions ≥ 2 Bohr (4σ) from box walls; VTI
  cadence WRITE_EVERY=10 (p3-identical, 88 frames — p3 comparability overrides
  the 300-frame target for the null branch).
- **Abort conditions:** NaN / complex energy / GS checkpoint missing. GPU is the
  default; NVML/nvidia-smi is broken but compute works — never CPU-fall-back on
  an nvidia-smi error; probe with cudaMemGetInfo where needed.
- **Budget enforcement:** the orchestrator records t₀ and never starts a pilot
  whose floor-size estimate exceeds the remaining 14 h wall; deferred pilots are
  emailed as deferred, not silently dropped.
- **Idempotent resume:** any run whose `run_summary.txt` shows
  `run_completed = true` is skipped on orchestrator restart.
- **Verdict hygiene:** no "NG validated/violated" claim from the null branch
  alone; a null-branch splitting is reported as an artefact-channel finding with
  the spreading check attached. The positive theorem test is Phase 2.
</guard_rails>

<file_placement>
Per ADR-0007, sweep = `nazarov_gross`:
- Run machinery: `ResearchProject/systems/localised_jellium/scripts/nazarov_gross/`
  (`gs/run.cpp` env-spaced GS, `wp/run.cpp` env-driven mass/velocity WP clone,
  `orchestrate.py` Python orchestrator, per-run results under `wp/results/<run>`).
- New GS checkpoints: `shared_gs/slab_n234_L50_h0p35`, `shared_gs/slab_n234_L50_h0p40`.
- Analysis: `hypotheses/nazarov_gross/` (S(m) figure, ledger notebooks, pilot
  verdicts). Run names: `null_m0p5`, `null_m0p71`, `null_m1p41`, `null_m2p2`,
  `pilot_slow_m1`, `pilot_slow_m10`.
- Logs/results gitignored; provenance in run_summary.txt + handover.
</file_placement>

<rules>
- NEVER edit `inq/`; the mass fork is inq-study only (INQ_SOURCE=inq-study).
- Phases gate in order: GS → smoke → null runs → (budget-gated) pilots.
- The executor is the Python orchestrator (idempotent, per-stage try/except,
  full-traceback failure email, Gmail per stage via inqview.email).
- VTI is physical-order — NEVER fftshift a VTI; load via inqview.load_vti.
- Flip frontmatter `done` flags + `status` and update the handover as stages
  complete.
</rules>

<preflight>
- [x] Self-contained intent: falsifiable per-branch criteria; every task has a
      done-criterion.
- [x] Reproducible setup: geometry/GS/propagator/dt/steps/masses/velocities all
      locked with values + justification; observables enumerated; ADR-0007 paths.
- [x] New code pre-gated: NO new engine code (mass fork already validated by the
      muon-mass-fork campaign Phases 1–3, incl. bit-for-bit regression).
- [x] Validation & guard rails: smoke gate, cutoff guard, budget gate, abort
      conditions, light-projectile extraction rule, verdict hygiene.
- [x] Autonomous mechanics: Python orchestrator + per-stage Gmail + idempotent
      resume + handover pointer + frontmatter flips.
- [x] Grounding: NG source note in docs/sources; regime numbers derived from
      slab_n234_L50.hpp (r_s=3.996 ⇒ v_F≈0.48); mass-fork trust from the
      muon-mass-fork campaign.
</preflight>

<background_notes>
Original mind-dump (pre-grill, kept verbatim for provenance):

I will call this campaign nazarov-gross-validation. In this campaign, we are
going to test the fundamental theorem proposed by nazarov gross - The stopping
power of projectiles of the same charge but different masses, would also be
different. The difference is attributed to the size of the wavepacket in the
paper. Now I do not quite understand exactly what they mean by the the size of
the wavepacket. But, in this campaign, we are going to stick to one gaussian
width (for both classical and quantum wavepackets).

In this campaign, we consider the quantum wavepacket runs that have successfully
been executed (concentrated wavepacket runs in which we've observed that the
total excess energy plateaus). This excess energy plateuing helped make an
estimate of the quantum stopping power. Now, we need to find 3 such runs that
have been successful with different total energy of the wavepacket. Then, we are
going to make mass of the electron heavier. Consider 2 new masses (3 in all
considering the electron runs as one of the masses). Then, we are going to make
the same simulation.

However, before we finalise this simulation, we need to carefully think about
the energy of the wavepacket. We need to be in the slow energy regime for the
mass to play a massive effect. in the high energy limit, the different particles
should behave similarly. So, the above runs that are mentioned, should not show
appreciable differences in stopping power. However, we need to consider a low
velocity run for an electron and the two other masses. Run the stopping power
simulation, find the stopping power can them make comparisons.

This was tell us if Nazarov Gross has been successfully validated.

There are a few other questions nazarov gross paper prompts
1. What do they mean by "mass is linked to the size of the particle"? What
should be the effective gaussian width to simulate a classical particle? (I've
answered this question by running classical ehrenfest simulations, and fidnign
the stopping power. I found that around 0.5 bohr of gaussian width for the
radial potential, I found that they classical tddft results matched the
expectation in the high velocity regime). Hence, I have been trying 0.5 bohr as
the keystone to run these simulations. However, I want to understand if my
assumption of usin gthe classical gaussian width as a benchmark to set the
classical wavepacket is right. This begs the question - "Is ther an intrinsic
size of an electron?" The coronene simulations might help address this.
</background_notes>
