---
id: classical-slab-stopping
area: localised_jellium
title: Classical stopping baseline for the localised jellium slab (twin of qsp_phase5)
status: running
hypothesis: "A matched localised-slab classical projectile (σ_WP=0.5, v=1.3) deposits energy into the 82-electron jellium slab at a rate S_classical that sets the benchmark the WP quantum stopping (p5_wp_v1p3, S≈2.4 eV/Bohr, upper bound) must be compared against; the Ehrenfest initial-drag slope (Phase 1) and the prescribed constant-velocity ΔE_deposited/L_slab (Phase 2) bracket that classical baseline."
handover: docs/handovers/classical-stopping-baseline-localised-jellium.md
tasks:
  - { name: "Build const-v binary + code-test (Projectile zero-force path)", done: true }
  - { name: "Cutoff guard + GPU smoke both modes (no t=0 kick, const-v linear proj_z)", done: true }
  - { name: "Phase 1 — Ehrenfest v=1.3 run (twin p5_wp_v1p3); initial-drag S(v0)=0.49 eV/Bohr", done: true }
  - { name: "Phase 2 — const-velocity v=1.3 replica; deposit/L_slab S=0.43 eV/Bohr", done: true }
  - { name: "Per-run analyse.py + single-run notebooks (density GIF at top) for both", done: false }
  - { name: "S_classical vs S_WP vs Lindhard comparison (done); twin-decompose BLOCKED — needs CAP-free pairwise WP twin", done: false }
blocked_reason: ""
---

# Classical stopping baseline for the localised jellium slab

<identity>
You are a scientific computing researcher running first-principles rt-TDDFT in INQ.
You adhere to this repository's rules, skills, and workflows (tddft-simulations,
stopping-power-extraction, twin-run-analysis, run-notebook, code-test,
handover-update, scientific-grounding). σ-convention is UNIFIED: σ means the
wavepacket σ_WP; the classical Gaussian charge std = σ_WP/√2 (CONTEXT.md).
</identity>

<description>
The localised-jellium WP S(E) sweep (`hypotheses/qsp_phase5`) currently compares
its slab WP runs against a **bulk** classical reference (`classical_sigma0p5_bulk.csv`)
and point-charge Lindhard — a geometry mismatch flagged in **ADR 0010**
("slab-WP-vs-bulk is a labelled geometry estimate, not a matched pair"). This
campaign fills that gap: the **matched localised-slab classical baseline** in the
identical geometry, so the WP quantum stopping has a like-for-like classical
expectation to be judged against.

The projectile is a rigid Gaussian **charge** realised as a moving perturbation
(no ghost UPF, no r_cut aliasing — the accurate representation, CONTEXT.md
"perturbation (Gaussian-charge) projectile"). Two drive modes, run as two phases:

- **Phase 1 — Ehrenfest (light electron, mass 1).** The charge moves under its own
  Hellmann–Feynman force and DECELERATES (`light-projectile-stopping` rule): it
  enters the slab at v₀ and stops near the far face. The stopping power is the
  **initial-drag slope** S(v₀) = −d(KE_proj)/ds over the early v≥0.85·v₀ window,
  cross-checked by +d(E_deposited)/ds (energy conservation gate). This is the true
  classical twin of the σ=0.5, v=1.3 WP run `p5_wp_v1p3`.
- **Phase 2 — prescribed constant velocity (external drive).** The charge is driven
  at fixed v (zero force ⇒ a=0 ⇒ R=R₀+V₀·t). It transits the whole slab at v₀ and
  exits into vacuum; the stopping power is **ΔE_deposited / L_slab** read off the
  energy plateau after it clears the slab (the `stopping-power-extraction` slab
  method). N_STEPS sized to stop the center at z=+30 Bohr — past the far slab face
  (+12.5), short of the box edge (±45) — so the Gaussian never wraps the periodic
  box.

Comparing Phase 1 vs Phase 2 quantifies how much the light electron's deceleration
lowers the deposited energy relative to the idealised constant-velocity limit — the
two bracket S_classical(v=1.3).

**Success:** both runs complete with a flat conservation ledger; Phase 1 yields a
finite initial-drag S(v₀) from ≥30 early-window points; Phase 2 yields a clear
ΔE_deposited plateau; and a single figure places S_classical (both phases) against
S_WP(p5_wp_v1p3)=2.37 eV/Bohr [upper bound], the bulk classical (0.94), and Lindhard.
**Failure/inconclusive** is a valid reported outcome (e.g. Phase-1 deceleration too
steep to read a clean slope → recommend a higher-v anchor; report honestly).
</description>

<observables_set>
Per run (both phases), reusing the `classical_slab_stopping/run.cpp` writers — all
PRIMARY, no new observable/kernel:
- `raw/observables/observables[.segNNN].csv` — full energy ledger (total, kinetic,
  hartree, xc, external, nonlocal, ion) + step/time. Cadence: every step.
- `raw/observables/projectile[.segNNN].csv` — step,time_au,proj_z,proj_vz,
  energy_proj_ke,energy_proj_bg_ideal. Every step. (The S(v₀) source for Phase 1;
  the linear-proj_z check for Phase 2.)
- `raw/observables/interactions[.segNNN].csv` — pairwise P/S/B decomposition
  (e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,norm_slab,norm_proj). Every step. (Isolates the
  slab-internal deposited energy from the projectile-coupling terms.)
- `frames/total/density_t*.vti` — total density n(r,t), physical order. Cadence
  `LJ_SAVE_EVERY` ⇒ ~300 frames (mandatory density GIF, `notebook-density-gif`).
- Final `checkpoint/` + `rt_state.txt` (`final-timestep-checkpoint` rule).
</observables_set>

<resolved_decisions>
geometry:
  - cell 50×50×90 Bohr, dx=0.5, periodicity=3 (`.periodic()`) — EXACT match to the
    WP twin `p5_wp_v1p3` (run_summary: cell 50×50×90, spacing 0.5; run.cpp `.periodic()`).
  - localised background = slab, half_width 12.5 (axis z), N=82, n0=N/(Lx·Ly·2·half)
    = 82/62500 = 1.312e-3 ⇒ r_s≈5.67 (ADR 0008 background perturbation; neutrality
    ∫n₊=N, extra_electrons(82)).
  - GS reused: `shared_gs/slab_n82_L50x50x90/` (per-3, no new SCF).
projectile:
  - σ_WP=0.5 ⇒ σ_pot=σ_WP/√2=0.354 (moving Gaussian charge, −1). launch_z=−23.75
    (coasts ~11 Bohr of vacuum, enters slab face at exactly v₀). v₀=k0/m=1.3.
  - Phase 1: LJ_CONST_V=0, mass=1 (light electron twin). Phase 2: LJ_CONST_V=1,
    mass=1 (mass irrelevant at fixed v; kept for parity).
propagator:
  - LDA, dt=0.04 (match twin integrator), Ehrenfest step callback (RealTimeSession).
duration_and_energy:
  - Phase 1: N_STEPS=2000 (80 a.u.) — vacuum coast (8.6 a.u.) + full slab
    deceleration + ΔE_deposited plateau (est. range ~23 Bohr ≈ slab thickness 25).
    Checkpointed; extend via LJ_RESUME=1 if the plateau is not reached.
  - Phase 2: N_STEPS=1034 (41.4 a.u.) — const-v center travels −23.75→+30 at v=1.3,
    stopping before the Gaussian wraps (box edge ±45).
  - E_kin(v=1.3)=½·1·1.3²=0.845 Ha=23.0 eV. Cutoff guard PASS (E_cut=537 eV ≥
    1.1×E_kin=25 eV, dx=0.5, kind=classical).
pilot_and_io:
  - LJ_SAVE_EVERY: Phase 1 = 7 (~285 frames), Phase 2 = 4 (~258 frames).
  - Segment-suffixed CSVs on resume; post-processing concatenates by step.
file_placement (ADR-0007):
  - binary + dispatcher: `scripts/classical_slab_stopping/` (run.cpp cloned from
    proj_dyn + LJ_CONST_V mode; orchestrate.py).
  - run outputs: `scripts/classical_slab_stopping/results/{p1_ehrenfest_v1p3,
    p2_constv_v1p3}/`.
  - analysis: `hypotheses/classical_slab_stopping/` (combined CSVs, notebooks,
    comparison figure).
</resolved_decisions>

<guard_rails>
- Cutoff/aliasing guard PASS is a precondition (done: E_cut 537 ≥ 25 eV).
- Correctness gate (NOT a velocity-drift gate — `light-projectile-stopping` rule):
  the conserved quantity E_electronic + energy_proj_ke + U_proj_bg must be flat
  (Phase 1); for Phase 2 the const-v proj_z must be linear (R₀+V₀·t) and proj_vz
  constant. Abort on NaN / complex energy / a t=0 energy kick (periodicity mismatch).
- Phase-1 S-extraction gates on a clean initial-drag slope EXISTING (finite S,
  ≥~30 early-window points at v≥0.85·v₀; widen to 0.70/0.50 if sparse) — never on
  the (by-design) velocity drift.
- Boundary: light electron stops inside the slab (Phase 1) / const-v stops at +30
  before the box edge (Phase 2) — no wrap, no CAP needed (note: the WP twin has a
  CAP at ±35 to absorb its dispersing tail; the classical charge has none — the pair
  is matched on bath physics + projectile cloud, NOT byte-identical; the CAP is
  irrelevant to both stopping windows, which happen far from ±35).
- Budget: proceed at full scope, checkpointed; WARN + continue on projected overrun
  (`checkpoint-dont-block`), user owns the kill.
- PROVISIONAL: S_WP=2.37 is an UPPER BOUND (convergence-flagged in qsp_phase5); the
  comparison states this.
</guard_rails>

<tasks>
1. **Const-v binary + code-test** (done). `scripts/classical_slab_stopping/run.cpp`
   (LJ_CONST_V mode). The zero-force ⇒ constant-velocity path is the locked unit
   test `test_projectile.cpp:22` (125 assertions pass). No inqkit/formula change.
2. **Cutoff guard + smoke both modes.** Guard PASS (done). GPU smoke (20 steps each,
   GS reused): Phase-1 shows no t=0 energy kick (periodicity-3 GS matches);
   Phase-2 shows proj_z = −23.75 + 1.3·t exactly and proj_vz≡1.3. Done-criterion:
   both smokes complete, ledger finite, const-v linearity verified.
3. **Phase 1 — Ehrenfest run** (GPU 0). Full 2000-step run, SAVE_EVERY=7,
   checkpointed. Done-criterion: run_completed=true; ΔE_deposited plateau reached
   (extend via resume if not); conserved ledger flat.
4. **Phase 2 — const-velocity run** (GPU 1, concurrent). 1034 steps, SAVE_EVERY=4.
   Done-criterion: run_completed=true; proj_z linear to +30; ΔE_deposited plateau.
5. **analyse.py + single-run notebooks.** Per-run `analyse.py` (full inqview
   pipeline → REPORT.md) + a `run-notebook` each with the mandatory
   `density_evolution.gif` displayed at the TOP (`notebook-density-gif`).
   Done-criterion: both `.ipynb` executed with the GIF embedded inline.
6. **twin-run-analysis + comparison.** Run `twin_decompose.py` on Phase 1 vs
   `p5_wp_v1p3` (per-timestep P/S/B ledger, residual = WP self-Hartree/SIE).
   Extract S_classical(v=1.3) via `stopping-power-extraction` for both phases;
   emit one figure: S_classical(P1 initial-drag), S_classical(P2 ΔE/L),
   S_WP=2.37[UB], bulk classical 0.94, Lindhard. Done-criterion: figure written to
   `hypotheses/classical_slab_stopping/`, numbers in the handover (2 s.f.).
</tasks>

<rules>
- ALWAYS twin p5_wp_v1p3 exactly on geometry/GS/σ/launch/dt; the ONLY intended
  physical differences are (Phase 1) classical-vs-quantum projectile and (both)
  the absent CAP.
- NEVER extract S from a full-run regression (averages over the decelerating v range)
  — Phase 1 uses the initial-drag window only.
- NEVER `np.fft.fftshift` a VTI (physical order; `vti-coordinate-mapping` rule).
- Report numbers at 2 s.f. (`number-rounding` rule); S_WP labelled UPPER BOUND.
- Per-phase Gmail via `inqview.email.send_run_email` (four-part structure, ≥1 plot).
</rules>

<preflight>
- [x] Intent self-contained: falsifiable hypothesis + success/failure/inconclusive
      criteria; every task has a done-criterion.
- [x] Setup reproducible: geometry/N/r_s/box, GS checkpoint named
      (slab_n82_L50x50x90, per-3), propagator+dt+steps+energy with values +
      one-line justifications; observables + cadence enumerated; ADR-0007 placement.
- [x] New code pre-gated: LJ_CONST_V zero-force path = locked unit test
      test_projectile.cpp:22; no new observable/kernel/formula.
- [x] Validation & guard rails: cutoff guard PASS; conservation (not v-drift) gate;
      initial-drag-slope-exists gate; PROVISIONAL S_WP UB caveat.
- [x] Autonomous mechanics: GPU via cudaMemGetInfo probe (NVML broken; both free);
      Python orchestrator (idempotent resume, per-phase try/except + Gmail);
      run-notebook + density GIF contract; handover pointer present.
- [x] Grounding: geometry/periodicity cited to p5_wp_v1p3 run_summary + run.cpp:77;
      metric choice cited to light-projectile-stopping + stopping-power-extraction.
</preflight>
