---
# Campaign frontmatter — single source of truth for the INDEX (build_index.py).
# id NEVER changes. The EXECUTING agent flips task `done` flags + bumps `status`.
id: muon-mass-fork
area: muon_mass_fork
title: "Per-state mass fork — validate, then all-muon jellium LDA vs muon-XC"
status: paused                            # draft → ready → running → blocked → paused → done (Phase-4 user pick)
hypothesis: "The inq-study per-state mass fork correctly simulates arbitrary-mass Kohn–Sham orbitals — provable by analytic dispersion/spreading oracles, conservation laws, and a bit-for-bit electron regression — and, once trusted, an all-muon localised jellium with a muon wavepacket incident reveals a measurable difference between standard electron-parameterised LDA and a mass-appropriate muon exchange–correlation prescription."
handover: docs/handovers/muon-mass-fork.md
tasks:
  - { name: "Phase 1 — code bug checks: GPU build + engine tests (kernel eigenvalue, expectation/ledger, wrong-slot, GPU-vs-CPU, MPI partition) all green + catalogue rows", done: true }
  - { name: "Phase 2 — physics tests: σ(t) spreading + v_group + KE/norm/time-reversal + particle-in-box + E(k) + mass-dial; xz density vs muon σ; Phase-2 notebook", done: true }
  - { name: "Phase 3 — bit-for-bit regression vs pristine inq for an electron-only system (GS + RT); Phase-3 notebook", done: true }
  - { name: "Phase 3b — a few simple real-time propagation sanity sims", done: true }
  - { name: "Phase 4 — research the muon XC functional (literature-review); present grounded candidates + PAUSE for user pick (checkpoint); docs/sources entries", done: false }
  - { name: "Phase 5 — recreate the r_s=5.69 bath all-muon + incident muon WP: LDA vs user-picked muon-XC; comparison run-set + notebook", done: false }
  - { name: "Index notebook — guided read-order over all phase + run notebooks", done: true }
blocked_reason: "Phase 4 checkpoint — awaiting user muon-XC pick (write muon_xc_pick.json)"
---

# Per-state mass fork — validate, then all-muon jellium LDA vs muon-XC

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository. σ always means the wavepacket width σ_WP. NEVER edit `inq/`; the mass
fork lives ONLY in `inq-study`.
</identity>

<description>
**Why this campaign exists.** A per-state inverse-mass fork was added to the
`inq-study` INQ replica so any Kohn–Sham orbital can carry an arbitrary mass
(muon, band-structure effective mass, …). The mechanism (design + call chain):
`docs/campaigns/muon_projectile/inq_study_engine_notes.md`; implementation +
validation matrix: `docs/plans/muon-mass-fork-implementation.md`. The code is
written and a first CPU kernel test (Tier-1) passes; this campaign takes it from
"compiles + one test" to "trusted + a physics result", autonomously and in order.

**The decision it informs.** (1) Is the mass fork correct enough to trust for
production muon / effective-mass physics? (2) For an all-muon jellium, does the
choice of exchange–correlation functional (naive electron-LDA vs a
mass-appropriate muon prescription) change the stopping/response of an incident
muon wavepacket measurably?

**Falsifiable success/failure.**
- The fork is TRUSTED iff every Phase-1/2/3 check passes its numeric criterion
  (below). ANY failure ⇒ status `blocked`, email, stop — do not proceed to physics.
- The XC comparison is CONCLUSIVE iff Phase-5 shows a difference in a primary
  observable (S, wake, ⟨T⟩ ledger) between LDA and the muon-XC that EXCEEDS the
  run-to-run/SIE noise floor bounded by a vacuum-WP control; otherwise the honest
  result is "no measurable XC sensitivity at this r_s", which is still reportable.

**Structure (ordered phases; each gates the next).**
1. Code bug checks & detection (no physics claims).
2. Physics tests via runs/sims — ALL analytic oracles + conservation laws here.
3. Bit-for-bit electron regression vs pristine `inq`.
   3b. A few simple RT sanity sims.
4. Muon XC functional research.
5. All-muon localised jellium + incident muon WP: LDA vs muon-XC.
Plus: run-notebooks (per consequential run), phase-notebooks (per phase), and one
INDEX notebook that dictates read order.
</description>

<observables_set>
Reuse the ADR-0006 minimal/maximal set (`inqkit/observables/minimum_observable_set.hpp`)
and the current cadence. Per run-type:
- **Vacuum WP (Phase 2):** `wp_real_space_stats` (σ_ρ(t) — the spreading oracle),
  `wp_momentum_stats`, centre-of-density (v_group), ⟨T⟩ (KE conservation), norm.
- **Jellium GS (Phase 3/5):** E_GS ledger (E_total, E_kinetic, E_H, E_xc), density
  VTI (xz slices), eigenvalues.
- **Jellium RT + WP (Phase 5):** the full WP suite — n(k,t) coherent-peak → S,
  density_wp/density_system/density_delta VTI, E_total ledger, classical ion track
  where a classical control is run.
**NO new observable/kernel** is introduced by the physics phases (they reuse the
existing suite). The ONLY new engine code is the mass fork itself, which is
pre-gated in Phase 1 (code-test + the plane-wave/expectation oracles +
catalogue rows) BEFORE any expensive run.
</observables_set>

<resolved_decisions>
**Engine fork (LOCKED — implemented; see `docs/plans/muon-mass-fork-implementation.md`):**
- Per-state `inverse_mass_` on `electrons` (mirrors `occupations_`), default 1.0.
- Opted into `ks_hamiltonian` via `set_inverse_mass()`; empty-factor guard routes
  all-mass-1 through the ORIGINAL scalar path (bit-for-bit). Applied by
  `laplacian_states`/`laplacian_add_states`/`laplacian_expectation_value_states`.
  RT opt-in `real_time/propagate.hpp:79`; GS opt-in `ground_state/calculator.hpp`.
- **Muon mass = 206.77 m_e** (PDG μ/e mass ratio).
- Gamma-only (kpin 0); spinor_dim==1; zero vector potential (WP momentum in the
  orbital phase). Vector-potential A²/2m term NOT mass-scaled → valid in this regime.

**Phase-2 vacuum σ(t) sim (LOCKED):** empty cubic box, single injected Gaussian WP,
free Ehrenfest (ETRS). σ_WP=0.5 Bohr, L=48 Bohr, spacing 0.4, dt=0.02 a.u.
- Panel A: k₀=0, t→12 a.u., masses {1, 206.77}. Oracle `σ_ρ(t)²=σ_ρ0²+(t/2mσ_ρ0)²`.
- Panel B: muon-only, t→120 a.u. (muon reaches σ_ρ~0.6).
- Panel C: k₀=0.5 Bohr⁻¹, centroid slope = k₀/m.

**Phase-2 xz-density-vs-σ visualisation (LOCKED set):** muon WP GS/early-RT xz
density slices at σ_WP ∈ {0.5, 1.0, 2.0, 4.0} Bohr → into the Phase-2 notebook
(canonical theme; VTI in physical order — NEVER fftshift a VTI, `inqview.load_vti`).

**Phase-3 regression system (LOCKED):** the smallest converged electron jellium
already in the repo (or an electron GS+short-RT), run with the fork present but ALL
mass=1, compared bit-for-bit to the SAME input built against pristine `inq`.

**Phase-5 target (LOCKED geometry, PROVISIONAL functional):** RECREATE the
canonical r_s=5.69 "Li-like" jellium bath — density 1.296e-3 e/Bohr³, N=162,
L=50 Bohr cubic, dx=0.40 (GS template
`ResearchProject/systems/jellium/save_gs/gs_L50_cubic_N162_dx0p40`) — but with ALL
bath particles muon-mass (global inverse_mass=1/206.77) and a muon WP incident.
Rationale (user, 2026-07-06): at r_s=5.69 CORRELATION dominates over exchange, and
making the bath muon-mass lowers the kinetic energy (∝1/m) so the interaction (XC)
dominates further — PRECISELY the regime where the electron-parameterised LDA is
most suspect, making the LDA-vs-muon-XC comparison most informative. The all-muon
GS must be recomputed (kinetic ∝1/m ⇒ a DIFFERENT GS from the electron bath — the
GS-mass path `calculator.hpp` handles this); the original electron r_s=5.69 run is
the reference. Caveat: at the same PHYSICAL density the muon effective coupling is
r_s,eff ≈ 5.69×206.77 (deep-correlation) — LDA is stretched, which is the point of
the comparison. The muon-XC prescription is chosen in Phase 4 (user-picked
checkpoint) and is PROVISIONAL until then.

**File placement (ADR-0007).** All notebooks under
`ResearchProject/systems/localised_jellium/hypotheses/muon_mass_fork/`:
`index.ipynb` (read-order guide), `phase1_code_checks.ipynb`,
`phase2_physics.ipynb`, `phase3_regression.ipynb`, `phase4_xc_research.ipynb`,
`phase5_lda_vs_muonxc.ipynb`, plus per-run deep-dive notebooks (run-notebook skill).
Vacuum σ(t) runs under a `vacuum_wp/` subfolder there. Engine tests in
`inq-study/tests/`. Run outputs: logs gitignored, provenance only.
</resolved_decisions>

<guard_rails>
- **Phase gating is strict.** Phase N+1 does not start until Phase N's numeric
  criteria all pass. A failure ⇒ set `status: blocked`, `blocked_reason`, send the
  failure email (full traceback), STOP. One phase's failure never silently
  proceeds.
- **Phase 1 is the pre-gate for ALL physics.** No GPU physics run launches until
  the GPU build is green AND every Phase-1 engine test passes (kernel eigenvalue
  k²/2m; expectation ⟨T⟩=k²/2m + ledger identity; wrong-slot; GPU-vs-CPU 6 s.f.;
  MPI 1-vs-N-rank). Catalogue rows added BEFORE physics.
- **Bit-for-bit gate (Phase 3) is a hard trust gate.** If the all-mass-1 fork does
  NOT reproduce pristine `inq` to machine precision on E_total and density, the
  fork is BROKEN — stop, do not report any muon physics.
- **Boundary + cadence:** 4σ/1σ launch-stop; 300-frame VTI cadence for every WP
  run. Vacuum σ(t): WP tail must stay < box boundary (abort if 4σ_ρ(t) ≥ L/2).
- **Abort conditions:** NaN / complex energy / GPU occupied by another user
  (cudaMemGetInfo probe; NVML/nvidia-smi is broken but compute works — never CPU-
  fall-back on an nvidia-smi error).
- **SIE bound:** before any Phase-5 "XC-difference" claim, bound the muon
  self-interaction with a vacuum-WP control; report the difference only if it
  exceeds that floor.
- **PROVISIONAL:** the muon-XC functional (Phase 4 output) and the strongly-
  correlated-LDA caveat (Phase 5) are provisional until Phase 4 grounds them.
</guard_rails>

<tasks>
**Phase 1 — code bug checks & detection** *(code-test; no physics)*
GPU-build `inq-study` (mirror `inq/build`: cuda-12.5 nvcc, arch 80); confirm the
`_states` kernels compile under nvcc; rerun the Tier-1 CPU test on GPU. Add engine
tests: (a) `kinetic_expectation_value` on an orbital_set with a muon slot →
⟨T⟩=k²/2m per state; (b) ledger identity (apply factor == expectation factor);
(c) wrong-slot guard; (d) GPU-vs-CPU 6 s.f.; (e) MPI `-np 2` 1-vs-N-rank identical
(the partition-alignment guard). *Done when:* all green + catalogue rows in
`docs/validation/test-catalogue.md`.

**Phase 2 — physics tests** *(tddft-simulations, simulation-validation, notebook-making)*
Run the Phase-2 vacuum sim (Panels A/B/C) and the analytic-oracle checks: σ(t)
spreading (fit m to <5%; rate ratio 206.77), v_group=k₀/m (<2%), KE conservation
(<0.1% drift), norm conservation (<1e-4), time-reversal (recover WP), particle-in-
a-box Eₙ∝1/m (<1%), dispersion E(k)=k²/2m, mass-dial continuity (rate ∝ 1/m linear
through origin). Produce the **xz density vs muon σ ∈ {0.5,1,2,4}** slices.
*Done when:* every oracle passes its criterion + `phase2_physics.ipynb` (with the
xz-density-vs-σ figure) executed.

**Phase 3 — bit-for-bit electron regression** *(code-test)*
Run an electron GS + short RT with the fork present (all mass=1) and against
pristine `inq`; diff E_total, E_kinetic, and the density field. *Done when:*
identical to machine precision + `phase3_regression.ipynb`. FAILURE = fork broken.

**Phase 3b — simple RT sanity sims** *(tddft-simulations)*
A few short real-time propagations (e.g. a muon WP in a weak external potential /
a small muon jellium) to confirm stable, physical dynamics end-to-end. *Done when:*
runs complete without NaN, energy bounded, qualitatively sensible; noted in a
run-notebook.

**Phase 4 — muon XC functional research + USER CHECKPOINT** *(literature-review)*
Research exchange–correlation for a heavy/muon Fermi gas: mass-rescaled LDA (HEG
effective-units scaling), two-component / multicomponent DFT, muon-in-matter
functionals, Car–Parrinello fictitious-mass (DISTINGUISH — not physical XC).
Present the grounded, implementable candidate prescriptions (how each differs from
naive electron-LDA, implementation cost, validity at r_s,eff), then **PAUSE and
EMAIL the user to PICK** (user decision, 2026-07-06 — NOT agent-chosen). Set
`status: paused`, `blocked_reason: "awaiting user muon-XC pick"`; resume Phase 5
only after the user selects. *Done when:* candidates grounded in `docs/sources/` +
`phase4_xc_research.ipynb` + the user's pick recorded in the handover.

**Phase 5 — all-muon r_s=5.69 jellium: LDA vs muon-XC** *(tddft-simulations, notebook-making)*
Recompute the all-muon GS of the canonical r_s=5.69 bath (N=162, L=50 cubic,
dx=0.40; calculator.hpp mass path), then run an incident muon WP TWICE: standard
electron-LDA and the user-picked Phase-4 muon-XC. Extract the primary observables
(S via coherent-peak, wake, E_total ledger); compare LDA vs muon-XC AND vs the
original electron r_s=5.69 run, with the SIE floor bounded by a vacuum-WP control.
*Done when:* comparison run-set + `phase5_lda_vs_muonxc.ipynb` with a verdict
(difference exceeds floor → XC-sensitive at r_s=5.69; else → not sensitive here).

**Index notebook** *(notebook-making)*
`index.ipynb`: a guided read-order over all phase + run notebooks (Phase 1→5),
one line per notebook stating what it establishes and its trust-gate role. *Done
when:* index executes and links every notebook.
</tasks>

<rules>
- NEVER edit `inq/`. The fork lives in `inq-study` only; keep `inq/` as the
  bit-for-bit regression reference (Phase 3).
- Phases run IN ORDER; a phase's gate must pass before the next starts.
- The autonomous executor MUST be a Python orchestrator (idempotent resume,
  per-phase try/except + full-traceback failure email, one-shot retry), NOT bash.
  Reference: `ResearchProject/systems/localised_jellium/scripts/campaign_autorun/orchestrate.py`.
- GPU is the default; schedule via `cudaMemGetInfo` probe; warn if a GPU is
  occupied by another user; never CPU-fall-back on an nvidia-smi/NVML error.
- σ always means σ_WP; report numbers at 2 s.f. (3 s.f. for near-equalities).
- VTI is physical-order — NEVER fftshift a VTI; load via `inqview.load_vti`.
- No muon-physics claim before the Phase-1 gate AND the Phase-3 bit-for-bit gate.
- Per-phase Gmail (email-notifications skill); update the handover + flip
  frontmatter `done`/`status` as each phase completes.
</rules>

<preflight>
- [x] Intent self-contained: falsifiable hypothesis + per-phase numeric criteria;
      each task has a done-criterion.
- [ ] Setup reproducible: engine fork LOCKED (line-refs in the plan); vacuum-sim
      params LOCKED; Phase-3 regression system + Phase-5 geometry LOCKED; **OPEN:**
      Phase-4 muon-XC prescription is agent-resolved (research gate), Phase-5
      functional PROVISIONAL until then; confirm GS sources for Phase-3/5.
- [ ] New code pre-gated: mass fork → code-test + oracles + catalogue rows in
      Phase 1 BEFORE any physics run (GPU build still pending).
- [x] Validation & guard rails: strict phase gating; bit-for-bit trust gate;
      pilot/boundary/NaN aborts; SIE floor; PROVISIONAL caveats named.
- [ ] Autonomous mechanics: Python orchestrator + cudaMemGetInfo probe + per-phase
      Gmail + notebook contracts (phase + run + index) + handover pointer +
      frontmatter flips — ORCHESTRATOR NOT YET WRITTEN.
- [x] Grounding: engine claims carry line-refs (the plan); Phase-4 functional
      choice grounded via literature-review before Phase 5.
</preflight>
