---
id: localised-jellium-energy-book-keeping
area: localised_jellium_parameter_study_2
title: Energy book-keeping analysis
status: running
hypothesis: "The WP-classical total-energy difference in localised-jellium runs is fully accounted for by (i) the WP quantum self-energy (zero-point KE 3/(4 sigma^2), self-XC, self-Hartree) and (ii) the classical run's missing projectile-background Coulomb term E_proj_bg; after adding (ii) the ledger closes (residual <= 3 eV per row) at every radius, and for launched pairs dKin_WP-CL = KE_projectile + 3/(4 sigma^2)."
handover: docs/handovers/localised-jellium-energy-book-keeping.md
tasks:
  - { name: "A1 gate — periodicity 2 vs 3 ledger verdict (user)", done: true }
  - { name: "A2 gate — launched-pair 100 eV KE audit (3 sub-claims)", done: true }
  - { name: "A3 — semi-empirical far-field forensics (5 sub-checks)", done: true }
  - { name: "A4 — localisation-energy 3/(4 sigma^2) derivation", done: true }
  - { name: "A5 — WP effective-potential cutoff model + sigma scaling", done: true }
  - { name: "A6 — long-range-effect synthesis (advisor-converted; user resolves on reading)", done: true }
  - { name: "B1 — E_proj_bg (advisor: post-hoc, no re-runs) + exact d(H+E) decomposition", done: true }
  - { name: "B2 — SCF-with-projectile screening pair (ghost SCF r={4,12,28} + 83e)", done: true }
  - { name: "B3 — timestep-by-timestep ledger diff (qsp_phase3 pair + E_pb(t))", done: true }
blocked_reason: ""
---

# Energy book-keeping analysis

**Interactive gate-stopped campaign** (user decision 2026-07-10): the agent runs
until a task marked ⏸ completes, presents evidence NEUTRALLY (verdicts are the
user's, never the agent's), stops, and proceeds only on the user's decision.
The full autonomy checklist is relaxed accordingly; reproducibility, grounding,
and guard rails still apply.

## Question

Is the WP-vs-classical energy gap in localised-jellium runs fully explained by
known, calculable terms — or does a residual survive that constitutes a genuine
quantum effect? Decision informed: whether quantum-vs-classical stopping-power
differences are physical or book-keeping artefacts.

- **Success**: post-B1 residual |WP−CL − known terms| ≤ 3 eV per ledger row at
  every radius, and dKin_WP−CL ≈ KE_proj + 3/(4σ²) (≈ 100 + 82 ≈ 180 eV for the
  σ = 0.5, 100 eV pair) confirmed on a launched pair.
- **Failure**: a residual survives — then B3 localises which component and which
  timesteps carry it (that residual is the "quantum effect").

## Locked decisions

| Decision | Value | Why |
|---|---|---|
| Shape | Interactive, gate-stopped; Phase A (existing data / understanding) before Phase B (new runs) | user 2026-07-10 |
| Shape AMENDED | A3–A6 + B1–B3 run AUTONOMOUSLY; decisions that would need the user are made by a Fable 5 advisor agent (each ruling logged); deliverable = one executed notebook (what/why/results per task); the USER interprets the results afterwards | user 2026-07-11 |
| Default periodicity | 2, if any run must proceed before the A1 verdict | user draft |
| Closure tolerance | ≤ 3 eV residual per ledger row | user "Agree" 2026-07-11 to "a few eV" proposal |
| T3 scope | Supporting investigation; pass/fail rides on the ledger only | user "Agree" 2026-07-11 |
| Output placement | Analysis artefacts → `ResearchProject/systems/localised_jellium/hypotheses/campaign_autorun_study/`; `docs/notes/localised-jellium-parameter-study-2.md` is the USER's thinking file (agent writes only A6's resolution there, in the user's words) | user lock 2026-07-11 |
| Ledger comparison columns | dKin, dXC, summed d(H+E), WP−CL total only; raw dHartree/dexternal are charged-cell-convention-poisoned (−274 eV p2 vs −29 eV p3 at r=40) | notebook cell 39 |

## Key verified anchors (checked 2026-07-10/11, absolute paths)

- Ledger notebook: `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/campaign_autorun_study/theoretical_slab_model.ipynb`
  (p2 ledger cell 22; p3 mirror cells 31–39; both decompositions exact to ~1e-13 Ha).
- Insertion runs are **t=0, projectile at rest** (k0 = 0): dKin = 81.7 eV already
  matches 3/(4σ²) = 81.6 eV to 0.1% (cell 39). The +100 eV prediction therefore
  concerns LAUNCHED runs only (A2).
- Meeting deck: `/local/data/public/skcb2/tddft/docs/reports/09-07-2026-meetng-emilio/Meeting 7.pptx`.
- Semi-empirical chain (A3 target): notebook cells 11–19 →
  `docs/reports/09-07-2026-meetng-emilio/assets/make_s1_3_field_potential.py` →
  `ResearchProject/systems/localised_jellium/hypotheses/plate_model/build_plate_model.py` + `VALIDATION.md`.
- "w parameter" = erfc edge width `EDGE_WIDTH_BOHR` (shared configs; = 1.0 in
  production slabs; the h0 GSs used edge_width **0**, sharp). Existing w-sweep:
  `scripts/campaign_autorun/runs/h1/gs_w{0,0.25,0.5,0.75,1,1.5,2,3}`.
- GS references (both slabs identical: 50×50×120 Bohr, half-width 12.5, N=82,
  spacing 0.5, LDA, edge_width 0):
  p2 `runs/h2/gs_p2_lz120` E_GS = 60.38307052445239 Ha;
  p3 `scripts/h0_base_difference/gs` (= `runs/h2/gs_lz120`) E_GS = −108.5336851082701 Ha
  (checkpoint `shared_gs/slab_n82_L50x50x120`, per `orchestrate.py:27-28`).

## Tasks (⏸ = gate: stop for the user)

### Phase A — existing data / understanding (no new runs)

- **A1 ⏸ Periodicity verdict.** p2 vs p3 ledgers side-by-side, convention-free
  columns, r ∈ {4,12,20,28,36,40}. Done: table presented, user verdict recorded,
  periodicity locked downstream.
  → Artefact: `hypotheses/campaign_autorun_study/a1_periodicity_ledger_comparison.md`.
  **DONE 2026-07-11 — verdict "use p2 for now": periodicity 2 locked downstream.**
- **A2 ⏸ Launched-pair 100 eV audit.** Locate matched launched classical/WP pair
  (100 eV, σ = 0.5); per-step CSVs; test (a) classical ionic 100 eV absent from
  E_total, (b) WP 100 eV present in E_kinetic, (c) dKin_WP−CL ≈ 180 eV.
  Done: three sub-claims answered with numbers; user verdict.
- **A3 ⏸ Semi-empirical far-field forensics.** (i) line-by-line recheck of the
  chain above → bullet summary; (ii) enclosed charge Q(<z>) vs z from the GS VTI;
  (iii) density-spill check at box edges; (iv) plate-thickness (dz) sensitivity;
  (v) sharp vs erfc-softened (w = 1) background comparison (reuse h1 w-sweep).
  Done: sub-checks with numbers, causes ranked, user verdict on the cause.
- **A4 Localisation-energy derivation.** Derive 3/(4σ²) from the Gaussian WP's
  ⟨T̂⟩ (literature-review grounded); confirm 82 eV at σ = 0.5. Done: written note.
- **A5 WP effective-potential cutoff model.** Closed-form Gaussian-charge
  potential erf(r/(√2 σ_ρ))/r; define r_cut criterion; tabulate σ-scaling.
  Done: formula + plot + table.
- **A6 ⏸ Long-range-effect conversation.** Interactive; draws on A3 + A5 +
  Meeting 7 results; until the user declares it resolved. Resolution recorded in
  the user's notes file IN THE USER'S WORDS.

### Phase B — new runs (matrix numbers locked at the B1 gate)

- **B1 ⏸ E_proj_bg column + re-run sweep.** Surgical tracker of the classical
  projectile↔background Coulomb energy (dynamics untouched; new CSV column).
  PRE-GATED: code-test + formula-validation vs the analytic Gaussian-charge/slab
  integral + catalogue row BEFORE any re-run. Then re-run the classical vs WP
  insertion sweep at locked periodicity/radii. Done: closure table vs 3 eV; user
  hypothesis verdict.
- **B2 ⏸ SCF-with-projectile screening pair.** Converge SCF with WP present and
  with classical projectile present; compare screening via the established plots.
  Done: matched figures; user verdict.
- **B3 ⏸ Timestep ledger diff** (conditional on B1). Matched launched pair,
  per-component dE at every step, earliest steps finest. Done: attribution plot +
  table; user verdict on what is "the quantum effect".

## Rules

- ALWAYS present evidence neutrally at gates; the user owns every verdict.
- NEVER compare raw dHartree/dexternal across periodicities or cite them as
  self-energies (charged-cell G=0 convention).
- σ always means σ_WP; classical UPFs are generated at σ_pot = σ_WP/√2 but
  labelled by σ_WP (sigma-wp-convention rule).
- Round reported numbers to 2 s.f. (3 s.f. only for near-equal differences).
- New observables/kernels are pre-gated (code-test + formula-validation +
  catalogue row) before any expensive run.
- GPU is the default for all Phase-B runs (cudaMemGetInfo probe; NVML mismatch
  is not a blocker).

---

## Original draft notes (user, verbatim — 2026-07-10)

In this note, I am going to think carefully about the synthesised results. Then, I am going to come up with new experiments with the knowledge gained to experiment the hypotheses. 

Let's call this campaign - "Energy booking analysis"

## Energy book keeping analysis
For the following analysis, I am going to look at the results from "docs/reports/09-07-2026-meetng-emilio/Meeting 7.pptx"
Also, the previous analysis upon which these tasks have been designed have their results in theoretical_slab_model.ipynb. 

### Test the hypothesis that there is no (or minimal) difference between the periodicity 2 and 3 cases in enregy book keeping. 
In this task, I essentially compare the tables for periodicity 2 and periodicity 3 I've made here that has as the columns dE_WP (eV), dE_CL (eV), WP-CL (eV), dKin (eV), dXC (eV), d(U_H+U_ext) (eV). I can visually examine the two tables and determine on my own if these results are the same. 

This has to be done first, and then, I examine the results. If, it ends up happening that you are having to run the simulations without me deciding on this, I want the default for the next runs (and tasks) to be periodicity 2. 

### 
The enegy book keeping I am referring to is essentially a talbe I've made here that has as the columns dE_WP (eV), dE_CL (eV), WP-CL (eV), dKin (eV), dXC (eV), d(U_H+U_ext) (eV), U_proj_bg. In this plot, I am not sure why the delta K.E column showcases this value. We claim that this is the localisation energy. But, I was thinking, the classical projectile's energy of 100 eV does not show up in the total energy (test this claim). In the wavepacket case, I would assume the 100 eV K.E of the wavepacket would show up. So the difference between the classical and the kinetic energy values must be the energy of the projectile + the localisation energy [completed 2026-07-10].


### Semi empirical model for the total system
In this plot, we considered the density distribution of the electronic system and the positive background. Here, using the understanding of the analytical expression for an infinte plate, we found the total effective potential of the jellium plus the background. However, there is an important feature here. Consider the net charge distribution in the system. After a certain distance away from the center of the jellium slab, at about 15 bohr on either side, the total enclosed charge within this region is 0. Hence, by Gauss' law, the total field out of this region would be 0. Howevver, we see in the plot that this is not quite the case in teh semi empirical field. The field comes to a constant that is not 0. This might represent that there is some charge that is spilling out of this region. This needs to be investigated carefully. 

For this, I want you to recheck everything that has been done in this analysis, line of code by line of code. I want you to give me a summary of all the steps, clearly mentioned, and succicntly explained in bullet points. Then, perhaps, we should investigate if the effect of the thickness of the thin plate we consider is too big for this analysis. Then, I want you to integrate the simulation cell in plates, and calculate the net charge as a function fo the z axis. Then, we should make a sanity check as to where there is almost 0 density. Perhaps, its genuinely the case that there is some non-zero density flowing out of the box that extends until the end of the simulation box (up to the boundaries on either side from the jellium slab).

Perhaps, the impact of the w parameter must be considered.


### Fixing the energy booking by including the E_projectile_bg term
The enegy book keeping I am referring to is essentially a talbe I've made here that has as the columns dE_WP (eV), dE_CL (eV), WP-CL (eV), dKin (eV), dXC (eV), d(U_H+U_ext) (eV), U_proj_bg. The only energy that is not accounted for in the classical simulation is the coulombic energy term that is attractive in nature between the classical projectile and the positive background charge. So, in the simulations, I want to track this additional term and broadcast and save it in the energies. This way, the book keeping would be complete. This term would be added as an additional column in the table. 

I think, we do not need to change the actual energy that is being used. I only need to surgically include this term so that it can be saved in the csv file that we stored. Then, I need to re-run the classical vs wavepacket cases at different radii all over again. 


### Understanding the long range effects on different quantities in the results
In this task, we aim to understand why there is a long range effect that shows up in the electric field, electric potential and the energy book keeping d E_total (where we compare the wavepacket and the classical cases). I suspect this has to do something with the radial cutoff of the pseudopotential. However, I need to understand the physics behind it. So, I would want to engage in a conversation with a grounded model in all the results that we have form the presentation attached. Then, we are going to discuss physics, and test ideas, until I understand why this effect happens. 


### Investigate, what is the effective radius cutoff of the radial potential produced by the gaussian wavepacket. Does this change with sigma?
In this task, we are going to think about modelling the electrostatic potential of the wavepacket as a pseudopotential's coulombic potential. This way we understand the difference between the different width wavepackets. This would help make a mental model better. 

### Screening effects of wavepacket and classical porjectile
In this task, I want to compare the impacts of a classical projectile and a wavepacket on the GS jellium system. To do this, before making the energy and density measurements, I want to run SCF again. This ensures that the jellium system and the wavepacket in the quantum case, self adjust to each other's impact. The same can be done for the classical projectile case. Then, the screening effect using the plots that were made in the presentation (and the theoretical_slab_model.ipynb). This way, I understand exactly how the identical localised jellium system reacted to both these projectiles, and I would be able to make a comparison. 

### Understanding why 81.7 eV is the localisation energy of the 0.5 bohr wavepacket
In this task, I want to understand how one could come up with the value of localisation energy of 3/(4(sigma^2)). I want to understand the derivation of this to get an intuitive understanding for what this is. 

### (Future task) 
This task is conditional on the previous tasks of fixing the energy book keeping and the subsequent tasks. In this task, we keep a track of the energy book keeping each run. Then, we are going to find a run of a wavepacket (concentrated and with a well defined quantum stopping power) which also has a classical analogue. Then, we are going to examine the book keeping for the classical and the wavepacket runs at each timestep. Whenever, we observe some energy difference in the runs that we did not account for, then, we are going to analyse this. We are going to plot this, understand which components of total energy contribute to this energy difference. We do this timestep by timestep. We carefully analyse the first few timesteps, as these would be the most comparable. These inform us about specifically what's the quantum effects in these runs (meaning difference between the classical and the wavepacket runs).
