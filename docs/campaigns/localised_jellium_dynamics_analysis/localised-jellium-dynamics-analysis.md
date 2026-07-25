---
id: localised-jellium-dynamics-analysis
area: localised_jellium_dynamics_analysis
title: Localised jellium — complete the energy ledger (E_proj_bg) + projectile r_cut sweep
status: draft
hypothesis: "Adding the projectile–background Coulomb term E_proj_bg (= ∫ n_proj·v_bg) completes the classical energy ledger so d(U_H+U_ext+U_proj_bg) reproduces the measured WP−CL electrostatic gap; and sweeping the classical projectile pseudopotential radial cutoff r_cut isolates its effect in E_external (E_proj_bg,ideal is r_cut-invariant), revealing the WP's effective r_cut."
handover: docs/handovers/localised-jellium-energy-book-keeping.md
tasks:
  - { name: "inqkit projectile_background_energy capability + two observables.csv columns (code-test + formula-validation gate)", done: true }
  - { name: "Phase 1 — 6 classical insertion re-runs (r∈{4,12,20,28,36,40}) carrying E_proj_bg; ledger rebuilt with U_proj_bg closing WP−CL", done: true }
  - { name: "Phase 2 — 5 classical r_cut runs at fixed r=20 (r_cut∈{10,20,30,40,50}); energy-component + effective-r_cut table", done: true }
  - { name: "Phase 3 — open-z (p2) vs PBC (p3) at w=0: leaking-charge comparison across Lz∈{90,120,160,240}; reuse existing p2/p3 GS runs + add p3_lz120", done: true }
  - { name: "Phase 4 — semi-empirical E(z)/φ(z) for the w-sweep (analysis-only; already computed in semiempirical_spillout, consolidate here)", done: false }
  - { name: "Phase 5 — screening: RT-propagate WP + classical (at rest, r=12) from the same GS; total + induced (bath-only) density-difference GIFs", done: true }
  - { name: "campaign notebook (.ipynb): completed ledger + Phase-2 r_cut + Phase-3 boundary + Phase-4 w-field + Phase-5 screening analysis", done: false }
blocked_reason: ""
---

# Localised jellium — complete the energy ledger (E_proj_bg) + projectile r_cut sweep

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository.
</identity>

<description>
The classical-vs-WP energy ledger (task "Energy book-keeping", A1) is INCOMPLETE
for the classical projectile: INQ's fixed 8-term total has no slot for the Coulomb
interaction between the classical ghost (an ion/pseudopotential) and the jellium
positive background (a `localised_background_perturbation`, not an INQ ion). The WP
run has no such hole — the WP is electron density, so its background coupling is
already inside `E_external`. This asymmetry is the missing `U_proj_bg` column.

**Phase 1** adds a reusable inqkit capability that computes this term from the
perturbation's cached `v_bg`, exposes it as two `observables.csv` columns, and
re-runs the 6 classical insertion runs so the completed ledger
`d(U_H+U_ext+U_proj_bg)` reproduces the WP−CL electrostatic gap.

**Phase 2** fixes the projectile at r=20 Bohr and sweeps the pseudopotential radial
cutoff `r_cut` (the radius beyond which electrons stop feeling the ghost). Because
`E_proj_bg,ideal` uses the projectile's TRUE Gaussian charge it is r_cut-invariant,
so the r_cut effect is cleanly isolated in `E_external`. By matching the classical
`E_ext+E_H(+U_proj_bg)(r_cut)` to the WP reference we infer the WP's effective r_cut.

Success: (P1) the completed classical ledger closes the WP−CL gap to within the
campaign's 3 eV row tolerance at every radius; (P2) `E_proj_bg,ideal` is flat across
r_cut (< 0.1 eV) while `E_external` varies monotonically, and an effective WP r_cut
is read off. Failure: the completed ledger still leaves an unexplained > 3 eV gap,
or the two E_proj_bg formulations disagree at r_cut=50 (full UPF) beyond the
reciprocity tolerance.

**Phase 3** (density-leakage / boundary-condition thread — user's Q3). In the
semiempirical_spillout study the electron density leaks beyond the positive-background
region, and that leak *decreases with the background edge width w*. Open question:
"even at w=0, PBC vs open-z might have an impact on the leaking charge." Test it by
comparing, at **w=0**, the leaking charge between **open-z (p2)** and **fully-periodic
(p3)** cells at matched Lz ∈ {90,120,160,240}. Reuse the existing p2 {90,120,160,240}
+ p3 {90,160,240} GS runs and add the one missing pair member **p3_lz120**. Success:
a clean p2-vs-p3 leaking-charge table/overlay per Lz; failure = boundary condition
makes no measurable difference (which would itself answer the question).

**Phase 4** (semi-empirical fields for the w-sweep — user's Q4). Build E(z) and φ(z)
from each w-run's EMPIRICAL density via the sheet-stack method, to verify that a small
w removes the spurious far field (not just the density pile-up). **Already computed**
in the semiempirical_spillout notebook (far-field |E|: 0.098 eV/Bohr at w=0 → 0.000 at
w=1,2,4); this phase consolidates that result into the campaign notebook. Analysis-only,
no new runs.
</description>

<observables_set>
Per-run (t=0 insertion, static ghost — 2 steps, dt=0.01): the full energy
decomposition already streamed by campaign_autorun/classical (energy_total,
kinetic, hartree, xc, external, nonlocal, ion, ion_kinetic, exact_exchange, nvxc,
eigenvalues) PLUS the two NEW columns:
  - `energy_proj_bg_ideal` = ∫ n_proj·v_bg  (n_proj = Gaussian σ_pot=σ_WP/√2, ∫=1, at z_proj)
  - `energy_proj_bg_impl`  = −∫ n₊·v_ion     (as-implemented pseudopotential; = the B1 ∫n₊·v_ghost)
NEW code — pre-gated by code-test + formula-validation + a test-catalogue row
BEFORE any expensive run (see <tasks> task 0). No VTI cadence needed (static, t=0).
</observables_set>

<resolved_decisions>
geometry: p2 (periodic x,y, open z), Lx=Ly=50, Lz=120, half_width=12.5, N=82,
  spacing=0.5, LDA, edge_width=0 (sharp) — matches the A1 ledger exactly so dE_WP is
  reused unchanged. n0 = 82/(50·50·25) = 1.312e-3.
gs_source: shared checkpoint `scripts/campaign_autorun/runs/h2/gs_p2_lz120/checkpoint`
  (E_GS = 60.38307052445239 Ha, verified run_summary anchor).
projectile: classical Gaussian ghost, UPF `jellium/shared/pseudopotentials/
  electron_gaussian_wpsigma0p5.upf` (σ_pot=0.354=σ_WP/√2, z_valence=0, V(r)=+erf(r/0.5)/r,
  mesh r_max=50). σ_WP label = 0.5 (sigma-matching convention).
E_proj_bg formula: `E_proj_bg,ideal = ∫ n_proj·v_bg`, v_bg = the cached
  `localised_background_perturbation` potential (background_perturbation.hpp:65-68,
  the SAME field validated as infinite-plate-like); n_proj = normalised Gaussian of
  std σ_pot at the projectile z. `E_proj_bg,impl = −∫ n₊·v_ion`, n₊ from
  perturbation.background_density(), v_ion = poisson(atomic_pot.ionic_density) +
  atomic_pot.local_potential (self_consistency.hpp:102 recipe; electrons.atomic_pot()
  accessible at electrons.hpp:496). NEITHER enters energy_total (diagnostic only).
phase1_radii: r ∈ {4,12,20,28,36,40} (A1 ledger radii). launch_z = −(half+r).
  Reuse existing wp_r{r}_p2 for dE_WP (WP energy unaffected by the new term).
phase2: fixed r=20 (launch_z=−32.5); r_cut ∈ {10,20,30,40,50}. UPFs:
  `campaign_autorun/cutoff_test/upfs/electron_gaussian_wpsigma0p5_rc{10,20,30,40}.upf`
  + the full UPF (= rc50). Reuse wp_r20_p2 as the WP reference. Truncation = potential
  absent beyond r_cut (make_cutoff_upfs.py), z_valence=0 → no long-range part.
phase3 (open-z vs PBC, w=0): reuse GS runs under
  `scripts/semiempirical_spillout/runs/` — p2: {lz90,lz120,lz160,lz240}; p3:
  {p3_lz90,p3_lz160,p3_lz240}. ADD p3_lz120 (semiempirical_spillout/gs/run, LJ_PERIODICITY=3,
  LJ_LZ=120, N=82, EDGE_W=0, es=20) to complete the matched {90,120,160,240}×{p2,p3} set.
  Leaking-charge metrics (from the semiempirical_spillout notebook helpers): spill beyond
  |z|>12.5, near-edge/inter-slab density within 4 Bohr of |z|=Lz/2, enclosed-charge deficit
  Q(|z|<25). p2 = open z (Rozzi slab-truncated Poisson, solvers/poisson.hpp:188-208); p3 =
  periodic (image slabs at ±Lz, no free z boundary).
phase4 (w-sweep E/φ, analysis-only): reuse the p2 w-sweep GS runs {lz160(w0),w1,w2,w4} and
  the ALREADY-COMPUTED semi-empirical E(z)/φ(z) (sheet-stack, each run's empirical density +
  own background). Result to consolidate: far-field |E| = 0.098 (w0) → 0.000 (w1,2,4) eV/Bohr;
  peak vacuum |E| 0.30→0.15→0.076→0.045. NO new runs.
phase5 (screening, RT): SAME GS gs_p2_lz120 for both. Two RT runs from it — WP (injected orbital,
  k0=0) and classical (ghost as static external potential, v=0), BOTH AT REST at r=12
  (launch_z=−24.5), p2, Lz=120, σ_WP=0.5. Propagator ETRS (CN renormalises the WP each step,
  breaking the picture — reference_inq_propagator_mask_absorber). Default kinematics locked:
  N_STEPS=500, dt=0.01 (5 a.u. ≈ 0.1 plasmon period ω_p=0.128 a.u.), save total-density VTI every
  25 steps → ~20 GIF frames. Diagnostics per frame: (a) TOTAL diff = n_total(WP)−n_total(CL)
  (includes the WP orbital blob); (b) INDUCED diff = [n_bath(WP)−n_GS] − [n_CL−n_GS], where
  n_bath(WP)=n_total(WP)−n_WP_orbital (canonical bath density, reference_canonical_bath_density) —
  isolates the BATH screening response. Both rendered as GIFs (shared colorbar for compared panels,
  feedback_shared_colorbar_rule; load via inqview.load_vti, NEVER fftshift a VTI). N_STEPS is the
  only real compute knob — raise later for a longer GIF if the onset is too brief.
file_placement (ADR-0007): run machinery in
  `scripts/localised_jellium_dynamics/{phase1,phase2}/`; new inqkit header in
  `inq-stack/include/inqkit/jellium/projectile_background_energy.hpp`; wrapper test in
  `inq-stack/tests/include/inqkit/jellium/`; analysis notebook + builder in
  `hypotheses/localised_jellium_dynamics/`. Phase 3/4 reuse semiempirical_spillout runs +
  notebook helpers (do NOT duplicate the GS matrix).
</resolved_decisions>

<guard_rails>
- ABORT (correctness): NaN / complex energy; missing GS checkpoint; a GPU occupied
  by another user (warn via cudaMemGetInfo probe — NVML is broken, compute still works).
- CODE GATE (blocks all runs): the projectile_background_energy known-case test MUST
  pass first — `_ideal` matches the analytical infinite-plate estimate for the σ=0.5
  ghost at r=12 to a stated tolerance, AND `_impl` matches the B1 post-hoc
  `∫n₊·v_ghost` (b1_decomposition.py) to < 1 eV; at r_cut=50 the two formulations
  agree by reciprocity to the same tolerance. Fail ⇒ do NOT launch Phase 1/2.
- These are cheap static t=0 insertion runs (2 steps) — no wake/traversal sizing, no
  transient exclusion, no CAP. Not light-projectile stopping (the ghost is at rest).
- Cost: 11 classical runs total (~a few min each on GPU). No budget concern.
</guard_rails>

<tasks>
0. **inqkit capability + observables columns (CODE GATE).** Add
   `inqkit::jellium::projectile_background_energy` computing `_ideal` and `_impl`
   from the cached v_bg / background_density + a constructed Gaussian n_proj. Extend
   `StepContext` + `ObservableSelection` + `ObservablesWriter` with two columns
   (default 0, backward-compatible) and a hook for the run to inject the computed
   (static) values. code-test: write→known-case-test→confirm; formula-validation
   agent on the ∫n_proj·v_bg formula vs its source; add a test-catalogue row.
   Done = test passes, columns appear, no existing-run regression.
1. **Phase 1 — completed ledger.** Re-run 6 classical insertions (r∈{4,12,20,28,36,40})
   from gs_p2_lz120 with the two new columns. Reuse wp_r*_p2. Rebuild the ledger table
   (dE_WP, dE_CL, WP−CL, dKin, dXC, d(U_H+U_ext), U_proj_bg columns). Done = the
   completed classical channel closes WP−CL to ≤ 3 eV at each radius (or the residual
   is quantified and explained).
2. **Phase 2 — r_cut sweep at r=20.** 5 classical runs, r_cut∈{10,20,30,40,50}, fixed
   r=20, same geometry, all components + both E_proj_bg columns. Reuse wp_r20_p2. Done =
   a table of energy components vs r_cut, `_ideal` flat (<0.1 eV) while E_external
   varies, and an inferred effective WP r_cut (the r_cut where classical E_ext+E_H
   matches the WP).
3. **Phase 3 — open-z vs PBC at w=0.** Reuse existing GS runs
   `scripts/semiempirical_spillout/runs/{lz90,lz120,lz160,lz240}` (p2) and
   `{p3_lz90,p3_lz160,p3_lz240}` (p3); add the one missing pair member **p3_lz120**
   (same GS variant, LJ_PERIODICITY=3, Lz=120, N=82, w=0, es=20). Compare the leaking
   charge — spill beyond the slab (|z|>12.5), near-edge/inter-slab density, enclosed-
   charge deficit — between p2 and p3 at each matched Lz. Done = a p2-vs-p3 leaking-
   charge table + n_e(z) overlay per Lz, with a stated whether-PBC-changes-the-leak verdict.
4. **Phase 4 — semi-empirical E(z)/φ(z) for the w-sweep (analysis-only).** Already
   computed in `hypotheses/campaign_autorun_study/semiempirical_spillout.ipynb` (the Q4
   cells, using each w-run's empirical density + sheet stack). Consolidate the E(z)/φ(z)
   panels + the far-field-|E|-vs-w table into the campaign notebook. Done = the panels +
   table appear with the w=0→w=1 far-field collapse shown. No new runs.
5. **Phase 5 — screening (RT, at rest, r=12).** Two ETRS RT runs from gs_p2_lz120 — WP
   (injected orbital) + classical (ghost) — N_STEPS=500, dt=0.01, total-density VTI every
   25 steps. Compute TOTAL diff (n_WP−n_CL) and INDUCED diff (bath-only, subtract the WP
   orbital) per frame; render both as GIFs (shared colorbar). Done = both runs complete
   (energy finite, no NaN), the two GIFs render, and a t=0-vs-final screening panel is saved.
6. **Campaign notebook (.ipynb).** Per notebook-making: completed ledger + Phase-2
   r_cut + Phase-3 boundary + Phase-4 w-field + Phase-5 screening (GIF frames + panels),
   executed auto after the runs. Done = 0-error execution with all tables + figures.
</tasks>

<rules>
- ALWAYS: E_proj_bg is DIAGNOSTIC — never added to energy_total or the SCF/dynamics.
- ALWAYS: label the classical run and its WP partner by σ_WP=0.5 (sigma-matching);
  σ_pot=σ_WP/√2 appears only in the n_proj construction / UPF methods note.
- NEVER: modify inq/ or inq-study/ — the capability is a wrapper-only inqkit header.
- Autonomous executor = a PYTHON orchestrator (idempotent resume, per-phase try/except
  + failure email), NOT bash (2026-06-27 rule); reference campaign_autorun/orchestrate.py.
</rules>

<preflight>
- [ ] Intent self-contained: hypothesis + P1/P2 success/failure criteria above; each task has a done-criterion.
- [ ] Setup reproducible: geometry/N/box/GS/UPFs/radii/r_cut all locked with values in <resolved_decisions>.
- [ ] New code pre-gated: task 0 (code-test + formula-validation + catalogue row) BLOCKS Phase 1/2.
- [ ] Guard rails: correctness aborts + the code gate numeric criteria; cheap static runs, no pilot sizing needed.
- [ ] Autonomous mechanics: cudaMemGetInfo GPU probe (warn if occupied); Python orchestrator; per-phase Gmail; notebook auto-built; this handover updated + frontmatter flipped as tasks complete.
- [ ] Grounding: engine claims carry source line-refs (background_perturbation.hpp:65-68, self_consistency.hpp:102, electrons.hpp:496); formula grounded via formula-validation.
</preflight>
