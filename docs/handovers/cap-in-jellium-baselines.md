# Handover: CAP-in-jellium baselines

Plan: `/local/data/public/skcb2/tddft/docs/plans/cap-in-jellium-baselines.md` (authoritative build order + test matrix).
Spec (resolved decisions in XML tags): `/local/data/public/skcb2/tddft/docs/prompts/cap_in_jellium/baseline_runs.md`.
Glossary: "CAP in jellium baselines" in `/local/data/public/skcb2/tddft/CONTEXT.md`.

> Running AUTONOMOUSLY from 2026-06-17 (user stepped away). Email per baseline to
> chiddukanna@gmail.com. ALL absorption numbers PROVISIONAL until Task #7.

## Milestone: 2026-06-22 — all 7 qvc runs DONE; `run-notebook` skill created + run-notebooks built

**All 7 runs complete** (Study A R4/R1/R2 + Phase 2 benchmark E150/300/450/600), all
`run_completed=true`, `nan_seen=false`. CAP runs drain ~37%; no-CAP benchmark conserves
N to ~1e-7; classical projectile decelerates (stopping signal present).

**NEW skill `run-notebook`** (grilled, then built): a deep SINGLE-RUN analysis notebook
(sibling of the study notebook). Skill-local + shippable:
`.claude/skills/run-notebook/{SKILL.md, run_notebook_builder.py}`. It is an ASSEMBLER over
`inqview.pipeline` (runs the phases → embeds figures) + adds density-matrix carpets, the
lead xz GIF, and the loss-function note. Auto-gated by run-type (WP/classical/baseline).
Glossary term in CONTEXT.md ("run-notebook"). Grill decisions: (a) 2D k_z–k_⊥ scattering
+ bath/WP real-space density deferred to a FUTURE observable upgrade (runs save 1D |k| and
total density only); carpets + ONE GIF format; loss function always produced w/ low-res
note; WP orbital shows TOTAL ⟨ψ|H|ψ⟩ only; the full standard battery accepted.

**Run-notebooks BUILT (all 0 errors):**
- `hypotheses/qvc_cap_sigma3/`: run_qvc_wp_s3_E300 (66 cells, 51 imgs), run_qvc_cl_s3_E300
  (54, 42), run_qvc_b1_s3 (53, 41).
- `hypotheses/qvc_nocap_sigma3/`: run_qvc_bench_cl_s3_E{150,300,450,600} (53 cells each).
Build cmd: `PYTHONPATH=.../inq-stack/python venv/bin/python3 run_notebook_builder.py
<results_dir> <out.ipynb> [--baseline <dir>] [--run-cpp <path>]`.

**STILL TODO:** (1) coherent-peak momentum tracker (the clean WP stopping observable;
½⟨k²⟩ unusable); (2) the SWEEP study notebooks for Study A and the S(v) benchmark (vs
Lindhard); (3) the future observable upgrade (2D momentum + |ψ_WP|²); (4) Study B (no-CAP
twin) + R3 vacuum SIE control on user command.

## Milestone: 2026-06-21 (night) — Study A (WITH-CAP σ=3/300eV) LAUNCHED; Study B (no-CAP) planned

Plan: `docs/plans/quantum-classical-sigma3-stopping.md` (both studies).

**Study A WITH-CAP — RUNNING autonomously** (PID 824203, GPU 1, detached nohup):
`scripts/cap_baselines/qvc_dispatch_phase1.sh` runs R4 (b1 baseline) → R1 (b3 WP) →
R2 (b2 matched classical), sequential, ~3 h, emails per run + all-done (family
`[qvc-cap-s3-E300]`). Config: 50³, N=162, reuse GS; CAP L=20/η=−1.0; launch z0=−6.5;
V0=4.696 (300 eV); N_STEPS=336 (τ=6.72 a.u.=rigid full-exit); WRITE_EVERY=5. Matched
UPF `electron_gaussian_wpsigma3p0.upf`. Outputs `results/qvc_{b1_s3,wp_s3_E300,cl_s3_E300}/`.
**Binary smoked OK** (b1/b2/b3, exit 0, no NaN, WP injected idx=100 max_overlap 1.5e-7).
run.cpp gained env hooks `CAP_WP_SIGMA`, `CAP_PROJ_PSEUDO`; rebuilt vs inq-study.
Generic emailer `email_run.py`.

**R3 vacuum-WP SIE control DEFERRED** (needs new single-electron/no-bath code; periodic
neutralizing-background subtlety). Bounds the ~7 eV SIE; required before any quantitative
quantum-component claim.

**Phase 2 (no-CAP classical S(v) benchmark) QUEUED** (watcher PID 2055543, GPU 1,
`qvc_dispatch_phase2.sh`): waits for Phase 1, then runs matched-width (charge std 2.121,
no CAP, CAP_ETA=0) classical at 150/300/450/600 eV, launch z0=−16, ~1.5 box periods each,
emails per run + final ALL-DONE. (Was deferred last turn; set up now per user.)

**Study B (NO-CAP 300 eV twin) — written as a PROMPT, NOT run** (user: ready-to-run on
command): `docs/prompts/quantum_classical_nocap/run_nocap_sigma3.md`. Config = Study A
with CAP_ETA=0, launch **z0=−16** (3σ_wp from box edge — user chose 3σ; unified σ=σ_wp),
N_STEPS=532 (1 box period, ~32-Bohr clean traversal). Runs BW (WP) + BC (classical), full
suite + analysis + notebook `hypotheses/qvc_nocap_sigma3/`. R3 vacuum control: NOT built
(user said no).

**NEXT:** (1) build the coherent-peak momentum post-processing kernel (+test) — the
stopping observable (½⟨k²⟩ unusable); (2) notebooks for A (+ Phase 2) once runs land;
(3) Study B + R3 on the user's command.

## Milestone: 2026-06-21 (late) — σ-convention UNIFIED in code; quantum-vs-classical design decisions LOCKED

User decisions (grill): (1) reuse the **50³ cubic box, N=162, existing GS** (no new GS);
(2) σ-convention **unified, WP as truth**; (3) WP **σ_wp=3** paired with a width-matched
classical (option (a), made exact by the code change); (4) document everything, keep
legacy runs as-is (option (a) relabel scope).

**CODE CHANGED + TESTED (`inqview.io.gaussian_psp`):** `generate_gaussian_psp(sigma_wp)`
now builds its erf charge at std σ_wp/√2 internally → classical & WP at the SAME σ present
the identical cloud. `GaussianPspResult` carries `sigma_wp`+`sigma_charge`. Tests updated,
**8/8 pass**. Generated + validated the matched UPF
`shared/pseudopotentials/electron_gaussian_wpsigma3p0.upf` (σ_wp=3, charge std 2.121).
Full convention doc + **legacy registry** (every old charge-std UPF/run → unified σ_wp×√2)
in CONTEXT.md "σ-convention unification".

**Fixed existing mislabel:** cap_baselines **B2-vs-B3 is NOT width-matched** (B2 classical
= unified σ_wp 0.707; B3 WP = σ_wp 0.5; √2 apart). Caveat added to the notebook builder
B2-vs-B3 section; notebook rebuilding. Do not read B2/B3 as a matched pair.

**Agent verdict on the new run (GO-WITH-CHANGES), still to fold into the spec:** CAP
≥10/side at η≈−1.0 (NOT 7.5; leak-through, not reflection, is the risk for the fast
projectile); energy a **200–400 eV ladder** (500 eV minimizes the quantum effect); the
**vacuum-WP SIE control is the measurement** (~7 eV SIE ≫ ~0.4–2 eV signal) — mandatory
primary run per (E,σ); stopping observable = WP **coherent-peak** momentum (sub-bin
centroid fit, ≥30–50 momentum dumps over the transit).

**NEXT:** fix the projectile energy (geometry in the 50³ box: CAP 10/side → free |z|<15,
launch 4σ_density≈8.5 from edge → transit ≈21.5 Bohr; recompute spread vs E), then the
classical S(v) benchmark high-E points at the matched width, then run + notebook.

## Milestone: 2026-06-21 (eve) — quantum-vs-classical stopping DESIGN grilled (2 fresh-agent validations)

Design discussion (grill-with-docs) to extend the baselines toward a **quantum
component of stopping** (S_WP − S_classical). Two independent fresh-context
validations run. Key resolved facts (all in CONTEXT.md "CAP in jellium baselines"):

- **ΔE total-energy subtraction REJECTED as a stopping observable.** During the
  Hermitian transit, stopping energy is *redistributed* into bath excitations that
  stay in E_total → ΔE(t) is flat-in-the-clean-limit; its drop is differential CAP
  drainage, not stopping. Use the projectile's mechanical (drift) KE instead:
  classical = −dKE/dz (ion); WP = **coherent-peak** of n_wp(k,t), NOT the 2nd moment
  ½⟨k²⟩ (scattering-inflated; verified non-monotonic on B3: 182→185→140 eV).
- **σ-matching convention LOCKED (user):** keep WP σ_WP; classical UPF at
  σ_pot=σ_WP/√2 (charge std == WP density std). Recorded in tddft-simulations skill
  Phase 2d′, CONTEXT.md, and memory `reference_sigma_matching_convention`. The
  shipped `electron_gaussian_sigma0p5.upf` is a real erf Gaussian charge (σ_pot=0.5),
  NOT a point charge (agent corrected an earlier mis-reading).
- **Spreading-suppression at σ≈0.5 is IMPOSSIBLE** (uncertainty-conjugate; needs
  ~0.2–0.8 MeV / relativistic). Point-like (σ_WP≲0.7) and non-spreading (σ_WP≳2.1 at
  10%) regimes are DISJOINT → an electron WP cannot reach the classical point limit.
  Salvage: **matched-pair σ-scan + extrapolate ΔS(σ)→0** (see plan).

### SELF-INTERACTION ERROR — decision + FUTURE TODO (user 2026-06-21)
- The WP carries a **~7 eV residual SIE** (self-Hartree 21.7 eV − LDA-x 14.7 eV) the
  classical ion lacks. **DECISION: do NOT correct it in the simulations for now.**
- **FUTURE TODO (record, don't drop):** add an SIE-correction path (SIC functional
  or vacuum-WP subtraction) before any *quantitative* quantum-component claim.
- **ANALYSIS USE:** every S_WP − S_classical "quantum component" MUST be reported
  with the SIE flagged + bounded via a **vacuum-WP control** (same WP, no bath).

## Milestone: 2026-06-21 — ALL 4 baselines complete + in the notebook (B2 added)

**All runs DONE** (`results/*/run_summary.txt`, all `run_completed=true`, `nan_seen=false`):
| run | η | N₀→N_final | absorbed | wall |
|---|---|---|---|---|
| `b1_eta0p50` | −0.5 | 162→4.61 | 97.2% | 12.4 h |
| `b1_eta0p10` | −0.10 | 162→7.28 | 95.5% | 10.6 h |
| `b2_classical_E100` | −0.5 | 162→4.58 | 97.2% | 30.2 h |
| `b3_wp_E100` | −0.5 | 163→4.59 | 97.2% | 12.0 h |

**B2 (classical) ADDED to the notebook** (`cap_baselines_study.ipynb`, now **96 cells,
0 errors**). New B2 block placed B1→**B2**→B3 (natural baseline order):
- **Stopping (headline):** projectile decelerates v_z 2.711→2.540 (KE 3.675→3.226 Ha);
  lost **12.2 eV total**, **6.6 eV over the clean first traversal** (27.6 Bohr) ⇒
  **S ≈ 0.240 eV/Bohr** at 100 eV. KE(0)=100.00 eV exactly (sanity check).
- **CAP does NOT absorb the classical ion** (it's in `ions`, not density) → flies
  through the periodic box 6.9× to z_final=346 Bohr. Clean measurement = transit only.
- B2 sections mirror B3: stopping v_z/KE, trajectory (abs+folded), ΔE=E_B1−E_B2,
  energetics, ∫J_z, density/wake/E_z GIFs — full + transit (T*=10.33). Plus a
  **B2-vs-B3 ΔE comparison** cell before the takeaway.
- `precompute_b2.py` (mirror of precompute_b3.py): reuses **cap_b3_clim.json** so
  B1/B2/B3 GIFs share one colour scale; overlays the EXACT projectile z(t) (folded)
  as a green tick. Artefacts `fig_b2_*.{png,gif}` (+ `_transit`).
- **HONESTY FIX:** ΔE=E_B1−E_B2 is confounded by CAP drainage (both baths bleed tens
  of e⁻; projectile rearranges the bath → drainage-trajectory difference ~7.6 Ha
  swamps the eV-scale deposition). Reframed B2.3 + TB2.2 + B3 Result-6 markdown to say
  ΔE is NOT clean "energy deposited"; the robust loss is the projectile KE drop.

Notebook path: `ResearchProject/systems/jellium/hypotheses/cap_baselines/cap_baselines_study.ipynb`.

**STILL DEFERRED (unchanged):** current-density FIELD VTI + ∮J_z flux screens + the
9 plane screens; total-system n(k). A clean single-pass classical stopping number
wants a non-periodic launch or stopping the ion at the first CAP crossing.

---

## Milestone: 2026-06-18 — B1 analysed (free region survives wake window); B2+B3 LAUNCHED

**B1 COMPLETE + analysed.** Region-resolved (`hypotheses/cap_baselines/build_b1_drainage.py`,
`cap_b1_region_drainage.csv`, `fig_b1_region_drainage.png`, `fig_b1_density_carpet.png`):
free region [−15,15] holds 97.2 e⁻ at t=0; **survives the ~10 a.u. projectile
transit** — η=−0.5: 98%/89%/73% at t=5/10/15; η=−0.10: 99%/94%/85%. The 95–97%
TOTAL drainage is slabs + late-time collapse; the wake window is usable. **B2/B3
viable with B1 subtraction.** Emailed user w/ figures (msg 178175968597…).

**run.cpp now does b1/b2/b3** (one binary, env CAP_MODE). b2=classical σ=0.5 e⁻
(electron_gaussian_sigma0p5.upf, "H"/m_e), insert z₀=−13, v=2.711, `.ehrenfest()`,
electron_track.csv. b3=σ=0.5 WP inject z₀=−13, k0=2.711, + MomentumDistribution.
Both: CAP η=−0.5 + full suite. **Both SMOKED (3 steps, exit 0):** b2 projectile at
z=−13 v=2.711 moves +z, track OK; b3 N0=163 (162+WP), WP injected, no NaN.

**B2 + B3 LAUNCHED ~06:26 2026-06-18, η=−0.5, 7000 steps:**
- GPU0 `b2_classical_E100` — PID 3405344 — ~16.9 s/step → **~33 h**.
- GPU1 `b3_wp_E100` — PID 3405345 — ~5.2 s/step → **~10 h**.
Launcher `run_baseline_launch.sh` (nohup) emails on completion via
`email_on_done.py`. Logs `b{2,3}_*.log`.

**B3 (WP) DONE + ADDED TO NOTEBOOK** (2026-06-18): N0=163, absorbed 97%, no NaN,
305 frames. Notebook now 46 cells, 0 errors. New B3 section (Results 6–13):
energy difference ΔE=E_B1−E_B3 with t*=10.3 a.u. transit line; B3 energetics;
integrated current ∫J_z (with explicit note that the current-density FIELD is
deferred); momentum n(k) before/after (total + WP); WP centroid track; + 3 GIFs
(total density, wake=B3−B1, E_z field) via `precompute_b3.py`. Builder extended.
B2 still running.
TRANSIT-WINDOW section added (user-requested): all B3 analyses recomputed on
t<=t*=10.3 a.u. (WP reaches free-zone edge z=+15). Cells T6–T13 mirror Results
6–13 with the mask + transit GIFs (`precompute_b3.py` is now window-aware via
`T_MAX`/`SUFFIX` env; `_transit` artefacts). Notebook now 64 cells, 0 errors.

SHARED-CLIM FIX (user-requested 2026-06-18): all GIFs now use ONE fixed colour
scale per quantity, computed once from the full window (`cap_b3_clim.json`: density
[0,1.38e-3], wake ±7.9e-5, E_z ±0.063) and reused for the transit GIFs + the B1
density GIF — so full vs transit and B1 vs B3 are directly comparable. WP is
clipped (bath/response dynamics stay visible; the wake GIF isolates the WP-induced
response). `precompute_b3.py` writes/loads the clim; `gif()` takes explicit clim.
Codified as **report-figures skill production rule 7** (shared scale across
compared figures AND animation frames; set clim once at imshow, only set_data per
frame). Notebook 64 cells, 0 errors.

**B1 STUDY NOTEBOOK BUILT** (user-requested 2026-06-18):
`hypotheses/cap_baselines/cap_baselines_study.ipynb` (executed, 0 errors, 25 cells,
4 figs) via `build_cap_baselines_report.py` (house narrative: question → jellium
scales one-per-cell → setup → sources → 4 results → takeaway; PROVISIONAL flagged,
De Giovannini 2014 cited). Covers B0/B1 + the drainage finding. When B2/B3 land,
EXTEND this builder for them and rebuild (add the auto-build tail to
`run_baseline_launch.sh` — not yet wired).

**NEXT (resume when BOTH done — ~33 h):**
1. Verify `results/b{2,3}_*/run_summary.txt` (completed, nan_seen).
2. Post-process: projectile track → stopping (dKE/dt); B3 COD track; wake =
   (B2|B3 density) − (b1_eta0p50 density) frame-aligned; E-field on wake frames;
   momentum n(k) before/after (B3); per-orbital energies / occupations.
3. **Build `hypotheses/cap_baselines/cap_baselines_study.ipynb`** (notebook-making
   skill: context → formulas → setup → linked sources → results → takeaway),
   covering B0(GS ref)/B1/B2/B3 + the drainage finding. Auto-build script
   `build_cap_baselines_report.py` (NOT yet written).
4. Final email (all 4 baselines).

**STILL DEFERRED (new code + tests):** current-density FIELD VTI (covariant→
Cartesian) + ∮J_z flux screens + the 9 plane screens; total-system n(k). E-field
kernel + per-orbital energy + occupations + momentum-dist(WP) already DONE.

---

## Milestone: 2026-06-17 (eve) — PILOT PASSED + KEY FINDING (bath over-drains at η=−0.5)

**Pilot (100-step B1, η=−0.5) PASSED the engine gate:** N₀=162 → N_final=127.25,
no NaN, energy real+finite throughout (−11.9 → −1.86 Ha), "ended normally".
**The inq-study CAP works in interacting LDA jellium.** Pilot output:
`scripts/cap_baselines/results/pilot_b1/` (pilot.csv + run_summary.txt).

**CRITICAL PHYSICS FINDING:** the CAP **absorbed 21.5% of the WHOLE bath in 2.0 a.u.**
Checked physical, not a bug: CAP slabs = 40% of box volume (~65 bath e⁻ sit
permanently inside the absorber); η=−0.5 removal timescale ~1–2 a.u. → ~half the
slab density gone in 2 a.u. **The vacuum-optimal CAP (η=−0.5, L=20) is too
aggressive for a jellium bath** (in vacuum the CAP only saw a transient packet).
OPEN decisive question (needs region-resolved density): is it mostly the SLAB
bath draining (fine — it's inside the absorber, never observed) or is the FREE
region [−15,+15] collapsing too (would confound B2/B3 wake)? → answered by the
full-B1 density VTI (region-N free vs slab, post-processed).

**In flight now (2 GPUs):**
- η-comparison sweep (reuses protected binary `scripts/cap_baselines/run_pilot_validated`):
  b1, 100 steps, η ∈ {−0.05,−0.10,−0.20,−0.30} → drainage(η) to pick a usable
  jellium η. Logs `scripts/cap_baselines/eta_*.log`, out `results/eta_compare/`.
  Driver `eta_compare.sh`, `eta_compare.out` (look for ETA_COMPARE_DONE).
- Full-B1 binary REBUILD (`build_b1.log`): `run.cpp` rewritten to the full suite
  (density VTI system series @WRITE_EVERY, StateEnergyWriter=⟨ψ|H|ψ⟩+var,
  OccupationsWriter, DensityDelta, observables.csv, per-step N(t)
  `electron_number.csv`) + the CAP. Based on the proven full-suite template
  `run_wp_n162_L50_E100/run.cpp` (which ALREADY wires nearly every observable I
  thought was new — state energies, occupations, momentum dist, density VTI).

### UPDATE (later eve) — drainage(η) mapped, full-B1 binary built+smoked, TWO B1 runs LAUNCHED

**drainage(η)** (100 steps = 2.0 a.u., reuses protected binary): −0.05→3.7%,
−0.10→6.9%, −0.20→12.1%, −0.30→16.1%, −0.50→21.5%. Sub-linear; slab bath drains
at every η. `results/eta_compare/eta_*/`.

**Full-B1 binary** (`scripts/cap_baselines/run.cpp`, the enhanced full-suite
version; binary `./run`) BUILT + runtime-SMOKED (3 steps, exit 0, all writers OK:
density VTI frames, state_energies.csv with E_expect_ha+variance, occupations,
observables.csv, electron_number.csv, eigenvalues). Smoke: `results/smoke_b1/`.
NOTE minor: electron_number.csv has a duplicate step-0 row (init + first callback);
dedup in post.

**TWO full Baseline-1 runs LAUNCHED ~17:43, ~5.7 s/step → ~11 h each, finish
overnight 2026-06-18:**
- GPU0 `b1_eta0p50` (η=−0.5, chosen config) — PID 3312782
- GPU1 `b1_eta0p10` (η=−0.10, gentler alt) — PID 3312783
Each: CAP_N_STEPS=7000 CAP_WRITE_EVERY=23 (~300 VTI frames). Launchers
`run_b1_launch.sh` (nohup) EMAIL the user on completion via `email_on_done.py`
(`inqview.email.send_run_email`, family `[cap-jellium-baseline]`). Logs
`b1_eta0p{50,10}.log`. Status email to user already sent (msg
178171470218…). **NOTE: these launcher PIDs are detached nohup — NOT
harness-tracked, so they will NOT auto-re-invoke the agent.**

**NEXT (resume after B1 lands — ~11 h):**
1. Read both `results/b1_eta0p{50,10}/run_summary.txt` (absorbed_frac, completed).
2. Post-process: region-N (free [−15,15] vs slabs) + n(z,t) from density VTI;
   E-field via `inqview.analysis.electric_field` on a few frames; COD. THE
   decisive output: does the FREE region survive at η=−0.5, or only at −0.10?
   → choose the production η for B2/B3.
3. Implement b2 (classical: copy `run_classical_n162_L50_E100/run.cpp`, insert e⁻
   pseudo-ion at z₀=−13, `.ehrenfest()`, add CAP) and b3 (WP: copy the
   `run_wp_n162_L50_E100/run.cpp` template, σ=0.5, z₀=−13, add CAP). Both reuse
   the full writer suite already proven here. Run at the chosen η (both GPUs).
4. dispatcher + auto-build `hypotheses/cap_baselines/cap_baselines_study.ipynb`
   (notebook-making skill). Email per baseline.
5. Deferred new observables (need careful new code + tests, NOT yet done):
   current-density FIELD VTI (covariant→Cartesian via cell metric) + ∮J_z flux
   screens; total-system n(k). per-orbital energy + occupations ALREADY done
   (StateEnergyWriter/OccupationsWriter). The 9 plane/flux screens NOT yet wired.

---

## Milestone: 2026-06-17 — grill complete, E-field kernel locked, pilot launched

### Design (LOCKED via grill — see prompt XML tags + CONTEXT.md)
- Baselines: **B0** pure-jellium GS (reuse `checkpoints/gs_L50_cubic_N162_dx0p40`,
  verified pure jellium — no projectile in the GS build); **B1** CAP, no
  projectile (drainage ref); **B2** CAP + classical σ=0.5 e⁻ @100 eV; **B3** CAP +
  σ=0.5 WP @100 eV. B1–B3 share ONE 140-a.u. (~7000-step, dt=0.02, ETRS) window.
- Geometry: 50-Bohr cubic, two-sided sin² CAP 10 Bohr/side, η=−0.5 → slabs |z|∈
  [15,25]; free [−15,+15]. Fractional: `absorbing(η, ±0.4, 0.2)` via `+` (sum).
  Launch z₀≈−13 (B2/B3), exit through far CAP. Propagator ETRS (INQ default; NEVER CN).
- Aim of B1 = CHARACTERISE bath drainage (not minimise). Continuity/CAP-sink
  check: dN_free/dt vs ∮J·dA at the CAP-edge flux screens.

### DONE
- **E-field kernel (idea 1) — LOCKED.** `inq-stack/python/inqview/analysis/efield.py`
  (FFT Poisson, native atomic units, optional SI; G=0→0 = neutralizing background
  → field of δn). Registered in `analysis/__init__.py`. 4/4 known-case tests PASS
  (`inq-stack/tests/python/inqview/analysis/test_efield.py`: uniform→0; sinusoid
  analytic E_z=−(4πA/G₀)sin(G₀z); Gaussian erf field; SI rescale). deps-clean holds.
  **formula-validation agent CONFIRMED all 5 points** (4π factor, −iGφ̃ gradient
  sign, ρ=−n direction, G=0 background, analytic single-mode). Both electron-only
  and electron+projectile variants = caller's choice of input n_grid.
- **Plan + handover + glossary** written. inq-study verified = inq + ONLY the two
  sanctioned CAP edits (self_consistency complexify + absorbing_monomial.hpp).
- **Pilot run.cpp written + LAUNCHED**: `ResearchProject/systems/jellium/scripts/
  cap_baselines/run.cpp` (env-driven; only CAP_MODE=b1 implemented). Builds vs
  inq-study. 100-step B1 pilot building+running in background, log at
  `scripts/cap_baselines/pilot_build_run.log`, output `results/pilot_b1/`
  (pilot.csv + run_summary.txt). PASS = energy finite/real, N(t) decreases
  smoothly (absorption), nan_seen=false.

### IN FLIGHT / NEXT (resume here)
1. **Check pilot result** (`results/pilot_b1/run_summary.txt`: nan_seen, absorbed_frac;
   `pilot.csv`: N(t) monotone-ish down, energy finite). If FAIL → diagnose CAP in
   interacting RT before anything else.
2. If PASS: implement the new observables (each code-test'd):
   - current-density field VTI writer (inqkit/io/) from `observables::current_density`
     (inq/src/observables/current.hpp:26 — native, incl. nonlocal [r,V_nl]).
   - plane/flux screens: PlaneScreen (exists) + NEW ∮J_z flux reducer (inqkit/screens/).
   - region N(t) free + per-slab reducer.
   - per-orbital energy ⟨ψᵢ|H|ψᵢ⟩(t) @300; total-system n(k)=Σfᵢ|FFT ψᵢ|² (reduced
     per frame, full 3D t=0/t=end).
3. Extend run.cpp to b2 (classical, insert projectile ion + `.ehrenfest()`, launch
   z₀=−13) and b3 (WP via inqkit WavePacket inject). B3 track = COD post.
4. `dispatch.py` (2-GPU, cudaMemGetInfo probe — NVML broken; email per baseline).
5. Long runs B1→B2→B3 (~7000 steps each, 300-frame VTI). Auto-build
   `hypotheses/cap_baselines/cap_baselines_study.ipynb`. Emails.

### Key facts / gotchas
- N(t) = `operations::integral_partial_sum(observables::density::calculate(electrons), min(2,set_size))`.
  Needs the captured NON-const `electrons` (calculate takes non-const ref).
- propagate call: `real_time::propagate(ions, electrons, cb, theory, rt, cap)` (cap LAST).
- Build recipe: `INQ_SOURCE=…/inq-study INQ_SHARE_PATH=…/inq/install/share
  PSEUDOPOD_SHARE_PATH=…/inq/install/share/pseudopod inq-run --reconfig` from the
  run.cpp dir. ~14 s/step self-reported (verify from pilot).
- NVML/nvidia-smi broken (driver mismatch) — compute fine; schedule via cudaMemGetInfo.
- GS load needs matching electrons config (spacing/extra_electrons/extra_states/
  temperature) — Common_E100_L50_cubic matches the saved GS (sv runs prove it).
