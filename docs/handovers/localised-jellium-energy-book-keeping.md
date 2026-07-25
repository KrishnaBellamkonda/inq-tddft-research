# Handover: Energy book-keeping analysis (localised jellium parameter study 2)

Campaign file: `/local/data/public/skcb2/tddft/docs/campaigns/localised_jellium_parameter_study_2/localised-jellium-parameter-study-2.md`
(id `localised-jellium-energy-book-keeping`, status `running`, interactive gate-stopped).

---

## Milestone: 2026-07-15 — OVERNIGHT twin campaign launched (classical-vs-WP diversity)

Autonomous overnight campaign to extract maximal classical-vs-WP (quantum) differences.
Orchestrator: `…/localised_jellium_dynamics/orchestrate_overnight.py` (detached, GPU0, ~9h
soft budget, idempotent-resume, checkpointed, per-pair Gmail + pairwise analysis). Log:
`orchestrate_overnight.log`; per-pair `analysis_<pair>.txt`. Design informed by a FABLE
agent (advice captured in this session) + `docs/plans/wide-wp-and-notebook-enhancements.md`.

**Matrix (m=1, fable-informed, budget-sized ~300–400 steps/run so runs reach the slab):**
- `p5_null_s2_k4` (σ2,k4.2,z−24.5): null control — Δ→0 except SIE; falsifies the framework.
- `p1_reflect_s2_k04` (σ2,k0.4,z−17): quantum reflection at the attractive surface.
- `p4_capture_s2_k11` (σ2,k1.1,z−17): capture-vs-escape + plasmon ringing.
- `p2_tunnel_s2_k05` (σ2,k0.5,z0 INSIDE slab): tunnelling leakage (riskiest — WP injection
  inside slab; 4th so failure doesn't block earlier pairs).
- `p6_ladder_s1_k11` (σ1,k1.1,z−17): σ-ladder ZPE/SIE scaling (vs σ2, σ0.5).

All pairs pass the aliasing guard (kmax<0.9·Nyquist). Both run.cpp emit interactions.csv
(pairwise P/S/B), density frames (proj_dyn frame-saving ADDED this session), + final
checkpoint. Waits for the σ=0.5 200-step run (pdyn_k1_200 done; wp_k1_200 finishing) to
free GPU0, then builds (recompiles proj_dyn frames + interaction_energies) and runs.

**Deferred (need WP fictitious mass, muon-mass-fork):** the m=10 Bragg-peak stopping (v=v_F)
and adiabatic-image pairs — the fable's headline S_WP-vs-S_cl number. Also deferred:
notebook upgrades (n(r)+Δn maps, pairwise GIF, WP−cl bar plot) + per-run notebooks
(`docs/plans/wide-wp-and-notebook-enhancements.md`).

**Campaign COMPLETE 2026-07-15 09:05** (6.8h, all 5 pairs, 0 failures incl. the inside-slab
tunnelling pair). Gauge test passes in EVERY pair (ΔE_SS=ΔE_SB=ΔE_BB=0.0000 → no gauge, all
Δ physical). σ-ladder power laws confirmed across σ=0.5/1/2: ZPE(dKin_loc) 81.6/20.4/5.1 eV
(1/σ²), self-Hartree(R) 20.8/9.8/4.4 (1/σ), SIE 4.34/1.13/−0.25 (collapses with σ). Per-pair
`analysis_<pair>.txt`. **5 per-pair analysis notebooks built** (0 errors, 8 figs each):
`hypotheses/twin_dynamics/<pair>/study.ipynb` — each covers BOTH twins with pairwise
decomposition + gauge test, n(z,t) density carpets + Δn, WP−classical bar plot, and a
pairwise-energy GIF. Builder `twin_notebook_builder.py` extended with these upgrade cells.

**Launch fix (2026-07-15 02:17):** the first `wait_for_gpu` used a bare `ps grep /run$`
that wrongly matched an UNRELATED GPU1 job (another session's `sigma1_masspair` fictitious-
mass run) → orchestrator hung "GPU0 busy" though GPU0 was free. Fixed to a device-specific
`cudaMemGetInfo` probe on physical GPU0 (throwaway subprocess, no held context; CUDA compute
works despite the NVML mismatch). Killed + relaunched; P5 classical propagating from 02:17.
Orchestrator pid 3150322 (detached, ppid=1).

**Verify in the morning:** `orchestrate_overnight.log`, per-pair emails, `analysis_*.txt`;
completed pairs' `interactions.csv`/frames/checkpoints; re-run orchestrator to resume any
unfinished pair (the fixed GPU0 probe now also guards the resume).

---

## Milestone: 2026-07-14 — Rung 2 P-dyn: moving projectile built + first dynamic pair analysed

Built the `Projectile` Ehrenfest infrastructure (wrapper-only, NO inq/ edit) and ran the
first dynamic twin pair. Design: `docs/plans/twin-run-rung2-dynamic-spec.md`.

New inqkit (all `inq-stack/include/inqkit/dynamics/`):
- `projectile.hpp` — `Projectile{mass,charge,R,V}`, velocity-Verlet (KDK, symplectic).
  PURE (detail::Vec3 only); unit-tested `tests/include/inqkit/dynamics/test_projectile.cpp`
  (4 cases / 125 assertions PASS: zero-force, constant-force exactness, F/m scaling,
  harmonic energy conservation). Registered in tests/include/CMakeLists.txt (pure tier).
- `projectile_force.hpp` — finite-difference Hellmann-Feynman force `F_z = -d/dR ∫n_proj·φ_drag`,
  `φ_drag = poisson(n_e - n_+)`; self-force zero by symmetry. NB `operations::subtract` does
  NOT exist — use linearity `F(φ_e-φ_+)=F(φ_e)-F(φ_+)` (two force calls), no field subtraction.
- `moving_gaussian_projectile_perturbation.hpp` — tracks a live `Projectile*` (survives INQ
  copy-by-value), recaches φ on move (mask-absorber pattern).

Run: `proj_dyn/run.cpp` (classical, Ehrenfest in the callback) emits `projectile.csv`
(step,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal — Hartree). WP twin = existing
`phase5_wp/run.cpp` (already supports LJ_K0). Engine `twin_decompose.load_run` now MERGES the
auxiliary `projectile.csv`/`wp_centroid.csv` by step; `check_twin --dynamic` requires
projectile.csv (incl. U_proj_bg from it).

**P-dyn pair (GPU0, k0=1.0, 50 steps, dt=0.05, σ=0.5, same GS h2/gs_p2_lz120):**
- Classical `proj_dyn/results/pdyn_k1`, WP `phase5_wp/results/wp_k1_dyn`.
- GATES PASS: energy conservation classical **0.0003 eV** (force validated), t=0 collapses
  EXACTLY to golden (20.81/81.7/-16.47/4.34), check_twin --dynamic OK.
- QUANTUM EFFECT: residual R collapses **20.81 → 1.03 eV** as the WP DISPERSES
  (σ_z 0.35→3.66 Bohr, matches analytic free-dispersion σ(t)=σ0√(1+(t/2mσ0²)²)→3.56).
  Centroids TRACK (both +2.5 Bohr) so it's spreading, not trajectory divergence. R∝1/σ so
  dispersion crushes the WP self-Hartree — the classical rigid Gaussian cannot spread.
- Classical stopping: proj KE 13.61→13.50 eV → **S~0.044 eV/Bohr** over 2.5 Bohr.
  Quantum stopping proxy E_deposited_wp = -0.047 eV (negligible; WP still approaching + short run).
- Caveat: 50 steps at k0=1 shows DISPERSION more than stopping (WP disperses before reaching slab).

Analysis notebook: `hypotheses/twin_dynamics/pdyn_k1_study.ipynb` (7 cells, 0 errors; figs
twin_decomposition/conservation/dynamics/wp_spreading.png). Builder `twin_notebook_builder.py`
extended with dynamic sections (conservation gate, trajectory/stopping/residual-collapse, WP
dispersion vs analytic). NOT git-committed. Next options: wire per-step wp_centroid.csv;
longer/slower/heavier pair for stopping; per-run notebooks.

---

## Milestone: 2026-07-14 — TWIN-RUN skills built (energy decomposition, classical vs WP)

Built two composable skills that operationalise the whole energy book-keeping into a
reusable, tested classical-vs-wavepacket comparison. Design doc:
`docs/plans/twin-run-energy-decomposition-skills.md`.

- **`.claude/skills/twin-run-generation/`** — makes a *twin pair* (same GS, identical
  config, full energy decomposition on in both; only the projectile differs) and gates
  it with `check_twin.py` (parity + observable-completeness + `U_proj_bg` availability →
  `twin_manifest.json`). Motivated by a real inconsistency: existing pairs
  (`h0_base_difference`) emit only total/kinetic/hartree/xc — no external, no proj_bg.
- **`.claude/skills/twin-run-analysis/`** — deterministic engine `twin_decompose.py`
  (parses both runs, asserts parity, computes per-step `d(·)=WP−classical`, residual
  `d(E_H+E_ext)−U_proj_bg`, `SIE=R+dXC`, known attributions → findings table) +
  `twin_notebook_builder.py` (executed analysis `.ipynb`) + SKILL.md carrying the
  interpretation rules (gauge caveat, localisation KE, WP self-Hartree, SIE). Agent
  narrates; Python does the arithmetic.

**Ladder:** Rung 1 (static, known-answer) DONE + VALIDATED; Rung 2 (dynamic) FULL SPEC
`docs/plans/twin-run-rung2-dynamic-spec.md`; engine made representation-aware + dynamic-capable
(10/10 tests). **G-static built + verified 2026-07-14** (GPU0): new ghost-UPF run
`proj_ghost/run.cpp` (p2/open-z, matched to the WP twin) → engine residual **8.85 eV**.
KEY FINDING: for the ghost/pseudopotential representation `U_proj_bg` is **ADDED**, not
subtracted (`R = d(E_H+E_ext) + U_proj_bg`) — INQ omits the z_valence=0 background-comp.
term; using −U_proj_bg gave a spurious −260 eV. Ghost residual (8.85) sits ~12 eV below the
clean perturbation (20.81) = the ghost tail aliasing (documented ~14 eV); SIE not clean for
the ghost. Codified in `reference_ghost_u_proj_bg_sign` + `twin_decompose.py` (`u_sign`).
Remaining Rung-2: `Projectile` Ehrenfest class + moving perturbation (P-dyn), ghost-ion
dynamics (G-dyn) — not built.

**Golden pair (already on disk, reused as the regression):**
`…/localised_jellium_dynamics/proj_perturbation/results/proj_pert_dx0p5` (classical) +
`…/proj_perturbation/stress_scratch/s0p5_r12_lz120_p2/results/wp` (WP) — same GS
(`campaign_autorun/runs/h2/gs_p2_lz120/checkpoint`), σ_WP=0.5, r=12, Lz120, p2.
Engine reproduces **dKin 81.74 (loc 81.63) · dXC −16.47 · residual 20.81 · SIE 4.34 eV**
exactly. Tests: 7/7 pass (`tests/test_twin_decompose.py`, synthetic + golden). Note the
WP observables lack `energy_external`/`energy_ion` in some older runs — the golden
classical run has the full ledger; `U_proj_bg=134.69 eV` comes from its
`run_summary.txt` (constant, static case).

**Verified:** engine CLI + all 7 tests pass; `check_twin.py` PASS on golden / FAIL on a
broken pair; notebook executes end-to-end (0 errors). No `inq/` edits (wrapper/skill
only). Skills self-contained per the shippable rule.

---

## Milestone: 2026-07-14 — EMPIRICAL boundary-matched self-Hartree: residual fully accounted, gauge gap CLOSED

User's insight: replace the analytic free-space self-Hartree (21.71 eV, which needed a 0.9 eV
"charged-cell gauge" subtraction) with the WP self-Hartree computed EMPIRICALLY using INQ's OWN
Poisson solver in the actual run cell. INQ picks the boundary-matched kernel from cell periodicity
(inq/src/solvers/poisson.hpp:190): periodicity 3 → fully-periodic FFT (Makov–Payne); periodicity 2 →
Rozzi et al. (2006) 2D Coulomb-cutoff = OPEN-Z. So E_self in a p2 cell IS the open-z self-Hartree —
no analytic correction. (My earlier Python "in-cell 21.49" was wrong for p2: it used a full-3D-periodic
FFT = the p3 kernel.)

New code (wrapper-only, GPU0): `proj_perturbation/self_hartree.cpp` (single-point: build WP Gaussian
density std σ_ρ, phi=poisson in the cell, E_self=0.5∫n·phi; no GS/dynamics) + `self_hartree_sweep.py`
(driver over per/dx/Lz/σ) → `hypotheses/perturbation_method/self_hartree_empirical.csv`.

RESULT — E_self reproduces the MEASURED residual to ~0.01 eV, grid-for-grid, EVERY axis:
- p2 open-z: 20.82 (dx0.5)=meas 20.81 ; 20.65 (dx0.3)=meas 20.65.
- p3 full-PBC: 21.50 (dx0.5)=meas 21.49.
- Lz {90,160,240} p2: E_self flat 20.82 == measured flat residual (open-z has no z-images).
- σ {0.7,1.0} p2: 14.45/9.79 == measured exactly. σ=0.35 lone outlier (38.1 vs 34.9): σ_ρ=0.25 only
  0.5 grid pts/σ at dx0.5, under-resolved — discretisation artifact, not physics.
- max|R − E_self| over σ≥0.5 both BC grid-matched = 0.013 eV.
CONSEQUENCE: the 0.9 eV "gauge" was just the open-z-vs-free-space reference mismatch; the residual is
now FULLY the WP self-Hartree, no fudge. Only genuinely unaccounted energy = LDA SIE = R+dXC ≈ 4.3 eV.
Notebook section "Stress test 6" added (fig `empirical_self_hartree.png`: parity plot + gap-closure bar).

BUILD NOTE: /tmp was 99% full (140M) → nvcc "Invalid argument" writing intermediates. Fix: export
TMPDIR=/local/data/public/skcb2/tddft/.buildtmp (1.4T free) for any inq-run/nvcc build. Persistent
env-level issue; set TMPDIR for all future builds until /tmp is cleared.

---

## Milestone: 2026-07-14 — PERTURBATION-METHOD stress campaign COMPLETE + notebook built

The Gaussian-charge perturbation method (prior milestone) was stress-tested across five axes on GPU0.
All runs done; notebook executes 0-error.

Notebook: `ResearchProject/systems/localised_jellium/hypotheses/perturbation_method/perturbation_method_study.ipynb`
(builder `build_perturbation_report.py`; figures sigma_sweep / r_independence / lz_gauge / grid_sweep .png).
Data CSVs: `hypotheses/perturbation_method/stress_{baseline,sigma,r,lz,p3vp2}.csv` +
`scripts/localised_jellium_dynamics/proj_perturbation/grid_sweep.csv`.

RESULTS (all confirm the residual R = d(E_H+E_ext) − U_proj_bg is the WP Hartree self-energy):
- **Baseline** (σ=0.5,r=12,Lz120,p2): R=20.81 eV; analytic free 21.71, FFT in-cell 21.49 → gauge 0.91 eV.
  dKin=+81.74 (loc 81.63), dXC=−16.47, SIE=R+dXC=**+4.34 eV** (the only unaccounted energy = LDA SIE).
- **σ-sweep** {0.35,0.5,0.7,1.0}: R tracks analytic self-Hartree, scales 1/σ; max|R−SH_incell|=0.68 eV
  (σ=0.35 worst — σ_ρ=0.25 under-resolved at dx=0.5).
- **r-independence** {4,12,20,28}: R flat to 0.01 eV (20.80→20.81) — confirms self-energy (position-free).
- **Lz open-z** {90,120,160,240}: R invariant to 0.001 eV — open-z has no periodic-image gauge growth.
- **p3 vs p2** (fresh WP decomp for p3): R p2=20.81, p3=21.49, Δ=+0.68 eV = Makov–Payne periodic-image
  gauge; p3 lands exactly on FFT in-cell periodic self-Hartree (21.49). Individual E_H/E_ext flip sign
  (charged-cell G=0) but R robust. dKin/dXC unchanged (gauge-free).
- **grid** dx {0.5,0.4,0.3}: R = 20.81/20.65/20.65 eV — grid-stable (vs pseudopotential 7.4 eV artifact).

Notebook-build bug fixed: `s.loc` returned the pandas .loc indexer (len-1), not the 'loc' column →
use `s['loc']`. No other cells affected.

Added section "The projectile potential — Gaussian-charge vs pseudopotential vs Coulomb"
(fig `proj_potential.png`): the UPF ghost's radial V(r) = erf(r/0.5)/r Ha is IDENTICAL to our
Gaussian-charge Poisson potential (max|V_upf − erf/r| = 3e-11 Ha); both cap the core at the finite
plateau 2/(√π·σ_WP) = 2.26 Ha and merge with 1/r Coulomb beyond ~0.9 Bohr; INQ's periodic Poisson
reproduces the analytic form up to the G=0 constant (~0.03 Ha). PUNCHLINE: the potential SHAPE is
correct in both representations — the 7.4 eV artifact is a treatment effect (z_valence=0 →
truncated short-range local), NOT a wrong potential.

---

## Milestone: 2026-07-14 — CLEAN residual via Gaussian-charge PERTURBATION = 20.8 eV (self-Hartree)

User's method (their idea): represent the classical projectile as a STATIONARY GAUSSIAN CHARGE
perturbation (Poisson potential added to the KS potential, like the background) instead of a UPF
ghost — no ion, no pseudopotential, no r_cut, no aliasing. New code (wrapper-only, no inq/ edit):
- `inq-stack/include/inqkit/jellium/gaussian_projectile_perturbation.hpp` — adds v_proj=+poisson(n_proj)
  (−1 projectile ⇒ repulsive), composable via `perturbations::sum(background, projectile)`.
- `scripts/localised_jellium_dynamics/proj_perturbation/run.cpp` — loads bare GS, applies
  sum(background,projectile), 2-step, tabulates E_*, computes clean U_proj_bg=−∫n_proj·φ₊ (= ideal).

RESULT (dx=0.5, r=12, GPU0, energy-conserved 55.243 Ha):
- U_proj_bg = +134.69 eV (clean ideal, r_cut-free) — vs pseudopotential impl −524.5 eV (aliased).
- (E_H+E_ext)_pert = 60.9725 Ha; d(E_H+E_ext) vs WP = +155.5 eV; **residual = d(E_H+E_ext) − U_proj_bg
  = 20.81 eV** = the WP Hartree self-energy (in-cell p2; free-space 21.71; gauge = 0.9 eV).
- Confirms the pseudopotential's 7.4 eV (rc120) was ENTIRELY the aliasing artifact: 20.81 − 13.45 = 7.36.
- Full clean WP−classical ledger: dKin +81.7 (localisation 3/(4σ_WP²)), residual +20.8 (self-Hartree),
  dXC −16.5 → net self-interaction (H+XC) = **+4.3 eV = LDA SIE** (the only genuinely unaccounted energy).

r_cut sweep of the pseudopotential residual (eval_projpot with GS, `rcut_residual_sweep.py`): the ACTUAL
observable drifts monotonically 15.6 eV (r_cut→0, projectile absent) → 7.3 eV (r_cut=120) — never
reaches the clean 20.8 (the truncated ghost never faithfully represents the projectile). impl(r_cut)
grows linearly unbounded past r_cut≈37 (`rcut_impl_sweep.py/.png`): the erf/r tail integrated against a
laterally-infinite (periodicity=2) background is conditionally convergent — real-space truncation
diverges, Poisson/Ewald (the ideal) converges. Root cause: ghost UPF declares Z_valence=0 → INQ places
the WHOLE 1/r tail as "short-range" local potential truncated at r_max=r_cut.

IN PROGRESS: grid sweep (`proj_perturbation/grid_sweep.py`, GPU0) — GS+WP+perturbation at dx=0.4,0.3 to
confirm 20.8 is grid-stable (unlike the pseudopotential 7.4). Artifacts: `rcut_impl_sweep.{csv,png}`,
`rcut_residual_sweep.{csv,png}`, `grid_sweep.{csv,png}` (pending). NOT yet folded into ledger_rcut.ipynb.

---

## Milestone: 2026-07-13 (analysis) — SOLVED: residual 7.4 eV is a pseudopotential artifact, NOT a gauge

Interactive analysis of the electrostatic residual `d(E_H+E_ext) − E_proj_bg` (the notes'
`docs/notes/energy-decomposition-skill.md`). User asked to explain why it reads 7.4 eV (r_cut=120)
instead of the expected WP self-Hartree. FULLY RESOLVED and verified:

- **Physical residual = WP Hartree self-energy** ≈ 20.8 eV (INQ p2 open-z) / 21.49 (fully-periodic
  FFT) / 21.71 (free-space analytic 1/(2·σ_ρ·√π), σ_ρ=0.354). The classical projectile is an
  external potential with NO self-energy in the ledger; the quantum WP repels itself → the entire
  r-independent residual is that self-Hartree. `consistent_ideal_residual.py`: residual=21.49 exactly
  (distortion term J[n_slab−n_bg, n_WP−n_proj]=0.000; WP density std 0.3532 == σ_pot 0.3536 → WP IS
  the rigid Gaussian at t=0).
- **The 7.4 eV = self-Hartree − pseudopotential representation error.** Exact decomposition (INQ's own
  convention, via new `eval_projpot/run.cpp`): residual = self_Hartree − ∫(n_slab−n_bg)·(v_ion−V_proj_ideal).
  Reproduced to ±0.05 eV: r_cut=50 → 20.83−6.82=14.0 (INQ 13.99); r_cut=120 → 20.83−13.45=7.4 (INQ 7.36).
- **Three proofs it is representation, not gauge:** (1) error term doubles 6.8→13.5 eV as r_cut 50→120;
  (2) grid sweep: r_cut=120 impl SWINGS SIGN with dx (−524→−269→+224→+34 for dx=0.5,0.4,0.3,0.25) — a
  numerical pathology, while ideal is grid-stable ~135 eV and r_cut=50 impl stable ~−140; (3) ideal
  (true-Gaussian) recomputation gives the clean self-Hartree.
- **Mechanism:** UPF ghost v(r)=erf(r/0.5)/r ≈ 1/r is LONG-RANGE (Z_valence=0 → whole tail placed as
  "local potential", no reciprocal long/short split). Tabulated to r_cut and gridded, the tail aliases;
  r_cut=120 (≫ Lx=Ly=50) wraps 2.4× → WORSE than r_cut=50. Bigger r_cut is worse, not better — the
  clean route is the analytic ideal term, never a larger cutoff.
- **Charged-cell gauge is <1 eV**, only in the self-Hartree channel (21.71 free → 21.49 periodic → 20.83
  p2-open-z). The user's earlier "whatever remains is the gauge" is falsified; my earlier "wrap adds
  −14 eV additively" was also wrong (wrap cancels in a consistent convention — it's the imperfect
  cancellation vs the slab spillout).
- Also confirmed: dKin=81.7 eV = 3/(4·σ_WP²) localisation; dE_xc=−16.47 eV vs directly-evaluated bare-GS
  E_xc=−230.87 (eval_gs_xc); LDA XC is local → no WP-slab XC → dE_xc r-independent & pure WP self-XC.

Artifacts (all in `hypotheses/localised_jellium_dynamics/`): `gauge_selfhartree_probe.py`,
`consistent_ideal_residual.py`, `distortion_vs_time.py`, `wrap_mechanism.py`, `projpot_wrap_findings.md`;
new eval `scripts/localised_jellium_dynamics/eval_projpot/run.cpp` (+ eval_gs_xc). GLOBAL RULE added:
no LaTeX in chat (`feedback-no-latex-in-conversation` memory). NOT yet folded into ledger_rcut.ipynb.

---

## Milestone: 2026-07-13 (08:00–09:12) — Campaign RAN autonomously to completion (5/7 tasks)

Orchestrator bug fixed (inq-run has NO --build flag; build() now calls plain inq-run which
auto-configures GPU when build/ absent). Re-launched 07:49; GPUs freed 07:44 after the
semiempirical extra runs. All phases ran, NO failures/aborts, per-phase + COMPLETE emails sent.
- Gate PASS (a smoke reciprocity, b phase12 columns emit n_proj_norm=1.0003). P1 6/6, P2 5/5,
  P3 p3_lz120 done, P5 wp+cl 500 steps (21 frames each, ~59 min/run). Outputs in
  hypotheses/localised_jellium_dynamics/: ledger.{csv,png}, rcut.{csv,png}, ledger_rcut.ipynb,
  screening_{total,induced}.gif + 42 frame PNGs.
- PHASE 2 = clean win: proj_bg_ideal FLAT 178.67 eV at ALL r_cut∈{10..50}; E_external 3994→4091,
  proj_bg_impl 0→−95.9 as r_cut grows. r_cut effect isolated in E_external exactly as designed.
- PHASE 1 ledger reproduces A1 EXACTLY (dE_WP 81.2..79.6, dKin 81.7, dXC −16.5) + adds U_proj_bg
  columns (ideal 90.7→288.6 linear in r = gauge/plate; impl −183.8→−10.5). CAVEAT: my auto
  "closure residual" = WP−CL−dKin−dXC−dHE is ≡0 BY CONSTRUCTION (INQ energy identity), so
  resid_plus_upb just echoes U_proj_bg — NOT a meaningful closure test. The real interpretation
  of how U_proj_bg completes the classical ledger (it reshapes dHE, doesn't close a leftover) is
  the USER's; gauge cancels in WP−CL only for matched single-e pairs (formula-validation caveat).
- REMAINING (2/7): Task 4 (consolidate w-field E/φ) + Task 6 (single unified campaign notebook) —
  analysis exists but spread across ledger_rcut.ipynb (P1/P2) + semiempirical_spillout.ipynb
  (P3 refreshed with p3_lz120, P4 w-fields) + the GIFs. Not merged into one notebook yet.
- set_task_done() index logic drifts once multiple flip; frontmatter corrected by hand to 5/7.

---

## Milestone: 2026-07-12 (later⁵) — Campaign FULLY BUILT + autonomous orchestrator LAUNCHED

All campaign code written, compiled, and the autonomous orchestrator is running (waiting for GPUs).
- Observables writer EXTENDED: `inqkit/io/observables_writer.hpp` + StepContext untouched — added
  ObservableSelection flags energy_proj_bg_{ideal,impl} + ObservablesWriter::set_proj_bg(ideal,impl)
  + header/append columns (default off → backward-compatible; per-run constant, not per-step).
- Runs written+compiled (CPU-verified compile; orchestrator clean-builds GPU): `phase12/run.cpp`
  (Phase 1/2 classical, LJ_LAUNCH_Z + LJ_PROJ_UPF, emits the 2 columns + n_proj_norm guard);
  `phase5_wp/run.cpp` + `phase5_cl/run.cpp` (RT screening, per-step density VTIs, ETRS).
- Orchestrator `scripts/localised_jellium_dynamics/orchestrate.py` (PYTHON, syntax-clean, LAUNCHED
  → orchestrate.log, bg task b3nhobk22): wait_for_gpus (semiempirical runs) → clean GPU builds (.gpu_built
  marker; wipes prior --cpu builds) → GPU smoke gate (a) smoke_eprojbg PASS (b) phase12 emits columns +
  n_proj_norm≈1 → ABORT+email if gate fails → P1(6 cl, 2-GPU parallel) → P2(5 cl r=20 rc10-50) →
  P3(p3_lz120, reuse semi gs binary) → P5(wp+cl RT) → notebooks. Idempotent resume (_done), per-phase
  try/except + Gmail, set_task_done flips campaign frontmatter.
- Analysis builders: `build_ledger_notebook.py` (completed ledger + U_proj_bg closure test →
  ledger.png/csv; r_cut sweep → rcut.png/csv; ledger_rcut.ipynb) and `build_screening_gifs.py`
  (total + induced/bath-only density-diff GIFs, load_vti no-fftshift, shared clim). Both syntax-clean,
  robust to partial runs. Output → hypotheses/localised_jellium_dynamics/.
- SAFETY: correctness gates only (build fail / smoke fail → abort, no garbage runs); GPU binaries
  unverifiable pre-launch (GPUs busy) but the gate is the end-to-end GPU test before any expensive run.
  User comes back to: all phases done + per-phase emails, OR a clear gate-failure/abort email.
TO RESUME/MONITOR: tail orchestrate.log; re-run orchestrate.py to resume (completed runs skipped).

---

## Milestone: 2026-07-12 (later⁴) — Campaign Task-0 code implemented + BOTH gates passed

DONE (code gate for the whole campaign):
- New inqkit header `inq-stack/include/inqkit/jellium/projectile_background_energy.hpp`:
  `projectile_background_energy(pert, electrons, ions, proj_center, sigma_pot)` → {ideal, impl,
  n_proj_norm}. ideal = ∫n_proj·v_bg (v_bg=−poisson(n₊), n_proj = normalized Gaussian σ_pot,
  ∫=1); impl = −∫n₊·v_ion (v_ion = poisson(atomic_pot.ionic_density)+atomic_pot.local_potential).
  gpu::run (NOT inq::gpu::run) — namespace fix. Wrapper-only, no inq/ edit.
- CPU smoke `scripts/localised_jellium_dynamics/smoke_eprojbg/run.cpp` PASSES: n_proj_norm=1.0003;
  ideal (+92.18 eV) == same-kernel reciprocity ∫n₊·poisson(ρ_proj) to 0.0000 eV (validates the
  poisson+integral machinery). impl=−186.78 eV (UPF gauge, r_cut-diagnostic — legitimately differs
  from ideal by gauge, NOT a bug). GATE checks: norm, reciprocity, finite. Catalogue row added.
- formula-validation agent VERDICT: CONFIRM (Jackson §1.11 cross-term; ledger-consistent with
  E_external=∫n·v_bg self_consistency.hpp:189-191). KEY CAVEAT to enforce: absolute value is
  gauge-dependent (came out +92 eV, positive, due to open-z G=0 gauge — do NOT report absolute as
  physical); gauge −c·∫n_proj cancels in WP−CL ONLY for matched single-electron pairs (N_WP−1=0),
  same box/BC/solver; guard n_proj_norm≈1. Our σ=0.5 pairs satisfy this.

REMAINING to reach full autonomy (NOT yet done — GPU-blocked by current semiempirical_spillout runs):
1. Extend shared inqkit StepContext + ObservableSelection + ObservablesWriter with 2 columns
   (energy_proj_bg_ideal, energy_proj_bg_impl; default 0, backward-compat) + injection hook.
2. Classical run.cpp for Phase 1/2 (compute+inject the 2 columns; parametrised r + UPF/r_cut).
3. Phase 5 RT run.cpp (WP + classical, ETRS, N_STEPS=500 dt=0.01, save total-density VTI every 25).
4. PYTHON orchestrator scripts/localised_jellium_dynamics/orchestrate.py: wait for current GPUs →
   GPU smoke gate (emit columns, sane) → Phase1 → Phase2 → Phase3(p3_lz120, reuse gs binary) →
   Phase5 → notebooks; per-phase Gmail; idempotent resume; per-phase try/except. Bash only for smoke.
5. Notebook builders (ledger + r_cut + screening GIFs).
Compile-testable now (CPU); actual runs + GPU smoke wait for GPUs (~2 h out).

---

## Milestone: 2026-07-12 (later³) — Campaign authored: localised-jellium-dynamics-analysis

Grill-with-docs + campaigns skill session. New campaign
`docs/campaigns/localised_jellium_dynamics_analysis/localised-jellium-dynamics-analysis.md`
(id localised-jellium-dynamics-analysis, status draft, 0/6 tasks, in INDEX). FOUR phases
(P1 E_proj_bg ledger, P2 r_cut sweep, P3 open-z vs PBC leaking charge at w=0, P4 w-sweep E/φ).
P3+P4 formalise the semiempirical_spillout Q3/Q4 work (REUSE those runs; P3 adds only p3_lz120;
P4 is analysis-only, already computed: far |E| 0.098→0.000 as w 0→1).

LOCKED decisions (all via user AskUserQuestion locks):
- E_proj_bg = new inqkit capability `projectile_background_energy`, TWO observables.csv
  columns (user chose real column, as an inqkit wrapper addition, using the cached v_bg):
  `energy_proj_bg_ideal` = ∫ n_proj·v_bg (n_proj = Gaussian σ_pot=σ_WP/√2, ∫=1) — r_cut-INVARIANT;
  `energy_proj_bg_impl` = −∫ n₊·v_ion (as-implemented pseudopotential) — r_cut-dependent.
  Both DIAGNOSTIC (never in energy_total). Recipe: v_bg from background_perturbation.hpp:65-68;
  v_ion = poisson(atomic_pot.ionic_density)+atomic_pot.local_potential (self_consistency.hpp:102);
  electrons.atomic_pot() at electrons.hpp:496.
- Phase 1: re-run 6 CLASSICAL insertions r∈{4,12,20,28,36,40}, p2/Lz120/σ0.5, GS gs_p2_lz120;
  REUSE wp_r*_p2 (WP energy unchanged). Rebuild ledger with U_proj_bg closing WP−CL (≤3 eV).
- Phase 2: fixed r=20, r_cut∈{10,20,30,40,50}; UPFs cutoff_test/upfs/…_rc{10,20,30,40}.upf + full
  (=rc50). Reuse wp_r20_p2. E_proj_bg,ideal flat vs r_cut → r_cut effect isolated in E_external →
  infer WP effective r_cut. cutoff UPFs truncate potential beyond r_cut (make_cutoff_upfs.py).
- CODE GATE: projectile_background_energy known-case test (ideal vs analytic infinite-plate;
  impl vs B1 ∫n₊·v_ghost <1 eV; reciprocity at rc50) MUST pass before Phase 1/2.

NEXT (implementation, code-test gated): (0) inqkit header + StepContext/ObservableSelection/
ObservablesWriter 2-col extension (backward-compatible, default 0) + classical run injection;
build+validate; (auto-start) PYTHON orchestrator waits for current semiempirical_spillout runs to
free both GPUs, runs the validation smoke (abort on fail), then Phase 1 → Phase 2 → notebook.
Code compiles without GPU (do now); validation+runs wait for GPU (~3-4 h out).

---

## Milestone: 2026-07-12 (later²) — Semi-empirical spill-out study (GS matrix + notebook)

Task: "Semi empirical model for the total system" — why the semi-empirical field doesn't
vanish far from the slab. Grill-with-docs interview locked: both GPUs, p2 open-z, 10-run
matrix (box/Lz, softer-w, confinement n0, solver es). New GS variant
`scripts/semiempirical_spillout/gs/run.cpp` (adds LJ_EXTRA_STATES, LJ_TEMP_EV). Orchestrators
`orchestrate.sh` (10 p2) + `orchestrate_extra.sh` (finish 4 pending p2 + 3 p3/PBC).

KEY FINDINGS (data-backed, notebook `hypotheses/campaign_autorun_study/semiempirical_spillout.ipynb`,
built by `build_semiempirical_spillout.py`, resilient to partial matrix):
- The far-field plateau A3 saw is NOT a uniform floor: interior vacuum tail DECAYS and is
  Lz-independent (n_e@20≈1.3e-6 all boxes) = physical; the plateau is driven by a NEAR-EDGE
  DENSITY PILE-UP at the open-z box boundary (~5e-5 e/Bohr³ for w=0 at Lz=160).
- w REMOVES it: at Lz=160, near-edge density 5.1e-5 (w=0) -> 5.4e-15 (w=1,w=2) ~ machine zero.
- Q4 (semi-empirical E/φ for w-sweep, user-requested verification): far-field |E| = 0.098
  eV/Bohr (w=0) -> 0.000 (w=1,2,4); peak vacuum |E| 0.30->0.15->0.076->0.046 as w=0->1->2->4.
  A w=1 Bohr smoothing kills the spurious far field entirely.
- First-metric bug caught+fixed: original "floor" window scaled with Lz (conflated distance);
  now interior tail at FIXED distances + near-edge pile-up measured separately.
STILL RUNNING (orchestrate_extra.sh, 2026-07-12 ~20:19): N328, lz240, N164, es60 (finish p2),
  p3_lz90/p3_lz160/p3_lz240 (Q3: does the pile-up vanish under PBC?). Re-execute the notebook
  when done to fill Sweep A/C/D + the Q3 open-z-vs-PBC overlay. N164/N328 empirically test the
  user's untested "denser electrons -> less spill" expectation.

---

## Milestone: 2026-07-12 (later) — KE-bookkeeping experiment: INQ does NOT save projectile KE

### Question resolved (source + run)
Does INQ record the projectile's kinetic energy as an ionic energy, and do the WP-drift
and classical-ion KEs cancel?

- **INQ does NOT populate `energy_ion_kinetic`.** `ions::kinetic_energy()=Σ½mv²` exists
  (`systems/ions.hpp:267`) and the energy slot + CSV column exist, but the ONLY
  `ion_kinetic(value)` setters in the whole tree (inq AND inq-study) are unit tests
  (`energy.hpp:300`, `results.hpp:108`). The propagator sets `energy.ion(...)` but never
  `ion_kinetic`. So the column is always 0.
- **Experiment** (`scripts/ke_check/`, new): 4 runs, p2 GS `gs_p2_lz120/checkpoint`, 2 steps,
  GPU (WP binary is CUDA-only; step-0 accounting only — no dynamics needed). New binary
  `ke_check/classical/run.cpp` = campaign classical + `LJ_VZ` velocity + `LJ_EHRENFEST` knobs.
  - WP drift: kinetic(k0=1)−kinetic(k0=0) = 13.59 eV ≈ ½k0²=13.61 eV, flows into E_total
    (Δtotal=13.59). WP KE IS inside E_total (electronic kinetic).
  - Classical cl_v1 (MOVING, ehrenfest=1, v=1, species.mass()=1 m_e): `ion_kinetic`=0.000,
    Δtotal(cl_v1−cl_v0)=0.000. The ½Mv²=13.61 eV is entirely ABSENT.
- **Conclusion:** ½k0²(WP)=13.59 ≈ ½Mv²(CL)=13.61 (mass+velocity matched) are equal by
  construction, but INQ keeps only the WP's. So `WP−CL` at finite velocity carries a spurious
  un-cancelled +½k0²; the classical total needs a HAND-ADDED ½Mv² (second dropped term, exactly
  analogous to E_proj_bg — both because the ghost is a potential/ion, not charge/wavefunction).
  At v=0 (A1/A2 runs) this vanishes, so it never bit before; any MOVING-projectile stopping
  comparison must add it. Sanity anchor: dKin(v=0)=156.72−74.97=81.75 eV = localisation 3/(4σ²).
- Files: `scripts/ke_check/classical/run.cpp` (new, GPU-built), `scripts/ke_check/runs/` (outputs).
- **Validation notebook** `hypotheses/campaign_autorun_study/ke_bookkeeping.ipynb` (built by
  `build_ke_bookkeeping.py`, executed 0 errors / 3 figures): component table + 3 validations
  (WP KE into E_total 13.59 eV; ion_kinetic=0 even moving; matched to 0.016 eV) + v=0 anchor
  dKin=81.75 eV = localisation. This is the data-backed table/plot form of the finding above.

---

## Milestone: 2026-07-12 — USER interpretation session (A1 clarified) + v_bg sanity-check section added

### What was done (interactive, with the user going through A1 results)
- Clarified A1 ledger column semantics from source (`build_energy_book_keeping_notebook.py:88-97`):
  `dE_WP`/`dE_CL` = E_total(run,t=0) − E_GS (two separate insertion energies vs the
  bare GS); `WP−CL` = their difference; `dKin`/`dXC`/`d(H+E)` = per-component **WP−CL**
  differences. The `d` prefix means TWO reference points (vs GS for dE_*, vs CL for the
  components).
- Confirmed from run code: A1 WP has `k0=0` (`campaign_autorun/wp/run.cpp:56`) and the
  classical ghost is inserted at rest (`classical/run.cpp:58`), so `dKin=81.7 eV` is PURE
  localisation zero-point `3/(4σ²)` with no drift KE; `energy_ion_kinetic` IS recorded
  (both runs, `sel.energy_ion_kinetic=true`) but =0 here.
- Confirmed INQ does NOT register E_projectile_bg: the projectile is an INQ ion, the
  jellium background is a `localised_background_perturbation` (not an ion), so no
  Ewald/`ion` term couples them — structurally absent (`classical/run.cpp:8,85` say so).
- Established (from `inq/src/hamiltonian/self_consistency.hpp:164-197`): the KS scalar
  potential is REBUILT each step from `vion_` (line 176), the perturbation re-adds `v_bg`
  every step (line 189), but `φ` is Poisson-solved ONCE and cached
  (`background_perturbation.hpp:65-68`, `mutable optional phi_`). `v_bg` is time-INDEPENDENT.
- Verified from source that INQ `periodicity(2)` changes ONLY the Poisson kernel
  (`solvers/poisson.hpp:188-208`, Rozzi slab truncation); the kinetic/Laplacian operator
  has no periodicity branch (`operations/laplacian.hpp`, FFT grid), so **p2 does NOT stop
  wavefunction wrapping** — that needs a CAP. Recommendation to user: p2 (correct isolated-
  slab electrostatics; the ~6 eV p3−p2 offset is a charged-cell G=0 convention artefact).

### New deliverable — v_bg vs infinite-plate sanity check (theoretical_slab_model.ipynb)
- New INQ dumper `scripts/campaign_autorun/dump_vbg/run.cpp` (built+run, inq-study, GPU):
  builds the p2 slab background, `φ=poisson(n₊)`, writes `v_bg=−φ` z-lineout to
  `dump_vbg/results/vbg/vbg_lineout.csv`.
- Added section "v_bg from the Poisson solver vs the infinite plate" to
  `hypotheses/campaign_autorun_study/theoretical_slab_model.ipynb` (cells 11-12, via
  `build_theoretical_model.py` const `VBG`). Notebook rebuilt + executed in-place, 0 errors.
- RESULT (verified): `v_bg` is parabolic inside, linear outside, outside slope
  numeric 5.50 vs analytic 4πn₀a = 5.38 eV/Bohr (2% = discretised half-width 12.0 vs 12.5);
  far field is a constant ≠0 (charged-plate hallmark). Behaves like an infinite plate.
  Caveat surfaced: zero-far-field-by-Gauss applies to the NET density, not the background alone.

---

## Milestone: 2026-07-11 (06:30) — AUTONOMOUS PHASE COMPLETE: all 9 tasks done; notebook executed (0 errors, 4 figures)

### Final state
- All B2 runs complete: ghost SCF E = 67.086014 (r=4), 65.508394 (r=12),
  62.404374 (r=28) Ha; screening gains −3.0/−2.0/−0.4 eV (monotone → 0, far
  control PASS); n83_r12 (83-e unconstrained SCF) = 60.028121 Ha — the extra
  electron BINDS 9.7 eV below the bare GS, ~90 eV below the σ=0.5 packet
  (the "WP SCF is ill-posed" illustration).
- Notebook `hypotheses/campaign_autorun_study/energy_book_keeping_campaign.ipynb`
  executed in-place: 11/11 code cells, 0 errors, 4 figures (builder fix: no
  Agg override — inline backend needed for figure capture).
- Frontmatter 9/9 tasks done; status stays `running` pending the USER's
  interpretation (their explicit intent: "I am going to interpret the results").
- Campaign INDEX regenerated.

### For the user's interpretation session
Read the notebook top-to-bottom; the takeaway table separates what closed
(A3 mechanism, A4 derivation, A5 scaling, B1 decomposition ±4 eV, B2 trend)
from what is theirs to judge (naive hypothesis NOT supported as worded; the
UPF-truncation asymmetry; the ±4 eV bounded residual; B3 CAP-era differences).
`docs/notes/localised-jellium-parameter-study-2.md` is reserved for their notes.

---

## Milestone: 2026-07-11 (superseded) — AUTONOMOUS PHASE (user-ordered): A3–A6 + B1 + B3 done; B2 runs in flight

### Shape change (user, 2026-07-11)
All remaining tasks run autonomously; decisions a human would make are ruled by a
Fable 5 ADVISOR agent (rulings logged verbatim below in compact form); deliverable
= one executed notebook; the USER interprets all results afterwards.

### Advisor rulings (Fable 5 agents a48d253211f6e91ed, a9b9ca1f36038153c)
1. B1 = hybrid: post-hoc E_proj_bg, NO re-runs (stationary 2-step insertions are
   deterministic — explicitly overrode the pre-autonomy "re-run all over again");
   C++ per-step tracker filed as follow-up for moving-projectile runs.
2. Dual-route validation sufficient (independent implementations, σ_pot entered
   separately, limiting cases, convention stated, docs/validation record).
3. B2 = classical-ghost SCF r∈{4,12,28} + ONE 83-e SCF at r=12 (labelled
   illustration; frozen-WP ≡ ghost caveat up to WP self-XC −16.5 eV).
4. B1 closure failure → exact decomposition, gated: parse UPF by data FIRST;
   known-case gate (reproduce measured d(H+E)); then ≤2 new save radii (approved).
5. B3 includes post-hoc E_pb(t) alongside raw diff; finest window [0, 6.4] au.

### Results so far (all evidence neutral; user interprets)
- A3: far-field plateau = 2πQ_enc/A EXACTLY; Q_enc = +0.39 e = vacuum-floor density
  (8.4e-6 e/Bohr³) pooled at box edges (user's spill suspicion = the mechanism;
  floor labelled Inference: SCF numerical floor). w and dz (native) EXCLUDED.
  Artefact: hypotheses/campaign_autorun_study/a3_far_field_forensics.md
- A4: 3/(4σ²) derived (σ = ψ-width ⇒ density std σ/√2, consistent with σ-matching);
  production-grid numeric ⟨T⟩ = 3.0038 Ha = 81.74 eV == measured 81.7 (incl. grid +0.1%).
- A5: V(r) = −erf(r/σ_WP)/r; r_cut(1%) = 1.82 σ_WP; σ-independent long range.
- B1: E_pb dual-route PASSED (0.20 eV; docs/validation/e-proj-bg-dual-route.md).
  Naive closure FAILS (+2.5 eV/Bohr residual slope). UPF parsed by data:
  V = +erf(r/0.5)/r, pure +1/r tail ENDING at mesh r_max = 50 Bohr. Exact 4-term
  t=0 identity (b1_decomposition.py, per-G p2 Poisson) reproduces measured d(H+E)
  at r={4,12,28,40} to ±4 eV; ablation DECISIVE: only ghost-truncated-at-50-with-
  lateral-images matches (alternatives off by 96–510 eV). Two new screening saves
  r={28,40} (advisor-approved).
- A6: three long-range mechanisms synthesised (A3 floor artefact; UPF truncation
  CONFIRMS user's cutoff suspicion; σ core-only) — in notebook, user resolves.
- B3: per-step ledger diff computed + E_pb(t) (Lz=90 re-validated ≤0.23 eV);
  classical twin decelerates below 0.85·v0 by t≈2 au (expected physics).
- B2: gs_ghost/run.cpp written+built (only new code; full energy decomposition in
  summary). Runs: ghost_r28_p2 (GPU 1) + ghost_r12_p2 (GPU 0) IN FLIGHT;
  ghost_r4_p2 + n83_r12 queued next. NOTE: r=28 run outputs land in
  scripts/campaign_autorun/gs_ghost/results/ (launched via inq-run from the
  binary dir — MOVE into runs/ghost_r28_p2/ before analysis); r=12 correctly
  writes under runs/ghost_r12_p2/.
- Notebook: hypotheses/campaign_autorun_study/energy_book_keeping_campaign.ipynb
  built (builder: build_energy_book_keeping_notebook.py); EXECUTE after B2 lands
  (venv python3 builder --execute). Catalogue rows appended (test-catalogue.md).

### Exact next steps
1. `gs_ghost/finish_b2.py` is RUNNING in background: waits for ghost_r12_p2 +
   ghost_r4_p2 (slow SCF, ~67 s/iter, de ~1e-1 at iter 48 — possibly hours),
   then launches n83_r12 (GPU 1), then rebuilds+executes the campaign notebook
   (nbclient; builder fixed for the pyenv-shim nbconvert issue), then writes
   `gs_ghost/B2_ALL_DONE.marker`. Idempotent — safe to re-run.
2. ghost_r28_p2 DONE: E = 62.404374 Ha (dE = 55.0 eV vs insertion 54.6 → far-field
   control passes, screening ≈ 0.4 eV at r=28). Outputs already moved into
   gs_ghost/runs/ghost_r28_p2/results/.
3. After the marker appears: verify notebook has 0 error outputs, flip B2 task
   done in campaign frontmatter, regenerate INDEX (campaigns build_index.py),
   final summary to the user. Status stays `running` until the user reads +
   interprets (their explicit intent).
4. Frontmatter already flipped for A3–A6, B1, B3; INDEX regenerated.

---

## Update: 2026-07-11 — A2 computed, STOPPED AT A2 GATE

Status: A2 evidence complete, awaiting user verdict. Pair =
`scripts/qsp_phase3/{wp,classical}/results/p3_{wp,classical}` (matched twins,
100 eV, σ=0.5, 3D-periodic Lz=90 box, E_GS = −70.22568216820937 Ha; "p3" = qsp
phase 3, not periodicity). Results: (a) classical 100 eV NOT in E_total —
E_tot(elec)+KE_ion conserved to 2.9 eV while terms swing >100 eV; (b) WP 100 eV
IS in E_kinetic; (c) dKin(0) = 180.8 eV vs predicted drift+ZP = 181.6 eV
(−0.8 eV, 0.5%). Incidental: run_summary `ke_ion_initial_ha` field is
mislabelled (holds FINAL ½vz²); classical m_e projectile never entered the slab
(bounces outside the surface). Artefact:
`hypotheses/campaign_autorun_study/a2_launched_pair_100ev_audit.md`.
Next: user verdict on A2 → flip task, then A3 forensics.

---

## Update: 2026-07-11 — A1 gate PASSED: verdict "use p2 for now"

Status: A1 done; periodicity 2 locked for all downstream tasks (pragmatic lock;
the 4–7.5 eV p3−p2 WP-channel offset stands unadjudicated). A2 (launched-pair
100 eV audit) started. Candidate pair identified from the (independently
authored) debugging-quantum-stopping-power campaign file: `qsp_phase3` runs
`p3_wp` + `p3_classical` (100 eV, σ = 0.5, box 50×50×90, N = 82, CAP,
E_GS = −70.22568216820937 Ha) — match to be verified from run_summary before use.

---

## Milestone: 2026-07-11 — campaign designed, locked, and started; A1 table built, STOPPED AT A1 GATE

### Current status
Campaign designed via the campaigns-skill grill (Decisions 1–4 locked by the user
2026-07-10/11) and flipped to `running`. Task A1 (periodicity 2 vs 3 ledger
comparison) is COMPUTED and its artefact written; the session is **stopped at the
A1 gate awaiting the user's verdict**. Nothing else has started. No new
simulations were run — A1 used existing insertion-run CSVs only.

### Locked decisions (user)
- Interactive gate-stopped campaign; agent stops wherever the user's judgement is
  needed. Phase A (existing data / understanding) before Phase B (new runs).
- Task ladder A1–A6, B1–B3 with done-criteria (see campaign file). Derive/
  understand tasks (A4, A5, A6) at the END of Phase A (user reorder 2026-07-10).
- Hypothesis: WP−CL difference fully accounted for by WP self-energy
  (3/(4σ²) zero-point + self-XC + self-Hartree) + missing classical E_proj_bg;
  closure tolerance ≤ 3 eV per ledger row; launched-pair prediction
  dKin_WP−CL ≈ 100 + 82 ≈ 180 eV. T3 (far-field anomaly) is supporting, not
  pass/fail.
- Default periodicity 2 if runs must proceed before the A1 verdict.
- Artefacts → `ResearchProject/systems/localised_jellium/hypotheses/campaign_autorun_study/`;
  `docs/notes/localised-jellium-parameter-study-2.md` is the user's own thinking
  file (agent writes only A6's resolution there, in the user's words).

### What changed
- `/local/data/public/skcb2/tddft/docs/campaigns/localised_jellium_parameter_study_2/localised-jellium-parameter-study-2.md`
  — full campaign body written (question, locked decisions, verified anchors,
  task ladder with gates, rules); frontmatter: status running, 9 tasks, this
  handover pointer; user's original notes preserved verbatim at the bottom.
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/campaign_autorun_study/a1_periodicity_ledger_comparison.md`
  — NEW: the A1 side-by-side ledger (p2, p3, p3−p2) with provenance and neutral
  observations. Verdict line deliberately open.

### A1 result (awaiting user verdict — do NOT judge for them)
- dKin (81.7 eV) and dXC (−16.5 eV) identical across periodicities to < 0.05 eV.
- dE_CL agrees to ≤ 2.1 eV. dE_WP differs by +6.0…+6.4 eV (p3 higher),
  all in the d(H+E) channel, drifting 4.2 → 7.5 eV over r = 4 → 40.
- Inference (labelled, unverified): net −1 charged cell under WP insertion ⇒
  p2 vs p3 Poisson G=0 conventions differ; classical projectile lives in
  external potential ⇒ dE_CL nearly convention-free.

### Verified facts this design rests on (checked in-session, 2026-07-10/11)
- Insertion ledgers are t=0, projectile AT REST (k0 = 0, `run_summary.txt`):
  the +100 eV test belongs to LAUNCHED runs (task A2).
- p2 ledger: notebook cell 22, runs `scripts/campaign_autorun/runs/h0_p2/`,
  E_GS = 60.38307052445239 Ha (`runs/h2/gs_p2_lz120`).
- p3 ledger: cells 31–39, runs `runs/h0_p3/`, E_GS = −108.5336851082701 Ha
  (`scripts/h0_base_difference/gs` == `runs/h2/gs_lz120`; checkpoint
  `shared_gs/slab_n82_L50x50x120`; anchors `orchestrate.py:27-28`).
- Both GS slabs identical apart from periodicity: 50×50×120 Bohr, half-width
  12.5, N = 82, spacing 0.5, LDA, edge_width 0 (sharp).
- Raw dHartree/dexternal are charged-cell-convention-poisoned
  (−274 eV p2 vs −29 eV p3 at r = 40, cell 39) — never compare them.
- "w parameter" = `EDGE_WIDTH_BOHR` (erfc edge softening; 1.0 in production
  configs, 0 in the h0 GSs). Existing w-sweep for A3(v):
  `runs/h1/gs_w{0,0.25,0.5,0.75,1,1.5,2,3}` (GS energies span −69.87 to
  −71.83 Ha).
- Semi-empirical chain (A3 target): notebook cells 11–19 →
  `docs/reports/09-07-2026-meetng-emilio/assets/make_s1_3_field_potential.py` →
  `hypotheses/plate_model/build_plate_model.py` + `VALIDATION.md`. The s1_3
  builder uses the Lz=160 p3 GS (`runs/extend_r160/gs_lz160_p3`), sharp top-hat
  background, gauge fixed at ±(edge−5 Bohr), Gaussian-convolved with
  σ_pot = σ_WP/√2.

### Commands run
```bash
# A1 computation (venv python inline): read observables.csv row 1 for
# runs/h0_p2/{wp,cl}_r{4,12,20,28,36,40}_p2 and runs/h0_p3/..._p3,
# columns dE_WP, dE_CL, WP-CL, dKin, dXC, d(H+E) vs the two GS energies above.
```

### Tests and validation
- Proposed: none needed for A1 (pure re-tabulation of existing CSVs; the
  notebook's own sum(parts)==total exactness checks pass at ~1e-13 Ha).
- Run: cross-checked cell-22 printed table reproduces exactly from raw CSVs.
- Remaining gaps: A2 launched-pair identification not yet done; B1 formula
  validation pending (pre-gated).

### Known issues / blockers
- BLOCKED (by design) on the user's A1 verdict.

### Assumptions still in play
- The ~6 eV p3−p2 WP offset is convention, not physics — labelled Inference,
  untested; the user may direct A-phase work at it.
- A2 assumes a matched launched 100 eV σ=0.5 classical/WP pair exists with
  per-step energy CSVs (candidates in run sets h0..h5 / quantum-stopping runs;
  not yet located).

### Exact next steps
1. Present the A1 gate to the user (three tables + neutral observations);
   record their verdict in the campaign file (flip task A1 done) and in
   `a1_periodicity_ledger_comparison.md`; lock periodicity downstream.
2. Start A2: locate the matched launched pair (100 eV, σ = 0.5, classical + WP,
   per-step energies); test the three sub-claims; stop at the A2 gate.
3. Then A3 forensics (5 sub-checks, reuse the h1 w-sweep for (v)).
4. Regenerate `docs/campaigns/INDEX.md` after any frontmatter change
   (`venv/bin/python3 .claude/skills/campaigns/build_index.py docs/campaigns`).
```
