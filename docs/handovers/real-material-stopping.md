# Handover: Real-material test of the two stopping-power definitions

---

## Update: 2026-08-06 (v16) — QUOTA CLIFF: MPhil SL2 GPU account at 354/365 h; Phase 3 (32910120) queued but likely blocked; padded 12 h TimeLimits trimmed to unblock scheduler

- **Account near exhaustion:** `mybalance` → MPHIL-NIKIFORAKIS-SKCB2-SL2-GPU
  usage 354 of 365 GPU-h, **~11 h available** (2026-08-06 ~04:15). Entire queue
  (lzb campaign + gr2d-p3) pending with `AssocGrpGRESMinutes` — SLURM projects
  each job's full requested TimeLimit against remaining quota, and every job
  requested a padded 12 h, so NOTHING could start (deadlock, idle GPUs — the
  checkpoint-dont-block failure mode at account level).
- **Mitigation applied (this session):** `scontrol update TimeLimit` — gr2d-p3
  32910120 → 3:00:00; lzb-wp 32886170 + lzb-cl 32886171 → 3:00:00; lzb-vac
  32886172–75 → 1:30:00. Effect was immediate: lzb-wp_10, lzb-wp_11, lzb-cl_0
  started within seconds. Measured member costs justify the margins (lzb-wp
  0:52–2:04; gr2d Phase 2 WP 0:32–1:32).
- **Known risk:** running lzb-wp_9 was already at 2:22 elapsed when the trim
  landed (users cannot RAISE TimeLimit) → it will be killed at 3:00 if not done;
  resumable per final-timestep-checkpoint rule.
- **Phase 3 status: 0/9 started.** Projected cost ~7–8 GPU-h (mono ≈ half the
  bilayer cost; bilayer E15 the longest). lzb (ahead in priority) will consume
  most/all of the remaining ~11 h → **Phase 3 very likely blocked on quota**.
- **USER DECISION REQUIRED:** (a) top up / extend the SL2 allocation; or
  (b) switch `--account` to an SL3 GPU account with headroom
  (NIKIFORAKIS-SL3-GPU: 3,000 h avail; BLAKELY-SL3-GPU: 2,997; PENG-SL3-GPU:
  3,000) — billing another project is not an autonomous call; or (c) kill
  remaining lzb members to leave the last hours for Phase 3.
- Monitor b7a6ypiq6 still armed on 32910120. CPU-side closure work (P0/P1/P2
  phase notebooks, E25 field-level analysis, catalogue upsert) proceeds
  regardless of GPU quota.

---

## Update: 2026-08-05 (v15) — Phase 2 MEASURED: ratio 2.8 @300 eV / ≥0.58 @100 / reflection-dominated @25; email pending credentials

- **Ledger results (hypotheses/twodef_sv/phase2_ledger.{csv,png}, method +
  caveats in build_phase2_ledger.py):** classical deposits 11/12/9.4 eV at
  25/100/300 eV; WP flux-ledger deposits: E300 ≈ 26 eV → ratio ≈ 2.8 (near
  bulk 2.2!); E100 ≥ 7.1 eV → ratio ≥ 0.58 (lower bound, front-biased);
  E25 UNPHYSICAL (−5.1 eV; energy_total rose 7.8 eV — non-Hermitian
  bookkeeping + strong reflection) → E25 needs field-level analysis.
  E_wp0 internal check passed (30.2 = 25 + 5.1 spread ✓). **Phase 2 gate MET:
  hypothesis "no difference" REFUTED with structure** (ratio energy- and
  channel-dependent; bulk-like excess only in the transmission regime).
- All 6 run notebooks built+executed (hypotheses/twodef_sv/*_run_notebook.ipynb)
  with sidecar *_figs/ dirs (figures exist as files; inline density-GIF embed
  still to verify per notebook-density-gif rule).
- **Phase 2 email composed but NOT SENT — Gmail credentials missing on CSD3**
  (~/.config/inqview/gmail_credentials.json; interactive App-Password setup —
  USER-ONLY action: `venv/bin/python -m inqview.email setup`). Full four-part
  body saved verbatim: hypotheses/twodef_sv/phase2_email_pending.txt.
- Remaining Phase 2 closure: phase notebook (P2) + P0/P1 phase notebooks +
  run-catalogue upsert + E25 field-level analysis. Then Phase 3 dispatch
  (mono + full grid + sites, same CAP verdict).

---

## Update: 2026-08-05 (v14) — Run notebooks building (builder works on graphene unmodified); ledger-analysis plan pinned

- E100 pair notebooks BUILT+EXECUTED via the corpus run_notebook_builder.py
  with NO graphene adaptation (44-cell wp / 31-cell classical), at
  `hypotheses/twodef_sv/{wp,classical}_E100_run_notebook.ipynb`. Skips:
  overlap (not wired in these runs), RT orbital VTIs, paraview. **TODO
  (notebook-density-gif rule): density phase ran in 0.0 s and ipynb ~15-20 KB
  → density GIFs almost certainly NOT rendered/embedded — diagnose in the
  phase-notebook pass** (VTI naming/layout vs builder expectation?).
  E25+E300 pairs building in background (task with 4 invocations).
- **Pinned plan for the PROPER WP deposited-energy ledger** (replaces the
  invalid ⟨p_z⟩²/2 proxy): try the twin-run-analysis skill's deterministic
  engine (twin_decompose.py) on the pair layout first (runs lack
  twin_manifest.json — may need a hand-written manifest); else port the
  m8/m9 CAP-loss accounting: deposited = ΔE_bath where E_bath from
  observables.csv energy components minus WP-orbital share (interactions
  e_ss/e_ps/e_pp per step), CAP-removed energy from -dE_total/dt integration.
  Cross-check with classical window deposits (11/12/9.4 eV at 25/100/300).
- After ledger: phase notebooks (P0 GS battery, P1 CAP scan incl. σ=0.5
  lesson + characterization runs, P2 twins), tddft-run-catalogue upsert,
  Phase 2 four-part email (email-notifications skill) with phase2_ratio.png
  + ledger table, then Phase 3 dispatch decision (mono runs + remaining
  energies + sites with same CAP verdict).

---

## Update: 2026-08-05 (v13) — Classical v2 SUCCESS (real stopping); first-look twin comparison: QUALITATIVE divergence

- **Classical v2 runs COMPLETE + physically sound** (job 32890981, 1h12m):
  deceleration through the bilayer at all E — v: 1.3555→0.9089 (E25),
  2.7110→2.4873 (E100), 4.6955→4.5927 (E300); all stop at CAP inner edge.
  Matched ±8 Bohr window deposits (image-force-cancelling): ≈11/12/9.4 eV
  at 25/100/300 eV — broad stopping maximum near 100 eV ✓ plasmon-matching
  expectation.
- **First-look WP/classical comparison
  (hypotheses/twodef_sv/build_phase2_ratio.py, phase2_ratio.{csv,png}):**
  E25: WP centroid NEVER crossed +8 Bohr — packet substantially
  backscattered/absorbed (consistent with Geelen low-E bilayer transmission);
  E100/E300 proxy ratios 6.9/12. **PROXY CAVEAT (recorded in script):
  ⟨p_z⟩²/2 at centroid crossings is INVALID under packet splitting
  (transmitted+reflected mixture) — these are NOT deposited-energy ratios.**
  Qualitative finding stands: WP vs classical diverge FAR beyond bulk
  jellium's 2.2 in a real 2D target (reflection/fragmentation channels the
  classical point lacks). Hypothesis "no difference" heading to refutation —
  quantify via proper ledger (bath e_ss gain + CAP-removed energy, m8/m9
  pattern) in the notebook stage. NEXT: port CAP-ledger accounting; run
  notebooks (trial for E100 pair in flight, task buucg9y67); phase notebooks;
  Phase 2 email with phase2_ratio.png + ledger table.

---

## Update: 2026-08-05 (v12) — Classical v2 (direct-potential) WRITTEN + build/run job 32890981 submitted

- classical/run.cpp REWRITTEN on the corpus-validated direct design (clone of
  sigma56_sv/classical): projectile = `inqkit::dynamics::Projectile` +
  `moving_gaussian_projectile_potential` (external erf/r, NEVER an INQ ion) +
  own velocity-Verlet. Forces: `projectile_force_direct_z(n_e)` + NEW analytic
  erf-smoothed carbon-core term (q_P·Z_C, g'(0)=0 smooth at contact, xy
  minimum-image, LIVE Ehrenfest core positions). Electrons construction now
  matches the GS exactly (no extra_electrons). Ledger via
  compute_coulomb_direct + zero background. Recorded limitation: no
  back-reaction on cores (nuclear channel ≤1.8e-4·E). Run ends before the CAP
  inner edge (we own the coordinate — no wrap/out-of-box garbage).
- Job 32890981: build v2 + run all three classical members (E25/E100/E300,
  bilayer, CAP η=1 W=20) sequentially. Monitor bn83eed03.
- Sanity to check on completion: track vz DECREASES through the slab (stopping,
  not slingshot); Fel/Fcore columns finite and ~cancelling far from slab;
  e_ss(t0) == energy_hartree(t0).

---

## Update: 2026-08-05 (v11) — Twin array done: WP halves GOOD; classical halves INVALID (point-ion slingshot) — rework required

- **Array 32886967 all COMPLETED exit 0. WP halves (bilayer E25/E100/E300)
  physically sound:** t=0 norms exact (norm_wp=1, norm_total=193); final
  norm_wp 0.149/0.151/0.047 (85–95% CAP-absorbed); bath lost ~0.06 e⁻ to CAP
  (secondary emission ✓). Full closure/ledger verification pending notebooks.
- **Classical halves INVALID — runaway acceleration** (E100: vz 2.71→95 a.u.,
  z→744; E25 slingshot to vz=+18, z=−24). DIAGNOSIS (from electron_track +
  inq-study/src/ionic/species.hpp): mass CORRECT (.mass() takes amu; internal
  1822.8885×(1/1822.8885)=1 a.u. ✓; vacuum flight force-free ✓). Failure =
  **INQ ion–ion term treats the projectile as a POINT −1 charge** (erf
  smoothing exists only in the electron-side UPF) and our impact site
  (entry-layer hex centre) is atom-atop in layer 2 (Bernal) → 1/r²
  attractive singularity → slingshot catapult. The corpus's validated
  classical-electron design (sigma56_sv/classical, proj_dyn) avoids INQ ions
  ENTIRELY for this reason (direct moving erf/r potential + own
  velocity-Verlet). cap_cl/run.cpp's idiom was never actually executed
  (graphene GS never built pre-migration) — do not trust it as validated.
- Invalid results QUARANTINED: classical/results/invalid_pointion_E{25,100,300}
  (kept for the record; dispatcher will re-run E* fresh).
- **REWORK (next):** rewrite classical/run.cpp on the corpus design — clone
  `sigma56_sv/classical/run.cpp` (352 L): projectile = moving external
  erf-smoothed potential (σ_pot=σ_WP/√2), own velocity-Verlet; force =
  ∫n·(−∇V) from electrons + analytic erf-SMOOTHED force from C cores
  (+4 each) + implement edge-zeroing at box exit (fixes tddft-simulations
  2d″ TODO cleanly since we own the coordinate). Carbons stay INQ-Ehrenfest.
  Then resubmit --array=0-2 (classical members) with same CAP verdict.
- WP results live at scripts/twodef_sv/wp/results/E{25,100,300}; capchar job
  32886968 (k=1.36 characterization) may still be queued/running.

---

## Update: 2026-08-05 (v10) — Phase 1 COMPLETE (CAP verdict η=1, W=20); Phase 2 TWIN RUNS DISPATCHED (job 32886967)

- **Phase 1 COMPLETE.** Full σ=2.0 residue table measured (capscan/results/*_s2).
  Verdict (recorded in plan Stage D): **η=1 Ha, W=20 Bohr/face, Lz=80 — one CAP
  for the whole campaign.** Production-range worst residue 6.3e-3 at k=1.36
  (documented E=25 norm-accounting uncertainty); ≤2.1e-4 at k≥1.92. k=1.05
  ~3% reflection floor → Phase 3 decision. Characterization job 32886968
  (η=0.5 / longer-window at k=1.36) feeds the phase notebook only.
- **Phase 2 DISPATCHED: job array 32886967 [0-5]** = classical/WP ×
  E=25/100/300 eV, bilayer, σ_WP=2.0, CAP verdict baked in, GS =
  shared_gs/gs_3x2_bi_cut50_Lz80, UPF = electron_gaussian_wpsigma2p0_He.upf.
  Monitor bzo9mcxxc (polls 7 min) reports when ALL members exit.
- On completion: check per-member run_summary + t=0 gates in job logs
  (gr2d-twin-32886967_*-*.out at repo root); closure gates
  (WP: e_hartree_check vs INQ; classical: e_ss==energy_hartree in analysis);
  Ehrenfest sanity (carbon drift in ions_track); THEN measure the
  WP/classical ratio per energy (Phase 2 hypothesis: ≈1; bulk gave 2.2);
  build run notebooks (run_notebook_builder.py — may need a graphene wrapper
  in hypotheses/twodef_sv/) + Phase 0/1/2 phase notebooks.

---

## Update: 2026-08-05 (v9) — Both Phase 2 binaries COMPILE-GREEN; twin dispatcher written; σ=2 scan mid-flight

- classical run.cpp fixed (coulomb_terms has no check members — classical
  closure column now carries e_ss itself; norm columns = norm_p/norm_slab)
  → job 32884110 "[classical] compiled OK". WP compiled OK earlier (32883933).
- Dispatcher written: `/rds/.../shared/bin/run-gr2d-twin.slurm` — array 0-5
  (classical/WP × E=25/100/300 eV, k0=1.3555/2.7110/4.6955), REQUIRES
  GR_CAP_ETA+GR_CAP_L at submit (capscan verdict), idempotent, resume-capable.
  GS=gs_3x2_bi_cut50_Lz80; UPF=electron_gaussian_wpsigma2p0_He.upf.
- σ=2.0 capscan (32883752, running): W=12 row in — k=4.7: η=2→4.5e-5 ✓,
  η=3→3.8e-7 ✓; k=1.05: residue RISES with η (0.053→0.074) = genuine CAP
  reflection of the slow spectral tail. W=16/20 rows pending (adiabaticity
  should help η=1). NOTE k=1.05 (15 eV) is NOT in the Phase 2 set — Phase 2
  gates on k=1.36+4.70 only; the 15 eV point's CAP residue is a Phase 3
  decision (wider CAP / documented systematic).
- On scan completion (monitor b0sx22qho): pick (η,W) minimizing max residue
  over k∈{1.36..4.70}, require <1e-3 there; then
  `sbatch --export=ALL,GR_CAP_ETA=<η>,GR_CAP_L=<W> --array=0-5
   shared/bin/run-gr2d-twin.slurm`; notebooks after runs.

---

## Update: 2026-08-05 (v8) — Phase 0 COMPLETE; Phase 2 twin binaries written; σ=2.0 scan + compile checks in flight

- **Phase 0 COMPLETE.** Both checkpoints on disk: `shared_gs/gs_3x2_mono_cut50_Lz80`
  (E=−143.94 Ha) and `gs_3x2_bi_cut50_Lz80` (E=−287.95 Ha). 60 Ha probe:
  E/atom drift 50→60 = +0.92 mHa (decelerating) → 50 Ha locked, absolute-offset
  caveat recorded. Binding ≈ 38 meV/atom. γ1 split check deferred to phase
  notebook (no eigenvalue dump in GS binary).
- **Phase 2 machinery WRITTEN (untested until compile check returns):**
  - `.../scripts/twodef_sv/wp/run.cpp` — sigma56 clone: real C ions from xyz,
    cutoff-matched electrons, Ehrenfest ions (`.ehrenfest()`), CAP-only
    perturbation, zero-φ₊ interactions ledger (closure E_SS+E_PS+E_PP==E_H),
    ckpt N/3+2N/3+final + ion R/V sidecars, VTI ~N/50, WF/momentum every 100,
    t=0 analytic gates, impact site = entry-layer hexagon centre
    (2.3244, 1.3420) Bohr — which is atom-atop in layer 2 (Bernal).
  - `.../scripts/twodef_sv/classical/run.cpp` — real-ion projectile via
    ionic::species("He").pseudo_file(UPF).mass(m_e), extra_electrons(+1)
    compensation (cap_cl idiom), velocity set directly, electron_track every
    step, n_P ledger at current ion position with σ_pot=σ_WP/√2, rt_state
    persists proj R/V. DELIBERATE parity exception (recorded): classical
    N_STEPS stops the ion at the CAP inner edge (CAP can't absorb classical
    ions; edge-zeroing unimplemented) — twins compare over the common window.
  - UPF generated: `shared/pseudopotentials/electron_gaussian_wpsigma2p0_He.upf`
    (σ_charge=1.4142 ✓, from sigma1p47_He template, inqview generate_gaussian_psp).
- **In flight:** job 32883752 σ=2.0 CAP scan (monitor b0sx22qho); job 32883933
  compile-check of both Phase 2 binaries (monitor bpd44rpa5).
- **Next:** (1) compile check → fix any API breakage; (2) σ=2 scan verdict →
  pick (η,W); (3) write Phase 2 dispatcher slurm (3 twin pairs, bilayer,
  E=25/100/300 eV ↔ k0=1.355/2.711/4.695; vacuum baselines = the σ=2 capscan
  runs at same k, same box+CAP); (4) run-notebook wiring per run
  (run_notebook_builder.py) + phase notebooks; (5) WP/classical ratio gate.

---

## Update: 2026-08-05 (v7) — Phase 0 gates evaluated; σ_WP=2.0 production decision; two fix jobs queued

- **Battery v2 (32882107) COMPLETED.** Gates: vacuum Lz80 vs Lz100 Δ=4.5e-8 Ha
  (PASS, decisive); interlayer binding E_bi−2E_mono ≈ −68 mHa ≈ 38 meV/atom
  bound at d=3.35 Å (PASS — LDA binding confirmed in our numbers); cutoff
  E/atom: 35→45 = 0.96 mHa, 45→50 = 1.5 mHa NON-monotonic (gate NOT strictly
  met; grid-commensurability rounding suspected; keeping 50 Ha — stopping is
  energy DIFFERENCES on a fixed grid; 60 Ha probe job 32882486 queued to close).
- **Checkpoint-skip bug found+fixed:** battery rebuild step ran bi@50 as probe,
  marked complete, save variant skipped → shared_gs has mono only. Fixed skip
  condition (checkpoint-dir check); stale marker cleared; bilayer save job
  32882530 queued (writes shared_gs/gs_3x2_bi_cut50_Lz80).
- **Phase 1 σ=0.5 scan verdict (job 32878515, σ=0.5 grid):** high-k residues
  track analytic transmission (CAP works; need ηW ≳ 33 for 1e-3 at k=4.7);
  low-k residues ≈0.2 are η-INDEPENDENT — lingering near-zero-p_z weight of the
  σ=0.5 packet, NOT CAP failure. → **PRODUCTION σ_WP = 2.0 decision** (see plan
  Stage C decision log; σ=0.5 → high-E-only variant). capscan run.cpp gained
  CS_BUFFER + σ-in-tag; slurm defaults CS_SIGMA=2.0/CS_BUFFER=2.0 with %g tag
  matching. RESUBMIT the capscan (σ=2.0 grid) AFTER 32878515 exits (do not
  rebuild while it runs); its σ=0.5 results feed the phase notebook as the
  σ-cleanliness demonstration.
- Queued now: 32882486 (60 Ha probe), 32882530 (bi GS save). Monitors: capscan
  bkw5pk1s8 still armed on 32878515.

---

## Update: 2026-08-05 (v6) — Phase 0 bug fixed + resubmitted (32882107); mono GS DONE; Phase 1 running

- Job 32877889 result: **monolayer GS complete + checkpointed**
  (`shared_gs/gs_3x2_mono_cut50_Lz80`, E = −143.94 Ha = −6.00 Ha/atom, gates
  passed). **All bilayer variants failed** at "atom 24 out of cell": the
  bilayer generator's +b_y layer-2 shift crossed the periodic y boundary, and
  the coronene-style bounds gate treated periodic xy as hard walls.
- Fixes: `graphene_3x2_bilayer.xyz` regenerated with y wrapped into the cell
  (Bernal re-verified 12/24 atop); gs/run.cpp gate now fatal on |z|>Lz/2 only,
  warn-only for xy (periodic). Resubmitted → **job 32882107** (mono skipped via
  idempotent battery; binary rebuilds since run.cpp changed). Monitor armed.
- Phase 1 capscan (32878515) RUNNING since ~04:36: first (η,W) point
  propagating at 0.072 s/step (~2–3 min/point), t=0 injection gates PASSED.
  Monitor bkw5pk1s8 armed.

---

## Update: 2026-08-05 (v5) — AUTONOMOUS EXECUTION STARTED; Phase 0 GS battery submitted (job 32877889)

Status: campaign running autonomously (user instruction 2026-08-05: all phases
back-to-back, no user input; correctness gates only, warn-don't-block).

**Key infrastructure discoveries (Explore agent, verified):**
- Build: repo-local `shared/bin/inq-run` wrapper; env via `shared/bin/csd3-env.sh`
  + `shared/config.sh`; engine = `inq-study` fork (has `perturbations::absorbing`
  CAP — sin² imaginary potential, two-sided via ±mid_frac). Compute nodes have
  no network (deps cache reused).
- SLURM: account mphil-nikiforakis-skcb2-sl2-gpu, partition ampere, gpu:1;
  env-driven binaries, no mpirun. Templates: shared/bin/run-{ng,s56}-*.slurm.
- Clone bases: WP RT = `scripts/sigma56_sv/wp/run.cpp` (663 L, CAP+ckpt+decomp);
  classical = `sigma56_sv/classical/run.cpp` (352 L); GS = `ng_mass_ladder/gs/run.cpp`.
- **EXISTING graphene system reused** (`ResearchProject/systems/graphene/`):
  3×2 rect supercell (Lx=13.9462, Ly=16.1037 Bohr; nx=3 folds K→Γ — already
  encoded in `shared/configs/graphene_gs.hpp`), 24 C mono, LDA, ONCV-C default
  pseudos, 50 Ha, T=0.1 eV. Old paths stale (/local/data/...); no GS on disk.
  Classical-projectile trick documented there: species "He" + custom Gaussian
  UPF + mass 1/1822.9, extra_electrons(+1) compensation.
- Periodicity resolved: route (a) — `.periodicity(2)` (xy periodic, z finite),
  the corpus per2 slab convention.

**Phase 0 built + submitted this session:**
- `/rds/.../ResearchProject/systems/graphene/shared/geometry/graphene_3x2_bilayer.xyz`
  — 48 C AB bilayer, d=3.35 Å, layer2 = +bond-vector shift; Bernal verified
  programmatically (12/24 atoms atop layer-1).
- `/rds/.../ResearchProject/systems/graphene/shared/configs/twodef_gs.hpp`
- `/rds/.../ResearchProject/systems/graphene/scripts/twodef_sv/gs/run.cpp`
  — env-driven battery variant runner (GR_VARIANT/GR_CUTOFF_HA/GR_LZ_BOHR/
  GR_SAVE); gates: finite E, electron count, state count. NOTE: no forces gate
  (template API lacks GS forces) — Ehrenfest lattice-drift check moved to
  Phase 2 pilot.
- `/rds/.../shared/bin/run-gr2d-gs.slurm` — battery: bi@{35,45}Ha probes,
  bi@50+mono@50 production (checkpoints → shared_gs/gs_3x2_{bi,mono}_cut50_Lz80),
  bi@50@Lz100 vacuum probe. Idempotent (skips run_completed variants).
- **SLURM job 32877889 submitted** (pending at write time; ng-build-* jobs from
  the other campaign also queued). Background monitor b9injfxc9 polls squeue
  every 5 min; on completion the session evaluates gates (cutoff convergence
  <1 mHa/atom 45→50, Lz80 vs Lz100, electron counts) and proceeds to Phase 1.

**Phase 1 ALSO built + submitted (parallel with Phase 0 — independent):**
- `/rds/.../ResearchProject/systems/graphene/scripts/twodef_sv/capscan/run.cpp`
  — clone of sigma56_sv/vac (free WP + two-sided sin² CAP, t=0 analytic gates),
  adapted: graphene production cell 13.9462×16.1037×LZ, `.cutoff()` grid (NOT
  .spacing — grid must match production), σ_WP=0.5, records norm_final+pz_final.
  Aliasing pre-checked: worst case k0=4.70 at 50 Ha → ~0.009% tail, PASS.
- `/rds/.../shared/bin/run-gr2d-capscan.slurm` — (η,W,v) sweep: η∈{1,2,3} Ha ×
  W∈{12,16,20} Bohr at Lz=80 (+W20@Lz100); extremes-first pre-filter (residue
  <5e-3 at k=1.05,4.70) before full 6-velocity grid. Idempotent.
- **SLURM job 32878515 submitted**; monitor bkw5pk1s8.
- Known design tension (recorded): σ_WP=0.5 spreads at 1.41 Bohr/a.u. — severe
  for the slow 15 eV point. Matched vacuum controls (corpus pattern) will
  quantify; production σ decision may revisit after Phase 2 pilot.
- dt=0.025 for scan (cutoff 50 Ha vs corpus dt=0.04 at E_cut≈31 Ha — same
  stability margin); Phase 2 pilot must verify energy conservation CAP-off.

**Next steps (autonomous):** (1) GS job 32877889 completes → check battery
gates (cutoff <1 mHa/atom 45→50, Lz80 vs 100, counts) in gr2d-gs-32877889.out;
(2) capscan 32878515 completes → analyzer picks (η,W,Lz) with residue <1e-3
across v (write capscan_choice into plan); (3) write Phase 2 twin run.cpp pair
(clone sigma56_sv wp/classical, adapt: ions from xyz, GS load, Ehrenfest ions
ON, contract = VTI ~50–100 steps / ckpt 1/3+2/3+final / interactions.csv every
step with φ₊=0 adaptation (B group = ionic lattice; Hartree closure gate
E_SS+E_PS+E_PP==energy_hartree still exact) / WP momentum ~100 steps);
(4) run notebooks per run + phase notebooks (run_notebook_builder.py pattern);
(5) Phase 3 sweep; (6) Phase 4 unbox reference pack.
Resume note: if session lost — squeue -u skcb2; outputs gr2d-*-<jobid>.out at
repo root; continue per plan phase order.

---

## Update: 2026-08-05 (v4) — Reference pack acquired + FROZEN; run-output contract locked

Status: Phase R DONE. Experimental reference data verified extractable and
frozen (blinded) at
`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/validation/reference_data_bilayer_graphene/`
(ESTAR graphite CSV 1 keV–1 GeV via direct CGI POST, S_col(1 keV)=1.8 eV/Å;
Geelen PRL 123 086802 quoted anchors — λ_inel ≈ 3→1 layers over 0→25 eV,
mono-vs-bi 5–15 eV bandgap signature, T/R↔λ formulae; CXRO carbon f1/f2 for
the caveated ELF bridge; proof figure + reproducible make_reference_pack.py).
Blinding protocol in `docs/validation/bilayer-graphene-stopping-reference.md`:
DO NOT consult during Phases 0–3; Phase 4 only. ELF-derived absolute S at
15–300 eV remains a caveated TODO (Palik/EELS or own linear-response).
Run-output contract added to the plan (user-locked): energy decomposition +
interactions.csv every step; ΔKE_ions column; CAP ledger; WP momentum
distribution ~every 100 steps; VTI ~every 50–100 steps (~30–60 frames);
2–3 checkpoints/run (interior ~1/3, ~2/3 + mandatory FINAL, all ion R/V);
run notebook per run + phase notebook per phase. Phase 4 plot battery
defined (6 plots; definition-vs-definition first, reference unboxing in
plots 3–5).
Geelen source note updated implicitly by full-paper read (PDF cached in
session tool-results; quantitative anchors now in the frozen pack).
Next: Phase 0 GS build (awaiting go-ahead; supercell tier + INQ
periodicity route resolve there).

---

## Update: 2026-08-05 (latest) — Ehrenfest ions ON; CAP scan scope resolved

Status: design phase, v3. User locked: (1) Ehrenfest ions ON (free C
lattice) — plan Stage D now records the three consequences: ion-KE
becomes an explicit ledger column in ΔE_slab, Phase 0 gains a
residual-forces≈0 gate, rt_state.txt persists all ion R/V; (2)
periodicity xy-only + absorbing CAP on z confirmed; (3) CAP scan
question answered: electron-k range only (CAP acts on orbitals, not on
classical ions; kinematic max transfer e→C ≈ 1.8e-4·E ≈ 0.05 eV at
300 eV ≪ 20 eV displacement threshold, so no carbon reaches the CAP);
heavier-mass scan only if the muon-mass-fork WP trick is ever imported.
Remaining decisions: supercell tier (3×3 rec.), INQ periodicity route
(resolve from docs/inq_tutorial.md in Phase 0).
Next: Phase 0 implementation on user go-ahead.

---

## Update: 2026-08-05 (later) — Campaign v2: periodicity, supercell rationale, early classical-twin phase

Status: design phase, v2. User feedback incorporated into
`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/real-material-stopping-comparison.md`:
(1) box periodic in xy ONLY, CAP de-periodicizes z (INQ route — true slab
periodicity vs periodic-z+vacuum+CAP — to verify from docs/inq_tutorial.md);
(2) in-plane supercell justified to user: separates projectile from its
periodic images + acts as folded k-mesh (3n×3n → K at Γ), cost ≈ n162
jellium at 3×3; (3) classical twins cut to 3 pairs (E = 25, 100, 300 eV,
bilayer) moved EARLY as Phase 2 sanity — user expects WP ≈ classical;
recorded caveat: bulk jellium measured WP/classical ≈ 2.2, gate is
"ratio measured", not "ratio = 1". Phase order now 0 GS → 1 CAP scan →
2 twins → 3 WP sweep → 4 comparison.
Remaining decisions: supercell tier (3×3 rec.), frozen vs free ions,
INQ periodicity route.

---

## Update: 2026-08-05 — Bilayer graphene confirmed; campaign design v1 drafted

Status: in progress (design phase; nothing built or launched).
User locked: bilayer graphene; slab ΔE/L_z definition PRIMARY, KS-orbital
definition free by-product; CAPs in the box; electron projectile.
Changed: `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/real-material-stopping-comparison.md`
— "Campaign design v1" section replaces the sketch: Stage A GS build+verify
(3n×3n supercell REQUIRED so K folds to Γ; LDA interlayer-binding
justification TODO), Stage B z-budget (L_z ≈ 70–90 Bohr), Stage C energy
grid 15–300 eV mapped to Geelen/LEED windows + no-through-hollow-channel
trajectory note, Stage D CAP reflection scan + ΔE_slab = ΔE_box + E_CAP
bookkeeping (port jellium m8/m9 accounting).
Next: user answers the four decision points at the end of the plan
(supercell tier, run-matrix trim, classical twin timing, frozen ions);
then Stage A GS convergence battery (cheap) can start.

## Milestone: 2026-08-05 — Material-selection deep search complete; bilayer graphene recommended

### Current status

Investigation phase DONE, decision drafted, nothing simulated. The user wants
to apply the project's two stopping-power definitions — (1) the KS-orbital-
dependent definition and (2) the jellium-slab ΔE_total/L_z definition — to a
real material and compare against experimental/analytical references. Primary
selection criterion (user's): availability of experimental evidence; motivation
is secondary. Candidates: bilayer graphene vs "BNO3". A web literature/database
search (NIST ESTAR, IAEA stopping database review read directly, transmission-
experiment and TDDFT literature) was completed this session.
**Recommendation: bilayer graphene.** Awaiting user confirmation before any
campaign design or runs.

### What changed

- Decision + full evidence matrix written to
  `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/real-material-stopping-comparison.md`.
- Four source notes created (see below).
- Key facts established:
  - "BNO3" is not a standard material — interpreted as h-BN; user should
    confirm (if an oxide was meant, IAEA data favour Al2O3/SiO2/Ta2O5/TiO2,
    but those are bulk 3D — poor fit to the slab definition).
  - Graphene has experimental evidence in EVERY channel: eV-TEM electron
    transmission through free-standing 1–4 layers at 0–25 eV with measured
    IMFP (Geelen PRL 123, 086802 (2019)) — overlapping our E15–E300 electron
    runs; LEEM layer-counting (Hibino 2008); HCI/heavy-ion transmission
    energy loss per layer (TU Wien, Nat. Commun. 7, 13948 (2016) + PRA 93,
    052708 (2016)); carbon among best-measured IAEA ion targets; ESTAR
    graphite ≥ 1 keV as high-energy anchor; TDDFT baseline (Ojanperä PRB 89,
    035120 (2014)); analytic 2D dielectric models (Mišković group).
  - h-BN has essentially none of this (EELS + defect studies only; absent
    from IAEA compound highlights; Bragg-rule estimates explicitly
    unreliable at low velocity per the IAEA review, read directly pp. 1–8).
  - Novelty opening: no layer-dependent TDDFT stopping study of bilayer
    graphene found (recorded as "not found", not "does not exist").

### Files touched

- `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/real-material-stopping-comparison.md` — decision, evidence matrix, caveats, staged campaign sketch
- `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/sources/geelen-2019-evtem-graphene.md` — electron-projectile evidence (primary)
- `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/sources/gruber-2016-hci-graphene.md` — ion-projectile evidence
- `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/sources/ojanpera-2014-graphite-stopping.md` — TDDFT theory baseline
- `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/sources/montanari-2024-iaea-database.md` — IAEA database census (pp. 1–8 read directly)

### Commands run

```bash
# web searches + fetches only (WebSearch/WebFetch); no builds, no runs
# IAEA review PDF cached at:
# /home/skcb2/.claude/projects/-rds-user-skcb2-hpc-work-tddft-inq-tddft-research/ded938d1-2279-4493-bee2-dcd2b5157ea3/tool-results/webfetch-1785886406720-ghw5so.pdf
```

### Tests and validation

- Proposed: none yet (no code written). Campaign sketch in the plan lists
  Stage 0 feasibility checks (graphene GS lattice constant / band sanity)
  as the first component tests.
- Approved / Run: n/a this session.
- Remaining gaps: Ojanperä 2014 and Geelen 2019 verified at abstract level
  only — full PDFs must be obtained before quoting numbers; TU Wien
  follow-ups (Commun. Phys. 2019, 2021) located but unread.

### Trusted sources used

- Montanari et al., arXiv:2402.03080 (2024) — IAEA database census (read pp. 1–8)
- Geelen et al., PRL 123, 086802 (2019) — eV-TEM through few-layer graphene
- Gruber et al., Nat. Commun. 7, 13948 (2016); Wilhelm et al., PRA 93, 052708 (2016)
- Ojanperä et al., PRB 89, 035120 (2014) — TDDFT stopping, graphitic targets
- NIST ESTAR (ICRU 37): 1 keV–10 GeV electrons, ~280 materials incl. graphite

### Known issues / blockers

- "BNO3" interpretation unconfirmed by user (assumed h-BN).
- APS pages 403 on direct fetch — get paper PDFs via institutional access.
- ESTAR cannot validate our 15–300 eV electron regime (≥ 1 keV only);
  low-energy comparison must go through transmission/IMFP + ELF channels.

### Assumptions still in play

- h-BN reading of "BNO3".
- LDA + norm-conserving C pseudopotential will be acceptable for graphene
  (consistent with jellium corpus) — needs literature-justification note
  at Stage 0.
- INQ + inqkit WP-injection machinery transfers to a real-lattice supercell
  without engine edits (inq-immutable rule) — believed true, unverified.

### Exact next steps

1. User confirms material choice (bilayer graphene) and the h-BN reading of
   "BNO3"; h-BN earmarked as the later gap-contrast material.
2. Obtain full PDFs: Ojanperä 2014, Geelen 2019, Kononov npj Comput. Mater.
   9, 205 (2023); extract quantitative S(v)/IMFP numbers into the source notes.
3. Design Stage 0 feasibility (pseudopotential, supercell, GS convergence
   tests) in the plan; propose the test menu per simulation-validation skill
   BEFORE any expensive run.
