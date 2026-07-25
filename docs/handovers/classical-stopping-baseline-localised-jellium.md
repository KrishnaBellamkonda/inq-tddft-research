# Handover — Classical stopping baseline for the localised jellium slab

Campaign: `docs/campaigns/localised_jellium/classical-stopping-baseline.md`
(id `classical-slab-stopping`, area `localised_jellium`).
Branch: `overnight-gaussian-classical`. Started 2026-07-15 (grill-with-docs).

## Goal (one line)

Produce the **matched localised-slab classical stopping baseline** (twin of the WP
`qsp_phase5` S(E) sweep) so the slab WP quantum stopping has a like-for-like
classical expectation instead of the ADR-0010 bulk/point-charge reference.

## Design (locked via grill)

- Projectile = moving Gaussian **charge** (perturbation, not ghost UPF), −1,
  σ_WP=0.5 ⇒ σ_pot=0.354. Twin of `p5_wp_v1p3` (σ_WP=0.5, v=1.3).
- **Two phases, run concurrently (one per GPU):**
  - **Phase 1 (GPU 0) — Ehrenfest, light electron (mass 1).** Decelerates
    (`light-projectile-stopping`); enters slab at v₀, stops near far face. Metric =
    **initial-drag slope** S(v₀) = −d(KE_proj)/ds over v≥0.85·v₀, cross-checked by
    +d(E_deposited)/ds. `LJ_CONST_V=0`, N_STEPS=2000, SAVE_EVERY=7.
  - **Phase 2 (GPU 1) — prescribed constant velocity.** Zero force ⇒ R=R₀+V₀·t;
    transits slab, exits. Metric = **ΔE_deposited/L_slab** off the plateau.
    `LJ_CONST_V=1`, N_STEPS=1034 (center −23.75→+30, no periodic wrap), SAVE_EVERY=4.
- **Geometry (EXACT twin of p5_wp_v1p3):** cell 50×50×90, dx=0.5, **periodicity 3**
  (`.periodic()`), slab half-width 12.5 (axis z), N=82, launch_z=−23.75, dt=0.04.
  GS reused: `shared_gs/slab_n82_L50x50x90/` (per-3). **No CAP** (note: the WP twin
  has a CAP at ±35 for its dispersing tail; the classical charge needs none — pair
  matched on bath physics + cloud, not byte-identical).
- **Deliverables (user chose full scope, Q6b):** both runs → per-run `analyse.py` +
  single-run `run-notebook` with `density_evolution.gif` **at the top** (2026-07-15
  rule) → `twin-run-analysis` pairwise P/S/B decomposition vs `p5_wp_v1p3` → one
  comparison figure: S_classical(P1 initial-drag), S_classical(P2 ΔE/L),
  S_WP=2.37 eV/Bohr [UPPER BOUND], bulk classical 0.94, Lindhard.

## Files

- Binary: `scripts/classical_slab_stopping/run.cpp` (cloned from
  `localised_jellium_dynamics/proj_dyn`, added `LJ_CONST_V` zero-force mode +
  const-v run_summary field). Build dir `.../classical_slab_stopping/build/`.
- Orchestrator: `scripts/classical_slab_stopping/orchestrate.py` (concurrent 2-GPU,
  idempotent resume, per-phase try/except + Gmail).
- Outputs: `scripts/classical_slab_stopping/results/{p1_ehrenfest_v1p3,p2_constv_v1p3}/`.
- Analysis (to be created): `hypotheses/classical_slab_stopping/`.
- Campaign + CONTEXT glossary ("classical stopping baseline", "drive mode") updated;
  `docs/campaigns/INDEX.md` regenerated.

## Status (2026-07-15)

DONE / verified:
- [x] `run.cpp` const-v mode written. **code-test:** the zero-force ⇒ constant-velocity
  path is the locked unit test `test_projectile.cpp:22` (ctest green, 125 assertions).
- [x] **Cutoff guard PASS** — classical, dx=0.5, E_kin=23 eV: E_cut=537 eV ≥ 25 eV.
- [x] GPUs: both free (~23.8 GB), no other users, compute OK (NVML mismatch is cosmetic).
- [x] Campaign file, CONTEXT.md, INDEX.md, orchestrator, this handover.

- [x] **GPU smoke PASSED** (both). Binary builds as `./run` (NOT `build/run`).
  Phase-1: smooth energy ledger, no t=0 kick (per-3 GS matches), proj decelerates
  1.3→1.2984. Phase-2 const-v: proj_z EXACTLY linear (−23.75+1.3·t, step Δ=0.052),
  proj_vz≡1.3, KE≡0.845 Ha. (First build hit `/tmp` 100%-full — FIXED via
  `TMPDIR=/local/data/.../.build_tmp`; `/tmp` is a tiny 9.8 G volume, ALWAYS set TMPDIR.)
- [x] **PRODUCTION LAUNCHED 2026-07-15 18:55** — `orchestrate.py` (nohup, PID 3497162
  at launch). Phase 1 GPU 0 (2000 steps), Phase 2 GPU 1 (1034 steps), concurrent,
  checkpointed. Both GPUs ~19 GB used. Logs: `run_p1_ehrenfest_v1p3.log`,
  `run_p2_constv_v1p3.log`, `orchestrate.log`. Per-phase Gmail on completion.

PENDING (post-processing, runs AFTER sims complete ~1–2.5 h):
- [ ] Per-run `analyse.py` (full inqview pipeline → REPORT.md) for both runs.
- [ ] Single-run `run-notebook` each with `density_evolution.gif` at TOP.
- [ ] `twin-run-analysis` `twin_decompose.py` on Phase 1 vs `p5_wp_v1p3` (P/S/B ledger,
  SIE residual).
- [ ] Comparison figure: S_classical(P1 initial-drag), S_classical(P2 ΔE/L),
  S_WP=2.37[UB], bulk 0.94, Lindhard → `hypotheses/classical_slab_stopping/`.
  (These are NOT yet wired into the orchestrator — build them against the real
  outputs once `run_completed=true`, then flip campaign tasks 5–6 + status→done.)

## RESULTS (2026-07-15, both runs complete)

- **P1 (Ehrenfest, light electron):** proj decelerated v 1.30→0.52, lost 19.3 eV KE,
  transited the slab (did NOT stop — proj_z_final=38.8). **S(v₀=1.3) = 0.49 ± 0.004
  eV/Bohr** — initial-drag −dKE/ds over the in-slab v≥1.1 window (n=282, mean v=1.22).
  (Raw E_total is coupling-contaminated: peaks +68 eV mid-slab; NOT the clean deposit.)
- **P2 (const-v):** transited at v≡1.3, center stopped at +30 (no wrap). **S = 0.43 ±
  0.18 eV/Bohr** — coupling-subtracted deposit/L_slab (raw E_total/L=0.40, consistent).
- **Both agree ⇒ S_classical(v=1.3) ≈ 0.43–0.49 eV/Bohr.** N-drift 0.00% (no CAP).
- **Headline:** S_WP(p5_wp_v1p3)=2.37 [UPPER BOUND] is ~5× the classical baseline; below
  bulk σ=0.5 (0.94) and Lindhard-point (0.57). Inference: the WP UB is likely inflated
  (SIE / spreading / 2nd-moment), not a 5× real enhancement — benchmark quantum vs ~0.45.
- Artefacts in `hypotheses/classical_slab_stopping/`: `analyse_classical_baseline.py`,
  `S_comparison.png`, `stopping_extraction.png`, `S_summary.csv`, `density_evolution.gif`
  (65-frame n(z,t), P1 vs P2, physical order).
- Result email sent (`[classical-slab-stopping] RESULT …`).

## REMAINING (campaign tasks 5–6, not blocking the headline)

- [ ] Full single-run `run-notebook` per phase — **needs builder adaptation** (attempted
  2026-07-17, timed out at 300 s + convention mismatch): `run_notebook_builder.py` keys
  classical detection + the projectile/stopping panels off `raw/observables/electron_track.csv`
  with ion-track column names, but these proj_dyn-derived runs emit `projectile.csv`
  (cols step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal). To wire it up:
  either (a) teach the builder to also read `projectile.csv` (map proj_z→ion_z,
  proj_vz→v, energy_proj_ke→KE), or (b) emit an `electron_track.csv` alias in the run
  layout. The `inqview.pipeline` runner also overran 300 s here — raise the timeout /
  restrict phases. A partial unexecuted `p1_ehrenfest_v1p3.ipynb` (18 KB) was left; treat
  as scratch. The core figures (`S_comparison.png`, `stopping_extraction.png`,
  `density_evolution.gif`) already stand alone in the hypotheses folder.
- [ ] `twin-run-analysis` `twin_decompose.py`: **BLOCKED** — `p5_wp_v1p3` is NOT a valid
  twin for the engine (verified 2026-07-17): (1) it has **no `interactions.csv`** and its
  `observables.csv` lacks the pairwise/`energy_proj_bg_ideal` columns the engine needs;
  (2) **CAP parity mismatch** — the WP run has a two-sided CAP (drains N), the classical
  runs don't, so the skill's N-conservation guard + gauge test (ΔE_SS≈0) fail; (3) no
  validated `twin_manifest.json`. A clean residual/SIE needs a **CAP-free, pairwise-
  emitting WP twin** at σ=0.5, v=1.3 (new follow-on run via `twin-run-generation`).
  What we CAN say analytically now (interpretation rules): expected `dKin` ≈ 3/(4σ²)+k0²/2
  = 81.6+23 ≈ 105 eV WP-localisation kinetic; expected σ=0.5 `SIE ≈ 4.3 eV` (LDA
  one-electron self-interaction) — a real but small (~4 eV) contribution, NOT enough to
  explain the ~1.9 eV/Bohr·25 Bohr ≈ 48 eV S_WP−S_classical gap on its own ⇒ the WP
  UPPER-BOUND S is dominated by extraction inflation (2nd-moment / spreading), consistent
  with the qsp_phase5 convergence flag.

## Resume recipe

- Check progress: `tail run_*.log`; `grep step run_p1_ehrenfest_v1p3.log | tail`.
- If a phase died: re-run `nohup venv/bin/python3 orchestrate.py &` — it is idempotent
  (skips completed runs; else `LJ_RESUME=1` from the last checkpoint).
- To extend a completed run: raise `LJ_N_STEPS` for that phase in `orchestrate.py`,
  re-run (resumes from checkpoint).

## Gotchas / rationale to survive compaction

- **`/tmp` is 9.8 G and full** → every INQ build here MUST set
  `TMPDIR=/local/data/public/skcb2/tddft/.build_tmp`, else nvcc dies "No space left".
- **periodicity 3, not 2.** proj_dyn defaulted to per-2; the twin `p5_wp_v1p3` uses
  `.periodic()` (per-3, run.cpp:77). A per-2 RT on a per-3 GS gives a spurious t=0 kick.
- **Do NOT gate on velocity drift** (light-projectile rule): Phase-1 is SUPPOSED to
  decelerate. Gate on a clean initial-drag slope existing (≥30 early points), not v-drift.
- **S_WP=2.37 is an UPPER BOUND** (convergence-flagged in qsp_phase5) — label it so.
- Engine: classical run is real-valued → default `inq/` (not inq-study; that was only
  for the WP's complex CAP). Correct and faster.
- Est. wall-clock: Phase 1 ~2000 steps, Phase 2 ~1034 steps; twin WP took ~4.2 s/step
  ⇒ P1 ≈ 2–2.5 h, P2 ≈ 1–1.5 h (concurrent).
