# Handover: energy-oscillation debugging (post-Arm-B thread)

Rolling handover for the thread that follows the cap_fix + pbc_open_z campaigns:
identify the CLOCK of the drain-then-rise E_total oscillation. Predecessors (both
done): `/local/data/public/skcb2/tddft/docs/handovers/cap-fix-experimentation.md`,
`/local/data/public/skcb2/tddft/docs/handovers/pbc-open-z-oscillation.md`.
Advisor spec-diff analysis (documented, committed f6a85a6):
`/local/data/public/skcb2/tddft/docs/notes/oscillating-vs-clean-run-spec-comparison.md`.

## 2026-07-14 (19:45) — m=1 engine-drift rerun of the clean p3_wp run LAUNCHED

**User request:** recreate the qsp_phase3 clean run (E_total decays to a fixed
value) exactly, m=1, "using the new library (inq-study) instead of inq".

**Premise correction (verified):** the original p3_wp run was ALREADY built
against inq-study — `CMakeCache.txt` in
`/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/qsp_phase3/wp/build/`
has `inq_SOURCE_DIR=/local/data/public/skcb2/tddft/inq-study`; `run_summary.txt`
says `engine = inq-study`. So the rerun is NOT an inq-vs-inq-study twin. It IS an
**engine-drift regression test**: six inq-study headers changed after the original
binary was built 2026-06-25 (`ks_hamiltonian.hpp`, `propagate.hpp`,
`electrons.hpp`, `laplacian.hpp`, `calculator.hpp`, `initial_guess.hpp` — the
mass-fork surgery). Question: does TODAY's inq-study at m=1 still reproduce the
clean decay?

**Config (exact original recipe):** `scripts/qsp_phase3/wp/run.cpp` unmodified;
`LJ_CAP=1 LJ_DT=0.04 LJ_N_STEPS=2500 LJ_LAUNCH_Z=-23.75 LJ_WRITE_EVERY=8
LJ_WF_EVERY=8`; GS `shared_gs/slab_n82_L50x50x90` (exists, untouched). Only
operational deltas: GPU 1 (GPU 0 occupied by another user), output
`results/p3_wp_m1_rerun` (original `results/p3_wp` preserved), and
`TMPDIR=$PWD/build/tmp` because the 9.8 GB `/tmp` partition is 100% full with
other users' files (killed the first build attempt with "No space left on
device"; do NOT delete others' /tmp files — always set TMPDIR to the big disk).

**Done + verified:**
- Rebuild against current inq-study: OK (full dep rebuild — libxc/catch2/spdlog).
- 50-step smoke on GPU 1: exit 0, all energies finite, propagation ended normally
  (`wp/smoke_m1_rerun.log`, results in `wp/results/smoke_m1_rerun/`).
- **KEY EARLY RESULT — bit-for-bit reproduction:** smoke energies match the
  original `prod_wp.log` to all 12 printed decimals at steps 1, 2, 10, 25, 50
  (e.g. step 50: e = -63.414464508084 both). The post-June-25 mass-fork engine
  edits are numerically inert at m=1. Expect the full run to reproduce the clean
  decay exactly.

**In flight:** production rerun, 2500 steps, launched ~19:45 on GPU 1,
log `wp/prod_m1_rerun.log`, output `wp/results/p3_wp_m1_rerun/`. ~4.2 s/step +
VTI-write spikes → ETA ~3.5–4.5 h (original wall was 5.7 h on GPU 0). No
checkpointing in this run.cpp (exactness prioritised over the checkpoint rule —
a kill loses the run; user aware).

**Not done:** final verdict (full-trace comparison vs original observables.csv),
notebook/doc update, commit of this milestone. S1/S2/S3 discriminators from the
advisor note remain designed-but-not-launched; S3 (absorbed-fraction analysis)
is zero-GPU and can run any time.

**Verify on completion:** exit 0 + `run_completed = true` in
`results/p3_wp_m1_rerun/run_summary.txt`; diff
`results/p3_wp_m1_rerun/raw/observables/observables.csv` against
`results/p3_wp/raw/observables/observables.csv` (expect identical or
noise-level); confirm E_total(t) decays to a fixed value with no late rise.

## 2026-07-14 (23:45) — RERUN COMPLETE, VERDICT: engine drift EXONERATED at m=1

Run finished clean: exit 0, `run_completed = true`, wall 14372 s (4.0 h, GPU 1).

**Primary result — BIT-IDENTICAL reproduction:** all 313 rows × all columns of
`observables.csv` match the original `p3_wp` exactly (max |dE_total| = 0.0).
Today's inq-study (post-mass-fork headers) at m=1 reproduces the 2026-06-25
clean run bit-for-bit over the full 2500 steps. The mass-fork engine edits are
numerically inert at m=1 — engine drift is ELIMINATED as a cause of the
oscillation. The oscillating campaign runs differ from the clean run by CONFIG
ONLY (advisor note suspects S1 CAP standoff / S2 perturbation strength).

**Secondary observation (present in BOTH runs, since bit-identical):** the
"clean" run's E_total minimum sits at t = 97.9 (of 100) with a terminal rise of
0.11 eV — marginally above the 0.1 eV noise floor. Under the period-lengthening
reading this hints even the clean box may turn, with a very long period; a
longer clean run would decide. Treat as a watch-item, not a claim.

**Shape numbers (2 s.f.):** drain 130 eV over the run; E_min at t = 98;
rise 0.11 eV; last-quarter slope −8.6e-2 eV/a.u. (still draining).

**Next:** S1/S2 discriminators (designed in the advisor note) and/or the
zero-GPU S3 absorbed-fraction analysis. Rerun outputs kept at
`wp/results/p3_wp_m1_rerun/`; smoke at `wp/results/smoke_m1_rerun/`.

## 2026-07-15 — sigma1_masspair campaign planned (grill interview) + launching

**User request:** new runs in the CLEAN geometry (energy must decay), same slab/
CAP/everything, only the WP changed: σ=1, higher mass, aliasing-safe, spreading
<4% at the slab, checkpoints, full + pairwise energy decomposition (twin
contract observables, no classical twin). Plus (mid-implementation): "something
must run for the full ~9 h; a fable agent intervenes if runs stop midway."

**Plan (approved, full interview trail):**
`/local/data/public/skcb2/tddft/docs/plans/sigma1-masspair-decay-runs.md`.
Key resolution: spreading-at-arrival is MASS-FREE — σ(d)=σ0√(1+(d/(2k0σ0²))²);
only k0 controls it and the grid caps k0 (CONTEXT.md "Spreading-at-arrival
law"). User relaxed 4%→10%; locked: launch −16.5 (4σ from slab face), k0=4.5
(guard WARN, tail 0.59% — verified), mass PAIR m=2 (v=2.25, 138 eV) + m=3
(v=1.5, 92 eV), all else byte-identical to p3_wp_m1_rerun. Gate: guard +
50-step smoke incl. checkpoint/resume round-trip.

**Built artefacts (uncommitted at this milestone):**
- `scripts/sigma1_masspair/wp/run.cpp` — qsp_phase3 full suite + mass fork
  (`inverse_mass()[0][wp_idx]`, re-applied on resume) + interior ckpt every
  200 steps + `LJ_RESUME=1` segment resume (from phase5_wp/muon references) +
  per-step `interactions.csv` (compute_coulomb_wp; closure checks in-file) +
  ledger extended with energy_external/nonlocal.
- `scripts/sigma1_masspair/orchestrate.sh` — serial m2→m3 chain, GPU 1
  (SMP_GPU), idempotent (skips completed, auto-LJ_RESUME=1 partials), stall
  watchdog (no log growth 15 min → kill, exit 42). Requires pre-built ./run.
- CONTEXT.md glossary entry "Spreading-at-arrival law (k0-only)".

**Smoke gates (all passed, 2026-07-15 02:15–02:32):**
- First smoke at `sigma(1.0)` CAUGHT a σ-convention error: the sigma() param is
  the WAVEFUNCTION width (density std = σ/√2 = 0.707), so the interview's
  spreading numbers (computed for density std 1) were unreachable — min-unc
  core ~35 % at the face, measured second-moment ~73 % (orthogonalisation
  broadening σ_pz 1.10 vs 0.707 min-unc at 4 Bohr standoff; actual wrapped
  weight only ~0.01 %, the 4.5 % Gaussian estimate was moment-inflated; the
  p3 run shows the same effect at a mild 8 %).
- **Resolution (autonomous):** production runs at `sigma(√2)` → density std
  1.0 Bohr (house label σ_WP=1.41). Rationale: satisfies the user's
  twice-negotiated ≤10 % spreading bound and ALL approved plan numbers; guard
  improves WARN→strict PASS (tail 0.02 %). Documented in the plan ("σ
  convention correction") — FLAG FOR MORNING REVIEW.
- Corrected smoke: σ_z(0)=1.0000, σ_pz=0.525 (5 % birth broadening),
  pz_mean=4.497, closure ≤4e-10 Ha, projected spreading at face 10.4 %.
- Resume round-trip: ckpt 50→70 extension OK, seam ~2e-6 Ha, segment velocity
  2.253 → inverse_mass correctly re-applied; `.from50` segment CSVs written.

**PRODUCTION IN FLIGHT (launched 02:32:44, detached via setsid — survives
session loss):** `orchestrate.sh` on GPU 1: wp_m2_k4p5 (2500 steps, ~4–5 h)
then wp_m3_k4p5. Logs: `scripts/sigma1_masspair/orchestrate.log` +
`wp/wp_m2_k4p5.log`. GPU 0 carries the OTHER campaign
(localised_jellium_dynamics `proj_dyn/p5_null_s2_k4_cl`,
CUDA_VISIBLE_DEVICES=0) — no contention, verified from /proc environs.
A watcher background task notifies the session when the orchestrator ends
(success or watchdog exit 42/crash) → then spawn a fable agent to diagnose +
relaunch (idempotent auto-resume, max 200 steps lost).
**On completion:** verify both run_summaries `run_completed=true`, spreading
at face vs 10.4 % projection, E_total decay shape vs the clean rerun, pairwise
ledger narrative (twin-run-analysis skill), study notebook, commit.

## 2026-07-15 (15:25) — CHAIN COMPLETE; ARTIFACT RETURNS IN CLEAN GEOMETRY

Both runs completed clean mechanically (rc=0, run_completed=true; m=2 3.3 h,
m=3 3.6 h; chain 02:32→09:26; watcher had a pgrep self-match bug → no
notification, found at 15:22 status check; GPU idle 09:26→15:22).

**HEADLINE: the drain-then-rise artifact is BACK, at full strength, in the
CLEAN geometry** (same box/CAP/GS that decays monotonically at σ=0.5/m=1):
- wp_m2_k4p5: drain 60 eV → E_min at t=32.6 → RISE +180 eV by τ=100.
- wp_m3_k4p5: drain 44 eV → E_min at t=50.2 → rise +125 eV.
- t_min ratio 50.2/32.6 = 1.54 ≈ velocity ratio 1.5, and t_min ≈ 1.44× the
  WP→far-CAP arrival time in BOTH → **the clock tracks the PROJECTILE's CAP
  arrival, not slow-spill** (new, sharpest clock datum).
- H1 (CAP standoff geometry protects) REFUTED as sufficient. Mass correlation
  returns: every oscillating run has m_eff>1 (2.10 family; now 2 and 3);
  every clean run has m=1. Suspect: CAP-filtered norm-divided kinetic ledger
  applied to the massive WP state itself.
- Spreading at face 14.0/15.1 % (free-flight projection 10.4 %; slab
  interaction adds the rest). Absorbed ≈1.01/1.02 e (WP fully absorbed).
- Per-step pairwise interactions.csv exists for BOTH runs — the decomposition
  can now finger WHICH term rises. NOT yet analysed.

**Next:** pairwise-ledger analysis (twin-run-analysis interpretation rules),
study notebook, commit of scripts/sigma1_masspair + plan + handover.
