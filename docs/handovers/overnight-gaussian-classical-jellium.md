# Handover: Overnight Gaussian classical projectile in jellium — S(v)

Rolling handover. Task plan: `docs/plans/overnight-gaussian-classical-jellium.md`.
Branch: `overnight-gaussian-classical` (off `main`, which now includes the merged
`rejuvenation/claude-ecosystem`). Started 2026-06-11T23:59Z, 10 h budget.

## MORNING SUMMARY (read first)

What this run does: a classical −1 / σ=0.5 Bohr erf-smoothed Gaussian electron
(mass = m_e, NOT a fictitious mass) decelerating by free Ehrenfest dynamics
through the validated r_s=5.69 jellium bath (N=162, L=50, dx=0.40), to map the
electronic stopping power S(v) and compare to a corrected Lindhard reference.

**THE DOMINANT CONSTRAINT (read this):** the 125³-grid / 101-state propagation
costs **~15–16 s/step** on this GPU. This makes the plan's full Stage-4 kick
spectroscopy (needs ~20k steps ≈ 78 h) and full-traversal Stage-6 runs
physically impossible in 10 h. The work was **re-scoped** to what fits:
short decel-segment S(v) runs + the analytic Lindhard reference + the kick
*capability* (cost-blocked for production). This is honest graceful degradation,
not silent truncation.

### Results — COMPLETE (all 6 S(v) runs done)
| v | measured S (Ha/Bohr) | Lindhard S_LR | ratio | regime |
|---|---|---|---|---|
| 2.98 | 0.0082 | 0.0067 | 1.22 | above LR |
| 1.94 | 0.0173 | 0.0131 | 1.33 | above LR |
| 1.13 | 0.0336 | 0.0301 | 1.12 | above LR |
| 0.58 | 0.0387 | 0.0543 | 0.71 | BELOW LR (Barkas) |
| 0.40 | 0.0360 | ~0.032 | 1.12 | friction tail (noisy) |
| 0.77 (σ=0.4) | 0.0382 | 0.0480 | 0.79 | σ-sens = 1.05× σ=0.5 |

**Three headline numbers (morning ask):**
1. Bragg peak v≈1.0, S≈0.046 Ha/Bohr — shifted higher + broadened vs the Lindhard
   peak (v≈0.6/0.054): nonlinear-screening signature.
2. Z=−1 Barkas crossover confirmed (above LR for v≳1, below near LR peak v≈0.6).
3. Low-v friction Q NOT cleanly determined (extraction noisy below v≈0.5) — honest
   limitation; follow-up = longer averaging / force-based S.
- f-sum: corrected Lindhard ELF = 1.000 at all q (old pipeline.lindhard FAILS 0.01–0.13).
- **Executed `report.ipynb`** (0 errors, 4 figs) + REPORT.md + money plot all current.
- Catalogue: `docs/runs_catalogue.csv` has a run_sv_sigma0p5 row (dir-level; per-velocity
  detail in the report table — one-dir-one-run model mismatch noted).
- **Executed `docs/reports/.../report.ipynb`** (0 errors, 4 figures) + REPORT.md done.
- Refresh: `run_sv_sigma0p5/analyse_sv.py` (money plot) and
  `docs/reports/.../build_notebook.py` (notebook) — re-run as runs complete.

## Stage status

| Stage | Status | Notes |
|---|---|---|
| 0 git foundation | DONE | rejuv merged → main (dc488dc); branch overnight-gaussian-classical; plan committed |
| 1 baseline+GS+timing | DONE | GS valid (162 e⁻, SCF conv); baseline tests 114p/5xf/1xp; **16 s/step** measured |
| 2 erf psp + tests | DONE | inqview.io.gaussian_psp; 6 tests pass; σ0.5/0.4 .upf generated. Static-impurity GS run DEFERRED (psp validated analytically + smoke-loaded in DFT) |
| 3 dt check | DEFERRED | dt=0.020 inherited from validated runs; dedicated dt/2 check skipped (time budget) |
| 4 cosine kick | CODE DONE | inqkit/perturbations/cosine_kick.hpp + test (norm conservation analytic). Build+run DEFERRED. Production loss-function **COST-BLOCKED** (kick needs ~20k steps) |
| 5 Lindhard ref | DONE | inqview.analysis.lindhard_elf; **10 tests pass**, f-sum exact at all q; r_s corrected 5.74→5.69 |
| 6 S(v) ladder | RUNNING | 2 GPUs, see below |
| 7 k-points | NOT DONE | skipped (no GPU budget); log "blocked, deferred" |
| 8 report | IN PROGRESS | this handover + report/ipynb pending |

## Stage 6 — live runs (poll these)

Run dir: `/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/run_sv_sigma0p5/`
Binary `./run` (built, smoke-tested 3 steps OK). Runtime env: PROJ_V0, SV_N_STEPS,
SV_WRITE_EVERY, SV_OUT_SUBDIR, SV_PSEUDO. Dispatcher: `dispatch_ladder.sh <gpu> "<queue>"`.

- **GPU0 queue** (`dispatch_gpu0.log`): v0=3.0/300/v3p0 → v0=1.3/700/v1p3 → v0=0.6/700/v0p6
- **GPU1 queue** (`dispatch_gpu1.log`): v0=2.0/450/v2p0 → v0=0.8/700/v0p8 → σ0.4 v0=1.0/700/sig0p4_v1p0

Outputs per velocity: `results/<subdir>/electron_track.csv` (step,time,pos,vel every
step) + `observables.csv` (energy) + `run_summary.txt` (run_completed=true when done).
ETA: v0=3.0 ~02:47, v0=2.0 ~03:15, 700-step runs ~06:00, extras ~09:15.

**To assemble S(v) from completed runs:**
```python
from inqview.analysis.stopping_extract import load_track, stopping_vs_v
from inqview.analysis import lindhard_elf as E
tr = load_track(".../results/v3p0/electron_track.csv")   # m=1, axis z
v, S = stopping_vs_v(tr, transient_bohr=3.0, window=21)   # local S(v) binned by v(t)
# overlay E.stopping_power_sigma(v, E.kF_from_rs(5.69), 0.5) Lindhard curve
```

## Deliverables (committed on branch)
- inqview.io.gaussian_psp (+test) — erf psp generator
- inqview.analysis.lindhard_elf (+10 tests) — corrected ELF + σ-stopping
- inqview.analysis.stopping_extract (+2 tests) — S(v) from track
- inqkit/perturbations/cosine_kick.hpp (+deferred test)
- ResearchProject/.../shared/configs/sv_ladder_L50_sigma0p5.hpp + run_sv_sigma0p5/

## Verified vs unverified
- VERIFIED: psp form/FT/V(0); Lindhard f-sum/TF/convergence; stopping_extract on
  synthetic data; GS validity; per-step cost; smoke run loads new psp + propagates.
- UNVERIFIED: production S(v) values (runs in flight); cosine_kick compile/linearity
  (deferred); Friedel screening (static run deferred); ALDA-vs-RPA at r_s=5.69.

## Subtask-4 (ecosystem finetune) observations
- 8 plan-vs-convention conflicts caught in grilling (see plan §Subtask-4).
- **Hidden M=1836 fictitious mass** in the plan silently redefined S(v); user
  caught it → reverted to m_e + free Ehrenfest. Headline finding.
- Commit-message hook **fired correctly** (blocked a 74-char subject >72).
- r_s slip 5.74→5.69 caught against existing module's canonical value.
- **Per-step cost (14–16 s) is the binding real-world constraint** the plan's
  10 h budget ignored — argues for a standing "cost calibration before scope
  lock" step in tddft-simulations.

## Next actions (resume here)
1. Poll `run_sv_sigma0p5/results/*/run_summary.txt` for `run_completed=true`.
2. As each finishes, extract S(v), upsert tddft-run-catalogue, journal entry.
3. Assemble S(v) money plot over the Lindhard σ=0.5 curve (inqview theme).
4. Write `docs/reports/overnight-gaussian-classical-jellium/` md + executed ipynb.
5. docs/sources/ notes (Lindhard 1954; Lindhard-Winther 1964; Echenique 1981/86;
   Correa 2018). Final commit.
