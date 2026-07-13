# Handover: CAP energy-artifact removal (cap_fix experimentation)

Campaign: `/local/data/public/skcb2/tddft/docs/campaigns/localised_jellium/cap-fix-experimentation.md`
(`id: lj-cap-fix-experimentation`, status: running)
Loop state (source of truth for experiments):
`/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/cap_fix/`
— `autoresearch.md` (contract), `autoresearch.jsonl` (results), `experiments/worklog.md`
(narrative), `autoresearch-dashboard.md` (regenerated per run).

## 2026-07-13 — session 1: harness built, loop launched

**Goal (user, locked in-conversation):** find a CAP *setup* that removes the diagnosed
drain-then-rise energy artifact from ALL localised-jellium runs. Arms named by the user:
(1) unified wrap-around CAP; (2) tune η / CAP length / other params until the effect
disappears. Method: reputable online experiment-loop skill + small runs + iterate; start
immediately on GPUs.

**Skill:** installed `autoresearch` from github.com/drivelineresearch/autoresearch-claude-code
(MIT, v1.1.0, 327★) into `/local/data/public/skcb2/tddft/.claude/skills/autoresearch/`
(SKILL.md + ar-log.sh + LICENSE + PROVENANCE.md). Hooks NOT installed — loop is
agent-driven this session; state files follow the protocol exactly so the full plugin
can resume.

**DONE (verified):**
- `inq-stack/include/inqkit/perturbations/absorbing_wrap.hpp` — single cos² CAP bump
  peaking AT the periodic z-boundary (fractional coords; adapted from
  `inq/src/perturbations/absorbing.hpp`, MPL-2.0 credit in header; `inq/` untouched).
  Width 30 Bohr ⇒ same |z|>25 footprint and same ∫W dz = 15η as the default two-sided
  pair — topology is the ONLY difference. Discovery en route: the two-sided sin² bumps
  are ZERO exactly at the periodic boundary (verified analytically + in test).
- `ResearchProject/systems/localised_jellium/scripts/cap_fix/run.cpp` — ablation-binary
  clone + `EM_CAP_MODE=two|wrap` + `charge.csv` (∫n dV per write step; closes the
  diagnosis Part-IV gap). Built with inq-run (INQ_SOURCE=inq-study). 10-step wrap smoke
  PASSED on GPU 1 (`results/smoke_wrap`): N_total = 53.000 in-propagator (52 slab +
  1 WP → WP absorption IS measurable); stale pre-propagator step-0 row = 52.0 (WP not
  yet in density — known inqkit caveat, kept as diagnostic; run_metrics dedups
  keep-last).
- Locked harness in `scripts/cap_fix/`: `autoresearch.sh` (fixed 700-step witness
  benchmark, direct-binary exec), `run_metrics.py` (METRIC lines), `checks.sh`
  (absorbed_e ≥ 0.5 correctness gate). Primary metric `artifact_rise_eV` =
  E_total(final)−E_total(min), lower better; noise floor 0.1 eV (3× worst clean-run
  rise 0.033 eV); baseline signal 23.49 eV. Sanity gate passed from existing diagnosis
  data (known-bad 23.49 / known-clean ≤0.033).
- Tests recorded in `docs/validation/test-catalogue.md`: wrap-profile unit test PASS
  (`hypotheses/cap_fix/tests/test_wrap_profile.py`), smoke PASS. Follow-up noted: a
  wrapper-suite Catch2 engine test for absorbing_wrap is not yet written.
- Campaign doc + INDEX regenerated (34 campaigns).

**FAILED ATTEMPT (matters):** first launch ran runs 1+2 as concurrent `inq-run` in the
same build dir → cmake/nvlink race corrupted run 1's device link (run 2 won and
survived); the `| tail` pipe masked the non-zero exit. Fix: autoresearch.sh now
requires the pre-built `./run` (stale-check vs run.cpp) and executes it directly — no
per-run cmake; logs to `results/<name>.log`, no pipes. Harness locked after this fix.

**IN FLIGHT (ETA ~15:30–15:40 local):**
- run01_baseline_two_eta0p2 — GPU 0, two-sided η=−0.2, 700 steps. Reproduction gate:
  expect rise ≈ +23.5 eV on the new binary.
- run02_wrap_eta0p2_w30 — GPU 1, wrap topology equal-integral twin (user arm 1).

**NEXT (queued, one atomic change each):** run03 two-sided η=−2.0 (strong-CAP arm);
run04 two-sided η=−0.2 centre 35 width 10 (pushed-out footprint arm); then improve on
the best draft; finally re-init a confirmation segment (η=−1-strength, 950 steps,
phase-0 regime) for any config with rise ≤ 0.1 eV. Budget: maxRuns 24, ~55 min/run,
2 GPUs. On completion of each run: run checks.sh, append JSONL via ar-log.sh, update
dashboard + worklog, commit scoped campaign files only (no `git add -A` — repo carries
unrelated uncommitted work; commit-message rules apply).

**Mechanism predictions on record (test, don't assume):** topology alone does NOT kill
the rise (needs W·slow-density overlap either way); strong η and pushed-out footprint
do; if ALL setup arms fail → the fix is the absorbed-energy ledger accumulator
(separate campaign; valuable falsification).

**GPU state at session start:** both free (25.0/25.2 GB via cudaMemGetInfo; NVML
mismatch harmless as usual). No other users' jobs touched.
