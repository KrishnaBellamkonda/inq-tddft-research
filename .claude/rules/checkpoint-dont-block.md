# Rule: Budget overruns — checkpoint and run, never self-block

Apply to: every run orchestrator / dispatcher (`orchestrate.py`, `dispatch.sh`,
campaign executors) and every long TDDFT run definition (`run.cpp`) under
`ResearchProject/systems/**/scripts/` and `docs/campaigns/**`. Always on.

## The one rule

**When a run (or run-set) is projected to exceed its scheduled time budget, do
NOT block or refuse to launch. Ensure good checkpointing, launch at the full
planned scope, WARN the user (email + log) with the measured projection — and
leave the kill decision to the user.** (User decision, 2026-07-12, after the
Nazarov–Gross budget gate idled both GPUs for 8 hours awaiting a human.)

1. **Checkpointing is the precondition.** Every run long enough to matter saves
   an interior RT checkpoint every ~200 steps (`electrons.save(rt_ckpt_dir)` +
   `rt_state.txt` with `last_step`/`wp_idx`/`dt`), plus a final checkpoint, and
   supports `*_RESUME=1` (load ckpt, re-apply per-state masses/perturbations,
   `real_time::propagate(..., START)`, segment-suffixed CSVs). Reference
   implementations: `scripts/muon_mass_fork/sigma1_massonly/wp/run.cpp` and
   `scripts/nazarov_gross/wp/run.cpp`.
2. **A projected overrun is a WARN, not a gate.** The orchestrator measures the
   real cost (smoke), emails the projection ("~X h per run, exceeds remaining
   Y h — proceeding, checkpointed, kill+resume instructions inside"), and
   proceeds. A killed run loses at most one checkpoint interval.
3. **Hard gates remain only for correctness**, never for cost: cutoff-guard
   BLOCK, NaN/complex energy, missing GS, failed smoke physics (energy drift).
4. **Silence is still forbidden.** The overrun warning email is mandatory —
   proceeding quietly past a budget is as bad as blocking on it.

## Why

An autonomous executor that self-blocks on a cost projection converts a soft
constraint (wall-clock preference) into a hard stop that wastes idle GPUs until
a human notices (8 h lost, 2026-07-12). With interior checkpoints the cheap
recovery direction inverts: launching too much is recoverable (kill + resume),
launching nothing is not. Cost budgets belong to the user; correctness gates
belong to the orchestrator.

## How to apply

- New run.cpp for long runs: clone the checkpoint/resume block from a reference
  implementation (state parsers, RESUME branch, ckpt-in-callback, final ckpt,
  `run_summary.txt` always final-state with `start_step`/segment fields).
- Orchestrator: replace any `if est > remaining: stop` with the WARN email +
  proceed; include the kill (`kill <pid>`) and resume (`*_RESUME=1`) recipe in
  the email body.
- Post-processing: expect segment-suffixed CSVs (`observables.from600.csv`) and
  concatenate segments; note that density_delta L2 baselines reset per segment.
