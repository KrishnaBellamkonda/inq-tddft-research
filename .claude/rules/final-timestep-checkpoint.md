# Rule: Every run checkpoints its last timestep (resumable)

Apply to: every TDDFT run definition under
`ResearchProject/systems/**/scripts/**/run.cpp` and its dispatcher/orchestrator.
Always on. Composes with `checkpoint-dont-block.md` (interior checkpoints for long
runs) — this rule is the universal FINAL-state guarantee, required even for short runs.

## The one rule

**Every run saves a FINAL checkpoint at its last timestep and supports resuming
from it, so the run can be EXTENDED to more steps later without recomputing.**
(User decision, 2026-07-14: "checkpoint the last timestep for all the runs … this
allows us to continue progress once the last timestep is reached if more time is
required.")

1. **Final checkpoint (mandatory, every run).** After `real_time::propagate`
   returns, write:
   - `electrons.save(OUT + "/checkpoint")` — the collective RT state (all orbitals,
     incl. the WP for a wavepacket run);
   - `OUT/rt_state.txt` with at least `last_step=<N>`, `time_au=<t>`, `dt=<dt>`, and
     any per-run dynamical state NOT held in `electrons` — for a moving classical
     projectile: `proj_z`, `proj_vz`, `proj_mass`, `proj_charge`; for a WP:
     `wp_idx`. (`run_summary.txt` still records final state as usual.)

2. **Resume branch (mandatory).** On `LJ_RESUME=1` (or the run's env equivalent):
   read `rt_state.txt` → `START=last_step`; if `START >= N_STEPS`, exit 0 (clean
   no-op); else `electrons.load(OUT/checkpoint)` INSTEAD of the GS, restore the
   extra dynamical state (projectile R/V; do NOT re-inject a WP — it is already in
   the checkpoint), and call `real_time::propagate(…, START)` (the trailing START
   offset) to run `START → N_STEPS`.

3. **Segment outputs.** On resume, write observables to a segment-suffixed file
   (`observables.from<START>.csv`, `projectile.from<START>.csv`, VTI overwrite=false)
   so no data is lost. Post-processing CONCATENATES `observables*.csv` /
   `projectile*.csv` in step order (the twin-run engine's `load_run` does this).

4. **Extend, don't restart.** To add steps to a completed run, re-invoke the SAME
   run with a larger `LJ_N_STEPS` and `LJ_RESUME=1` — never recompute from step 0.
   A run whose last timestep was NOT checkpointed cannot be extended (must re-run):
   that is the failure this rule prevents.

## Why

Runs routinely turn out too short (a WP still approaching the slab at 50 steps; a
spectrum needing longer time). Without a final checkpoint the only option is to
recompute from scratch. A last-timestep checkpoint makes every run a resumable
prefix, so "run more steps" costs only the new steps.

## How to apply

- New long/dynamic run.cpp: clone the checkpoint/resume block from a reference
  implementation — `scripts/nazarov_gross/wp/run.cpp` (WP; `read_last_step`,
  RESUME branch, `electrons.save/load`, `propagate(…, START)`, segment CSVs) and
  `scripts/localised_jellium_dynamics/proj_dyn/run.cpp` (moving classical
  projectile; also persists `proj_z/proj_vz` in `rt_state.txt`).
- Dispatcher: expose `*_RESUME=1` + a larger `*_N_STEPS`; leave completed segments
  in place.
- Analysis: expect segment-suffixed CSVs and concatenate by step (density_delta
  L2 baselines reset per segment — see `checkpoint-dont-block.md`).
- See the `tddft-simulations` skill (run lifecycle) which references this rule.
