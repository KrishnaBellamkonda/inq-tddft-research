# inqkit / inqview rejuvenation — running ideas

> Raw, append-only working notes gathered during subtasks 1–3 of the
> unit-testing task (`docs/prompts/codebase_rejuvination/task_unit_testing.md`).
> Subtask 4 brainstorms these and promotes accepted items into the locked
> plan at `docs/plans/inqkit-rejuvenation.md`. Nothing here is a decision.

## Status legend
`[idea]` raw · `[under-review]` discussed · `[accepted]` → goes to the plan

## Restructuring / reformatting ideas

- `[idea]` **CLAUDE.md is drifted.** The "inqkit modules" table advertises
  `jellium/analytics.hpp` (0 bytes) and a populated `core/` (pipeline,
  session_context, task — all 0 bytes). Doc describes code that does not
  exist. Decide: implement the headers, or prune the doc to reality.
- `[idea]` **11 zero-byte placeholder headers** confirmed never implemented on
  any branch (checked across all of git history): `core/pipeline`,
  `core/session_context`, `core/task`, `jellium/analytics`,
  `config/simulation_config`, `ground_state/ground_state_tasks`,
  `detail/validation`, `detail/filesystem`, `detail/text_io`,
  `io/manifest_writer`, `io/text_summary_writer`. Decide per file: implement,
  delete, or keep as an intentional placeholder with a tracking comment.
- `[idea]` **`screens/plane_screen.hpp` carries an internal TODO** (167 L,
  real but in-progress). `leed_pattern_accumulator.hpp` (120 L) also flagged
  "in-progress" in CLAUDE.md. Decide what "done" means before testing the
  unfinished path.
- `[idea]` **`inqview.postprocess.test_lindhard.py`** lives *inside* the
  package rather than under `tests/`. Relocate to the test tree during
  restructuring.
- `[idea]` **One-off scripts** (`inqview/report1/**`, `inqview/scripts/**`)
  to be relocated out of the importable package during the final restructure
  (already excluded from the testing surface).

## Findings from subtask-1 mapping (knowledge graph + signature checks)

- `[idea]` **Pure/engine split is finer than the plan assumed.** The I/O
  writers and `center_of_density` / `density_delta` operate on POD
  `RealField3D`/`ComplexField3D`, so they are `pure` (no INQ) and testable in
  the cheap CI lane — see `docs/inqkit_map.md` refinement box. Net pure
  testing surface is larger than the plan's bottom-up assumption.
- `[idea]` **WP observable overlap.** `momentum_distribution`,
  `wp_momentum_stats`, and `wp_real_space_stats` all FFT/reduce the same WP
  orbital for closely-related moments. Candidate for a shared
  moment-extraction helper during restructuring (reduce duplication).
- `[under-review]` **Writer-trio tier unconfirmed.** Confirm whether
  `observables_writer` / `occupations_writer` / `state_energy_writer` take POD
  scalars (pure) or INQ objects (engine) at their subtask-2 turn.
- `[idea]` **`detail/` is nearly empty** — only `grid_layout.hpp` is real;
  `validation`/`filesystem`/`text_io` are `// TODO`. Consider collapsing.

## Test-surfaced codebase errors
(populated during subtask 3 — see `docs/validation/` for the formal record)

- _none yet_
