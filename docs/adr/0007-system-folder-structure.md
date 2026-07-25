# ADR 0007 — Canonical `systems/<name>/` folder structure

Date: 2026-06-13
Status: accepted (design; vacuum is the first system retrofitted)

## Context

Production research experiments live under `ResearchProject/systems/<name>/`
(jellium, coronene, vacuum). Each system accumulates four kinds of artefact that
were previously scattered with no enforced convention:

1. a reusable converged **ground state** (re-running SCF per job is wasteful);
2. **shared config** headers / `Common_`-derived structs;
3. the **production machinery** that generates runs (a build-once binary, a
   GPU dispatcher, per-run `analyse.py`);
4. the **interpretation** of a *set* of runs (combined CSVs, aggregation/plotting
   scripts, the study notebook) plus the task-specific tests that validate that
   the run implementation is correct.

Until now (4) landed in `docs/reports/<task>/` (e.g. the absorbing-boundary
notebooks), divorced from the runs it analyses, while (3) sat in ad-hoc per-task
subfolders (`vacuum/mfa_sweep/`, `vacuum/tests/`). Cross-run analyses had no
home, so combined CSVs and study notebooks drifted apart from their data.

## Decision

Every `systems/<name>/` uses exactly five top-level folders with fixed contracts:

Names align to the **established jellium/coronene convention** (decided in the
2026-06-13 grilling session: plural `hypotheses/` with numbered subfolders, FLAT
top-level `run_*`). `shared_gs/` is the one new name (existing systems used
inconsistent `checkpoints/`+`configurations/` vs `save_gs/`).

```
systems/<name>/
├── shared_gs/          # converged ground state(s), reused across runs (new
│                       #   unifying name for the prior checkpoints/save_gs)
├── shared/             # shared config headers / Common_-derived cfg structs
├── scripts/            # how runs are PRODUCED: build-once binary (run.cpp +
│                       #   build/), dispatcher, gpu_probe, per-run analyse.py
├── run_<type>_<params>/   # FLAT top-level, one dir per run; outputs (logs gitignored)
└── hypotheses/         # PLURAL — matches jellium/coronene
    └── <NN_purpose>/   # numbered: what a run-SET MEANS (e.g. 00_mfa_reflectivity)
        ├── *_combined.csv   # aggregated cross-run data
        ├── build_*.py       # aggregation + plotting scripts
        ├── *.ipynb          # the study notebook(s) + README.md + figures
        └── tests/           # task-specific implementation / mechanism checks
```

Two-tier test rule (interacts with ADR 0001, the inqkit test harness):

- **Library-generic feature tests** (a new `inqkit` capability) → the wrapper
  suite `inq-stack/tests/include/inqkit/<module>/`, Catch2, `_engine.cpp` suffix
  when a live `electrons` is required, pure host test otherwise.
- **Task-specific implementation / mechanism checks** (validating that *this
  run-set's* physics is wired correctly) → `hypotheses/<NN_purpose>/tests/`.

`hypotheses/` is the home for **system-local, run-tied** analysis. `docs/reports/`
is retained only for **cross-system / manuscript-level** writeups not owned by a
single system.

**Grandfathering:** jellium and coronene already follow this layout (flat `run_*`,
plural `hypotheses/`) and are NOT migrated — only their ground-state folders keep
their legacy names (`checkpoints/`, `save_gs/`). The standard applies to NEW
systems; `vacuum` is rearranged to it as the reference instance.

## Alternatives considered

- **Keep notebooks in `docs/reports/`** (option B in the grilling session). Simpler
  and matches the prior `file-placement.md` rule, but splits a study across two
  trees — the combined CSVs and tests live with the system while the notebook that
  reads them lives under `docs/`, so neither half is self-contained.
- **No `hypotheses/` layer; analysis scripts in `scripts/`.** Conflates "how runs
  were produced" with "what they mean", and gives cross-run combined data no home
  distinct from the per-run `analyse.py`.

We chose system-local `hypotheses/` so a run-set's data, scripts, notebook, and
validation tests travel together and are reproducible from one folder.

## Consequences

- **Hard to reverse** once tooling assumes the layout — hence recording it as an
  ADR. Existing systems are **grandfathered** (they already use flat `run_*` +
  plural `hypotheses/`); the standard binds NEW systems, and `vacuum` is rearranged
  to it as the reference instance.
- `file-placement.md` is updated: report/notebook artefacts for a run-set point into
  `hypotheses/<NN_purpose>/`, narrowing `docs/reports/` to cross-system writeups.
- The MFA feature exposed the gap this formalises: `MaskAbsorber` shipped with only
  task-specific checks (`vacuum/tests/{gate1_mask_absorber,mask_mechanism_check}`)
  and **no** library unit test. Under the two-tier rule those checks move to
  `vacuum/hypotheses/00_mfa_reflectivity/tests/`, and a new
  `inq-stack/tests/include/inqkit/absorbers/test_mask_absorber.cpp` (pure host test
  of `sin2_mask_value`) is owed.
- CONTEXT.md gains glossary entries for the folder terms.

## Amendment (2026-06-15) — runs grouped by sweep

Superseding the original "FLAT top-level `run_*`" rule (decided in the 2026-06-15
grilling session for `docs/plans/cap-thin-absorber-tuning.md`):

- Runs are **grouped by sweep**: `systems/<name>/<sweep_name>/<run_name>/` — one
  folder per sweep holding its run subdirectories — instead of flat `run_*` at the
  system root. A "sweep" is the run-set produced together for one hypothesis.
- `hypotheses/<sweep_name>/` uses **bare sweep names** (drop the `NN_` numeric
  prefix) so a sweep's runs (`<name>/<sweep_name>/`) and its analysis
  (`hypotheses/<sweep_name>/`) share an identical token.
- Production machinery stays in `scripts/<sweep_name>/` (build-once binary +
  dispatcher), distinct from the runs it produces.

**Motivation:** with several sweeps per system (vacuum already has MFA, the CAP
knobs study, and the new `cap_thin_L5` study), flat `run_*` at the root mixes
unrelated run-sets in one listing; grouping keeps each hypothesis's runs together
and one-to-one with its analysis folder.

**Scope:** binds NEW systems and new sweeps going forward. **jellium and coronene
stay grandfathered-flat** (not migrated). **vacuum is migrated** to the grouped
layout as the reference instance: `run_cap_* → cap_real/`, `runs/run_mfa_* →
mfa_sweep/`, top-level `mfa_sweep/` machinery → `scripts/mfa_sweep/`,
`hypotheses/01_cap_real/ → hypotheses/cap_real/`.

**Consistency:** `CONTEXT.md` "System folder structure" and
`.claude/rules/file-placement.md` are updated to match.
