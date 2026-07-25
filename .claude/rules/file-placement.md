# Rule: File Placement

Apply to: entire project

## Designated directories

| Type of file | Directory |
|---|---|
| Task plans | `docs/plans/` |
| Session handovers | `docs/handovers/` |
| Literature notes, source summaries, citations | `docs/sources/` |
| Test matrices, benchmark definitions, validation notes | `docs/validation/` |
| Cross-**system** / manuscript-level report drafts, figure captions | `docs/reports/` |
| Run-tied analysis of a **run-set** (combined CSVs, plotting scripts, study `.ipynb`) | `ResearchProject/systems/<name>/hypotheses/<sweep_name>/` (ADR 0007 + 2026-06-15 amendment) |
| Temporary working notes | `docs/notes/` |
| INQ run machinery (build-once `run.cpp`, dispatcher, per-run `analyse.py`) | `ResearchProject/systems/<name>/scripts/` (ADR 0007) |
| Tutorial examples (separate git) | `Tutorial/<name>/` |
| QBall reference calculations (separate git) | `QuantumKickExtension/<system>/` |
| inqkit C++ headers | `inq-stack/include/inqkit/<module>/` |
| inqview Python modules | `inq-stack/python/inqview/` |
| Always-on project rules | `.claude/rules/` |
| On-demand skills and reference material | `.claude/skills/` |

## Canonical `systems/<name>/` structure (ADR 0007 + 2026-06-15 amendment)

Every production system under `ResearchProject/systems/<name>/` uses these
folders with fixed contracts. Runs are **grouped by sweep** (2026-06-15
amendment), one `<sweep_name>/` folder per run-set, matching its analysis folder
`hypotheses/<sweep_name>/` (bare names, no `NN_` prefix). jellium/coronene stay
grandfathered-flat; vacuum is migrated as the reference instance.

| Folder | Holds |
|---|---|
| `shared_gs/` | converged ground state(s), reused across runs (new unifying name; legacy systems keep `checkpoints/` / `save_gs/`) |
| `shared/` | shared config headers / `Common_`-derived cfg structs |
| `scripts/<sweep_name>/` | how runs are PRODUCED: build-once binary (`run.cpp` + `build/`), dispatcher, `gpu_probe`, per-run `analyse.py` template |
| `<sweep_name>/<run_name>/` | runs grouped by sweep, one subdir per run; outputs only (logs gitignored). Supersedes flat top-level `run_*` |
| `hypotheses/<sweep_name>/` | what a run-SET MEANS: combined CSVs, `build_*.py` scripts, study `.ipynb`, `README.md` + figures, and a `tests/` subfolder |

Two-tier tests:
- **Library-generic** feature test (new `inqkit` capability) → wrapper suite
  `inq-stack/tests/include/inqkit/<module>/` (Catch2; `_engine.cpp` when a live
  `electrons` is needed, pure host test otherwise).
- **Task-specific** implementation / mechanism check → `hypotheses/<NN_purpose>/tests/`.

`hypotheses/` is for system-local, run-tied analysis; `docs/reports/` is only for
cross-system / manuscript-level writeups. **Grandfathered:** jellium/coronene
already match this and are not migrated.

## Rules

1. Do not create files outside the directories above without proposing the location first.

2. Do not scatter notes, scratch files, or reports into arbitrary locations (e.g. project root, `inq/`, `shared/`).

3. Do not create files inside `inq/` (the unmodified upstream INQ source, gitignored) unless modifying INQ itself is explicitly requested. New library code goes in `inq-stack/include/inqkit/` (C++) or `inq-stack/python/inqview/` (Python).

4. `Tutorial/` and `QuantumKickExtension/` are tracked by their own independent git repos (see branch `fixes/project-restructuring`). The main repo at `/local/data/public/skcb2/tddft/` ignores both directories. When working in either, run `git` commands from inside that directory so the correct repo is targeted.

5. If no suitable directory exists for a new file type, propose one or two sensible locations and ask the user to choose before creating files.

6. Prefer updating an existing document over creating a duplicate.

7. Do not create `README.md` or other documentation files unless explicitly requested.

8. Always save figures as `.png`. Never save figures as `.pdf` or `.svg` unless the user explicitly requests it.

## Naming conventions

- Plans: `docs/plans/<task-name>.md`
- Handovers: `docs/handovers/<task-name>.md` (rolling file with dated milestone sections)
- Source notes: `docs/sources/<author-year-keyword>.md`
- Validation notes: `docs/validation/<system-property>.md`
