# Plan: Report 1 submission package

A read-only, examiner-facing zip bundling the code and configuration files
that produced Report 1. Built on branch `report1/submission-package`.

## Goal

Produce `docs/reports/report1/code/` — a self-contained, human-readable
folder mirroring the source repo's layout — and zip it as the report's
supporting-code deliverable.

## Examiner profile (locked)

Domain-literate physicist (DFT/TDDFT familiar). Will **read** all
code to verify methodology. **Not** expected to rebuild or re-run any
simulation. Therefore the package is optimised for readability, not
reproducibility.

## Locked decisions (from grill session 2026-05-28)

| # | Decision |
|---|---|
| 1 | Examiner reads everything; does not re-run. |
| 2 | Scope: all 3 coronene run dirs (`run_save_gs_paper_replica`, `run_propagate_paper_replica`, `run_cc_bond`) + every jellium run dir that fed any report figure (incl. full σ-sweep behind A11 master plot). |
| 3 | Per-run **Lean tier**: `run.cpp`, the Cfg header it includes, `shared/cpp/run_template.hpp`, `shared/configs/boundary_rule.hpp`, `analyse.py`, `run_summary.txt`, `REPORT.md`. **Excluded:** CSV outputs, VTI series, GS checkpoints. |
| 4 | Staging path: `docs/reports/report1/code/` (mirror-repo layout inside). |
| 5 | Library payload: full `inq-stack/` (inqkit + inqview, incl. `report1/` figure scripts). |
| 6 | Draft5 figure scripts copied **into** the staging copy of `inq-stack/python/inqview/report1/` (not into the live source). |
| 7 | `docs/` excluded from the zip (handovers, journals, sources, plans, reports/* except the staging folder itself). |
| 8 | Comments: **L2** — top-of-file header (3–6 lines) + one-liner above every public class/struct/function. No parameter docstrings, no inline `// why` comments. |
| 9 | Architecture explanation lives in the single README, §5, with one ASCII data-flow diagram. |
| 10 | README outline: §1 What this is · §2 How to read · §3 Use of generative AI · §4 Project structure · §5 Library architecture · §6 Anatomy of a run dir · §7 Inventory of runs · §8 Provenance. |
| 11 | Worked-example run for README §6: `run_wp_n162_L50_E25_sigma1_v2`. |

## Workflow

1. ✅ Create branch `report1/submission-package` (no snapshot commit; packaging touches only new paths under `docs/reports/report1/code/`).
2. ✅ Write this plan + handover.
3. Build file manifest (`docs/plans/report1-submission-package-manifest.md`) listing every file to copy.
4. **Checkpoint: user audits manifest.**
5. Assemble staging folder by copying files per the manifest. No edits yet.
6. Copy draft5 figure scripts into staging `inqview/report1/`.
7. Apply L2 comments to inqkit C++ headers.
8. Apply L2 comments to inqview Python modules.
9. Write `README.md`.
10. Zip, record SHA256, final commit.

## What is deliberately omitted

- `inq/` (upstream library, gitignored, examiner clones separately if rebuilding).
- `venv/`, `ParaView-*/`, `Tutorial/`, `QuantumKickExtension/`.
- `docs/` (except the staging folder being built inside it).
- CSV outputs, VTI series, GS checkpoints inside each run dir.
- Legacy / exploratory run dirs not cited in the report (`legacy_jellium/`, `hypotheses/`, `_compare_*`, knudsen-sweep experimental dirs).

## AI-use disclosure (for README §3)

Generative AI used for:
1. Debugging C++/Python code, especially GPU-related (CUDA, INQ field handling).
2. Adding comments to existing code.
3. Styling and presentation of plots.
4. Writing bash and orchestration scripts that chained simulations.

All scientific decisions, parameter choices, interpretation, and validation
made by the author.
