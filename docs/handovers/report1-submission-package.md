# Handover: Report 1 submission package

## Current status

**Complete (uncommitted).** Branch `report1/submission-package`.

Deliverable: `docs/reports/report1/report1-submission-package.zip` (788 KB).
- SHA256: `e194801ceb7fead593d9e4cadc0db50129cdd8eb768e396d656aabf0193324d1`
- 307 files staged under `docs/reports/report1/code/` (3.0 MB unpacked).
- 29 jellium runs + 3 coronene runs (Lean tier: `run.cpp`, Cfg header,
  `analyse.py`, `results/run_summary.txt` where present).
- Full `inq-stack/` (inqkit 36 .hpp + inqview ~120 .py incl. report1 figure scripts).
- L2 comments applied across both library trees.
- README.md (8 sections per Q9 outline, AI-use disclosure in §3).

Awaiting user decision: commit on this branch + return to parent branch?

## What changed

- New branch `report1/submission-package`.
- New plan files:
  - `/local/data/public/skcb2/tddft/docs/plans/report1-submission-package.md`
  - `/local/data/public/skcb2/tddft/docs/plans/report1-submission-package-manifest.md`
  - `/local/data/public/skcb2/tddft/docs/plans/report1-submission-package-assemble.sh`
- New handover: `/local/data/public/skcb2/tddft/docs/handovers/report1-submission-package.md`
- Staging folder built at `/local/data/public/skcb2/tddft/docs/reports/report1/code/` (307 files, 3.0 MB) — pure copy from source, plus L2 comment normalisation.
- Zip: `/local/data/public/skcb2/tddft/docs/reports/report1/report1-submission-package.zip` (788 KB).

## Files touched

- Original sources at `/local/data/public/skcb2/tddft/inq-stack/` and `/local/data/public/skcb2/tddft/ResearchProject/systems/{jellium,coronene}/` — **untouched**. All comment edits applied only to the staging copies under `docs/reports/report1/code/`.

## Commands run

```
git checkout -b report1/submission-package
bash docs/plans/report1-submission-package-assemble.sh   # 307 files copied
# L2 comment passes via subagents on the staging copies (no source touched)
cd docs/reports/report1 && zip -rq report1-submission-package.zip code/ \
    -x '*/__pycache__/*' '*.pyc' '*.egg-info/*'
sha256sum report1-submission-package.zip
```

## Tests and validation

- Manifest cross-checked against `panels_plan.md` per-panel "Run:" lines and `stopping_power_data.py` aggregator.
- All 307 staged files counted; jellium = 29 runs, coronene = 3 runs, draft5 figure scripts = 21.
- Subagent comment-pass passes AST parse on every Python file (raw-string conflict in `paraview_3d.py` caught and reverted).
- One subagent attempted to invent contents for the empty 0-byte `inqkit/config/simulation_config.hpp`; reverted (file removed from staging since it has no includers).
- Zip integrity: SHA256 recorded above; 401 zip entries (307 files + dirs).

## Trusted sources used

- `/local/data/public/skcb2/tddft/docs/reports/report1/drafts/draft5/panels_plan.md` — authoritative list of which runs feed which report figures.
- Grilling session transcript with user 2026-05-28 — Q1–Q10 decisions.

## Attribution notes

Plan derived from grill-with-docs session 2026-05-28. AI-use disclosure
section to surface in the deliverable README per user request.

## Known issues / blockers

- None.
- Working tree on parent branch carries ~200 untracked dirs and ~12
  modifications — deliberately not committed (would produce an
  unreviewable snapshot commit; user confirmed option (a) — no
  snapshot, return-trip-safe because packaging touches only new paths).

## Assumptions still in play

- The "all jellium run dirs that feed any report figure" set will be
  enumerated by inspecting `panels_plan.md` §2 + Part C + per-panel
  Run: lines + `stopping_power_data.py` aggregation list. Final list
  produced in the manifest step for user audit.
- "Lean tier" includes `REPORT.md` and `run_summary.txt` only if they
  exist in the run dir — older runs may be missing them; manifest
  records actual presence.

## Exact next steps

1. User decides whether to commit on `report1/submission-package`. Recommended commit content: the plan/manifest/handover/assembly-script + the staging folder + the zip. Single commit (all research-side per the commit-message rule). Suggested subject: `report1: assemble examiner-facing submission package zip`.
2. After commit, user returns to `runs/electron-classical-wavepacket-jellium` (`git checkout runs/electron-classical-wavepacket-jellium`). All in-flight uncommitted work on that branch is preserved (it was never committed to the packaging branch; the working tree carries through).
3. If revisions are needed, edit files under `docs/reports/report1/code/` and re-run the assembly script + zip step (both are idempotent).
