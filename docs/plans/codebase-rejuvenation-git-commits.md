# Plan: Codebase rejuvenation — structured git commits

Source prompt: `docs/prompts/codebase_rejuvination/task_git_commits.md`
Branch: `report1/submission-package` → push to `origin` (sets upstream).
Convention: codified in `.claude/rules/commit-messages.md` (rules 4–8).

## Scope (locked)

Commit the **existing backlog** on this branch cleanly — no new tests/CI-CD
this pass (that is a separate rejuvenation task). The backlog = 15 modified
tracked files + ~45 jellium/coronene run dirs + library code + selected docs.

## Disposition decisions (locked)

### Gitignore (add to `.gitignore`)
- **Run logs / runtime artifacts:** `*.stdout`, `*.stderr`, `profile.dat`,
  `run.startup`, `*.log.mpi` (`*.log` already ignored).
- **Binaries (keep git pure text):** `*.pptx`, `*.xlsx`, `*.zip`, `*.epub`,
  `*.html`; LaTeX build junk `*.aux`, `*.toc`, `*.spl`, `*.wcdetail`.
- **Scratch / machine-local:** `Misc/`, `CLAUDE-backup.md`, `todo.txt`,
  `scripts/sync_to_laptop.sh`.
- **Large doc bodies kept local this pass:** `docs/reports/`,
  `docs/presentations/`, `docs/sources/`, `docs/inq-docs/`.

### Untrack (`git rm --cached`, keep on disk)
Because `docs/reports/` is now ignored, untrack the 4 already-tracked files:
- `docs/reports/report1/drafts/draft5/panels_plan.md`
- `docs/reports/plasmon-detection-verdict.md`
- `docs/reports/plasmon-vs-wrap-verdict.md`
- `docs/reports/qball-spectra-comparison.md`

`docs/sources/` is ignored but its **9 existing tracked notes stay tracked**
(future source notes must be force-added). Flagged: this is in slight tension
with `scientific-grounding.md`; user accepted.

### Relocations (move, then commit at new path)
- `ResearchProject/dispatch_v2_runs.py`,
  `ResearchProject/dispatch_additional_sims.py`,
  `Misc/run_all_wp_rt.sh` → `ResearchProject/systems/jellium/scripts/`
- `Misc/build_contribution_page.py`, `Misc/render_contribution_png.py`,
  `Misc/INQ-flow-chart.drawio` → `docs/diagrams/` (new tracked folder; the
  originally-planned `docs/reports/report1/figures/` is now ignored)

### Run-dir provenance policy
Commit `run.cpp` + all analysis `*.py` + run configs. Exclude logs
(`*.stdout`/`*.stderr`) — covered by gitignore. No `REPORT.md`/`run_summary.txt`
present in the new dirs. Heavy data (`*.raw`/`*.vti`/`*.npy`/`results/`/
`checkpoints/`/binaries) already ignored.

## Validation gate (locked)
**Light smoke + record unverified.** Per code commit: `python -c import` the
touched Python; confirm C++ still builds via one representative run; record
verified-vs-unverified in the handover. No full P4 formula audit this pass.

## Commit sequence (each lockable by user before committing)

**Infra**
1. `chore(repo)`: gitignore run logs, binaries, scratch + machine-local dirs
2. `chore(repo)`: untrack report drafts kept local (`git rm --cached` ×4)
3. `chore(repo)`: document commit-message conventions *(internal config: `.claude/rules/`)*
4. `chore(repo)`: track inq-run build wrapper + env config (`shared/bin/inq-run`, `shared/config.sh`)

**Relocations**
5. `chore(jellium)`: add run dispatchers + launcher under `scripts/`
6. `docs(diagrams)`: contribution architecture diagram tooling + source

**Library code** *(smoke-checked)*
7. `feature(inqkit)`: density.hpp helpers
8. `feature(inqview)`: thread `t_start_au` transient cutoff through FFT pipeline
9. `fix(inqview)`: write occupations.csv under `raw/observables/`
10. `feature(inqview)`: paraview rendering adjustments
11. `feature(jellium)`: run_template.hpp updates
12. `feature(coronene)`: run_template.hpp updates

**Simulation provenance — jellium (8 campaign buckets)**
13. `sim(jellium)`: base run defs
14. `sim(jellium)`: L50 energy-sweep run defs
15. `sim(jellium)`: high-density L30 run defs
16. `sim(jellium)`: very-high-velocity run defs *(was "knudsen", E700–E1100)*
17. `sim(jellium)`: plasmon-probe run defs
18. `sim(jellium)`: free-wp reference run defs
19. `sim(jellium)`: sigma/tilt/energy-variant run defs
20. `sim(jellium)`: bath-only (`_wf`) reanalysis run defs
21. `sim(jellium)`: cross-run analysis + hypotheses registry (`_compare_*`, `_final_rollup`, `hypotheses/`)

**Simulation provenance — coronene**
22. `sim(coronene)`: run defs + configurations (run dirs + `configurations/` + `shared/configs/*.hpp`)

**Docs**
23. `docs(qke)`: handover + journal updates (li_extensive_kick, qke journal + index)
24. `docs(researchproject)`: journal index update
25. `docs(plans+prompts)`: task plans + rejuvenation prompts
26. `docs(repo)`: TDDFT run catalogue + observables notes (`docs/runs_catalogue.csv`, `docs/observables/`)
27. `docs(literature)`: source notes (`ResearchProject/literature/**/*.md`, 5 files)

**Publish**
28. `git push -u origin report1/submission-package` (also publishes the 14
    prior unpushed commits). No PR.

## Commit-message conventions (summary; full rule in `.claude/rules/commit-messages.md`)
- Format: `action(scope): description` — lowercase imperative, no trailing
  period, subject ≤ 72.
- Action words (classify by FIRST match, top→bottom): `rename` > `cut` >
  `sim` > `docs` > `fix` > `feature` > `refactor` > `add` > `chore`.
- Scope = component: `inqview`, `inqkit`, `jellium`, `coronene`, `qke`,
  `repo`; multi-scope `a+b` allowed.
- Body (`-` bullets, ~72 wrap) required when commit spans >1 file or carries
  physics/run provenance; record concrete provenance, never invented values.
- Forbidden words `claude`/`anthropic`/`ai`; no co-author/generated trailers.

## Known issues / assumptions
- `docs/sources/` ignore vs `scientific-grounding.md`: accepted; future
  source notes need manual force-add.
- `*.html` global ignore assumed harmless (no authored html in scope).
- Action word for the diagram tooling (commit 6) provisionally `docs`; adjust
  if user prefers `add`/`chore`.
- The two `dispatch_*.py` were never tracked, so git records commit 5 as new
  files at the new path (not a rename in history); semantically a relocation.

## Exact next steps
1. User locks this plan.
2. Apply `.gitignore` edits + `git rm --cached` untracks + relocations.
3. Execute commits 1→27, presenting each commit's file list + drafted message
   for user lock; smoke-check code commits.
4. Push (step 28). Update handover.
