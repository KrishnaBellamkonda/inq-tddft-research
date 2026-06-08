# Handover: Codebase rejuvenation — git commits + repo hygiene (session 2026-06-08)

## Current status
**COMPLETE.** All active changes across the three repos are committed and
pushed. `main` now holds the full latest state (fast-forwarded from
`report1/submission-package`). A fresh `unit-tests/inq-stack` branch was cut
off `main` for the next task (code unit testing — see
`docs/prompts/codebase_rejuvination/task_unit_testing.md`).

Push state at session end:
- main repo `main` == `origin/main` (pushed by user).
- main repo `report1/submission-package` == `origin/...` (pushed, frozen snapshot).
- `Tutorial` `main` == origin (pushed).
- `QuantumKickExtension` `features/li-extensive` == origin (pushed).

## What changed (this session, in order)

### 1. Commit-message conventions (grilled + locked, then codified)
Interviewed the user to lock a house style, written into
`.claude/rules/commit-messages.md` (rules 4–8):
- Format `action(scope): description` — lowercase imperative, no trailing
  period, subject ≤ 72.
- **9-word action taxonomy, classified by FIRST match top→bottom:**
  `rename > cut > sim > docs > fix > feature > refactor > add > chore`.
- Scope = component (`inqview`, `inqkit`, `jellium`, `coronene`, `qke`,
  `repo`); multi-scope `a+b` allowed.
- Body (`-` bullets, ~72 wrap) required when a commit spans >1 file or carries
  physics/run provenance.
- Pre-existing rules kept: forbidden words `claude`/`anthropic`/`ai` in
  messages; no co-author/generated trailers; two-commit hygiene.

### 2. `.gitignore` overhaul (main repo)
Added: `*.stdout`/`*.stderr`/`run.startup`/`*.log.mpi`; office/archive binaries
(`*.pptx/xlsx/zip/epub/html`); LaTeX junk (`*.aux/toc/spl/wcdetail`);
`Misc/`, `CLAUDE-backup.md`, `todo.txt`, `scripts/sync_to_laptop.sh`;
`docs/{reports,presentations,sources,inq-docs}/`; compiled smoke binaries
(`smoke_C1/C2/C3`, `dryrun`, `syntax_check`).
**`.claude/` selective tracking** — the `.*` rule ignores all of `.claude/`;
re-included ONLY `.claude/rules/` + `.claude/skills/`. Everything else
(`.credentials.json`, `.claude.json`, `history.jsonl`, `settings*.json`,
`projects/` memory, `plugins/`, `tasks/`, `sessions/`, caches) stays local.

### 3. Untracked report drafts (kept on disk)
`git rm --cached` ×4: `docs/reports/report1/drafts/draft5/panels_plan.md`
and the 3 `docs/reports/*-verdict/comparison.md`. (`docs/sources/`'s 9 tracked
notes were left tracked though the folder is now ignored — see Known issues.)

### 4. Relocations
- `dispatch_v2_runs.py`, `dispatch_additional_sims.py`, `Misc/run_all_wp_rt.sh`
  → `ResearchProject/systems/jellium/scripts/`.
- `Misc/{build_contribution_page.py,render_contribution_png.py,INQ-flow-chart.drawio}`
  → new tracked `docs/diagrams/`.

### 5. Library code committed (inqkit + inqview)
- inqkit: `density.hpp` (modified) + new `observables/wp_momentum_stats.hpp`,
  `wp_real_space_stats.hpp`.
- inqview: 8 new `postprocess/` modules (`lindhard`, `wake`, `knudsen_ke`,
  `spectral_weight`, `spectral_weight_full`, `kl_divergence`, `energy_balance`,
  `test_lindhard`); FFT `t_start_au` transient cutoff (5 modified modules);
  occupations-csv path fix; paraview tweaks; `email.py` (Gmail SMTP);
  **`report1/` package — 52 figure modules**.

### 6. Simulation provenance committed (run.cpp + analysis .py + configs only)
- **jellium**: ~95 run dirs in 8 campaign buckets — base · L50 energy-sweep ·
  high-density L30 · very-high-velocity (was "knudsen", E700–E1100) ·
  plasmon-probe · free-wp · sigma/tilt/variant · bath-only `_wf`. Plus shared
  per-energy configs, `save_gs/`, `legacy_jellium/`, cross-run analysis
  (`_compare_*`, `_final_rollup`, `hypotheses/`), and `scripts/` tooling.
- **coronene**: broadening + cc-bond + paper-replica run defs + configs +
  `configurations/`.

### 7. Docs committed
handovers, journals (qke + researchproject), plans, prompts,
`runs_catalogue.csv` + `observables/`, literature `.md` notes (incl. 4
`*-claude.md` filenames — user chose commit-as-is).

### 8. Sub-repo `.gitignore` alignment (separate repos)
`Tutorial/` and `QuantumKickExtension/` aligned to the main repo format:
appended the missing categories + an "INQ source" block (`inq/`,
`inq-codebase/` — QKE's 11 GB INQ copy). `git rm --cached` the already-tracked
runtime artifacts (Tutorial 25 = 24 profile.dat + 1 png; QKE 33 = 18
profile.dat + 15 png) — files kept on disk, PNGs untracked for full
consistency (user's call).

### 9. Branch consolidation
`report1/submission-package` was 45 ahead / 0 behind `main` → **`--ff-only`
merged into `main`** (no merge commit). `main` advanced to `50cf3ed`. Created
`research/next-phase` (interim) and then `unit-tests/inq-stack` off `main`.

### 10. VSCode nested-repo display fix
`inq/` is its OWN git repo (`inq/.git`, branch `master`); VSCode auto-detected
it and showed its internal changes. Added **local-only** `.vscode/settings.json`
(`git.ignoredRepositories` + `git.repositoryScanIgnoredFolders` = `inq`).
NOTE: requires a VSCode "Reload Window" to take effect.

## Files touched
- `.gitignore`, `.claude/rules/commit-messages.md`, `.claude/skills/**` (now tracked).
- `docs/plans/codebase-rejuvenation-git-commits.md` (the plan).
- `docs/handovers/codebase-rejuvenation-git-commits.md` (this file).
- `docs/diagrams/**` (relocated), `.vscode/settings.json` (local-only).
- `Tutorial/.gitignore`, `QuantumKickExtension/.gitignore`.
- Full set: `git log --oneline origin/main~52..HEAD` on `main`.

## Commands run
`git add`/`commit` per the 38-commit plan sequence; `git rm --cached` for
report drafts + sub-repo artifacts; `git check-ignore` safety checks;
`git merge --ff-only`; `git checkout -b`; venv import smoke tests;
`pytest test_lindhard.py`. Pushes were run interactively by the user (ssh key
passphrase — not available to the assistant's shell).

## Tests and validation
- inqview: 15 core/postprocess modules import OK; **all 51 report1 figure
  modules import OK (0 fail)**; `test_lindhard.py` → **7 passed, 5 xfailed,
  1 xpassed** (only a `np.trapz` DeprecationWarning).
- inqkit C++ (density.hpp + 2 headers): **NOT rebuilt** — relied on prior run
  usage; formal `inq-run` rebuild deferred (user-approved "record unverified").
- Repo invariants verified at session end: `main` 0/0 vs `origin/main`;
  0 files tracked under `inq/`; `.credentials.json` ignored; no compiled
  binaries tracked; sub-repos 0 tracked runtime artifacts.

## Trusted sources used
Project rules (`.claude/rules/`), `CLAUDE.md`, jellium base-run spec, auto-memory.

## Attribution notes
- All commit messages avoid `claude`/`anthropic`/`ai`. Two bodies reference
  filenames `CLAUDE.md` / `CLAUDE-backup.md` — user accepted as factual
  filename refs (not branding).
- `.claude/` folder name + `*-claude.md` literature filenames committed as-is
  per explicit user decision.

## Known issues / blockers
- **SAFETY (resolved):** repo `.claude/` is the full Claude Code runtime dir
  (470 MB incl. `.credentials.json`); a blind `git add .claude/` would have
  leaked credentials. Only `rules/` + `skills/` are tracked; verified.
- `docs/sources/` is ignored but its 9 existing notes stay tracked — slight
  tension with `scientific-grounding.md`; future source notes need force-add.
- inqkit C++ headers not formally rebuilt (verified-by-prior-use).
- Sub-repos were committed/pushed but their own CI/remotes not otherwise touched.

## Assumptions still in play
- `*.html` global ignore harmless (no authored html in scope).
- `unit-tests/inq-stack` is the intended next-task branch name (user said
  "unit-tests/"; a valid suffix was added). Rename freely.

## Exact next steps (NEW TASK: code unit testing)
Task spec: `docs/prompts/codebase_rejuvination/task_unit_testing.md`. Subtasks:
1. **Map the inq-stack library** (first pass); record restructure/reformat
   ideas in a designated file. User intends to use the `understand-anything`
   plugin (`/plugin marketplace add Lum1104/Understand-Anything`).
2. Map components needing unit tests; brief per-test plans (with the maths) +
   a simple method; user reviews and locks.
3. Write the locked unit tests; a validation agent checks them; make code
   changes as needed.
4. (per spec) restructure the library if necessary; run integration tests.
5. Add a tests/ folder under inq-stack + a CI/CD pipeline.
Start on branch `unit-tests/inq-stack` (already checked out, off latest `main`).
