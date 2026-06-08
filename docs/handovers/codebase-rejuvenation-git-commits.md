# Handover: Codebase rejuvenation — structured git commits

## Current status
36 new commits made on `report1/submission-package` (now 50 ahead of
`origin/main`). Working tree clean except 2 deliberately-excluded literature
files. **Not yet pushed** — awaiting user approval at the push gate and a
decision on rewording 2 commit-message filename references.

## What changed
- Codified commit conventions in `.claude/rules/commit-messages.md` (rules 4–8):
  `action(scope): description`, 9-word action taxonomy (rename>cut>sim>docs>
  fix>feature>refactor>add>chore, first-match), component scope, body policy.
- `.gitignore`: ignore `*.stdout`/`*.stderr`/`run.startup`/`*.log.mpi`,
  office/archive binaries, LaTeX junk, `Misc/`, `CLAUDE-backup.md`, `todo.txt`,
  `scripts/sync_to_laptop.sh`, `docs/{reports,presentations,sources,inq-docs}/`,
  smoke binaries; **re-include only `.claude/rules/` + `.claude/skills/`**.
- Untracked (kept on disk) the 4 `docs/reports/` notes.
- Relocations: dispatchers + `run_all_wp_rt.sh` → `jellium/scripts/`; diagram
  tooling → new `docs/diagrams/`.
- Committed library code: inqkit (density.hpp + 2 WP-stats headers); inqview
  (8 new postprocess modules, FFT transient cutoff, occupations-path fix,
  paraview, email.py, 52-module report1 figure package).
- Committed ~95 jellium run defs in 8 campaign buckets + shared configs +
  save_gs + legacy + cross-run analysis + hypotheses + scripts tooling.
- Committed coronene run defs + configs; docs (handovers, journals, plans,
  prompts, runs catalogue, literature .md notes).

## Files touched
- `.gitignore`, `.claude/rules/commit-messages.md`, `.claude/skills/**` (now tracked).
- `docs/plans/codebase-rejuvenation-git-commits.md` (the plan).
- See `git log --oneline origin/main..HEAD` for the full set.

## Commands run
- `git add` / `git commit` per the plan's 36-commit sequence.
- `git rm --cached` ×4 for the report notes.
- `git check-ignore` safety checks for `.claude/.credentials.json` etc.
- venv python import smoke tests; `pytest test_lindhard.py`.

## Tests and validation
- inqview: all 15 core/postprocess modules import OK; all 51 report1 figure
  modules import OK (0 failures); `test_lindhard.py` → **7 passed, 5 xfailed,
  1 xpassed** (only a `np.trapz` DeprecationWarning).
- inqkit C++ (density.hpp + 2 headers): **NOT rebuilt this pass** — relied on
  prior run usage; formal `inq-run` rebuild deferred (approved "record
  unverified" gate).

## Trusted sources used
- Project rules (`.claude/rules/`), CLAUDE.md, jellium base-run spec, memory.

## Attribution notes
- Commit messages avoid `claude`/`anthropic`/`ai`. **Open:** 2 bodies reference
  filenames `CLAUDE.md` / `CLAUDE-backup.md` — flagged to user for reword vs
  accept.
- `.claude/` folder name + `*-claude.md` literature filenames committed as-is
  per explicit user decision.

## Known issues / blockers
- **SAFETY (resolved):** repo `.claude/` is the full Claude Code runtime dir
  (470 MB, incl. `.credentials.json`, history, plugins, project memory). Only
  `rules/` + `skills/` are tracked; everything else verified ignored.
- Push not yet done (step 28). Requires user go-ahead.

## Assumptions still in play
- `*.html` global ignore is harmless (no authored html in scope).
- inqkit headers compile (verified-by-prior-use, not re-built).

## Milestone: sub-repo .gitignore alignment (done)
Aligned `Tutorial/` and `QuantumKickExtension/` (separate git repos) to the
main repo's .gitignore format:
- Appended the missing categories (extended sim outputs, checkpoints/, media,
  run artifacts, office/archive, LaTeX) — one block per repo, preserving their
  per-repo exe conventions (`**/run`, `**/li_kick_*`).
- `git rm --cached` the already-tracked runtime artifacts: Tutorial 25
  (24 profile.dat + 1 png), QKE 33 (18 profile.dat + 15 png). Files kept on
  disk. Commits: Tutorial `951fad2`, QKE `15a312d`.
- Validated: 0 tracked artifacts remain; new `*.png`/`profile.dat` caught; no
  uncommitted changes. Decision: untrack PNGs too (full consistency).
- Sub-repos not pushed (their remotes/auth not addressed this session).

## Main-repo validation (passed)
51 commits ahead of origin/main; no uncommitted/staged changes; only the 2
excluded literature non-md files remain untracked; `.credentials.json` ignored;
no compiled binaries tracked.

## Exact next steps
1. Complete the main-repo push: `git push -u origin report1/submission-package`
   (user is running it interactively — needs the ssh key passphrase / a key
   registered on the GitHub account).
2. If desired, push the two sub-repos to their own remotes.
3. Accepted: the 2 `CLAUDE.md`/`CLAUDE-backup.md` body references stay as
   factual filename refs (no reword).
