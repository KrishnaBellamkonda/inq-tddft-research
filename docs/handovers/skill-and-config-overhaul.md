# Handover: skill and config overhaul

---

## Milestone: 2026-05-01 — Skills/rules migrated, env hard-coded, permission allowlist tightened

### Current status

Skills and rules infrastructure is now first-class:
- 5 legacy flat `.md` notes in `.claude/skills/` migrated to dispatcher-invocable `.claude/skills/<slug>/SKILL.md` files with YAML frontmatter and refreshed paths.
- 1 new skill added: `physics-correctness/SKILL.md`.
- All 6 skills are discoverable by the `Skill` tool (verified via system-reminder skill list after each write).
- 4 rules path-refreshed to current `inq-stack/` topology and to flag `Tutorial/` + `QuantumKickExtension/` as separate-repo directories.
- Pre-edit copies of all skills + rules archived to `docs/claude/skills/` and `docs/claude/rules/`.
- `.claude/settings.json` now hard-codes `INQ_SHARE_PATH`, `PSEUDOPOD_SHARE_PATH`, and a `PATH` that includes `shared/bin`, `inq/install/bin`, and the pyenv shims (so `inq-run`, `clang`, `clangd`, `pyright`, `python` all resolve in non-login Bash tool calls without `source ~/.bashrc`).
- `.claude/settings.json` permissions allowlist gained `Bash(nvidia-smi *)`. Most other high-frequency commands (grep/tail/ls/wc/head/find/ps/git status/log/diff/etc.) are already auto-allowed by the harness.

The repo restructure planned in `docs/plans/in-this-chat-we-binary-hearth.md` (split `Tutorial/` and `QuantumKickExtension/` into their own gits on branch `fixes/project-restructuring`) is **NOT YET STARTED**. Working tree is dirty (137 changes, mostly staged deletions of the old `ResearchProject/jellium/...` flat layout that pre-dates the move into `ResearchProject/systems/jellium/`). Restructure must wait until those changes are resolved.

### What changed

- `.claude/skills/build-run/SKILL.md` (new) — `inq-run` workflow + bashrc/env note
- `.claude/skills/handover-update/SKILL.md` (new) — handover template
- `.claude/skills/literature-review/SKILL.md` (new) — source-grounding protocol
- `.claude/skills/simulation-validation/SKILL.md` (new) — Tier A/B/C validation menu
- `.claude/skills/report-writing/SKILL.md` (new) — IMRaD + attribution rules
- `.claude/skills/physics-correctness/SKILL.md` (new) — write→known-case-test→fix→confirm gate
- `.claude/rules/file-placement.md` — added `QuantumKickExtension/` row, dropped stale `ResearchProject/jellium/<N_task>/` path, added rule 4 about separate-repo dirs
- `.claude/rules/testing.md` — `Apply to:` line refreshed (`inq/src/` → `inq-stack/`, added `QuantumKickExtension/`)
- `.claude/rules/scientific-grounding.md` — same refresh
- `.claude/rules/development-feedback-loop.md` — same refresh
- `.claude/settings.json` — added `env` block + `permissions.allow`
- `docs/claude/skills/*.md` (new copies) — archive of pre-migration flat skills
- `docs/claude/rules/*.md` (new copies) — pre-edit archive of all six rules
- `docs/plans/in-this-chat-we-binary-hearth.md` (new) — approved plan for this session
- `docs/handovers/claude-turbocharging.md` (this file)

### Files touched

All paths absolute under `/local/data/public/skcb2/tddft/.claude/` and `/local/data/public/skcb2/tddft/docs/` — see "What changed" above.

### Commands run

```bash
bash -lc 'which inq-run clang clangd pyright pvpython python'   # path discovery
mkdir -p docs/claude/skills docs/claude/rules
cp .claude/skills/*.md docs/claude/skills/
cp .claude/rules/*.md docs/claude/rules/
rm .claude/skills/*.md   # flat files removed; replaced by <slug>/SKILL.md
git status --short        # diagnosed 137 uncommitted changes
```

No destructive git ops were run.

### Tests and validation

- Proposed: skill-discovery check via post-write system-reminder skill list; env resolution via `which`; outdated-path scan via `grep`.
- Approved: implicitly via auto-mode.
- Run:
  - **Skill discovery** — after each `Write` of a `SKILL.md`, the system-reminder skill list confirmed the new skill appeared with its trigger description. ✓ All 6 skills (`build-run`, `handover-update`, `literature-review`, `simulation-validation`, `report-writing`, `physics-correctness`) are live.
  - **Env resolution (login shell)** — `inq-run`, `clang`, `clangd`, `pyright`, `python` all resolve; `INQ_SHARE_PATH` and `PSEUDOPOD_SHARE_PATH` are set. ✓
  - **Env resolution (non-login Bash tool)** — same five tools resolve and INQ vars are set even without `source ~/.bashrc`. ✓
  - **Outdated-path scan** — `grep -nH 'inq/src\|ResearchProject/jellium\|Tutorial/' .claude/rules/*.md` returns only intentional references after the edits. ✓ (Verified for the four edited rule files; remaining matches are in the path-table in `file-placement.md` and are correct.)
- Outcomes: all pass.
- Remaining gaps:
  - Skill *triggering* (not just discovery) hasn't been smoke-tested with a real prompt — defer to the next session.
  - Settings.json env block hasn't been verified after a `/reload-plugins` / fresh session — same.

### Trusted sources used

None (this session was infrastructure-only; no scientific claims made).

### Attribution notes

- Skill content largely preserved verbatim from the pre-migration archives in `docs/claude/skills/`. Diffs limited to (a) YAML frontmatter, (b) `inq/src/` → `inq-stack/include/inqkit/...` path refresh, (c) absolute paths for handover/source/report directories, (d) explicit cross-references to other skills.

### Known issues / blockers

- **Working tree dirty** (137 uncommitted changes). Most are staged deletions of `ResearchProject/jellium/01_ground_state/...`, `ResearchProject/jellium/jellium-wp-rt/...`, and `Tutorial/coronene-leed/run_*` files — pre-existing churn from the jellium relocation work that pre-dates this session. **Must be resolved (committed or stashed) before starting the `fixes/project-restructuring` branch**, otherwise the subtree-split for `Tutorial/` will inherit unrelated noise.
- The settings.json env-var hard-coding has not yet been validated after a session restart. If on next session start `inq-run` fails to resolve, fall back to `bash -lc 'inq-run …'`.
- The `docs/claude/rules/` archive currently mirrors all 6 rules even though only 4 were actually edited. Harmless but slightly over-archived.

### Assumptions still in play

- The user's existing `superpowers`, `context7`, `explanatory-output-style`, `skill-creator` plugins remain enabled — verified in `settings.json`.
- `inq-stack/` is now the canonical home for inqkit C++ headers and inqview Python (per CLAUDE.md and confirmed by directory listing).
- The jellium restructure (mirroring the coronene layout under `ResearchProject/systems/jellium/`) will happen in a separate conversation. Rules and skills already accommodate the `ResearchProject/systems/<material>/<task>/` convention.
- `Tutorial/` and `QuantumKickExtension/` will be split into independent git repos on `fixes/project-restructuring`. The rule text in `file-placement.md` already says so, but the actual split has not happened yet.

### Exact next steps

1. Resolve the 137 uncommitted changes — either commit on `main` (with a clean message describing the jellium relocation + this session's infrastructure changes), or split into two commits: one for the jellium deletions, one for the `.claude/` + `docs/claude/` + `docs/plans/` + `docs/handovers/` additions.
2. After the working tree is clean, in a fresh Claude session, execute the `fixes/project-restructuring` branch plan from `docs/plans/in-this-chat-we-binary-hearth.md` with explicit per-step confirmation:
   - Tarball backups of `Tutorial/` and `QuantumKickExtension/` first.
   - `git subtree split --prefix=<dir> -b extract/<dir>` for each.
   - Materialise standalone repos in `../<dir>-repo/` and verify commit-count parity.
   - `git rm -r --cached Tutorial/ QuantumKickExtension/`, add to `.gitignore`, commit on the branch.
   - Replace the working-tree contents of `Tutorial/` and `QuantumKickExtension/` with the standalone repos (clone or in-place re-init).
3. Smoke-test that the dispatcher fires the `physics-correctness` skill on a "the calculation is correct" type prompt.
4. (Optional) Run `/reload-plugins` and re-verify the skill list after a fresh session start so the env-block change is exercised.
5. The user reverted the `env` and `permissions` blocks in `.claude/settings.json` after this session — the file is back to plugins/effort/theme only. If env hard-coding is wanted again, re-add the block and confirm with the user.
