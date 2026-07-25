# ADR 0009 — Campaigns: prompt-file unit + agent-maintained frontmatter tracking

Date: 2026-06-22
Status: accepted (skill + 17-campaign backfill landed the same day)

## Context

Research work in this repo is driven by **campaign prompts** — `.md` files
(formerly `docs/prompts/`, now `docs/campaigns/`) that a fresh agent executes
autonomously, end-to-end, with no user in the loop. They had no shared authoring
process (a thin `template.md` plus bespoke mature examples like
`baseline_runs.md`) and, more pressingly, **no way to see at a glance which of
the ~15 in-flight campaigns were running, paused, blocked, or done.** Active work
also existed with *no* prompt file at all (the σ=0.5 S(v) sweep, the
σ-convergence sweep), so it was invisible to any portfolio view.

Two questions had to be settled before tooling could exist:

1. **What is the unit that carries a status?** A folder (`cap_in_jellium/`) or a
   single prompt file (`baseline_runs.md`)? The folders group several distinct
   hypotheses, so a folder-level status would average over unrelated work.
2. **Where does "x/N tasks complete" come from, and who maintains it?** Options
   ranged from auto-detecting completion by grepping run dirs/handovers, to a
   human-edited status field, to a machine-readable field the runner maintains.

## Decision

- **Campaign = one prompt `.md` file** under `docs/campaigns/<area>/`. The folder
  (`<area>`) is a purely organisational *theme*; the file is the unit. ("Each
  prompt is a campaign.")
- **Status lives in YAML frontmatter and is the single source of truth**, with a
  fixed schema: `id` (stable kebab slug, the index join key — never changes),
  `area`, `title`, `status`, `hypothesis`, `handover`, `tasks: [{name, done}]`,
  `blocked_reason`.
- **`status` enum:** `draft → ready → running → blocked → paused → done`.
  `blocked` = waiting on a dependency; `paused` = deliberately stopped. Both carry
  a `blocked_reason`.
- **The executing agent maintains the status** end-to-end — it flips `done` flags
  and bumps `status` as it runs, *with no user sign-off required*. The user does
  not hand-edit status.
- **The INDEX is derived, never authored.** `docs/campaigns/INDEX.md` is
  regenerated wholesale by `.claude/skills/campaigns/build_index.py` (stdlib-only,
  skill-local), scanning all frontmatter into one status-grouped table
  (running → blocked → paused → ready → draft → done) with a portfolio-count
  header and a do-not-hand-edit banner. Markdown only — run-level cataloguing is
  `tddft-run-catalogue`'s job, a finer layer.
- **The authoring process is a skill** (`campaigns`): five gated stages
  (Frame → Matrix → Research → Validate → Autonomy-ready), each composing existing
  skills by reference, ending in an autonomy-readiness checklist whose every box
  must be answerable from the prompt text alone.

## Alternatives considered

- **Campaign = folder.** Simpler listing, but a folder holds several unrelated
  hypotheses (e.g. `cap_in_jellium/` has baselines, a classical-vs-WP study, and a
  loss-function study); one status line per folder would be meaningless.
- **Auto-detect completion** by grepping run dirs / handovers. Most "magical" but
  fragile — a tracker that silently infers the wrong status is worse than none.
  Rejected in favour of an explicit, machine-readable field the runner already
  passes near (it updates the handover anyway).
- **User-maintained status.** Keeps the human in control but goes stale the moment
  an autonomous overnight run finishes unattended — defeating the purpose of a
  portfolio view of *autonomous* work.
- **CSV index** (like `runs_catalogue.csv`). Duplicates the frontmatter and the
  run-level catalogue; a regenerated markdown table is enough at the campaign tier.

## Consequences

- **Hard to reverse** once campaign files, the skill, and downstream habits assume
  the schema — hence this ADR. The `id` is deliberately immutable so renames don't
  orphan an index row.
- **`docs/prompts/` → `docs/campaigns/`** rename; live references updated, historical
  handovers/plans intentionally left pointing at the old path.
- **Backfill:** 14 existing prompts gained frontmatter and 3 retroactive files were
  written for previously-untracked work (17 campaigns; new area `jellium_stopping/`).
  Two superseded files (`codebase_rejuvination_complete.md` umbrella,
  `localised_jellium.md` stub) are excluded by carrying no `id`.
- **The autonomy-readiness stage now requires** every finished prompt to ship a
  frontmatter `tasks:` list (so `x/N` is computable) plus a compact `<preflight>`
  echo of the checklist for the executing agent to re-verify before burning GPU.
- **Skill artefacts stay skill-local** (`build_index.py` + its test inside
  `.claude/skills/campaigns/`) because the skills are intended to ship
  self-contained (memory `feedback_skills_self_contained_shippable`).
- `CONTEXT.md` gains a "Campaigns" glossary section; `file-placement.md` already
  routes campaign prompts to `docs/campaigns/`.
