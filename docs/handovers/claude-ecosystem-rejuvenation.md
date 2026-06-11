# Handover: claude ecosystem rejuvenation

Task: `docs/prompts/codebase_rejuvination/task_calude_ecosystem.md`
Authoritative design: `docs/plans/claude-ecosystem-rejuvenation.md` (read it first).
Method: grill-with-docs. Only the **user** marks `LOCKED`.

## Current status (2026-06-11)

- **Subtask 1 (eval set): COMPLETE — at review gate.** 14 eval specs under
  `.claude/evals/` (skills/4, rules/3, clusters/2, programmatic/5). Anti-circular.
- **Subtask 2 (modularisation plan): COMPLETE — all clusters locked.**
- **Subtasks 3–4 (implement, finetune): ready to start**, gated on user lock of
  the eval set. Implementation order in plan §7.
- **D5' revision:** all ecosystem evals live under `.claude/evals/` (NOT
  inq-stack) — inqkit ships as production code; keep the package clean.

## What changed (this session — design only, no ecosystem files modified yet)

Locked decisions (all in the plan, with IDs):
- **D1** trigger-routing template (de-dup the principle into 1 index rule;
  partition procedures into skills by trigger; deterministic slivers → hooks).
- **D2** Cluster T (testing): 1 index rule + `code-test` skill
  (physics-correctness ∪ development-feedback-loop) + `simulation-validation`
  skill + a commit/post-run hook. Boundary LOCKED: substantive code only;
  typos/comments/docs/renames exempt.
- **D3** full target architecture (commit-messages rule+hook; file-placement
  rule+backstop hook, **dir table stays always-on**; slim handovers/journal/
  grounding rules to thin trigger+invariant indexes; build-run env fix).
- **D6'** `jellium-base-run-spec` **DELETED** (stale N=138→now N=162; no "base
  jellium" concept). Pointer to live `shared/configs/`; guardrails survive in
  `boundary_rule.hpp` + run skill.
- **D4** evaluator mapping; **D5** split harness (programmatic →
  `inq-stack/tests/ecosystem/` (CI); behavioural → `.claude/evals/`).
- **Cluster O** (min observable set): canonical = `minimum_observable_set.hpp`;
  `tddft-simulations` Phase 3 stops restating tiers; drift eval enforces
  agreement of spec/skill/validator/catalogue with the emitted manifest.
- **Cluster R** (global figure standard): `inqview.visualisation.style` (ADR
  0004) = single executable standard EVERY figure path imports; `report-figures`
  rescoped to own the **global** figure standard + the 5 annotation rules;
  `tufte` reverts to general principles. Eval: theme-import enforcement (no rogue
  rcParams) + units drift.

## Files touched (created this session — all design artefacts, nothing wired)

- `docs/plans/claude-ecosystem-rejuvenation.md` (authoritative plan)
- `.claude/evals/skills/code-test.md`
- `.claude/evals/rules/cluster-t-index.md`
- `.claude/evals/rules/slimmed-pairs.md`
- `.claude/evals/clusters/cluster-o-min-observable-set.md`
- `.claude/evals/clusters/cluster-r-figure-standard.md`

## Commands run

- Read-only surveys of `.claude/{rules,skills,settings}`; grep for `TODO` in
  rules/skills; line-level read of report-figures/tufte/report-writing.
- One `Explore` subagent summarised all 12 skill bodies.

## Tests and validation

- None executed (design phase). All eval specs assert expected verdicts BEFORE
  the components exist (anti-circularity, per task).

## Trusted sources used

- The task spec; existing ADRs 0001–0006; `CONTEXT.md`; the 4 user `# TODO:`
  comments (jellium stale, reporting-overweight, min-obs-set link, min-obs-set
  standardisation) folded in directly.

## Attribution notes

- Subagent designs (formula-validation, test-validation) originate in
  `CONTEXT.md` (prior unit-testing grill); reused, not invented here.

## Known issues / blockers

- `build-run` skill claims env vars are in `settings.json` — they are NOT (no
  `env` block). Fix scheduled in D3.
- Remaining subtask-1 evals not yet written (see next steps).

## Assumptions still in play

- The 6 enabled plugins (outer ring) are out of scope except where a local skill
  overlaps one (report-writing ↔ claude-scientific-writer — noted, not acted on).
- "Global figure standard" = visual standard via the theme module, NOT routing
  auto-plots through the interactive grill (confirmed with user).

## Exact next steps

1. Author remaining subtask-1 evals (no interview needed): commit-hook +
   file-placement-hook accept/reject cases; formula-validation +
   test-validation planted-bug fixtures; build-run env-present check; the 3
   paired-skill functional rubrics. Place programmatic ones in
   `inq-stack/tests/ecosystem/`, behavioural in `.claude/evals/`.
2. Subtask-1 review gate: user locks the complete eval set.
3. Subtask 3: implement each component eval-first (write eval → component →
   run eval), one at a time, review gate each. Suggested order: 2 hooks →
   2 subagents → Cluster-T merge + index rule → slim rules → jellium delete →
   Cluster O de-dup → Cluster R theme enforcement → build-run env fix.
4. Final deliverable: `docs/claude-ecosystem-guide.md` (usage guide / good
   prompts for common use cases), written last to match streamlined triggers.
