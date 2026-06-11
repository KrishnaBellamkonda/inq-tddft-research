# Handover: ecosystem finetune (overnight) + rejuvenation follow-ups

Import this after clearing context. Self-contained kickoff for the next
overnight run. Authoritative deep history: `docs/handovers/claude-ecosystem-
rejuvenation.md` + `docs/plans/claude-ecosystem-rejuvenation.md`. Task spec:
`docs/prompts/codebase_rejuvination/task_calude_ecosystem.md` (read subtask 4 +
the principles/rules).

## Current status (2026-06-11)

Whole codebase-rejuvenation effort (3 tasks) is functionally complete & committed
**locally**; nothing since 2026-06-08 is pushed.

- **Task 1 (git-commits): COMPLETE** (main/report1/Tutorial/QKE pushed earlier).
- **Task 2 (unit-testing/restructure): core done** (26 inqkit tests, inqview
  4-layer restructure ADR 0003, CI/CD, minimum observable set — now rolled out).
  Engine-compile of inqkit (Vec3/density/orbital) **VERIFIED** via the Tier-3 GPU
  build. **Deferred (user-gated):** restructure TODOs T08/T14/T16/T23/T24/T26/T30;
  physics fixes **E02** (cached density at t=0) + **E04** (COD half-cell bias).
- **Task 3 (claude-ecosystem): subtasks 1–3 DONE**, tested end-to-end (all 4
  test tiers PASS incl. a real GPU run → manifest → validate_run PASS).
  **Subtask 4 (finetune) NOT started — THIS is the overnight job.**

Branch: **`rejuvenation/claude-ecosystem`** (local; ~57 commits ahead of
`origin/main`; **NOT pushed** — no GitHub SSH key in this env, `git push` fails
`Permission denied (publickey)`). Do NOT attempt to push.

## The overnight task — Task-3 subtask 4: finetune via the eval set

Run the **behavioural eval specs** (LLM-as-judge + trigger tests) that were
written but never executed, then iterate the skills/rules until they pass.
Per the task spec: **each iteration runs the full eval set and presents results
to the user before any further change; changes are `under-review`, only the user
locks; incremental only.**

Behavioural specs to run (each defines trigger ± and a functional rubric):
- `.claude/evals/skills/`: `code-test`, `handover-update`, `journal-writing`,
  `literature-review`.
- `.claude/evals/rules/`: `cluster-t-index`, `slimmed-pairs`, `jellium-removal`.
- `.claude/evals/clusters/`: `cluster-o-min-observable-set` (behavioural half),
  `cluster-r-figure-standard` (behavioural half).

**How to run each (LLM-as-judge, fresh-context subagents):**
1. **Trigger test:** give a judge subagent the positive/negative prompts from the
   spec + the live skill `description`s; it returns which skill fires per prompt.
   PASS = positives fire the target skill, negatives route elsewhere.
2. **Functional test:** give a judge subagent the spec's scenario + rubric + the
   skill body; have it (a) produce the skill's expected output for the scenario,
   then (b) score it against the rubric criteria. PASS = all hard criteria met.
3. Collect pass/fail + the judge's reasoning into a results table.
4. For each FAIL, propose a **minimal** skill/rule edit; mark `under-review`;
   present before changing. Re-run the affected eval after a locked fix.

Deterministic tier is already green (re-run to confirm, fast):
`for r in commit_hook file_placement cluster_o_drift cluster_r build_run_env; do
python3 .claude/evals/programmatic/run_${r}_eval.py; done` +
`venv/bin/python3 -m pytest inq-stack/tests/python/inqview/visualisation/test_theme.py -q`.

## Optional follow-ups (fold in if scope allows; all user-gated)

- **Cluster-O spec-key alignment** — align `docs/observables/minimum-set-spec.md`
  observable keys to `minimum_observable_set.hpp` so the drift-eval bridge drops.
- **literature-review TODO** — its source strategy changed (primary repo = the
  `literature/` folder + Drive, then internet); fold into the literature-review
  skill body.
- **test-catalogue rows** — add rows for the new ecosystem evals to
  `docs/validation/test-catalogue.md`.
- **Stale doc refs** — `docs/reflection.md`, `docs/claude-instructions/…` still
  name the deleted `testing`/`physics-correctness`/`development-feedback-loop`.

## Constraints (carry forward)

- Commit messages: `action(scope): desc`, ≤72-char subject, NO claude/ai/
  anthropic as prose words, no Co-Authored-By. The **commit hook enforces this
  live** — it will block bad messages (path tokens like `docs/claude/…` are OK).
- GPU is the default (NVML/`nvidia-smi` broken but compute works — verify via a
  `cudaGetDeviceCount` probe, not nvidia-smi).
- Use `venv/bin/python3` for Python; never edit results in `runs/`.
- Don't push (no auth) and don't merge unless the user sets up GitHub access.

## Exact first steps (fresh session)

1. Read this handover + `task_calude_ecosystem.md` (subtask 4 + principles).
2. Confirm branch `rejuvenation/claude-ecosystem`; re-run the deterministic evals
   (should be green).
3. Build the eval-run plan (one judge pass per behavioural spec); present it.
4. Run the behavioural eval set; present the results table.
5. Propose `under-review` fixes for any FAIL; iterate after the user locks each.
