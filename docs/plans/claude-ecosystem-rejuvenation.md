# Plan: claude ecosystem rejuvenation

Task spec: `docs/prompts/codebase_rejuvination/task_calude_ecosystem.md`
Method: grill-with-docs (interview → lock → implement → evaluate). Only the
**user** marks a decision `LOCKED`; everything else is `under-review`.

Mandate ordering: **evaluations before implementation** (subtask 1 → 2 → 3 → 4).

---

## 1. Ecosystem audit (the "general picture") — substrate, accepted

Current `.claude/` ecosystem: **9 always-on rules, 12 skills (+2 symlinks),
0 hooks, 0 subagents**, 6 enabled plugins (outer ring), `CLAUDE.md` at root.

### Layer/overlap findings

- **Spine:** `tddft-simulations` (752 lines, 8 phases) is the hub workflow; it
  calls into 6 satellite skills (build-run, simulation-validation,
  physics-correctness, tddft-run-catalogue, journal-writing, report-figures).
- **Cluster T (testing/validation, 4-way overlap):** `testing` (rule) +
  `development-feedback-loop` (rule) + `physics-correctness` (skill) +
  `simulation-validation` (skill) all restate "test before done."
- **Policy↔procedure pairs (rule + skill restate each other):**
  handovers↔handover-update, journal-entries↔journal-writing,
  scientific-grounding↔literature-review. Required-section/format lists
  duplicated across each pair.
- **Deterministic logic encoded as prose rules** (layer-mismatch, not just
  scatter): `commit-messages` (forbidden words / `action(scope):` / ≤72 char),
  `file-placement` (path allowlist), `tddft-run-catalogue` (`scan_runs.py --run`
  upsert), ADR-0006 observable-manifest validation. All hook candidates.
- **Context-cost flag:** `jellium-base-run-spec` is 4 KB of jellium-only physics
  injected into *every* session (even report-writing).
- **Stale content:** `build-run` claims env vars are hard-coded in
  `.claude/settings.json` — they are **not** (no `env` block exists there).

### Gaps (accepted)

- **Subagents already designed in `CONTEXT.md` but unbuilt:** formula-validation
  agent + test-validation agent (fresh-context, independence-enforcing).
  `agents/` is empty.
- **No hooks at all** (`hooks/` empty, no `hooks` block in settings).
- **No evaluation harness** for any skill/hook/rule (the point of this task).

---

## 2. Locked design decisions

### D1 — Trigger-routing template `LOCKED` (2026-06-11)

The unit of de-duplication is the **trigger point**, not the file. For every
overlapping cluster:

1. **One always-on rule = the index of triggers.** It states the shared
   principle once and routes each trigger to its owner. Holds *no* procedure.
2. **Skills partitioned by trigger point** (not by topic). A skill's
   `description` encodes exactly one trigger; that is how auto-invocation fires.
   Two skills never overlap because they fire at different moments on different
   objects.
3. **Deterministic slivers → hooks**, firing on the tool event (commit, post-run).

Rationale: a skill auto-fires only when the action matches its `description`.
Authors hedged by copying guidance into always-on rules as a safety net — that
*is* the duplication. Making each trigger explicit and owned removes the need
for the copies.

### D2 — Cluster T resolution `LOCKED` (2026-06-11)

Three distinct triggers, each owned once:

| Trigger | Owner (target) |
|---|---|
| Finishing a code function/feature (inqkit/inqview/utility) | **code-level known-case-test skill** = `physics-correctness` ∪ `development-feedback-loop` |
| About to launch an expensive simulation run | **run-level validation skill** = `simulation-validation` (unchanged) |
| About to commit / post-run | **deterministic hook** (test-file-exists, no-NaN, manifest valid) |

Plus **one index rule** ("nothing is 'done' on compile alone; every change ships
a known-case test; expensive runs get a pre-approved validation menu") that names
the three triggers and routes to the two skills + hook. Folds the `testing` rule
and `development-feedback-loop` rule into that single index rule.

---

## 3. Deliverables

- Streamlined rules / skills / hooks / subagents (per D1).
- **Evaluation set** built *before* each component (subtask 1): programmatic /
  LLM-as-judge / human-in-the-loop; trigger test + functional test per skill.
- Eval harness location + CI wiring (workspace note: a test folder in inq-stack +
  CI/CD — to be reconciled with the just-shipped `.github/workflows/ci.yml`).
- **Ecosystem usage guide** (`docs/claude-ecosystem-guide.md`, NEW) — written
  last: good prompts for the common use cases, matched to the streamlined
  triggers ("to run a sweep, say X; to validate, say Y; to report, say Z").

---

### D3 — Target architecture (subtask-2 modularisation plan) `LOCKED` (2026-06-11)

Applying D1 to every component:

| Component | Target |
|---|---|
| `commit-messages` | rule (9 action-words + scope + two-commit hygiene) **+ hook** (forbidden words, `action(scope):` regex, ≤72 char) |
| `file-placement` | rule (dir table + naming) **+ backstop hook** (warn on out-of-allowlist / `inq/` writes) |
| `testing` + `development-feedback-loop` | folded into the Cluster-T index rule (D2) |
| `handovers` / `journal-entries` / `scientific-grounding` | **slim each rule to a thin trigger+invariant index**; the procedure/lists live only in the paired skill (handover-update / journal-writing / literature-review). Duplicated lists deleted. |
| `jellium-base-run-spec` | **DELETED** (D6' — stale: N=138/L=30 superseded by N=162). No 'base jellium' concept anymore. Replace with a thin pointer in `tddft-simulations` to live `shared/configs/` as source of truth; guardrails survive in `boundary_rule.hpp` + the run skill. |
| `context-management` | unchanged (trim CLAUDE.md overlap) |
| `physics-correctness` ∪ `development-feedback-loop` | one **code-test skill** (D2) |
| `simulation-validation` and all other skills | unchanged |
| `build-run` | **fix stale env claim** — add `env` block to `settings.json` so the claim is true |

**Hooks to build (tool-event-bound only):** commit-message validator
(PreToolUse on `git commit`), file-placement allowlist (PreToolUse on
`Write`/`Edit`). **NOT hooks:** catalogue upsert + manifest validator — their
trigger is a sim process finishing, not a Claude tool event; they stay as
skill/CI-invoked scripts.

**Subagents to build:** `formula-validation` + `test-validation` agents
(already designed in `CONTEXT.md`, `agents/` currently empty).

---

### D4 — Evaluator-type mapping `LOCKED` (2026-06-11)

| Component class | Trigger eval | Functional eval | Type |
|---|---|---|---|
| Hook (commit, file-placement) | n/a | input → accept/reject cases | **programmatic** (CI) |
| Slimmed rule | "prompt routes to right skill?" | n/a | **LLM-judge** (trigger) |
| Skill | prompt → skill activates | output meets rubric | **LLM-judge + human** |
| Subagent (formula/test-validation) | n/a | planted-bug detection | **programmatic fixtures + human** |

### D5 — Eval harness location `LOCKED` (2026-06-11), **revised D5' (2026-06-11)**

**D5' (authoritative): ALL ecosystem evals live under `.claude/evals/`** — NOT
in `inq-stack/`. Reason: **inqkit ships as production scientific code**; the
claude-tooling evals must not pollute the shippable package. `inq-stack/` keeps
only library code + its own code-unit-tests (last phase).

```
.claude/evals/
├── programmatic/   CI-runnable: hook cases, subagent planted-bug fixtures,
│                   Cluster-O drift, Cluster-R theme-import, build-run env.
│                   Small runner pip-installs inqkit/inqview to import.
├── skills/   rules/   clusters/   behavioural judge specs
```

CI gets a NEW job triggered on `.claude/evals/**`, separate from the
`inq-stack/**` jobs — shippable package and tooling evals fully decoupled.
(Supersedes D5's `inq-stack/tests/ecosystem/` location and the task workspace
note, per explicit user instruction.)

---

### TODO-driven additions (from `# TODO:` comments in the rules/skills, 2026-06-11)

- **Cluster O — minimum observable set as canonical standard** `LOCKED`
  (2026-06-11). **Canonical = `minimum_observable_set.hpp`** (executable; writes
  the manifest at run start; already self-declares as source of truth).
  De-dup: `tddft-simulations` Phase 3 stops restating Tier 1/2/3 required lists
  — it references the min-obs-set per run-type and lists only run-specific
  *optional* extras; Phase 7 required analysis becomes manifest-driven (TODO 3).
  A **deterministic drift eval** asserts `minimum-set-spec.md`, skill Phase 3,
  `validate_run`, and `scan_runs.py` flags all agree with the manifest the
  `.hpp` emits. Eval lives in `inq-stack/tests/ecosystem/`.
- **Cluster R — reporting / global figure standard** `LOCKED` (2026-06-11).
  Deeper review: the 3 skills have distinct cores (production / Tufte critique /
  prose) — keep all 3 — but `tufte` accreted project-specific config. Resolution
  (**separate by nature**):
  - **Executable facts** (units fs/eV/Bohr, column widths, cmap roles, figure
    factories) → `inqview.visualisation.style` (ADR 0004) = the **single global
    standard EVERY figure path imports**: `analyse.py`, comparison `scripts/`,
    report figures. No ad-hoc `rcParams`/`plt.style.use` anywhere.
  - **Project production-rules** (5 annotation rules: no interpretive text, no
    leader lines, annotations in whitespace, contrast labels, units on axes) →
    `report-figures`, which is **rescoped from "report PDFs" to the project's
    global figure-production standard + workflow** and owns these rules.
  - **`tufte`** reverts to timeless general principles only (sheds the widths
    table + annotation rules + units).
  - Auto-generated analysis plots MUST use the theme but **skip** the interactive
    5-phase grill.
  - **Eval (programmatic):** (a) figure-producing code imports the theme — no
    rogue `rcParams`/`plt.style.use` (grep); (b) units stated in skills == the
    theme module (drift). Lives in `inq-stack/tests/ecosystem/`.

## 4. Open decisions (under-review — next grill targets)

- **Eval-set scope & authoring** (next): build evals for the *changing*
  components first — 2 hooks, 2 subagents, merged code-test skill + Cluster-T
  index rule, slimmed-rule routing, build-run env, jellium relocation.
  Programmatic ones authored directly; skill rubrics via interview.
- Implementation order once the eval set is locked (subtask 3).

---

## 5. Status

- Subtask 1 (eval set): **COMPLETE & LOCKED.** Specs under `.claude/evals/`;
  every eval fixes its expected verdict before the component exists.
- Subtask 2 (modularisation plan): **COMPLETE** — D1–D6', Cluster O, Cluster R,
  D4, D5' locked.
- Subtask 3 (implement): **COMPLETE.** C1–C11 done on branch
  `rejuvenation/claude-ecosystem` (eval-first, each committed). 5/5 portable
  evals green + theme pytest 9/9. Branch NOT pushed/merged (user-gated).
- Subtask 4 (finetune): not started — iterate when the user runs the behavioural
  (LLM-judge) specs.

## 6. Eval-set inventory (subtask 1 deliverable)

`.claude/evals/`
- `skills/` — code-test, handover-update, journal-writing, literature-review
- `rules/` — cluster-t-index, slimmed-pairs, jellium-removal
- `clusters/` — cluster-o-min-observable-set, cluster-r-figure-standard
- `programmatic/` — commit-hook, file-placement-hook, formula-validation-agent,
  test-validation-agent, build-run-env

## 7. Suggested implementation order (subtask 3, eval-first per component)

1. commit-message hook → 2. file-placement hook → 3. formula-validation +
test-validation subagents → 4. Cluster-T merge (`code-test` skill) + index rule →
5. slim the 3 policy↔procedure rules → 6. delete jellium rule (+ config pointer) →
7. Cluster O de-dup (Phase-3 ref + drift check) → 8. Cluster R theme enforcement →
9. build-run env fix → 10. CI job for `.claude/evals/programmatic/` →
11. `docs/claude-ecosystem-guide.md` (usage guide, last).
