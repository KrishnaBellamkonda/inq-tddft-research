# Driving the claude ecosystem — a usage guide

How to get clean, reproducible behaviour out of this project's `.claude/`
tooling. Four layers, each with a different job:

| Layer | When it acts | You... |
|---|---|---|
| **Rules** (`.claude/rules/`) | always-on, every session | don't invoke — they're ambient policy |
| **Skills** (`.claude/skills/`) | auto-fire when your request matches their trigger | phrase the request so the trigger matches (below) |
| **Hooks** (`.claude/hooks/`) | automatically, on a tool event | nothing — they just run (e.g. block a bad commit) |
| **Subagents** (`.claude/agents/`) | when a task needs fresh-context isolation | ask for them by purpose |

The design principle: **a fact lives in one place; everything else cites it.**
Deterministic checks are hooks; reasoning workflows are skills; always-on policy
is a thin rule that routes to the skill.

---

## Common tasks → what to say

Skills fire on their *trigger phrasing*. Say the task plainly; these prompts hit
the right trigger.

| You want to… | Say something like… | Fires |
|---|---|---|
| Run / sweep a TDDFT simulation | "run a jellium WP sweep at E=50,100,200 eV" · "set up a coronene LEED run" | `tddft-simulations` |
| Pick the validation before an expensive run | "what should I validate before launching this run?" | `simulation-validation` |
| Confirm new code is correct | "I finished the `vec3` subscript operator — is it done?" · "ready to commit this kernel?" | `code-test` (+ `validation-gates`) |
| Build/run a `.cpp` | "build and run this run.cpp" | `build-run` |
| Make ANY figure (analysis or report) | "plot the stopping curve" · "make the Fig 3 panel" | `report-figures` (global standard) |
| Critique a chart's integrity | "is this chart honest? reduce the chartjunk" | `tufte` |
| Draft a report / caption | "draft the methods section from the handover" | `report-writing` |
| Ground a physics claim | "justify this time-step choice with sources" | `literature-review` |
| Record a run in a journal | "log run_wp_n162_L50_E100, my observation: …" | `journal-writing` |
| Catalogue runs / query observables | "which runs have a loss_function?" | `tddft-run-catalogue` |
| Write/refresh a handover | "update the handover before we stop" | `handover-update` |
| Plan / stress-test a design | "grill me on this plan" | `grill-with-docs` |
| Debug something broken | "diagnose this NaN in the propagator" | `diagnose` |

---

## Worked examples (what fires, end to end)

**1 — Finishing a piece of code**
> *"I added `Vec3::operator[]` to inqkit — is it done?"*
- `validation-gates` (rule) routes the *finishing-code* trigger → **`code-test`**
  skill fires.
- It refuses "compiles ⇒ works", asks for a known-case test with the expected
  value fixed **up front** (e.g. `v[0]==v.x`), and — because `operator[]` is not
  formula-bearing — does not need a subagent.
- When you then `git commit`, the **commit-message hook** checks the message
  format; the **file-placement hook** is silent (inqkit is a designated dir).

**2 — Running a simulation and checking what it produced**
> *"set up and run a jellium WP run at E=50 on L=50"*
- **`tddft-simulations`** fires: validates the GS, plans observables from the
  **minimum observable set** (below), writes `results/observables_manifest.json`
  at startup, runs on GPU, post-processes, and upserts the run catalogue.
- Afterwards, `inqview`'s `validate_run(run_dir)` checks the produced
  observables against the manifest (4 tiers). See the contract below.

**3 — Making a figure**
> *"plot the stopping-power comparison"*
- **`report-figures`** fires and uses `inqview.visualisation.style` — canonical
  units via `style.axis_label("stopping_power")` → `"stopping power (eV/Bohr)"`,
  the semantic cmap role, the fixed-dimension factory. No ad-hoc styling.

**4 — Checking a formula is right**
> *"independently verify the `center_of_density` formula"*
- The **`formula-validation`** subagent re-derives `∫r·n/∫n` from its source,
  blind to the code, and returns CONFIRM/FLAG.

---

## What a run produces — the observable contract (Cluster O)

A run commits, at startup, to a **minimum observable set** for its run-type and
writes it to `results/observables_manifest.json`. The canonical definition is
`inq-stack/include/inqkit/observables/minimum_observable_set.hpp` — every other
view (the `tddft-simulations` Phase-3 tables, the spec doc, the validator) cites
it.

| Run-type | Required (beyond the universal core) |
|---|---|
| **universal core** (every run) | `energy_total/kinetic/hartree/xc`, `density_l2`, `gs_eigenvalues`, `gs_occupations`, `gs_system_density` (VTI), `run_summary` |
| **jellium-wp** | `wp_momentum_stats`, `wp_real_space_stats`, `momentum_distribution`, `state_energies`, `occupations_vs_time`, `density_system_rt`, `density_total_rt` (+ optional `density_wp_rt`) |
| **jellium-classical** | `electron_track`, `state_energies`, `density_system_rt` |
| **coronene** | `wp_momentum_stats`, `wp_real_space_stats`, `leed_screen_config`, `density_total_rt` |
| **free-wp** | `wp_momentum_stats`, `wp_real_space_stats` |

**Validating a run** (post-run, in the venv):
```python
from inqview.validation import validate_run
r = validate_run("ResearchProject/systems/jellium/run_wp_n162_L50_E100")
print(r.summary())          # 4 tiers: existence · schema · finite · invariant
assert r.passed
```
Tier 4 invariants (opt-in per observable): `energy_total` drift < 1 mHa,
`density_l2(0)=0`, `wp_real_space_stats` norm ∈ [0.97, 1.03].

---

## Verifying the ecosystem works

```bash
# 1. the deterministic tier (no GPU/INQ/LLM) — the same checks CI runs:
for r in commit_hook file_placement cluster_o_drift cluster_r build_run_env; do
  python3 .claude/evals/programmatic/run_${r}_eval.py
done

# 2. trip the commit hook (should be BLOCKED):
echo '{"tool_name":"Bash","tool_input":{"command":"git commit -m \"bad message\""}}' \
  | python3 .claude/hooks/commit_message_check.py ; echo "exit=$?  (2 = blocked)"

# 3. the inqview theme contract (needs matplotlib):
venv/bin/python3 -m pytest inq-stack/tests/python/inqview/visualisation/test_theme.py -q
```
- **Skills**: phrase a task from the table above and confirm the named skill
  engages (it announces itself).
- **Subagents**: register at the next session start; ask one to check a formula.

---

## What happens automatically (hooks)

- **Commit-message validation** — a `git commit` whose message has a forbidden
  word (claude/anthropic/ai as a standalone word), a bad `action(scope):` format,
  or a >72-char subject is **blocked** with the reason. Fix the message and
  retry. (`.claude/hooks/commit_message_check.py`)
- **File-placement nudge** — writing into `inq/` (upstream) or scattering a file
  at the repo root prints a non-blocking nudge toward a designated directory.
  (`.claude/hooks/file_placement_check.py`)

You don't invoke these — they ride the `git commit` / `Write` tool events.

---

## What's always-on (rules — ambient policy)

- **validation-gates** — nothing is "done" on compile alone; substantive code
  ships a known-case test (`code-test`); expensive runs get a validation menu
  (`simulation-validation`); GPU is the default.
- **commit-messages** — the `action(scope):` format + the 9 action words.
- **file-placement** — the directory map (where each file type goes).
- **scientific-grounding** — claims cite trustworthy sources; label inferences.
- **handovers / journal-entries** — when to write them + their hard invariants.
- **context-management** — token discipline; one subagent at a time by default.

---

## Subagents (ask for them when rigour matters)

- **formula-validation** — "independently check this `center_of_density` formula
  against its source." A fresh-context agent re-derives the math, blind to the
  code, and returns CONFIRM/FLAG.
- **test-validation** — "audit this test for circularity." Checks the expected
  value is independent of the code under test, with correct units/tolerance.

Use them for formula-bearing code (COD, momentum distribution, Lindhard,
stopping) before locking the formula.

---

## Running the evals

```bash
# portable ecosystem evals (no GPU/INQ/LLM) — also run in CI:
for r in commit_hook file_placement cluster_o_drift cluster_r build_run_env; do
  python3 .claude/evals/programmatic/run_${r}_eval.py
done

# the inqview theme units test (needs matplotlib):
venv/bin/python3 -m pytest inq-stack/tests/python/inqview/visualisation/test_theme.py -q
```

Behavioural specs (skill routing + functional rubrics, LLM-judge / human) live in
`.claude/evals/{skills,rules,clusters}/`. CI runs the portable tier on every push
touching `inq-stack/**` or `.claude/{evals,hooks}/**`
(`.github/workflows/ci.yml`).

---

## If a skill doesn't fire

Skills auto-fire on their `description` trigger. If the right one doesn't engage,
name it: "use the `code-test` skill on this." If you find yourself repeating an
instruction across sessions, that's a sign it belongs in a skill/rule — say so
and it can be added (with an eval first).
