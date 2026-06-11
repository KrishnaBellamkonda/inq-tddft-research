# Eval: Cluster-T index rule (LOCKED 2026-06-11)

Target rule (subtask 3): the always-on **index rule** replacing `testing` +
`development-feedback-loop`. States the test-before-done principle once and
routes 3 triggers to their owners. Holds NO procedure.

Evaluator: **LLM-as-judge** (routing) + structural check.

## Routing test (each trigger reaches its owner)

| Scenario | Correct owner |
|---|---|
| "finished a new inqkit function — is it done?" | `code-test` skill |
| "about to launch an E=50 jellium sweep" | `simulation-validation` skill (menu) |
| "about to commit code touching observables" | commit hook + `code-test` gate |

## Boundary (LOCKED) — what is EXEMPT from "ships a known-case test"

**Substantive code only** needs a known-case test: a new function, changed
numerical/logic behaviour, a new observable/kernel. **Exempt:** typos, comments,
docs, pure renames/moves, formatting, config-only edits.

Negative eval cases (must NOT demand a test):
- "fix a comment typo in density.hpp"
- "rename a variable, no logic change"
- "update a docstring"

## Functional rubric (all hard pass/fail)

1. Each of the 3 triggers routes to the correct owner.
2. No double-firing (a sim launch must NOT also trigger the code-test skill;
   a trivial edit triggers neither).
3. Structural: the rule body contains routing + principle only — **no
   procedure** (no known-case tables, no Tier A/B/C menu; those live in the
   two skills).

Evaluator: LLM-judge on 1–2; a grep/structural check on 3 (the rule file must
not contain the per-type test table or the Tier menu text).
