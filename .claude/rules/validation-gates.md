# Rule: Validation Gates

Apply to: `inq-stack/`, `ResearchProject/`, `Tutorial/`, `QuantumKickExtension/`, `docs/plans/`, `docs/handovers/`

The always-on index for "test before done". PROCEDURES live in skills; this rule
states the principle, names the triggers, and routes them. (Replaces the former
`testing` + `development-feedback-loop` rules and the `physics-correctness`
skill — their content moved to the `code-test` and `simulation-validation`
skills.)

## Principle

Nothing is "done", "correct", "working", or "ready" on compile alone. Never
claim correctness without evidence from at least one completed validation.

## Triggers → owner

| When | Invoke |
|---|---|
| Finishing **substantive code** (new function, changed logic/numerics, new observable/kernel) | `code-test` skill — write→known-case-test→fix→confirm |
| About to launch an **expensive simulation run** | `simulation-validation` skill — Tier A/B/C menu, user-approved |
| About to **commit** code | the commit-message hook + the `code-test` pre-commit checklist |

## Boundary — exempt from "ships a known-case test"

Trivial / non-logic changes: typos, comments, docs, pure renames/moves,
formatting, config-only edits. Everything **substantive** ships a test (or a
documented `xfail`/skip with a reason) in the SAME change — record its row in
`docs/validation/test-catalogue.md`.

## Always-on invariants

- **GPU is the default.** Run INQ on GPU (`inq-run`, not `--cpu`); for Python
  prefer CUDA backends. Fall back to CPU only if explicitly requested or no GPU
  is available; report if a GPU is currently occupied by another user.
- Every substantive change includes a validation update before it is complete.
- A handover states which tests were proposed, approved, run, and their outcomes
  (`handover-update` skill).
