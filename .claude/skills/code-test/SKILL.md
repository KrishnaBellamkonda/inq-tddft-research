---
name: code-test
description: Use when finishing — or about to call "done", "correct", "working", "ready" on — any code function, kernel, observable, utility, or data-processing change in `inq-stack/`, `ResearchProject/`, `Tutorial/`, or any utility/analysis code, and before committing such code. Enforces the write→known-case-test→fix→confirm loop and blocks "compiles ⇒ works" reasoning. (Merges the former physics-correctness skill + development-feedback-loop rule; the run-level benchmark menu is the simulation-validation skill.)
---

# Code-test gate

Wrong code → wrong results. This gate fires when code is about to be declared
done, committed, or fed into a longer program. Block on it; do not skip.

Triggered by the always-on **validation-gates** rule for the *finishing-code*
trigger. The sibling triggers: an expensive *run* → `simulation-validation`
skill; a *commit* → the commit-message hook + this gate.

## The hard rule

**Never claim code works because it compiles.** Test it with a known-answer
input first. Applies to ALL code — utility functions, data-processing scripts,
I/O helpers, callbacks, GPU kernels — not just physics. A two-line function that
computes the wrong index silently corrupts everything downstream.

## Required loop (every new function or non-trivial change)

1. **Write** — implement the function or change.
2. **Test with a known case** — at least one input with a known correct answer:
   - arithmetic/physics → analytic result, unit conversion, or published value;
   - I/O → the written file has the expected format and content;
   - GPU kernels → a trivial constant-value case before real computation;
   - data processing → a hand-crafted small input vs a manually computed output.
   The expected value is fixed **up front**, independent of the code's own output
   (never retrofit the assertion to what the code happens to produce).
3. **Fix discrepancies** — diagnose, change, re-run.
4. **Confirm inclusion** — only after the known-case test passes, integrate.
   State which test was run and what it showed (handover or comment).
5. **Repeat the loop, after the user's permission, until the goal is reached.**

## Known-case test menu

| Function type | Minimum known-case test |
|---|---|
| GPU orbital write (`inject_wp`, …) | constant 1.0 to all points; GPU reduce returns `n_pts`. |
| Coordinate computation | print one point; compare to a hand-calculated grid coordinate. |
| Norm / overlap | inject an analytical WP; norm ∈ [0.97, 1.03]. |
| 3D density save | write a 10-point toy array; read back and verify. |
| Overlap matrix | identity at t=0: diagonal ≈ 1, off-diagonal ≈ 0. |
| z-profile / density slice | integrated slice vs known total density at that plane. |
| Energy/momentum observable | vs the analytic free-particle expectation. |
| centre of density | uniform box → COD at box centre, |COD| < 1e-9. |
| SCF run | energy converged below tol, no NaN/Inf, GPU/CPU agree to 6+ sig figs. |
| RT propagation | energy drift < 0.1% over the run; t=0 dipole ≈ 0 for non-polar systems. |

## Exempt (the validation-gates boundary)

Trivial / non-logic changes need NO known-case test: typos, comments, docs,
pure renames/moves, formatting, config-only edits. The gate is for *substantive*
code: a new function, changed numerical/logic behaviour, a new observable/kernel.

## Pre-commit / pre-claim checklist

Before "this works", "correct", "fixed", "ready", or a `git commit` touching
substantive code:

- [ ] each new function has a recorded known-case test + the observed result (not "assumed correct");
- [ ] Tier-A validation status recorded for any simulation result claimed (`simulation-validation` skill);
- [ ] literature comparison cites its source (`literature-review` skill);
- [ ] GPU used (`inq-run`, not `--cpu`) unless explicitly requested or unavailable;
- [ ] the handover records what was tested and what remains unverified (`handover-update` skill);
- [ ] no correctness claim made about anything unverified.

If any box is unchecked, do not declare complete — run the missing test or tell
the user which evidence is missing.

## Recording what was tested

```
inject_wp:        GPU constant-write test, sum=1.75e6 == n_pts ✓
                  Gaussian WP norm = 1.002 ∈ [0.97, 1.03] ✓
save_orbital_3d:  not yet tested — must verify before next run.cpp launch
```
"Not yet tested" is acceptable. "Assumed correct" is not.

## Observable-producing code

Code that emits a run's primary observables must conform to the **minimum
observable set** — the canonical definition is
`inq-stack/include/inqkit/observables/minimum_observable_set.hpp` (ADR 0006,
Cluster O). Do not restate the required set here; reference it. The post-run
`validate_run` (inqview) checks a run against the manifest it emits.

## Formula-bearing code

For a function whose output is defined by a cited formula (center_of_density,
momentum_distribution, jellium shells, Lindhard, stopping), spawn the
**formula-validation** subagent (independent re-derivation) before locking the
formula, and the **test-validation** subagent to audit the test for circularity.

## Cross-references
- `validation-gates` rule — the always-on trigger index this skill serves.
- `simulation-validation` skill — the run-level tiered benchmark menu.
- `handover-update` skill — where test outcomes are recorded.
