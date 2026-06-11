---
name: physics-correctness
description: Use before claiming any physics-producing code or calculation is "correct", "working", "complete", or "ready" — and before committing changes that touch wavefunctions, densities, energies, forces, observables, propagators, or any GPU kernel under `inq-stack/include/inqkit/`, `ResearchProject/`, or `Tutorial/`. Enforces the write→known-case-test→fix→confirm loop and blocks "compiles ⇒ works" reasoning.
---

# TODO: I see that a lot of skills depend on this minimum observable set. So, we need to nail this down, and run a few runs to validate that what I wanted is happening. I believe, I have given information about this set in the previous tasks. But, I want to check that this is being standarding across the claude ecosystem too. 

# Physics Correctness Gate

Wrong code → wrong physics. This is the gate that fires whenever the work
is about to be declared done, or about to be committed, or about to feed
into a longer simulation. Block on it; do not skip.

## The hard rule

**Never claim a physics function works because it compiles.**
Test it with a known-answer input first.

## Required loop for every new function or non-trivial change

1. **Write** — implement the function or change.
2. **Test with a known case** — run on at least one input with a known correct answer.
   - Arithmetic / physics: compare against analytic result, unit conversion, or published value.
   - I/O: check the written file has the expected format and content.
   - GPU kernels: verify with a trivial constant-value case before testing real computation.
   - Data processing: hand-craft a small input and compare against a manually computed output.
3. **Fix discrepancies** — diagnose, change the code, re-run.
4. **Confirm inclusion** — only after the known-case test passes, integrate into the calling code. State explicitly (in handover or comment) which test was run and what it showed.

## Known-case test menu

| Function type | Minimum known-case test |
|---|---|
| GPU orbital write (`inject_wp`, etc.) | Constant 1.0 to all points; GPU reduce must return `n_pts`. |
| Coordinate computation | Print one point; compare to hand-calculated grid coordinate. |
| Norm / overlap | Inject analytical WP; check norm ∈ [0.97, 1.03]. |
| 3D density save | Write 10-point toy array; read back and verify values match. |
| Overlap matrix | Identity at t=0: diagonal ≈ 1, off-diagonal ≈ 0. |
| z-profile / density slice | Compare integrated slice against known total density at that plane. |
| Energy/momentum observable | Check against analytic free-particle expectation. |
| SCF run | Energy converged below tolerance, no NaN/Inf, GPU/CPU agreement to 6+ sig figs. |
| RT propagation | Energy drift < 0.1% over full run; t=0 dipole near zero for non-polar systems. |

## Pre-commit / pre-claim checklist

Before the words "this works", "this is correct", "fixed", "ready", or
before a `git commit` that touches physics-producing code, confirm:

- [ ] Each new function has a recorded known-case test with the observed result (not "assumed correct").
- [ ] Tier-A validation status is recorded for any simulation result being claimed (see `simulation-validation` skill).
- [ ] If the result is being compared to literature, the source is cited (see `literature-review` skill).
- [ ] GPU was used (`inq-run`, not `--cpu`) unless a CPU run was explicitly requested or no GPU is available.
- [ ] The handover for the task records what was tested and what remains unverified (see `handover-update` skill).
- [ ] No claim of correctness has been made about something that has not been verified with evidence.

If any box is unchecked, do not declare the work complete. Either run the
missing test or explicitly tell the user which evidence is missing.

## Recording what was tested

In every handover, for each new function, state:

```
inject_wp:        GPU constant-write test, sum=1.75e6 == n_pts ✓
                  Gaussian WP norm = 1.002 ∈ [0.97, 1.03] ✓
save_orbital_3d:  not yet tested — must verify before next run.cpp launch
```

"Not yet tested" is acceptable. "Assumed correct" is not.

## Cross-references

- `.claude/rules/development-feedback-loop.md` — full rule text
- `.claude/rules/testing.md` — test-tier requirements
- `simulation-validation` skill — tiered benchmark menu
- `handover-update` skill — where to record test outcomes
