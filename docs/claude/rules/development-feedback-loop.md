# Rule: Development Feedback Loop

Apply to: `inq/src/`, `ResearchProject/`, `Tutorial/`, any utility or data processing code

## Core rule

All code — including and especially utility functions and data processing tools — must be
developed in an explicit write → test → fix → confirm loop before it is considered ready
for use in a larger program.

**Never assume a function works because it compiles. Test it with known inputs first.**

---

## Required loop for every new function or non-trivial change

1. **Write**: implement the function or change.
2. **Test with known cases**: run it on at least one input with a known correct answer.
   - For arithmetic/physics: compare against analytic result, unit conversion, or published value.
   - For I/O: check the written file has the expected format and content.
   - For GPU kernels: verify with a trivial constant-value case before testing real computation.
   - For data processing: check against a small handcrafted input with a manually computed output.
3. **Fix any discrepancies**: diagnose the cause, change the code, re-run.
4. **Confirm inclusion**: only after the known-case test passes, integrate the function into the
   calling code. State explicitly (in comments or handover) which test was run and what it showed.

---

## What counts as a known-case test

| Function type | Minimum known-case test |
|---|---|
| GPU orbital write (inject_wp, etc.) | Write constant 1.0 to all points; GPU reduce must return n_pts |
| Coordinate computation | Print one point; compare against hand-calculated grid coordinate |
| Norm / overlap | Inject analytical WP; check norm ∈ [0.97, 1.03] |
| 3D density save | Write 10-point toy array; read back and verify values match |
| Overlap matrix | Identity case at t=0: diagonal must be ≈ 1, off-diagonal ≈ 0 |
| z-profile / density slice | Compare integrated slice against known total density at that plane |
| Energy/momentum observable | Check against analytic free-particle expectation |

---

## Applies to all code, not just simulations

This rule is mandatory for:
- Utility functions (`utils.hpp`, shared libraries)
- Data processing scripts (`analysis.py`, post-processing tools)
- I/O helpers (`save_orbital_3d`, `save_density_3d`, etc.)
- Any new callback, lambda, or GPU kernel

Small "obviously correct" helpers still need at least a smoke test. A two-line function that
computes the wrong index silently will corrupt all downstream results.

---

## Record of what was tested

Every handover must state — for each new function — which test was run and what result
was observed. Example:

```
inject_wp: GPU constant-write test showed sum=1.75e6 == n_pts ✓;
           Gaussian WP norm = 1.002 ∈ [0.97, 1.03] ✓
save_orbital_3d: not yet tested — must verify before run.cpp launch
```

"Not yet tested" is acceptable to record; "assumed correct" is not.
