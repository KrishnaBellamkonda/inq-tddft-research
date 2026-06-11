# Eval: test-validation subagent (LOCKED 2026-06-11)

Component (subtask 3): a fresh-context subagent given a written test + the
already-locked expected value; it audits for **circularity** (asserting code
output vs the independently-verified value), correct tolerances/units, and
isolation. Returns CONFIRM / FLAG. Designed in `CONTEXT.md`. Never sees the
formula-validation agent's derivation (independence is the point).

Evaluator: **planted-bug fixtures**, expected verdict fixed up front.

## Fixtures (test | planted flaw? | expected verdict)

| # | Test description | Expected |
|---|---|---|
| T1 | asserts `cod == (0,0,0)` for a uniform-box density (value from the COD *definition*, not the code) | **CONFIRM** |
| T2 | asserts `cod == run_function(uniform_box)` — i.e. expected captured from the same code under test | **FLAG** (circular) |
| T3 | asserts `norm ∈ [0.97,1.03]` for an analytic normalised WP | **CONFIRM** |
| T4 | asserts energy drift `< 1 Ha` (tolerance 1000× too loose; should be ~1 mHa) | **FLAG** (tolerance) |
| T5 | asserts `cod_x` in Ångström vs an expected value computed in Bohr | **FLAG** (unit mismatch) |
| T6 | golden compared with `np.allclose(atol=1e-9)` for a deterministic CSV parse | **CONFIRM** |
| T7 | "test" that only checks the function runs without raising (no value assertion) | **FLAG** (asserts nothing) |

## Pass criterion

≥ 6/7 match. A CONFIRM on T2 (circular) or T7 (no assertion) is a **hard fail** —
these are exactly the test-smells the agent exists to catch.

## Validation log

- 2026-06-11: agent prompt validated on a 2-fixture sample (T1, T2) via primed
  general-purpose agents. T1 (sound, analytic expected) → CONFIRM; T2 (circular,
  expected captured from the code under test) → FLAG. 2/2 expected. Full
  7-fixture suite: run once the agent is registered.
