# Eval: formula-validation subagent (LOCKED 2026-06-11)

Component (subtask 3): a fresh-context subagent given ONLY a formula-as-
implemented + its cited source; it re-derives/sanity-checks the math
independently and returns CONFIRM / FLAG. Designed in `CONTEXT.md`.

Evaluator: **planted-bug fixtures** — expected verdict fixed BEFORE the agent
exists (anti-circularity). The agent is correct iff it reproduces these verdicts.

## Fixtures (formula | source | planted? | expected verdict)

| # | Formula as implemented | Source | Expected |
|---|---|---|---|
| F1 | `cod = Σ r·n / Σ n` | COD definition ∫r·n/∫n | **CONFIRM** |
| F2 | `cod = Σ r·n` (missing `/Σn` normalisation) | same | **FLAG** (not normalised) |
| F3 | `momentum_dist = |FFT(ψ)|²` | \|ψ̃(k)\|² | **CONFIRM** |
| F4 | `momentum_dist = |FFT(ψ)|` (missing square) | same | **FLAG** |
| F5 | jellium magic numbers `2,8,18,20,34,40,…` closed-shell counts | free-electron-gas shell filling | **CONFIRM** |
| F6 | jellium magic numbers with `18`→`16` (planted) | same | **FLAG** |
| F7 | Lindhard ε(q,ω) static limit → Thomas–Fermi `k_TF² = 4k_F/π a0` | Ashcroft–Mermin | **CONFIRM** |
| F8 | stopping `S = -dE/dx` with sign flipped (`+dE/dx`) | stopping-power defn | **FLAG** |

## Pass criterion

≥ 7/8 verdicts match (the agent may legitimately hedge on F7's algebra but must
not FLAG a correct formula or CONFIRM a planted one). A CONFIRM on any planted-
bug fixture (F2/F4/F6/F8) is a **hard fail** — that's the failure mode the agent
exists to prevent.

## Independence guard

The agent must be given the formula + source ONLY — never the test, never the
main session's reasoning, never the test-validation agent's output.
