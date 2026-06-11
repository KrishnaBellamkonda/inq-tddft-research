# Eval: jellium-base-run-spec removal (LOCKED 2026-06-11, D6')

The always-on `jellium-base-run-spec` rule is DELETED (stale N=138; no "base
jellium" concept). Current jellium params come from live `shared/configs/`;
guardrails survive in `boundary_rule.hpp` + `tddft-simulations`.

## Removal proof (programmatic + LLM-judge)
- The file `.claude/rules/jellium-base-run-spec.md` no longer exists.
- In a non-run session (e.g. "draft the report intro"), no N=138/L=30/Nyquist
  jellium spec is injected into context (context-cost win).

## Load-on-demand / guardrail retention (LLM-judge + human)
- Trigger: "set up a jellium WP run at E=50" → `tddft-simulations` sources the
  current canonical params from `shared/configs/` (NOT the deleted N=138 rule)
  and applies the guardrails:
  - closed-shell electron count for the chosen box (current campaign: N=162);
  - Nyquist `dx ≤ π/(k₀+3σ_k)` — a violating `dx` is flagged (hard stop);
  - boundary 4σ/1σ launch-stop + ~300-frame VTI cadence (`boundary_rule.hpp`).
- Negative guard: the model must NOT cite the stale N=138/L=30 numbers as
  canonical.

## Pass criterion
File gone + no stale spec in non-run context + guardrails still fire on run
config + no N=138 citation. The physics guardrails must survive the deletion.
