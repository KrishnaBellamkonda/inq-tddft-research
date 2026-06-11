# Eval: literature-review skill (LOCKED 2026-06-11)

Skill carries the full procedure (the `scientific-grounding` rule slimmed to its
policy invariants). Evaluator: **LLM-judge + human gate**.

## Trigger test
- **Positive:** "justify this functional choice with sources" · "ground this
  claim about the time step" · "write a source note for Castro 2004".
- **Negative:** "draft the report" → `report-writing` · "propose a validation
  menu before the run" → `simulation-validation`.

## Functional rubric (all hard)
1. Writes/updates `docs/sources/<author-year-keyword>.md` with the template
   (full citation / relevance / key claims w/ page-section / limitations /
   cross-refs).
2. **Separates direct source statements from inference** — inferences explicitly
   labelled ("Inference:" / "This suggests"); never presents a guess as a
   sourced claim.
3. Prefers high-trust sources (peer-reviewed → textbooks → official INQ/libxc
   docs); flags low-trust sources as needing independent verification.
4. Credits adapted logic near the code (`// Adapted from <Author Year>`) and in
   the plan/handover.

PASS = 4/4. Negative guard: a claim with no citation must be labelled inference,
not asserted as fact.
