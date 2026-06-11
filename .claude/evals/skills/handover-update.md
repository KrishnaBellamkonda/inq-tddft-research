# Eval: handover-update skill (LOCKED 2026-06-11)

The skill now carries the FULL procedure (the `handovers` rule slimmed to a thin
trigger index). This eval verifies the procedure survived the slimming.
Evaluator: **LLM-judge + human gate**.

## Trigger test
- **Positive:** "record where we are before I compact" · "write a handover for
  this task" · (at a milestone) "we just finished the merge, update the handover".
- **Negative:** "log this run with my observation" → `journal-writing` ·
  "draft the methods section" → `report-writing`.

## Functional rubric (all hard)
1. **Prepends** a dated section to `docs/handovers/<task>.md` (does not append at
   the bottom; does not overwrite prior sections).
2. Contains the 10 required sections (status / what changed / files touched /
   commands / tests & validation / sources / attribution / known issues /
   assumptions / exact next steps).
3. All file paths **absolute**.
4. States what is done / partial / not done; records validation status
   (pass/fail/unverified) — no "assumed correct".
5. Resumable: another session could continue from it without guessing.

PASS = 5/5. (Minor-milestone `Update:` form is allowed but still needs 1, 3, 4.)
