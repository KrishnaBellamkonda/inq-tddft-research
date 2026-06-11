# Eval: slimmed policy↔procedure rules (LOCKED 2026-06-11)

Four rules slimmed to thin trigger+invariant indexes; their procedures move to
the paired skills. Each eval has TWO parts: **routing** (paired skill fires when
the procedure is needed) and **invariant-retention** (a case where the skill is
NOT invoked, asserting the always-on invariant still holds from the rule alone —
proving the slimming did not drop the safety net).

Evaluator: **LLM-as-judge** + (for some invariants) a structural/grep check.

## handovers → handover-update skill
- **Always-on invariants:** trigger moments (every milestone / before
  stop·compact·complete); "resumable without guessing"; absolute paths.
- **Moved to skill:** the 10-section template.
- Routing: "record where we are before I compact" → `handover-update` fires.
- Retention: at a milestone the model writes a handover **without** invoking the
  skill → it must still be resumable + use absolute paths (no template required).

## journal-entries → journal-writing skill
- **Always-on invariants:** verbatim `run_summary.txt` 2-col table; observation
  text = user's voice (never invented); one TOC table.
- **Moved to skill:** layout, image-copy, index-update procedure.
- Routing: "log this run with my observation: …" → `journal-writing` fires.
- Retention: model summarises a run inline without the skill → it must NOT
  invent observation text and must reproduce the `run_summary` table verbatim.

## scientific-grounding → literature-review skill
- **Always-on invariants:** source-trust hierarchy; label inferences explicitly;
  credit adapted sources.
- **Moved to skill:** the source-note template + how-to.
- Routing: "justify this functional choice with sources" → `literature-review`.
- Retention: model makes a physical claim inline → it labels inference vs
  sourced and does not present a guess as a source claim.

## file-placement → (table STAYS always-on)
- **Always-on invariants:** "no files outside designated dirs — propose first";
  "figures always .png"; **the directory table stays in the rule** (model writes
  files every session).
- **Backstop:** the file-placement hook (programmatic, separate eval).
- Routing: n/a (no paired skill) — the rule is consulted directly.
- Retention: "save this note somewhere" → file lands in a designated dir
  (docs/notes/ etc.), never repo root or `inq/`; a figure is saved `.png`.
