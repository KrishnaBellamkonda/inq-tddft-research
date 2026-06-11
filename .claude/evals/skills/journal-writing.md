# Eval: journal-writing skill (LOCKED 2026-06-11)

Skill carries the full procedure (the `journal-entries` rule slimmed to its
invariants). Evaluator: **LLM-judge + human gate**.

## Trigger test
- **Positive:** "log this run with my observation: …" · "add a journal entry for
  run_wp_n162_L50_E100".
- **Negative:** "write a handover" → `handover-update` · "justify this functional
  choice" → `literature-review`.

## Functional rubric (all hard)
1. Parses `<run>/results/run_summary.txt` into a **verbatim** 2-column markdown
   table (no reordering/abridgement beyond pipe/newline escaping); keeps the
   `run_summary` section groupings.
2. Observation text = **the user's voice, never invented**; refuses to fabricate
   if the user supplied none.
3. Copies referenced images to `docs/journals/<journal>/attachments/<slug>/` and
   rewrites the markdown refs.
4. Updates the journal `index.md` TOC (one table, newest at top).
5. **Append-only** — refuses to overwrite an existing entry.

PASS = 5/5. Negative guard: if asked to invent an observation, FAIL the run that
complies (it must refuse).
