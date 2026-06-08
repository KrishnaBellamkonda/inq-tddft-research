# Rule: Journal Entries

Apply to: `docs/journals/<journal>/`

## Rules

1. **Run-based entries must paste the full `run_summary.txt` as a
   two-column markdown table.** No editorial reordering, no abridgement
   beyond escaping pipes/newlines. Every run-based journal entry has a
   `## Run summary` (or `## Config snapshot`) section containing this
   table verbatim. The journal-writing skill (`.claude/skills/journal-
   writing/SKILL.md`, Step 2) already encodes this; this rule restates it
   so it cannot be silently dropped.

   **Why:** the table is the canonical config artefact that lets a
   reader (or a future session) reproduce the run without crawling the
   `run_summary.txt` file. Dropping it strips the entry of its scientific
   provenance.

   **How to apply:** if the source `run_summary.txt` has section headers
   (`1. Run identity`, `3. System configuration`, …) keep them as
   `### <heading>` subsections; the table inside each subsection lists
   the `key = value` pairs of that section. Do not collapse multiple
   sections into one giant table — keep the original grouping.

2. **Topical entries (no single run)** are exempt from rule 1. They have
   `Linked entries:` with a bullet list of the per-run entries they
   summarise. Run-summary tables stay in the per-run entries.

3. **Observation text is the user's voice.** Never invent observation
   text; ask the user. The `## Observations` section reproduces the
   user's words verbatim, with image references rewritten to
   `attachments/<slug>/...`.

4. **One TOC table.** `index.md` has exactly one TOC table, newest at the
   top.
