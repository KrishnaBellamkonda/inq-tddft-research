# Rule: Journal Entries

Apply to: `docs/journals/<journal>/`

PROCEDURE (parse `run_summary.txt`, copy images to `attachments/<slug>/`, slug +
`index.md` mechanics) lives in the `journal-writing` skill. This rule is the
always-on invariants — restated so they cannot be silently dropped.

## Invariants

1. **Run-based entries paste the full `run_summary.txt` VERBATIM** as a
   two-column markdown table, keeping its section groupings (no reordering or
   abridgement beyond escaping pipes/newlines) — the canonical config artefact
   that lets a future session reproduce the run.
2. **Topical entries** (no single run) are exempt from (1); they carry
   `Linked entries:` to the per-run entries they summarise.
3. **Observation text is the user's voice** — never invented; ask the user.
4. **One TOC table** in `index.md`, newest at top.
