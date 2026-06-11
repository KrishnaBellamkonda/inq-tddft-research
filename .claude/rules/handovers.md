# Rule: Handovers

Apply to: `docs/handovers/`, `docs/plans/`

PROCEDURE (the 10-section template, prepend mechanics, quality checklist) lives
in the `handover-update` skill — invoke it to write/update a handover. This rule
is the always-on trigger + invariants.

## When

Write or update the rolling handover at `docs/handovers/<task>.md`:
- every meaningful milestone
- before stopping, going idle, clearing context, or compacting
- before declaring a task complete

After compaction or session resume, read the latest handover before continuing
substantive work.

## Invariants (hold even without invoking the skill)

- One rolling file per substantive task.
- Resumable: sufficient for another session (or a different account) to continue
  without guessing.
- State what is done / partial / not done; record what was verified vs
  unverified, failed attempts that matter, and any rationale that must survive
  compaction.
- Absolute file paths, not relative. Concise and human-scannable.
