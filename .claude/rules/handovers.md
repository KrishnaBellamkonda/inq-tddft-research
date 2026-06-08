# Rule: Handovers

Apply to: `docs/handovers/`, `docs/plans/`

## Rules

1. Maintain one rolling handover file per substantive task in `docs/handovers/`.

2. Update the handover at:
   - every meaningful milestone
   - before stopping or going idle
   - before clearing context or compacting
   - before declaring the task complete

3. A handover must be sufficient for another session (or a different Claude account) to continue without guessing.

4. After compaction or session resume, read the latest handover before continuing substantive work.

## Required handover sections

```md
# Handover: <task name>

## Current status
## What changed
## Files touched
## Commands run
## Tests and validation
## Trusted sources used
## Attribution notes
## Known issues / blockers
## Assumptions still in play
## Exact next steps
```

## Handover quality rules

- State what is done, what is partially done, and what is not done.
- Record failed attempts if they matter for future sessions.
- Record what was verified and what remains unverified.
- Record any scientific or implementation rationale that must survive compaction.
- Keep handovers concise and human-scannable — avoid long narrative padding.
- Use absolute file paths, not relative paths.
