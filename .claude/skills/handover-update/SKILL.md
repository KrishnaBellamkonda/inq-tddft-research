---
name: handover-update
description: Use when writing or updating a session handover in `docs/handovers/<task>.md` — required at every milestone, before stopping/compacting/clearing context, and before declaring a task complete. Produces a standardised template another session (or a different Claude account) can resume from without guessing.
---

# Handover Update

Invoke this before stopping, before compacting, at each meaningful milestone,
and before declaring a task complete.

## 1. Collect current state

Before writing, check:
- What files were created or modified in this session? (`git status` or recent edits)
- What commands were run? (from memory / bash history if available)
- What tests were run? Did they pass?
- What remains unfinished?
- What decisions were made and why?
- What assumptions are still in play?

## 2. Locate or create the handover file

- File: `/local/data/public/skcb2/tddft/docs/handovers/<task-name>.md`
- If it already exists, **prepend** a new dated milestone section at the top.
- If it does not exist, create it with the full template.

## 3. Template

```md
# Handover: <task name>

---

## Milestone: <YYYY-MM-DD> — <one-line summary of what was achieved>

### Current status
<One paragraph: where things stand right now. Done / in progress / blocked.>

### What changed
- <file or module>: <what changed and why>

### Files touched
- `<absolute path>` — <purpose>

### Commands run
```bash
<exact commands that were run>
```

### Tests and validation
- Proposed: [list]
- Approved: [list]
- Run: [list]
- Outcomes: [pass/fail/partial]
- Remaining gaps: [list]

### Trusted sources used
- <Author Year> — <relevance>

### Attribution notes
- <file:line> adapted from <source>

### Known issues / blockers
- <issue>

### Assumptions still in play
- <assumption>

### Exact next steps
1. <first action the next session should take>
2. ...
```

## 4. Quality checks before finishing

- [ ] Can another session continue from this handover without asking the user?
- [ ] Are all file paths absolute (not relative)?
- [ ] Is the validation status recorded (even if "not yet run")?
- [ ] Are assumptions explicitly listed?
- [ ] Are next steps concrete and numbered?

## Quick update (for minor milestones)

If the change is small, prepend a brief section instead of the full template:

```md
## Update: <YYYY-MM-DD> — <one-line summary>

Status: <done/in progress>.
Changed: <files>.
Next: <one-liner next step>.
```
