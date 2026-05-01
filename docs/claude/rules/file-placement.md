# Rule: File Placement

Apply to: entire project

## Designated directories

| Type of file | Directory |
|---|---|
| Task plans | `docs/plans/` |
| Session handovers | `docs/handovers/` |
| Literature notes, source summaries, citations | `docs/sources/` |
| Test matrices, benchmark definitions, validation notes | `docs/validation/` |
| Report drafts, manuscript fragments, figure captions | `docs/reports/` |
| Temporary working notes | `docs/notes/` |
| INQ configuration files (user `.cpp`) | `ResearchProject/systems/<material>/<task>/` or `ResearchProject/jellium/<N_task>/` |
| Tutorial examples | `Tutorial/<name>/` |
| inqkit C++ headers | `inq-stack/include/inqkit/<module>/` |
| inqview Python modules | `inq-stack/python/inqview/` |
| Always-on project rules | `.claude/rules/` |
| On-demand skills and reference material | `.claude/skills/` |

## Rules

1. Do not create files outside the directories above without proposing the location first.

2. Do not scatter notes, scratch files, or reports into arbitrary locations (e.g. project root, `inq/`, `shared/`).

3. Do not create files inside `inq/src/` unless modifying the INQ source is explicitly requested.

4. If no suitable directory exists for a new file type, propose one or two sensible locations and ask the user to choose before creating files.

5. Prefer updating an existing document over creating a duplicate.

6. Do not create `README.md` or other documentation files unless explicitly requested.

7. Always save figures as `.png`. Never save figures as `.pdf` or `.svg` unless the user explicitly requests it.

## Naming conventions

- Plans: `docs/plans/<task-name>.md`
- Handovers: `docs/handovers/<task-name>.md` (rolling file with dated milestone sections)
- Source notes: `docs/sources/<author-year-keyword>.md`
- Validation notes: `docs/validation/<system-property>.md`
