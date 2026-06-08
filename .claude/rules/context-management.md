# Rule: Context and Token Management

Apply to: entire project

## Rules

1. Keep `CLAUDE.md` minimal (under ~200 lines). Move specialised content to skills and path-scoped rules.

2. Do not load more tools, agents, MCP servers, plugins, or files than the current task requires.

3. Do not spawn agent teams by default. Spawn an agent only if the user explicitly requests it, or the task is clearly decomposable and high-volume (e.g. log processing, large documentation retrieval).

4. Use at most one subagent at a time unless the user requests parallel work.

5. Read only the files needed for the current decision. Do not "read everything" unless genuinely necessary.

6. Avoid broad prompts like "understand the whole repository" when the task is narrow.

7. Before manual compaction, update the handover in `docs/handovers/`.

8. After compaction or resume, read the latest handover before continuing substantive work.

9. Use `/clear` only between unrelated tasks, not mid-task.

10. Prefer CLI tools over MCP servers when a trusted CLI already solves the problem.

## Token-efficient patterns for this project

- For build questions: read `docs/compilation.md` rather than scanning the full `inq/` source.
- For API questions: read `docs/inq_tutorial.md` first, then spot-check source headers if uncertain.
- For architecture questions: read `docs/inq_source_map.md` rather than reading all 178 source files.
- For session continuity: read `docs/handovers/<task>.md` at session start.
