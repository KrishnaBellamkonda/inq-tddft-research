# Claude Code Scientific Workflow Rules Pack

Version: 2026-04-04
Purpose: Portable hard rules, soft rules, folder conventions, and operating guidance for scientific coding projects using Claude Code.

## 1) Design goals

This rules pack is designed to achieve five goals:
1. Keep startup context and token use low.
2. Force explicit planning and periodic handovers.
3. Ground scientific claims, code, and reports in trustworthy sources.
4. Enforce validation and test-driven scientific coding.
5. Keep repository outputs organised, readable, and portable across sessions and projects.

---

## 2) Recommended repository layout

Use this structure unless the repository already has a better established equivalent.

```text
.claude/
  CLAUDE.md
  rules/
    testing.md
    scientific-grounding.md
    file-placement.md
    handovers.md
    context-management.md
  skills/
    literature-review.md
    simulation-validation.md
    report-writing.md
    handover-update.md

docs/
  plans/
  handovers/
  reports/
  notes/
  sources/
  validation/
```

Primary working memory should live in repo-local `docs/` files. Use `.claude/` for rules, hooks, and optional skills, not as the main store of evolving project memory.

### Folder intent
- `docs/plans/`: explicit task plans, one active plan per substantive task.
- `docs/handovers/`: milestone summaries and session continuation notes.
- `docs/sources/`: literature notes, source summaries, citation mapping, attribution notes.
- `docs/validation/`: test matrices, benchmark definitions, expected behaviours, comparison plots, validation notes.
- `docs/reports/`: report drafts, manuscript fragments, slide outlines, and figure-caption working files.
- `docs/notes/`: temporary but structured working notes that do not belong in plans, handovers, source summaries, or reports.
- `.claude/rules/`: always-on or path-scoped project rules.
- `.claude/skills/`: on-demand workflows and reference material to avoid bloating startup context.

---

## 3) Operating principles

### A. Keep startup context lean
- Keep `CLAUDE.md` short, stable, and under roughly 200 lines.
- Put only always-on project rules in `CLAUDE.md`.
- Move specialised workflows, literature summaries, validation protocols, and report-writing guides into skills.
- Use path-scoped rules in `.claude/rules/` so instructions load only when relevant files are touched.
- Disable or remove unused MCP servers, plugins, and agents.
- Do not spawn agent teams unless the user explicitly requests or approves them.
- Prefer a single main session.
- Use a subagent only when a task is high-volume and isolated, such as log processing, documentation retrieval, or large test output summarisation.
- Prefer CLI tools over MCP servers when a trusted CLI already exists.
- For scripted or one-off runs where deterministic minimal startup matters more than convenience, use `--bare` and pass only the needed settings explicitly.

### B. Explore first, then plan, then implement
- Before making substantive code changes, inspect the relevant files and create or update a plan in `docs/plans/`.
- Plans must state what is known, what is assumed, what files are relevant, how success will be checked, and what risks or rollback steps exist.
- Do not code against guessed architecture or guessed APIs.

### C. Ground scientific work
- All scientific explanations, implementation choices, and modelling recommendations must be grounded in trustworthy sources.
- Acceptable grounding sources include journal papers, textbooks, official documentation, authoritative institutional pages, lecture notes or slides from reputable universities, and project-specific notes explicitly approved by the user.
- Do not present scientific claims as fact without either direct evidence from sources or explicit labelling as an inference.
- When uncertain, say uncertain and verify.
- Record important sources in `docs/sources/` and cross-reference them from plans, handovers, reports, and code comments where relevant.

### D. Credit inspirations and prior work
- If code structure, formulas, algorithms, constants, workflows, or validation setups are adapted from a paper, package, repository, textbook, or prior internal note, document that influence.
- Credit should appear in the most appropriate place:
  - source note in `docs/sources/`
  - plan note in `docs/plans/`
  - handover note in `docs/handovers/`
  - code comment or module docstring near the adapted logic
  - report citation in the manuscript or slides
- Never imply original authorship for borrowed or adapted scientific ideas.

### E. Test-driven scientific coding
- No substantive code should be written without an explicit validation plan.
- Before implementation begins, propose:
  - unit or component tests for each core function
  - integration tests for the assembled workflow
  - scientific validation tests against known systems, analytic limits, conservation laws, benchmark datasets, or literature results
- Claude must suggest candidate tests to the user before running expensive, destructive, or long simulations.
- The user must decide which expensive tests or simulations to run.
- After implementation, Claude must report:
  - which tests were proposed
  - which tests were approved by the user
  - which tests were actually run
  - which tests passed or failed
  - what remains unverified

### F. Periodic handovers and continuity
- Maintain a current handover in `docs/handovers/` for every substantive task.
- Update the handover at each meaningful milestone, before stopping, before clearing context, before compaction, and before claiming the task is complete.
- A handover must be sufficient for another session to continue without guessing.

### G. File placement discipline
- Create files only in designated directories.
- Do not scatter notes, reports, scripts, or scratch files across arbitrary folders.
- If the repository lacks an appropriate folder, propose one or two sensible locations and ask the user to choose before creating many new files.
- Prefer updating existing documentation over creating duplicates.
- Keep created documents human-readable, succinct, and information-dense.

---

## 4) Hard rules

These are non-negotiable.

### Hard Rule 1: Minimal-context startup
- Never load or create more tools, agents, or external connections than the current task requires.
- Never use agent teams by default.
- Never enable or keep unused MCP servers for the session.
- Never move large specialised content into `CLAUDE.md` if it can live in a skill or path-scoped rule instead.

### Hard Rule 2: Explicit planning
- Before any substantive implementation, create or update a plan file in `docs/plans/`.
- No substantive implementation may begin until the plan includes objective, assumptions, target files, verification strategy, and next steps.

### Hard Rule 3: Mandatory handovers
- Before stopping, clearing context, compacting, or declaring a milestone complete, create or update a handover in `docs/handovers/`.
- If context becomes crowded, update the handover before any manual compact or clear operation.

### Hard Rule 4: Scientific grounding
- Do not make scientific claims, modelling recommendations, algorithmic choices, or validation claims without trustworthy grounding or an explicit uncertainty label.
- If a claim is inferred rather than directly stated by a source, label it as an inference.

### Hard Rule 5: Attribution
- If any implementation or explanation is adapted from prior work, record the source and the nature of the adaptation.
- Code inspired by external or internal prior work must contain appropriate credit near the relevant implementation or in the corresponding source note.

### Hard Rule 6: Testing and validation gating
- No code is complete until an explicit validation status is recorded.
- For each substantive code change, define component tests and at least one integration or benchmark-style validation route.
- Claude must not silently choose expensive scientific validation runs on the user’s behalf.
- Claude must ask the user which expensive tests or simulations they want to authorise.

### Hard Rule 7: No fabricated certainty
- Never invent paths, APIs, equations, constants, file contents, benchmark outcomes, or successful test results.
- When uncertain, say uncertain and verify.

### Hard Rule 8: Organised file creation
- All new files must be placed in designated directories.
- If no suitable directory exists, suggest a new folder structure before creating a spread of files.

### Hard Rule 9: Readability
- Plans, handovers, and documentation must be clear, concise, and easy for a human collaborator to scan quickly.
- Avoid long narrative padding.

---

## 5) Soft rules

These should usually be followed, but can be overridden when the user asks.

- Prefer one active plan file per task instead of many fragmented planning notes.
- Prefer one rolling handover file per task, with dated milestone sections, unless the project genuinely needs separate handovers.
- Prefer path-scoped rules over giant always-on rules.
- Prefer skills for literature packs, simulation recipes, reporting conventions, and benchmark procedures.
- Prefer subagents only for isolated, high-volume work.
- Prefer trusted CLI tools over MCP where both solve the same problem.
- Prefer updating existing docs over creating duplicates.
- Prefer exact quotes only when needed; otherwise summarise sources clearly.
- Prefer explicit assumptions lists when the scientific setup is under-specified.
- Prefer small, reviewable edits over broad refactors unless the user asks for broader change.

---

## 6) Required contents for plan files

Every substantive plan file in `docs/plans/` should contain these sections:

```md
# Task plan: <task name>

## Objective
## User request
## Relevant files / modules
## Trusted sources to use
## Known facts
## Assumptions
## Proposed approach
## Validation plan
### Component tests
### Integration tests
### Scientific benchmark tests
## Risks / failure modes
## Rollback or safety notes
## Open questions for the user
## Next action
```

### Plan rules
- Include source names for scientific grounding.
- State clearly which parts are evidence-backed and which are assumptions.
- For simulation work, include a benchmark or sanity-check section.
- If literature influenced the design, mention it explicitly.

---

## 7) Required contents for handover files

Every handover in `docs/handovers/` should contain these sections:

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

### Handover rules
- Be continuation-ready.
- State what is done, what is partially done, and what is not done.
- Record failed attempts if they matter for future sessions.
- Record what was verified and what remains unverified.
- Record any important scientific or implementation rationale that must survive compaction.

---

## 8) Scientific grounding protocol

Use this protocol for any scientific task.

1. Identify which claims require grounding.
2. Gather trustworthy sources before making strong recommendations.
3. Distinguish direct source statements from your inferences.
4. Record the key sources in `docs/sources/`.
5. Propagate source awareness into:
   - the task plan
   - the handover
   - any code comments for adapted logic
   - reports or manuscripts
6. If the task depends on unsettled assumptions, flag them explicitly.

### Minimum source-quality expectations
- Prefer peer-reviewed papers and authoritative textbooks.
- Use official documentation for software behaviour.
- Use university notes or slides when they are clearly authored and reputable.
- Treat forums and unreviewed snippets as low-trust unless independently verified.

---

## 9) Attribution protocol for code and reports

### In code
Add a short attribution comment when logic is directly adapted or strongly inspired.

Example:
```python
# Algorithm structure adapted from <paper/repo/textbook>, modified here for <project-specific reason>.
```

### In plans and handovers
Record:
- source name
- what was borrowed or adapted
- what was changed
- any known limitations of the borrowed approach

### In reports
Use formal citations and state when a method is adopted, adapted, or reproduced. Use the project’s required citation style. If none is specified, use a consistent author-year style in notes and plans, and the target journal style in reports.

---

## 10) Testing and validation protocol

### A. Component tests
For each meaningful function or module, define tests for:
- expected inputs
- edge cases
- failure modes
- units and dimensions where applicable
- deterministic behaviour where required

### B. Integration tests
For the assembled workflow, define tests for:
- correct module interaction
- expected artefact generation
- correct parameter passing
- correct I/O behaviour
- reproducibility or restartability where relevant

### C. Scientific validation tests
Where applicable, suggest tests against:
- analytically solvable cases
- conserved quantities
- limiting behaviours
- known benchmark systems
- published reference values
- symmetry checks
- mesh, timestep, cutoff, or convergence checks

### D. User approval rule
Before running expensive simulations, long test suites, or destructive workflows, Claude must present the proposed validation menu and ask the user which ones to authorise.

### E. Reporting rule
Every implementation handover must state:
- proposed tests
- approved tests
- executed tests
- observed outcomes
- remaining gaps

---

## 11) Context and token management rules

### Always-on guidance
- Keep the main session focused on the current task.
- Use `/clear` only between unrelated tasks.
- Use compaction for long related tasks.
- Before manual compaction, update the handover.
- After compaction or resume, read the latest handover before continuing substantive work.

### Strong defaults
- Do not ask Claude to broadly “read everything” unless that is truly necessary.
- Read only the files needed for the current decision.
- Avoid broad prompts like “understand the whole repo” when the task is narrow.
- When large background material exists, store it as a skill or source note and load it only when needed.
- For repeated scripted usage, prefer `--bare` plus explicit settings over a heavy global startup.

### Agent policy
- Default: no agent teams unless the user explicitly requests or approves them.
- Default: no more than one subagent unless the user requests parallel work or the task is clearly decomposable.
- Use subagents mainly for high-volume operations such as tests, logs, documentation fetching, or isolated codebase research.

---

## 12) Suggested hook policy

Use hooks to enforce the workflow instead of relying on memory alone.

### Recommended hooks
- `PreToolUse` on `Edit|Write`
  - block substantive edits if no active plan exists
  - block file creation outside approved directories
- `PostToolUse` on `Edit|Write`
  - update lightweight task-state metadata
  - refresh the active handover timestamp after milestone edits
- `PreCompact`
  - trigger handover update before compaction
- `PostCompact`
  - append the compact summary into the latest handover
- `SessionStart` with source `compact` or `resume`
  - inject the latest handover summary into context
- `Stop`
  - prevent stopping if the handover is stale or validation status is missing
- `TaskCompleted`
  - verify that tests and validation status have been recorded
- `SessionEnd`
  - write a final continuation-ready handover snapshot

### Test-gating hook idea
Use a prompt or agent hook so Claude checks whether:
- component tests are defined
- integration tests are defined
- expensive validation runs were explicitly approved by the user
before it is allowed to declare the task complete.

---

## 13) Suggested path-scoped rules

### `.claude/rules/testing.md`
Apply to code and test directories.

Rules:
- every substantive change must include a validation update
- define component and integration tests
- propose scientific benchmark tests when relevant
- do not claim success without reporting what was actually run

### `.claude/rules/scientific-grounding.md`
Apply to `src/`, `scripts/`, `docs/reports/`, `docs/plans/`, and `docs/handovers/`.

Rules:
- scientific claims must be grounded
- distinguish source statements from inferences
- track attributions for adapted code and methods
- propagate source awareness into code, plans, handovers, and reports

### `.claude/rules/file-placement.md`
Apply project-wide.

Rules:
- create files only in designated directories
- do not scatter scratch files
- propose folder additions before broad file creation

### `.claude/rules/handovers.md`
Apply to `docs/handovers/**` and `docs/plans/**`.

Rules:
- keep docs short, structured, and continuation-ready
- update handovers before stop, clear, compact, and completion

### `.claude/rules/context-management.md`
Apply project-wide.

Rules:
- keep CLAUDE.md minimal
- use skills for heavy reference content
- do not spawn agent teams unless explicitly requested
- use subagents only when their isolation benefit outweighs their setup cost

---

## 14) Suggested skills

### `literature-review`
Use for:
- source gathering
- evidence extraction
- uncertainty tracking
- source-note creation in `docs/sources/`

### `simulation-validation`
Use for:
- proposing benchmark systems
- defining convergence tests
- proposing component, integration, and physics validation tests

### `report-writing`
Use for:
- source-grounded writing
- attribution continuity
- converting handovers and source notes into report-ready prose

### `handover-update`
Use for:
- writing concise, standardised handover updates from current repo state

---

## 15) Suggested concise CLAUDE.md core

Use the following as the compact always-on core.

```md
# Core workflow rules

- Keep startup context lean. Do not use more tools, agents, MCP servers, plugins, or files than the task requires.
- Do not spawn agent teams unless the user explicitly requests or approves them.
- Keep `CLAUDE.md` minimal. Put heavy reference content into skills and path-scoped rules.
- Before substantive implementation, create or update a plan in `docs/plans/`.
- Before stopping, compacting, clearing context, or declaring completion, create or update a handover in `docs/handovers/`.
- Create files only in designated directories. If no suitable directory exists, propose one before creating many new files.
- Scientific claims, modelling choices, and validation claims must be grounded in trustworthy sources or explicitly labelled as uncertain inference.
- Record important sources and attribution notes in plans, handovers, reports, and code comments where relevant.
- No substantive code is complete without recorded validation status.
- For each substantive change, define component tests, integration tests, and scientific benchmark or sanity checks where applicable.
- Suggest test options to the user before running expensive simulations. The user decides which expensive tests to run.
- Never invent paths, APIs, equations, constants, file contents, or test results.
- When uncertain, say uncertain and verify.
- Prefer updating existing docs over creating duplicates.
- Keep plans, handovers, and notes clear, human-readable, and concise.
```

---

## 16) Practical user prompt additions

These short prompts help steer token use and scientific discipline.

### For low-token planning
- "Use only the minimum required tools and files. Do not spawn subagents or connect external tools unless needed."
- "Read only the files required for this decision. Do not scan the whole repository."
- "Keep startup context lean. Use existing skills or notes instead of loading large documents into always-on memory."

### For scientific tasks
- "Ground every scientific recommendation in trustworthy sources and label any inference clearly."
- "Record source and attribution notes in the plan and handover."

### For validation
- "Before implementing, propose unit, integration, and scientific benchmark tests. Ask me which expensive tests I want to authorise."

---

## 17) Final recommendation

For portability across projects:
1. Keep a very small shared `CLAUDE.md` core.
2. Keep reusable topic-specific rules in `.claude/rules/`.
3. Keep heavy scientific workflows and literature packs as skills.
4. Use hooks for enforcement.
5. Use a shared rules repository or plugin once the workflow stabilises.


