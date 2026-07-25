---
name: campaigns
description: Use when designing, authoring, refining, or status-tracking an autonomous research-campaign prompt under `docs/campaigns/` — a hypothesis-testing run-set a fresh agent must execute end-to-end without the user. Covers the gated authoring stages (concretise → frame → matrix → research → validate → autonomy-ready), the campaign frontmatter schema, the autonomy-readiness checklist, and regenerating `docs/campaigns/INDEX.md` (per-campaign `x/N` task status). Invoke also when the user asks "where am I across all campaigns?" or to refresh the campaign index.
---

# Campaigns

A **campaign** is one prompt `.md` file under `docs/campaigns/<area>/<campaign>.md`.
It defines a single hypothesis-testing run-set that a **fresh agent must run
autonomously, end-to-end, with no user in the loop**. The folder (`<area>`) is
just a theme grouping; the *file* is the unit that gets a status line.

This skill has **two modes**:

- **(A) Author / advance a campaign** — the gated stages below (a pre-framing
  **Stage 0** for fuzzy ideas + five core stages), run interactively in a grilling
  style (one decision at a time, each locked by the user before advancing).
  Produces an autonomy-ready prompt.
- **(B) Refresh the INDEX** — non-interactive: run `build_index.py`, which scans
  every campaign's frontmatter, recomputes `x/N`, and regenerates
  `docs/campaigns/INDEX.md`. Mode B **also runs automatically at the end of every
  Mode-A session** so the index is never stale.

This skill **composes** other skills by reference — it does not restate them.
The *executing* agent (later, autonomously) invokes them when it reads the prompt:
`tddft-simulations` (run lifecycle/dispatch/Gmail), `simulation-validation`
(tiered pre-run tests), `literature-review` (grounding), `notebook-making` (the
`.ipynb` output contract), `code-test` + `formula-validation` (new code/formulae).

This skill produces a **prompt document only**. It never launches runs, builds
notebooks, or dispatches — that is the executing agent's job.

---

## Campaign frontmatter schema (single source of truth for status)

Every campaign file begins with this YAML. The **executing agent** flips
`done` flags and bumps `status` as it runs; this skill (Mode B) only reads it.

```yaml
---
id: cap-jellium-loss-function           # stable unique kebab slug — NEVER changes (index join key)
area: cap_in_jellium                    # = folder name
title: Loss-function hypothesis check
status: draft                           # draft → ready → running → blocked → paused → done
hypothesis: "Plasmon loss L(q,ω)=|n_q(ω)|²/q² from the E15 run reproduces ..."
handover: docs/handovers/cap-in-jellium-loss-function.md   # pointer (may not exist yet)
tasks:
  - { name: "GS validation", done: false }
  - { name: "E15 2000-au production run", done: false }
  - { name: "loss-function analysis notebook", done: false }
blocked_reason: ""                      # filled ONLY when status: blocked
---
```

- `status` is **agent-set** end-to-end: the runner may flip `running`→`done` and
  set `blocked`/`paused` itself; no user sign-off is required for the status field.
  `blocked` = waiting on a dependency; `paused` = deliberately stopped (both carry
  a `blocked_reason`).
- `x/N` = (count of `tasks` with `done: true`) / (total `tasks`). Every finished
  campaign **must** carry a `tasks:` list or it cannot be tracked.

---

## Mode A — the gated authoring stages (Stage 0 + five core)

Run interactively, grilling style: **one decision at a time, each ending in an
explicit user "locked" before the next stage**. (Mirrors the project principle
"only locked changes are implemented".) Do not advance on your own initiative.

**Campaigns may start fuzzy.** Not every idea arrives frameable — often a
mind-dump needs literature + brainstorming before a falsifiable hypothesis can be
written. **Stage 0 (Concretise)** exists for exactly that. And **research /
grounding is a continuous thread available at *every* stage, not gated to Stage
3** — pull in `literature-review` or a codebase check the moment a gap blocks a
decision, whatever stage you are in.

| Stage | What happens | Deliverable | Compose |
|---|---|---|---|
| **0. Concretise** *(skip if already frameable)* | Turn a mind-dump into something frameable: literature + brainstorming to resolve core unknowns, surface prior art / existing drafts, and judge feasibility. Resolve the easily-answerable questions; defer fine numerics. Produces a **rough draft** a later pass sharpens. | rough campaign draft (`status: draft`) + open-questions list | `literature-review`, `brainstorming` |
| **1. Frame** | State the falsifiable hypothesis, the decision it informs, success/failure criteria, scope boundaries. | frontmatter `hypothesis` + the campaign "Question" | `brainstorming` |
| **2. Brainstorm the matrix** | Concrete run plan: systems, parameter sweeps, phases, observables; crystallise the `tasks:` list. **Build a complexity ladder** — decompose the target into increasingly-complex validated building blocks, simplest first (see below). | run matrix + `tasks:` + `<observables_set>` + complexity ladder | `tddft-simulations`, boundary-rule memory |
| **3. Research & ground** | Research is **continuous** (see the note above this table); this stage is where you **close any remaining gaps and verify** the design — two purposes, do not defer: (a) **close knowledge gaps** the campaign exposes (anything you or the user are unsure of: feasibility, whether the engine supports X, what a parameter should be, whether the approach has precedent), and (b) **verify** every claim, assumption, and number the design rests on (codebase facts via source line-refs; physical/numerical values via `literature-review`). A question that resolves neither a gap nor a verification (e.g. "should I research now or later?") does **not** belong in this stage — just do the research. Output: resolved gaps + verified facts written as justification + new `docs/sources/` entries. | resolved gaps + verified facts + `docs/sources/` entries | `literature-review` |
| **4. Validate** | Pre-run test menu (Tier A/B/C), GS validation, pilot-run criterion, transient-exclusion, guard rails, stop conditions. | `<guard_rails>` + validation plan | `simulation-validation`, `code-test` |
| **5. Autonomy-ready** | Run the autonomy checklist below; only on a full pass flip `status: draft → ready`, and echo the compact `<preflight>` into the prompt. | passed checklist + notebook output contract | `notebook-making` |

At the end of Mode A: run **Mode B** (regenerate `INDEX.md`).

### Explain as you grill (Stages 1–2, and wherever a term is introduced)

A campaign decision is only as good as the user's understanding of it. Whenever a
term, method, or concept enters the conversation that the user may not fully know
— or the user asks — **stop and give a clear, succinct, plain-language
definition before proceeding**. Define the word, say why it matters *for this
decision*, and only then ask them to choose. Contrast options in everyday terms
(e.g. "collinear = spins along one axis, two ↑/↓ channels; non-collinear = spins
in any 3D direction, heavier"). The goal is an **informed** lock, never a
deferential one — the user must understand what they are agreeing to. Prefer a
short definition list over a wall of prose; ground physics terms per
`literature-review` when precision matters.

### Build a complexity ladder (Stages 2 & 4)

A campaign that jumps straight to the full autonomous system is fragile.
Decompose the target into a **ladder of increasingly-complex building blocks**:
start from the simplest controllable sub-system, validate it with a **smoke test
or a single simple-hypothesis check**, then add **one** piece, validate, and
repeat — until the full system is reached and is ready to run autonomously. Each
rung's validation is cheap (a short pilot / smoke run) and gates the next rung;
the eventual autonomous run is trustworthy precisely because every block beneath
it was confirmed in isolation.

Equivalently, a campaign may be **two-phase**: a human-in-the-loop design loop
that climbs the ladder (one experiment at a time, user reviews between rungs),
then an autonomous production phase once the top rung is signed off.

---

## The autonomy-readiness checklist (Stage 5 gate)

**Rule: every box must be answerable from the prompt text alone** — a fresh
agent, no user in the loop, months later. If a box fails, the campaign is not
`ready`.

- [ ] **Self-contained intent** — falsifiable hypothesis + explicit
  success/failure criteria; every `tasks:` entry has an unambiguous done-criterion
  (no "looks reasonable").
- [ ] **Reproducible setup, zero guessing** — geometry/N/r_s/box; GS source
  (named validated checkpoint *or* a GS-validation task-0); propagator + dt +
  duration/steps + energy — all locked *with values and a one-line justification*;
  observable set enumerated per run + cadence; file placement spelled out per
  ADR-0007.
- [ ] **New code pre-gated** — any new observable/kernel routed through
  `code-test` + `formula-validation` + a catalogue row *before* the expensive
  runs, never after.
- [ ] **Validation & guard rails** — pilot-first gate with *numeric* pass
  criteria (for multi-piece systems, a validated complexity ladder underneath);
  abort conditions (NaN / complex energy / GPU occupied); boundary +
  cadence rules (4σ/1σ, 300-frame VTI); PROVISIONAL caveats and open-dependency
  tasks named.
- [ ] **Autonomous mechanics** — GPU scheduling via `cudaMemGetInfo` probe (NVML
  broken; GPU is the default; warn if a GPU is occupied by another user);
  dispatcher concurrency + per-phase Gmail; notebook output contract (per-phase
  and/or final, auto-built via dispatcher / `analyse.py` tail per
  `notebook-making`); handover pointer present; agent updates handover + flips
  frontmatter `done`/`status`.
  - **The autonomous executor MUST be a Python orchestrator, NOT a bash script**
    (user decision 2026-06-27). Bash autonomous dispatchers are brittle — no
    structured error handling, no resume, silent stalls. The Python orchestrator
    gives: structured logging, **idempotent resume** (skip runs whose
    `run_summary` shows `run_completed = true`, so a crash/restart continues),
    **per-phase `try/except` with full-traceback failure emails** (the chain
    continues; one phase's bug never kills the rest), one-shot retry on a sim
    failure, and direct reuse of the Python `analyse.py`/email code. Bash is fine
    only for a single build+run smoke; the *ladder* is orchestrated in Python.
    Reference implementation:
    `ResearchProject/systems/localised_jellium/scripts/campaign_autorun/orchestrate.py`.
- [ ] **Grounding** — every scientific/numerical choice cited
  (`literature-review`) or labelled "Inference:"; engine claims carry source
  line-refs (`inq/...:NN` / `inq-study/...:NN`).

The same checklist, compressed, is echoed into the finished prompt as a
`<preflight>` block (see template) so the executing agent re-verifies it before
burning GPU.

---

## Mode B — refresh the INDEX

```bash
/local/data/public/skcb2/tddft/venv/bin/python3 \
  .claude/skills/campaigns/build_index.py docs/campaigns
```

Writes `docs/campaigns/INDEX.md`: one status-grouped table
(running → blocked → paused → ready → draft → done; ties broken by area), a portfolio
count header, and a "do not hand-edit — regenerated by the campaigns skill"
banner. `template.md` and `INDEX.md` themselves are skipped (no campaign
frontmatter). Markdown only — no CSV (run-level cataloguing is
`tddft-run-catalogue`'s job, a coarser layer than campaigns).

---

## Output artefacts

| Artefact | Location |
|---|---|
| The campaign prompt | `docs/campaigns/<area>/<campaign>.md` (from `template.md`) |
| The generated index | `docs/campaigns/INDEX.md` |
| The canonical skeleton | `docs/campaigns/template.md` |
| The index generator | `.claude/skills/campaigns/build_index.py` (skill-local — skills ship self-contained) |

New terms resolved during a campaign grill go into `CONTEXT.md`; genuinely
hard-to-reverse, surprising trade-offs go into an ADR (offer sparingly).
