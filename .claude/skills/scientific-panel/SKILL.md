---
name: scientific-panel
description: Use when the user wants to deliberate on the INTERPRETATION of a hard scientific result — competing explanations for a puzzling simulation or experimental outcome, "what is really happening / why is S(v) flat / is this physics or an artifact" questions, or a stress-test of a physical interpretation that deserves more than one viewpoint. Convenes four INDEPENDENT domain-expert subagents (TDDFT + stopping-power methodologist; broad condensed-matter theorist; jellium-stopping specialist; run-data custodian) who each reason from the simplest reduced model upward, debate and rebut each other, then a judge synthesises a single grounded, provisional verdict — consensus, the live disagreements, and the decisive next test — for the user to push back on. Invoke whenever a result needs multiple expert perspectives, a genuine second opinion, or adversarial deliberation rather than one answer, and whenever the user says "convene the panel", "get the experts", "debate this", "what do the experts think", or asks you to deliberate on simulation results.
---

# Scientific panel

A structured way to get a **hard scientific question deliberated by four independent
experts and adjudicated by a judge**, instead of answered by a single model voice.
The point is *diversity of reasoning*: four subagents, each with a different lens,
reach their own conclusions before they see each other's, argue, and then a judge
weighs everything and hands the user a grounded but **provisional** verdict to push
back on.

Use it when an answer matters enough to be worth stress-testing — a puzzling
result, an interpretation that "feels right but might be an artifact", a fork in
the physics. Do **not** use it for routine work (running an analysis, extracting a
number, writing a plot) — that is what the ordinary skills are for. A panel that
deliberates over a one-step task is just expensive.

## Why a skill, not an existing command

No existing slash-command convenes named domain experts and a judge over a physics
question. This skill **reuses existing primitives** rather than inventing new ones:
the `Agent` tool spawns each expert as a real, context-isolated subagent (genuine
independence); the fan-out follows the `superpowers:dispatching-parallel-agents`
pattern; and the structure mirrors the `Workflow` tool's documented *judge-panel*
quality pattern (independent attempts → adjudication → synthesis). It runs in the
**foreground** so the user stays in the loop turn-by-turn, which a background
`Workflow` would not allow.

## The panel (default roster)

Four experts + a judge. The verbatim, copy-ready spawn prompts live in
`references/panel-roster.md` — **read that file when you spawn**, and paste the
shared evidence brief into each `<<BRIEF>>` slot. The charters below are just the
quick reference.

| Agent | Lens it brings | The distinctive challenge it presses |
|---|---|---|
| **E1 — TDDFT / stopping-power methodologist** | what the calculation actually computes; energy-ledger vs force vs momentum observables; CAP/absorber behaviour; self-interaction error; grid / dt / basis convergence | "is this number an *artifact of the method* — reference, absorber, grid, SIE — or physics?" |
| **E2 — condensed-matter generalist** | many-body screening, dielectric response, plasmons, e–h pairs, sum rules, Fermi-liquid intuition, EELS analogies | "what do linear response and the sum rules *demand*, and is the result consistent with them?" |
| **E3 — jellium-stopping specialist** | Lindhard/RPA stopping, the Bragg peak near v_F, the high-v Bethe ln(v)/v² tail, low-v friction Q, nonlinear screening (Echenique), Barkas/Z effects, ion-vs-electron projectile, finite-slab geometry | "how does this sit against the canonical jellium stopping curve and its known limits?" |
| **E4 — run-data custodian** | has read every run spec, handover, campaign doc, and results file; grounds every claim in the real numbers and flags partial / unconverged / aliased / mis-referenced data | "what do the *files actually say*, and which explanation is consistent with ALL the runs, not just the convenient one?" |
| **Judge (the bench)** | did not debate; weighs all four views | parsimony + rigor + decision-usefulness; actively hunts for groupthink |

The roster is the default for this project's TDDFT/jellium physics. For a question
in another domain, the convener may **swap the specific expertise** of E1–E3 to fit
(e.g. a surface-scattering specialist for a coronene LEED question) — but keep four
distinct lenses and **always keep E4** (the data anchor), since a debate ungrounded
in the actual runs drifts into hand-waving.

## Models, effort, and engine

Every expert and the judge run on the **opus** tier at **high** effort. The
deliberation is the part that most rewards a strong model and deep reasoning, and
offloading it to subagents keeps the heavy thinking **out of the main conversation**
— only the judge's verdict (and, if asked, short expert digests) return, so a long
back-and-forth doesn't bloat the user's context. (Per-subagent *minor* version
pinning — 4.7 vs 4.8 — is not exposed by the spawn tools; "opus" resolves to the
session's opus tier. Set the session model if an exact minor matters.)

**Default engine — the Workflow tool.** Run the openings → rebuttal → judge pipeline
as one Workflow, because `agent(prompt, {model:'opus', effort:'high'})` is the only
spawn path that pins BOTH model and effort per subagent, and it is the most
token-frugal (expert transcripts stay in the workflow; the script returns just the
judge verdict). The push-back loop = re-invoke with the user's new input folded into
the brief. Watch live with `/workflows`. A ready, parameterised script is in
`references/panel-workflow.md` — copy it, fill the brief/question/file-pointer slots,
and launch.

**Interactive fallback — the Agent tool.** When the user wants to watch the debate
unfold turn-by-turn or intervene mid-stream, spawn the experts in the foreground with
`Agent(..., {model:'opus'})` (the Agent tool exposes model but not effort, so effort
follows the session).

## Protocol

You (the convener / main agent) run this. Default is **opening → rebuttal →
judge**. Spawn the four experts **in parallel** (one message, four `Agent` calls).

**Round 0 — Build the evidence brief.** Assemble a single shared brief (see next
section). This is the common ground every expert receives. **Anchoring guard:** the
brief presents *evidence*, not your favoured conclusion. If you already have a
leading hypothesis (you usually will), keep it out of the brief — the value of the
panel evaporates if you seed all four with your prior.

**Round 1 — Opening positions (parallel, 4 experts).** Spawn E1–E4 with their
roster prompts + the brief. Each returns: the simplest reduced model + the number
it predicts, its leading explanation, the one discriminating test, its confidence,
and "what would change my mind." They do **not** see each other yet. Give E4 the
file pointers and tell it to actually read them.

**Round 2 — Rebuttal (4 experts).** Give each expert the other three openings and
ask it to: name the single strongest point against its view, concede what it
over-claimed, and sharpen or revise. Preserve each expert's reasoning by feeding it
**its own opening plus the other three** (Workflow: pass them into the round-2
prompt; Agent-tool fallback: `SendMessage` to continue the same agent — find it via
`ToolSearch` `select:SendMessage`). Forbid manufactured disagreement *and* agreeable
collapse: agreeing for a stated reason is fine; agreeing to be agreeable is not.

**Round 3 — Judge synthesis (1 judge).** Spawn the judge with the question, the
brief, and the full transcript (openings + rebuttals). It returns the structured
verdict (format fixed in the roster). Relay that verdict to the user **verbatim or
lightly framed** — do not overwrite it with your own opinion.

**Round 4 — User in the loop.** Present the verdict and **stop**. The user examines,
pushes back, adds thoughts. That is the deliverable — a provisional verdict the user
owns, not a closed answer.

## Building the evidence brief

A good brief is the highest-leverage part. Include, tightly:

- **The question**, stated sharply (one or two lines).
- **The cleaned data** — the actual numbers under discussion, as a small table, with
  units and which points are trusted vs excluded (and why).
- **The run specifics that matter** — system, method, grid/box/dt, σ, CAP, GS energy,
  convergence flags — only what bears on the question.
- **Known reference values / limits** — e.g. Lindhard S(v), v_F, sum rules — so the
  experts argue against anchors, not from memory.
- **File pointers** for E4 — handover, results CSVs, campaign docs, run dirs.
- **No conclusion.** State the puzzle, not the answer.

Keep it skimmable. Every expert reads it; bloat taxes all four.

## Resuming after the user pushes back

When the user replies with pushback or new thoughts, **re-deliberate** — don't just
answer yourself. Fold their points into an updated brief (a "Round N" addendum:
"the user observes X; reconsider in that light"), then continue the *same* experts
via `SendMessage` (so they remember the prior round) for a focused rebuttal, and
re-run the judge. If the pushback only touches one expert's domain, you may consult
just that expert and the judge. The loop continues until the user is satisfied or
calls it.

## Scaling knobs

- **Quick** (≈5 agents): openings + judge, skip rebuttal. For a fast second opinion.
- **Default** (≈9 agents): openings + rebuttal + judge. Use this unless told otherwise.
- **Deep** (≈12+ agents): add a second rebuttal round and/or a 2–3 member judge bench
  (split rigor / parsimony / decision-usefulness across judges, then merge). For
  high-stakes or genuinely unresolved questions.

Match the depth to the stakes. Say which mode you ran.

## Guardrails

- **Independence first.** Never paste your own conclusion into the brief or the
  expert prompts. The experts must be free to disagree with you and each other.
- **Reduced models first.** Every expert starts from the simplest model that could
  explain the observation and builds up only as it fails — this is the panel's house
  style and the user's stated preference. A clean simple model the whole panel can
  follow beats an elaborate one only its author understands.
- **Ground claims.** Cite known results and limits; separate evidence from inference
  ("Inference:"); carry units; round reported numbers to 2 significant figures
  (project rule). E4 grounds in the actual files, not memory.
- **Hunt groupthink.** If all four converge without stress-testing, the judge says so
  and treats the consensus with suspicion; the convener may add a devil's-advocate
  pass.
- **The user owns the verdict.** The judge gives a *provisional, best-grounded*
  recommendation framed for the user to ratify, refine, or reject — never a closed
  ruling. (Consistent with how "verify"-style requests are handled in this project:
  evidence and a recommendation, the human decides.)

## Spawn prompts

`references/panel-roster.md` holds the full, paste-ready prompts for E1–E4, the
rebuttal addendum, and the judge. Read it at spawn time and fill the `<<BRIEF>>`,
`<<OPENINGS>>`, `<<QUESTION>>`, and `<<TRANSCRIPT>>` slots. `references/panel-workflow.md`
holds a ready Workflow script that wires the whole openings → rebuttal → judge
pipeline with `model:'opus', effort:'high'` agents — the token-optimal default.
