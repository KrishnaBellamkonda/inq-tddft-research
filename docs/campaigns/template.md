---
# Campaign frontmatter — single source of truth for the INDEX (build_index.py).
# id NEVER changes once set (it is the index join key). The EXECUTING agent
# flips task `done` flags and bumps `status` as it runs; no user sign-off needed.
id: area-short-slug                       # stable unique kebab slug
area: area_folder_name                    # = the folder this file lives in
title: One-line human title
status: draft                             # draft → ready → running → blocked → paused → done
                                          # blocked = waiting on a dependency; paused = deliberately stopped
hypothesis: "The falsifiable claim this campaign tests, in one sentence."
handover: docs/handovers/<task>.md        # pointer (may not exist yet)
tasks:
  - { name: "Task 1 — short verb phrase with a done-criterion", done: false }
  - { name: "Task 2", done: false }
  - { name: "Task 3", done: false }
blocked_reason: ""                        # filled ONLY when status: blocked
---

# <Campaign title>

<identity>
You are a scientific computing researcher working on first-principles
simulations. You understand the first-principles domain, write scientific-standard
code, and adhere to the rules, principles, and workflows established in this
repository.
</identity>

<description>
<!-- The question and the plan in plain language. State: why this campaign
     exists, what decision its result informs, and the success/failure criteria
     that make the hypothesis falsifiable. Describe the run-set at a high level;
     the locked numbers live in <resolved_decisions>. -->
</description>

<observables_set>
<!-- WHICH observables, for WHICH run, at WHICH cadence. Reference the ADR-0006
     minimal/maximal set. Flag any NEW observable/kernel — it is pre-gated
     (code-test + formula-validation + catalogue row BEFORE the expensive runs). -->
</observables_set>

<resolved_decisions>
<!-- Every choice LOCKED, with a value and a one-line justification; engine
     claims carry source line-refs (inq/...:NN or inq-study/...:NN). Typical
     sub-blocks (see docs/campaigns/cap_in_jellium/baseline_runs.md for a worked
     example): geometry, propagator, duration_and_energy, pilot_and_io, screens,
     file_placement (ADR-0007), observable_enumeration. -->
</resolved_decisions>

<guard_rails>
<!-- Abort conditions (NaN / complex energy / GPU occupied), boundary + cadence
     rules (4σ/1σ, 300-frame VTI), pilot-first numeric gate, PROVISIONAL caveats
     and any open-dependency tasks named. -->
</guard_rails>

<tasks>
<!-- The same task list as the frontmatter, expanded into runnable detail.
     Each task: what to do, its done-criterion, which composed skill it uses
     (tddft-simulations / simulation-validation / literature-review /
     notebook-making / code-test). The agent flips the frontmatter `done` flag
     and updates the handover as each completes. -->
</tasks>

<rules>
<!-- Campaign-specific ALWAYS / NEVER on top of the repo's always-on rules. -->
</rules>

<preflight>
<!-- Compact echo of the autonomy-readiness checklist. The EXECUTING agent
     re-verifies every box from this prompt alone BEFORE burning GPU; if any box
     fails, stop and surface it rather than running. -->
- [ ] Intent self-contained: falsifiable hypothesis + success/failure criteria;
      every task has an unambiguous done-criterion.
- [ ] Setup reproducible, zero guessing: geometry/N/r_s/box; GS source (named
      checkpoint or a GS-validation task-0); propagator + dt + duration/steps +
      energy with values; observables per run + cadence; file placement (ADR-0007).
- [ ] New code pre-gated: any new observable/kernel → code-test +
      formula-validation + catalogue row BEFORE expensive runs.
- [ ] Validation & guard rails: pilot-first numeric gate; abort conditions;
      boundary + cadence rules; PROVISIONAL caveats + open deps named.
- [ ] Autonomous mechanics: GPU via cudaMemGetInfo probe (NVML broken; GPU is
      default; warn if occupied); dispatcher concurrency + per-phase Gmail;
      notebook output contract (auto-built via dispatcher / analyse.py tail);
      handover pointer present; agent updates handover + frontmatter done/status.
- [ ] Grounding: every scientific/numerical choice cited or labelled "Inference:";
      engine claims carry source line-refs.
</preflight>
