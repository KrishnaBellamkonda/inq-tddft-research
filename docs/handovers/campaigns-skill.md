# Handover: `campaigns` skill — autonomous-campaign authoring + status tracking

Rolling handover. Plan: `/home/raid/skcb2/skcb2/tddft/.claude/plans/snazzy-wiggling-naur.md`.
Glossary: "Campaigns (2026-06-22)" section of `/local/data/public/skcb2/tddft/CONTEXT.md`.

## Milestone 2026-06-22 — skill built, all campaigns backfilled, INDEX live

### What this is
A new `campaigns` skill: a front-door process for authoring **autonomous research
campaign prompts** (each a `.md` a fresh agent runs end-to-end, no user in the
loop) AND tracking every campaign's `x/N` task status in a regenerated INDEX.
Design was grilled + locked with the user (`/grill-with-docs`).

### DONE (verified)
- **Skill:** `/local/data/public/skcb2/tddft/.claude/skills/campaigns/SKILL.md`
  — dual-mode (A: 5 gated authoring stages Frame→Matrix→Research→Validate→
  Autonomy-ready; B: non-interactive INDEX refresh), frontmatter schema, the
  autonomy-readiness checklist.
- **INDEX generator (skill-local, shippable):**
  `.claude/skills/campaigns/build_index.py` — stdlib-only frontmatter parser
  (no PyYAML in venv), status-grouped table, portfolio header, do-not-edit banner.
  Skips `INDEX.md`/`template.md`/any file without an `id`.
- **Test:** `.claude/skills/campaigns/test_build_index.py` — standalone (no
  pytest). **VERIFIED PASSING** with the venv python. Covers flow+block task
  parsing, x/N, status order (running→blocked→paused→ready→draft→done), skip
  rules, handover-link rewrite, blocked/paused reason rendering.
- **Status enum** includes `paused` (deliberately stopped) distinct from
  `blocked` (waiting on a dependency). Propagated to SKILL.md, template.md,
  build_index.py, CONTEXT.md.
- **Rename:** `docs/prompts/` → `docs/campaigns/` (`git mv`). Live refs updated
  ONLY (CONTEXT.md ×3, `ResearchProject/systems/localised_jellium/shared/configs/slab_n234_L50.hpp`,
  `ResearchProject/systems/vacuum/hypotheses/twosided_cap_vs_mask/build_twosided_report.py`).
  Historical handovers/plans intentionally left pointing at old path.
- **Canonical skeleton:** `docs/campaigns/template.md` rewritten (frontmatter +
  tag anatomy from `baseline_runs.md` + `<preflight>` echo of the checklist).
- **Backfill (user-adjudicated statuses):** 14 existing campaign files got
  frontmatter prepended; 3 retroactive files created for previously-untracked
  work. **17 campaigns total.** `docs/campaigns/INDEX.md` generated & verified.
  - Retired (excluded, no frontmatter): `codebase_rejuvination/codebase_rejuvination_complete.md`
    (umbrella) and `localised_jellium/localised_jellium.md` (stub).
  - New area `jellium_stopping/`: `sv_sigma0p5_classical.md` (done),
    `sigma_convergence_sweep.md` (running). New `quantum_classical_nocap/run_withcap_sigma3.md`
    (Study A, running) beside its no-CAP twin.
- Portfolio at backfill: **5 running · 1 paused · 2 ready · 4 draft · 5 done**.

### Key correction captured during the grill
`ab-sigma0p5-baselines` is a **vacuum absorber** baseline (now UNBLOCKED — its
predecessor two-sided CAP-vs-mask sweep finished 2026-06-17 — → `ready`). The
"σ=0.5 WP+classical through jellium at energies" the user remembered is a
**different** campaign = `jstop-sv-sigma0p5-classical` (branch
`overnight-gaussian-classical`, **done**). The "classical blocker" was a
**cost/GPU-budget block** on the WP/loss-function *production* (~20k steps) +
deferred k-points + noisy low-v tail — NOT a crash.

### NOT done / open
- **4 target campaigns are `draft` with empty `tasks`** — to be authored next via
  the skill (one at a time, user's stated goal):
  `cap_in_jellium/loss_function_hypothesis_checking.md`,
  `cap_in_jellium/jellium_classical_vs_wavepacket_cap.md`,
  `td-hf/checking-wp-hf-orbital-approximation.md`,
  `ml-patterns/pattern-finding-in-wp-classical-runs.md`.
- **Untracked, needs user call:** `ResearchProject/systems/jellium/hypotheses/00_jellium_reference/`
  (2026-06-17) — may be the S(v) Lindhard reference or a standalone campaign; no
  prompt file created yet.
- Several `done`/task-level flags were the user's best recall (e.g.
  `rejuv-unit-testing` 1/5, `locjel-campaign` 2/4, graphene 2/6) — refine when
  the underlying handovers are next touched.
- Nothing committed (user has not asked). One-off backfill script at
  `/tmp/backfill_campaigns.py` (throwaway, not tracked).

### How to refresh the INDEX (anytime)
```bash
cd /local/data/public/skcb2/tddft
venv/bin/python3 .claude/skills/campaigns/build_index.py docs/campaigns
```

### Memory written
`feedback_skills_self_contained_shippable` — skills bundle ALL artefacts
skill-locally (drove `build_index.py` placement).

## Milestone 2026-06-22 (cont.) — first campaign authored via the skill: `tdhf-wp-orbital-approx`

Drove the 5-stage Mode-A flow for
`docs/campaigns/td-hf/checking-wp-hf-orbital-approximation.md` (now **fully
authored, status `blocked`, 0/7**). Highlights:
- **Hypothesis:** KS-WP orbital ≈ HF orbital; gap = SIE / physical-interpretability
  bound. Metrics LOCKED {fidelity F(t), density L1/L2, centroid + momentum-loss};
  **thresholds DEFERRED** (verdict descriptive until pinned).
- **Method:** 3 theory arms (full-LDA / exchange-only-LDA / HF) of the *identical*
  WP-injection setup; pairwise diffs decompose SIE vs correlation vs total.
  3 phases: A free (vacuum), B small high-density jellium slab (r_s≈4), C coronene.
- **Stage-3 verified (source/lit):** RT-TDHF supported in inq-study (spin-pol exact
  exchange ✓, WP injection theory-agnostic ✓, ACE-accelerated) — but **must use
  Crank-Nicolson** (ETRS asserts no exact exchange, etrs.hpp:26) and there is **no
  in-repo precedent** for HF+CN+WP ⇒ pilot is a hard gate. New memory
  `reference_inq_rt_tdhf_requires_crank_nicolson`. ε_x=−0.4582/r_s verified
  (−0.9163 is Rydberg); KS≈HF is known *statically* — the dynamic injected-orbital
  regime is the campaign's actual contribution.
- **Runs placed in home systems** (not a new system): Phase A→`systems/vacuum/
  tdhf_free/`, B→`systems/localised_jellium/tdhf_slab/`, C→`systems/coronene/
  tdhf_coronene/`; new kernel `orbital_fidelity`→`inqview/analysis/`; cross-system
  synthesis notebook→`docs/reports/td-hf-orbital/`.
- **Blocked because:** Phase B/C gated on `locjel` slab validation + the TD-HF
  pilot; thresholds deferred. Phase A + the kernel + pilot are runnable now.

### Skill refinements made this session (from user feedback)
- SKILL.md: added "**Explain as you grill**" guidance (give clear plain-language
  definitions of unfamiliar terms before asking the user to decide).
- SKILL.md: sharpened the **Stage 3** purpose — it exists ONLY to (a) close
  knowledge gaps and (b) verify claims/feasibility; process-logistics questions
  ("research now or later?") don't belong there.
- Added `paused` to the status enum (graphene); `build_index.py` renders
  blocked/paused reasons; tests updated + green.
- ADR `docs/adr/0009-campaign-tracking-convention.md` written.

### Next
User wants to author the remaining target campaigns:
`cap_in_jellium/loss_function_hypothesis_checking.md`,
`cap_in_jellium/jellium_classical_vs_wavepacket_cap.md`,
`ml-patterns/pattern-finding-in-wp-classical-runs.md` (all still `draft`, 0 tasks).
