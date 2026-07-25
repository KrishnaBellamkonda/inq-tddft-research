# Plan: organise the `campaign_autorun` run-set for review

**Task owner:** user (chiddukanna). **Created:** 2026-07-06.
**Scope decided via `/grill-with-docs` interview (2026-07-06).**

## Goal

Organise every notebook of the localised-jellium **`campaign_autorun`** run-set so
the user can read each result and independently confirm it. For each notebook state
**(1) the question the user was aiming to answer, (2) what was done, (3) the results**
— shown, not interpreted.

## Scope (locked)

- **Only** `campaign_autorun` — the localised-jellium GS ladder H0–H5
  (`ResearchProject/systems/localised_jellium/scripts/campaign_autorun/`, analysed in
  `.../hypotheses/campaign_autorun_study/`). NOT qsp_phase*, wide_wp, ml-patterns, etc.
- "campaign_autorun" is a **run-SET**, not a `docs/campaigns/` "campaign"
  (glossary clash recorded in `CONTEXT.md`, 2026-07-06).
- The run-set = ~90 runs: H0 (12 single-point wp/cl, r-sweep), H1 (8 GS, edge-w),
  H2 (11 GS, Lz + 2 open-z), H3 (9 GS, thickness), H4 (26 single-point WP p2/p3),
  H5 (26 single-point classical p2/p3). **No trajectories** — GS single points, or
  3-step frozen single points (a 4-row `observables.csv`, no density VTI).

## Hard constraints (user)

1. **No interpretation by the assistant.** Deriving learnings / further experiments
   is the user's job. Assistant presents question → method → results only.
2. Existing interpretive prose (per-hyp `Takeaway`, campaign verdict/conclusions/
   follow-ups) is **kept but quarantined** in a marked
   *"⚠ Provisional — author-generated, you own the verdict"* box, visually separated.
3. Every "question you were aiming to answer" is sourced from the user's **own
   wording** (the existing `hyp` field per hypothesis), never invented.
4. Figures travel **beside** the notebook (Jupyter won't serve outside its tree).
   `.png` only. VTIs read via `inqview.load_vti` (physical order — **no fftshift**).
5. venv python: `/local/data/public/skcb2/tddft/venv/bin/python3`. No INQ runs / GPU —
   analysis reads existing results only.

## Deliverables (all in `.../hypotheses/campaign_autorun_study/`)

### 1. Study notebooks H0–H5 — restructured (`build_notebooks.py`)
Fixed neutral order per notebook:
1. **Question you were aiming to answer** ← existing `hyp`
2. **What was done** ← existing `setup` + `method` (verified vs `run.cpp`)
3. **Results** ← existing plot + recomputed printed numbers (no verdict prose)
4. **⚠ Provisional — you own the verdict** ← existing `take`

### 2. `00_index.ipynb` — single canonical entry point (`build_notebooks.py`)
H0→H5 ladder as a table `[question | what was done | results/runs location | links to
study + run-evidence notebooks]`; one merged provisional box (old verdict + conclusions
+ follow-ups). Old `campaign_summary.ipynb` stays on disk, **dropped from the reading
path** (not deleted).

### 3. Six per-hypothesis run-evidence notebooks — new `build_run_evidence.py`
`runs/H0_runs.ipynb … H5_runs.ipynb`. One per hypothesis, listing **every run in the
sweep**: `run_summary.txt` config, converged/step-0 energy (E_tot, and E_tot(0)−E_GS
for H0/H4/H5), density slice/profile for GS runs — so each data point behind a
hypothesis plot is independently checkable. Built by a small aggregator (the single-run
assembler cannot span a sweep).

### 4. Four representative single-run deep-dives — the run-notebook assembler
Via `.claude/skills/run-notebook/run_notebook_builder.py` (extended, see below):
- `rep_H0_wp_vs_cl_r28` (single-point WP + classical at r28_p3 — the E_tot(0) gap)
- `rep_H2_gs_lz120` (reference full-PBC GS box)
- `rep_H4_wp_r28_p2_vs_p3` (open-z vs full-PBC WP)
- `rep_H3_gs_a15_N98` (a thickness GS point)

## Builder change (additive, shared skill-local file)
`run_notebook_builder.py` gains, **additively** (no change to existing trajectory
behaviour; also benefits jellium GS baselines):
- **GS-density fallback:** when `raw/vti/density_system/` is absent but
  `density_gs_system/*.vti` exists, render an xz density slice + n(z) profile via
  `inqview.load_vti` (physical order, no fftshift).
- **Single-point energy table:** when `observables.csv` is short (< ~8 rows), render a
  step-0 energy-decomposition table (E_tot, T, U_H, E_xc) and E_tot(0)−E_GS when an
  `--e-gs-ha` is supplied.

## Execution & validation
- Rebuild via the builders; **re-execute** with `python3 -m nbconvert --to notebook
  --execute --inplace` (module form — avoids the pyenv `jupyter` shim). 0 errors.
- Sanity: recomputed numbers match embedded plots (numbers read from run files, never
  re-converged); density orientation correct (slab at box centre, `expect_centered_axis`).
- No fabricated panels; sections auto-skip when an observable is absent.

## Out of scope / not doing
- No new INQ runs. No deletion of `campaign_summary.ipynb`. No assistant-authored
  learnings/next-steps (user owns those). No touching qsp_phase*/wide_wp/ml-patterns.
