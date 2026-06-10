# 4. inqview has one canonical, library-wide visualisation theme applied via semantic roles

Date: 2026-06-10
Status: accepted

## Context

inqview carried two disagreeing default systems:

- `config.py` / `defaults.py` — the nominal public theme: `cividis`
  sequential, `coolwarm` signed, `twilight_shifted` phase, figsize
  `(6.4, 4.2)`, dpi 160. Generic, never the designed report standard.
- `report1/_shared_style.py` — the carefully designed publication standard:
  `inferno` sequential, `RdBu_r` (zero-centred) diverging, a **fixed-column
  size scheme** (`ONE_COL_IN = (3.5, 3.0)`, `TWO_COL_W_IN = 7.0`) with a
  fixed axes rectangle so panels align, 10 pt fonts, save dpi 600, plus a
  `TufteCritic` save hook.

Meanwhile the pipeline phases used ad-hoc styling of their own (e.g.
`screens.py` hard-codes `viridis` and figsize `(5, 5)`). The result is three
visual languages in one library that is about to be released.

The user directed that the report1 standard (and the `report-figures` skill
it derives from) are "the defaults that are required to be made" — and that
this is a change to be **made and tested**, because the plot proportions are
deliberately designed.

## Decision

- Promote the report1 standard into **`inqview.visualisation.style`** as the
  one canonical theme; the generic `config.py`/`defaults.py` values are
  superseded.
- Apply it **library-wide via semantic cmap roles**, not literal cmaps:
  - `sequential → inferno` (positive-magnitude maps, diffraction intensity)
  - `diverging  → RdBu_r`, zero-centred (signed Δn, difference maps)
  - `phase      → twilight_shifted` (cyclic data)
  Every phase calls `cmap_for(role)`; no phase hard-codes a cmap name.
- Provide a **fixed-dimension figure factory**: `figure_one_col()` →
  3.5×3.0 in with a fixed axes rectangle, `figure_two_col()` → 7.0 in wide.
  The library emits **individual plots only**; panel composition stays a
  downstream LaTeX concern (per the `report-figures` skill).
- The theme is a **tested deliverable** (pure tier): assert `apply_theme()`
  sets the documented rcParams, `cmap_for(role)` returns the designed
  values, and the figure factory produces the exact geometry — including a
  regression guard that tight-bbox / constrained-layout do not silently
  override the fixed width (see memory
  `reference_fixed_dimension_plot_pitfalls`).

## Consequences

- One visual identity across publication figures, quicklook plots, and
  pipeline-phase GIFs.
- The cmap choice becomes centralised and testable; restyling the whole
  library is a one-line change to a role map.
- Migration cost: every phase that hard-codes a cmap/size is rewritten to
  request a role + use the figure factory. Guarded by the characterization
  tests during the restructure.
- Trade-off: less per-phase freedom. Accepted; a phase with a genuine domain
  reason still picks a *role*, and new roles can be added centrally rather
  than ad hoc.
- `TufteCritic` (the report1 save-time critique hook) is a candidate to move
  into `inqview.visualisation` as an optional save wrapper; deferred to the
  migration, not part of this decision.
