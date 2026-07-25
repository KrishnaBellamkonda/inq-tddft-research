---
name: scientific-figures
description: Use when producing ANY scientific figure for a presentation or report — system-design plots, S(v)/spectra, density maps, reflectivity curves, tables, workflow diagrams. Encodes the project's hard figure rules (density-based system design, minimal legends, colorbar geometry, titles/captions, .drawio workflows, table styling). Layers on report-figures (canonical theme) and is referenced by notebook-making.
---

# Scientific figures — the rule set for highest-quality figures

This is the **single source of truth** for figure rules. It LAYERS ON the
canonical theme + production workflow in `report-figures` (read that for
`inqview.visualisation.style`, the shared-colorbar rule, scientific-offset
colorbar ticks, and the linear+log density-panel rule). `notebook-making` and
`report-figures` both point here. Everything below is mandatory unless a rule
says "report-only" or "presentation-only".

Always: venv python (`/local/data/public/skcb2/tddft/venv/bin/python3`), canonical
theme (no ad-hoc rcParams), `.png` only at 600 DPI, individual figures (never
pre-composed panels — compose in LaTeX/slides). Never preview/Read the generated
PNG yourself — the user previews.

## 1. Presentation mode vs report mode

A figure is built for exactly one target; ask/know which.

| | **Presentation** | **Report** |
|---|---|---|
| On-canvas **title** | **PRESENT** (short, on the figure) | **CROPPED** (no title; the LaTeX `\caption` carries it) |
| **Caption** | lives in the **slide spec** (`caption:` line per slide), spoken/shown on the slide — NOT on canvas | LaTeX `\caption{}` |
| Font scaling | for on-slide downscale (see report-figures) | for column width |

Default for the Emilio deck = **presentation** (titles present; captions in
`03_nuts_and_bolts_plan.md`).

**Saving (FB-001 — the recurring clip bug).** The canonical `figure_one_col()`
pins a FIXED axes rect with ~zero headroom — a **report** factory. A presentation
figure that adds an on-canvas `set_title` or an external (`make_axes_locatable`)
colorbar will have them **clipped off-canvas** if saved with a bare
`savefig(dpi=…)`. Save every presentation figure with
**`style.save_presentation(fig, path)`** (`bbox_inches="tight"`); it sacrifices
the fixed on-page size, which is fine because a deck contain-fits each image.
Never save a titled / external-colorbar figure with the bare report save.

## 2. Nothing analytical on the canvas

- **No explanatory/analytical text** drawn on the plot ("Coulomb → ∞", "finite
  core", "converges to linear response", derivations). That is caption/slide
  content. The canvas shows data, axes, a short title (presentation), and at most
  a **short quantitative** label in dead whitespace (`ω_p = 3.47 eV`).
- **No arrowheads / leader lines** pointing at a data feature to "explain" it. If
  a feature needs calling out, it goes in the caption. (Markers/guide lines that
  ARE data — a dashed threshold, a slab edge — are fine; see §4, §5.)
  **Scope:** this bans leader/analytical arrows on **data plots only**. Flow /
  workflow diagrams (§7) are directional by nature — their connector arrows are
  REQUIRED, not banned (user, 2026-06-25).

## 3. Legends say *what the data is* — nothing else

- A legend entry names the **identity of a series** only: "linear response",
  "σ_WP = 0.5", "two-sided CAP". 
- **Forbidden in legends:** parameter dumps, method notes, definitions, run IDs,
  σ-convention asides, anything explanatory. Those go to the slide/caption.
- If a legend would need a title/subtitle to be understood, the information is
  in the wrong place — move it to the slide.

## 4. System-design plots are built from REAL density data — never cartoons

The system-design figure exists to give **confidence that everything is where it
should be**. Do NOT hand-draw boxes/atoms/arrows as a schematic. Instead:

- Plot the run's **actual density** — a **total-density** xz slice (preferred) or
  a **Δn** xz slice — loaded via `inqview.load_vti` (physical order; NEVER
  `np.fft.fftshift` a VTI; pass `expect_centered_axis="z"` for centred slabs).
- Overlay **dashed lines** marking the **slab extent** and the **CAP extent**
  (and box edges from the data axes). These are the only annotations.
- For a **bulk** run (no localised slab), omit the slab lines; show total density
  + dashed CAP extent + the projectile launch/trajectory line. The medium fills
  the box.
- **No legend of run parameters** (N, dx, r_s, v, T…). Those belong on the slide.
  The figure proves placement; the slide states the numbers.

This is the core philosophy: the geometry is read off the data, so it cannot be
mis-drawn, and the audience trusts the result.

## 5. Colorbars: outside the axes, same height as the panel

- The colorbar is placed **outside** the plotting axes and is the **same length
  (height) as the plotted square/panel** — never taller, never shorter, never
  overlapping the data. Use an axes-divider so the bar tracks the panel height:
  ```python
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  div = make_axes_locatable(ax); cax = div.append_axes("right", size="4%", pad=0.08)
  cb = fig.colorbar(im, cax=cax)
  ```
  (`make_axes_locatable` keeps `cax` at the panel height automatically; for a
  fixed-aspect `imshow` prefer `inset_axes`/`ImageGrid`.)
- Plus the report-figures rules: scientific-offset ticks (one shared `×10ⁿ`, ≤2
  sig-figs, ≤5 ticks), shared fixed clim across compared panels/frames.

### GIFs: colorbar AND gradations LOCKED across every frame (user, 2026-07-03)

For any GIF/animation the colorbar range (`clim`/`vmin`,`vmax`) **and its
gradations** (tick levels, the linear-vs-log scale, the shared `×10ⁿ` offset) must
be **computed ONCE over the whole frame stack and held FIXED for every frame** —
never recomputed per frame. A colorbar that rescales frame-to-frame makes frames
non-comparable: a feature that merely tracks the moving scale looks like it is
growing/shrinking. This is a **correctness** bug, not cosmetics.

- Compute global limits from the full time series up front
  (`vmin, vmax = stack.min(), stack.max()`; a **shared symmetric** clim for Δ
  GIFs), build **one** `Normalize`/`LogNorm` + tick locator, and pass the SAME
  objects to every frame. Draw the colorbar ONCE with fixed ticks — do not redraw
  it per frame.
- Density GIFs → one shared **log** clim across total/bath frames; Δ GIFs → one
  shared **symmetric diverging** clim. Fix once via `wake.shared_clim`, never
  per-frame.
- This extends the shared-clim rule to the **time axis**: identical clim + ticks
  across frames, not merely across side-by-side panels. Verify by eye that the
  colorbar labels are byte-identical on the first and last frame.

## 6. Tables: only the header row is coloured

- Header row: one accent fill + bold. **Body rows: no fill** (plain white, thin
  rules). No zebra striping, no per-cell colour. Keep it clean.

## 7. Workflow diagrams: Graphviz/Dot (diagram-as-code) — preferred

For a workflow/pipeline diagram, the **preferred, reproducible** technique is
**Graphviz/Dot** (FB-010, user decision 2026-06-25): write a text `.dot` source
(git-diffable) and render with the system `dot` (`dot -Tpng -Gdpi=300 …`). The
engine auto-routes clean orthogonal connectors, so the diagram re-lays-out from
the edited source — no pixel-pushing. Canonical builder:
`docs/reports/26-06-2026-meeting-emilio/build/build_workflows.py`.

- **Differentiate node types by role** with fill colour — input / process /
  output — never one uniform neutral box.
- **Directional arrows are kept** (the §2 no-arrowhead rule is data-plot-only).
- Keep labels plain-ASCII (Graphviz fonts don't render mathtext/exotic glyphs).
- Have the producing section script **NOT** also emit a matplotlib "cartoon"
  under the same filename (comment out that call) so the Graphviz PNG persists.

(Alternative when hand-layout is needed: the older `.drawio` source + matched
matplotlib PNG pattern — `docs/diagrams/build_contribution_page.py` emits the
`.drawio`; `render_contribution_png.py` renders the matched PNG; trimmed template
bundled here as `workflow_render_template.py`):

1. **`.drawio` source** — programmatically emit `<mxCell>` boxes/edges (or
   splice a page into an existing `.drawio`), grounded box labels = real
   module/function names from the codebase, MathJax equations (`math="1"`). This
   is the editable source of truth.
2. **Matched PNG** — reproduce the *same* node/edge layout in matplotlib at
   **1920×1080** in the report idiom: Carlito/Calibri sans, **muted restrained
   palette**, rounded `FancyBboxPatch` boxes, swimlane **containers**, **orthogonal**
   `FancyArrowPatch` edges, equations, balanced full-canvas layout. Large fonts
   (on-slide downscale ≈ ×0.69).

The PNG is NOT a cartoon precisely because it is laid out to this idiom. If a
diagram looks hand-doodled (freehand arrows, clip-art, uneven boxes, cramped
text), it is wrong — rebuild it.

## 8. Split crowded figures

If labels/numbers/series overlap or pile up, **split into separate individual
plots** (one quantity per figure) rather than cramming. Composition happens in
the slide/LaTeX, not by overloading one canvas.

## 9. Terminology

- The linear-response / Lindhard / RPA reference curve is **labelled "linear
  response"** on every figure and in prose (project terminology, see `CONTEXT.md`).
  The code module name (`lindhard_elf`, `lindhard`) is unchanged; only the
  displayed label is "linear response".

## 10. Numbers, labels, lineshape, provenance (Emilio-deck pass, 2026-06-26)

- **Scientific notation = `×10ⁿ`, never `1e-03`.** Any number that lands on a
  plot (offset text, tick, mathtext label) must render as `m×10ⁿ`. Canonical fix:
  `axes.formatter.use_mathtext=True` (in `style.apply_theme`) + `style.sci_notation(x)`
  for numbers embedded in mathtext labels. `print()` to stdout is exempt.
- **No exotic Unicode in saved labels.** Glyphs like `⟂`, `⟨⟩` fail to render in
  saved PNG/GIF labels. Use mathtext or a plain symbol (`n_WP`, not `WP density (⟂ sum)`).
- **Verify the functional SHAPE against the engine source, not just topology.** A
  plot can be topologically right yet physically misshapen (the CAP was two-sided
  but each band peaked at the wall instead of its centre). Match the formula to
  the engine (`absorbing.hpp`), not to a re-derived assumption.
- **Signal-vs-baseline:** plot the **delta from t=0** (`ΔE(t)=E(t)−E(0)`), not the
  absolute quantity, when the signal rides on a large offset.
- **Data provenance is part of correctness.** A figure must read the **canonical**
  run-set CSV. When two datasets exist (e.g. a small-L pilot vs the production
  sweep), confirm which is authoritative before plotting; the buildable plot shape
  is constrained by how the sweep was sampled.
- **Run parameters live on the slide, not the canvas.** Strip in-plot run-info
  boxes / threshold annotations; state them in the slide/notes (black text).
- **Tables in slides:** native pptx table, header-row fill + bold only, black text,
  no row banding (`tbl.horz_banding=False`).

## Pre-ship checklist (every figure)

- [ ] Built for the right target (presentation → title on; report → title cropped).
- [ ] Presentation figure saved with `style.save_presentation` (no clipped title/colorbar).
- [ ] No analytical text / arrowheads on canvas (arrows OK on workflow diagrams).
- [ ] Legend names series identity only.
- [ ] System-design = real density + dashed slab/CAP extents, no param legend.
- [ ] Colorbar outside, same height as panel; scientific-offset ticks; shared clim.
- [ ] All on-plot numbers as `×10ⁿ` (no `1e-03`); no exotic Unicode glyphs.
- [ ] Functional shape verified vs the engine source; signal plotted vs baseline.
- [ ] Reads the canonical run-set CSV (provenance checked).
- [ ] Tables: header row coloured only.
- [ ] Workflows: Graphviz `.dot` → PNG, role-coloured nodes, arrows kept.
- [ ] No overlap — split if crowded.
- [ ] "linear response" (not "Lindhard") on the curve.
- [ ] Caption written to the slide spec (presentation) / `\caption` (report).
- [ ] venv python, canonical theme, `.png` 600 DPI, individual files.
