---
name: report-figures
description: |
  Owns the PROJECT-WIDE figure standard (ADR 0004 / Cluster R) AND the
  interactive report-figure production workflow. Use (a) whenever producing ANY
  figure in this project — the standard below (canonical theme + annotation
  rules) applies to analysis plots and report panels alike; and (b) for the full
  spec → individual plots → LaTeX panels → compiled PDF workflow when the user
  has a plot specification. Covers:
  (1) Grilling the user to resolve all ambiguities in the spec
  (2) Making plots one at a time with iterative feedback
  (3) Composing panels in LaTeX with minipage layouts
  (4) Integrating into the report .tex file
  Invoke tufte-viz for general critique; uses inqview.visualisation.style.
---

# Report Figures — the global figure standard + production workflow

> **See also `scientific-figures` — the hard figure RULE SET** (density-based
> system-design plots, minimal legends, colorbar-outside-same-height, titles
> present-for-presentation/cropped-for-report, captions, `.drawio`+matplotlib
> workflows, header-only table colour). This skill owns the canonical THEME +
> production workflow; `scientific-figures` owns the rules. Read both.

## The project figure standard (applies to EVERY figure)

Every figure in the project — auto-generated `analyse.py` plots, comparison
`scripts/`, and interactive report panels — uses the canonical theme
`inqview.visualisation.style` (ADR 0004). No ad-hoc `plt.rcParams` /
`plt.style.use` / hardcoded `figsize` anywhere; go through the theme:

```python
from inqview.visualisation import style
style.apply_theme()                       # canonical rcParams
fig, ax = style.figure_one_col()          # fixed-dimension factory
cmap = style.cmap_for("diverging")        # semantic role, never a literal
ax.set_xlabel(style.axis_label("time"))   # canonical units ("time (fs)")
```

**Production rules (apply to every project figure; moved here from tufte):**
1. **No analytical/interpretive text on the figure canvas.** Explanatory
   paragraphs, physical arguments, derivation results belong in the LaTeX
   caption, never on the canvas.
2. **Short quantitative annotations** (e.g. `ω_p = 3.53 eV`) are permitted only
   in whitespace that obstructs no data ink; else they go in the caption.
3. **No long leader lines** from an annotation to a data feature — if it can't
   sit adjacent, move it to the caption.
4. **Curve labels must contrast their background** (offset into clear space or a
   minimal borderless white bbox).
5. **Standard units on every axis** — use `style.axis_label(...)`: energy eV,
   length Bohr, time **fs** (not a.u.), momentum Bohr⁻¹, stopping power eV/Bohr.
   No mixed-unit axes.
6. **Prefer the concrete value over an abstract label.** When a parameter has a
   known numeric value, put the value on the figure title/axis/caption and in
   the surrounding prose — `10 eV`, `L=20 Bohr`, `η=−0.5` — never a stand-in
   like "the anchor energy", "the reference width", or "η*". The reader should
   not have to look elsewhere to resolve what the label means; spell it out so
   no cognitive load is spent retrieving the number. (A variable like `ANCHOR_E`
   may still exist *in code*; what's forbidden is the abstract word leaking into
   the rendered title/text.) Reflectivity-style curves: also show a **linear
   companion panel (0–1)** beside any log-scale panel and mark threshold levels
   (e.g. a dashed 1 % line) explicitly.
7. **Shared, fixed colour scale across compared figures AND animation frames.**
   Any set of plots/GIFs a reader compares — the same quantity across time
   windows (full vs zoom/transit), runs (baseline vs perturbed), or conditions —
   must use **one identical colour scale**: same colormap, same `vmin/vmax`, same
   norm. Compute the clim **once** (from the union / the full window), persist it
   (e.g. a small `*_clim.json`), and pass it to every variant — never let each
   figure or each frame pick its own. For **animations specifically**: set
   `vmin/vmax` at the first `imshow` and only `set_data` per frame (never
   recompute the scale per frame, or the colour of a fixed density value drifts
   frame-to-frame and the GIF lies). Difference/diverging panels own a symmetric
   scale (`±max`) but it is still fixed across all frames and all compared
   variants. Sequential quantities with a huge localized spike (e.g. a wavepacket
   over a flat bath) are clipped at a robust percentile (e.g. 99.5th) so the
   background dynamics stay visible — the same clipped clim shared by every
   related GIF. (Codifies the long-standing shared-colorbar rule and extends it
   to time-series animations.)
8. **Readable colorbar ticks — no edge-clipping.** Small-magnitude fields
   (densities ~10⁻³–10⁻⁴, Δn ~10⁻⁴) produce long decimal tick labels that are
   clipped at the figure edge and unreadable. Force a **scientific offset**: one
   shared `×10ⁿ` multiplier shown once at the top of the colorbar, with **≤2
   significant-figure** tick values and ≤5 ticks. matplotlib:
   `fmt = ScalarFormatter(useMathText=True); fmt.set_powerlimits((0,0));
   cb.ax.yaxis.set_major_formatter(fmt)` (+ `MaxNLocator(5)`). Log/symlog axes keep
   their native `10ⁿ` decade ticks. Never leave raw long decimals on a colorbar.
9. **Every density map ships a LINEAR and a LOG panel, side by side.** For any
   density/field map (`n`, `Δn`, wake), render the linear and log versions in one
   figure (or one GIF) next to each other — linear shows magnitude, log reveals the
   low-amplitude structure (Friedel tails, wake fronts, the absorbed-WP remnant).
   Positive fields → `LogNorm` (floor at `vmax/10³`); signed/diverging fields (Δn) →
   `SymLogNorm` (`linthresh≈max/100`), keeping the diverging cmap. The fixed-clim
   rule (7) applies to BOTH panels. (Reference implementation:
   `ResearchProject/systems/localised_jellium/hypotheses/_density_views.py::_gif`.)

Auto-generated analysis plots must obey rules 1–5 + the theme, but do NOT need
the interactive Phase 0–5 workflow below (that is for report panels).

# Report Figures — Production Workflow

## Overview

This skill covers the full pipeline from a plot specification document to
publication-ready figures integrated into a LaTeX report. The workflow is
collaborative and iterative — the user drives visual decisions, the assistant
drives data processing and LaTeX composition.

## Phase 0: Grill the Spec

Before making any plots, resolve every ambiguity in the specification.
Use the `grill-with-docs` skill to interview the user about:

1. **Data sources**: Which run directory? Which CSV/VTI files? Confirm paths exist.
2. **Canonical run**: Which run is the "base" for all analysis? Pin it early.
3. **Units**: What units for each axis? (fs not a.u. for time, eV for energy, Bohr for length)
4. **Colour scales**: Linear, log, or symlog? Floor values? Symmetric or asymmetric norm?
5. **Panel composition**: Which plots go together in a single figure? What layout (2×2, 1×2, 2×1, etc.)?
6. **Appendix vs main text**: Which figures are primary, which are supplementary?
7. **Template constraints**: Column width, font package, existing figure macros.
8. **Free propagation / reference data**: Confirm availability of control runs.

Write the resolved spec to a `*_resolved.md` file alongside the original.
Never modify the user's original spec file.

## Phase 1: Style Setup

Before any plots, apply the canonical theme `inqview.visualisation.style`
(`style.apply_theme()`, `STYLE_CONFIG` lives there) — **not** the retired
report1-local `_shared_style.py`. For report panels, scale the theme's font/line
sizes for downscaling on top of the canonical base:

- **Font sizes must account for panel downscaling.** If plots will be placed
  2-per-row in LaTeX, they render at ~50% of their generated width. Font sizes
  must be scaled up accordingly:
  - Standalone plots: 10pt labels, 9pt ticks
  - 2-per-row panels: 14pt labels, 13pt ticks
  - 4-per-row panels: 18pt labels, 16pt ticks
- **Line widths scale similarly**: 1.2 → 1.8 for 2-per-row.
- **Tick geometry scales**: major size 3.5 → 5.0, width 0.6 → 0.9.
- Make all style values tuneable via a `STYLE_CONFIG` dict.
- Set `text.usetex: True` with preamble matching the report template.
- **Always produce individual plots, never pre-composed panels.**
  Panel composition is done exclusively in LaTeX via `minipage` layouts.
  Each script outputs one PNG per subplot. This allows the user to
  rearrange, resize, or drop individual subplots without regenerating.

## Phase 2: One Plot at a Time

Make plots sequentially, showing each to the user before proceeding.

### Per-plot checklist:

1. **Identify data source**: read the raw data, confirm format and range.
2. **Check existing scripts**: search for prior versions (`grep`, `find`).
   Adapt the data-loading logic but apply the new style.
3. **Generate at 600 DPI**: save to the designated output directory.
4. **Preview**: resize to ≤1200px wide for API display. If too large,
   use PIL to create a preview copy.
5. **Show the user**: display the preview and wait for feedback.
6. **Iterate**: apply feedback (colour scale, floor, norm, layout).
   Common iteration patterns:
   - Linear → log → symlog (for data with both signs)
   - Symmetric norm → asymmetric norm (when positive features are faint)
   - Noise floor tuning (start at 5%, adjust up/down)
   - Gaussian smoothing (σ=1 pixel) to suppress single-pixel noise
7. **Record in the resolved spec**: document the final choices
   (scale, floor, norm, data source, key learnings).

### Common pitfalls learned from experience:

- **FFT peak detection vs analytic positions**: for crystalline systems,
  use analytic Bragg peak positions from the lattice constant, not FFT
  peak detection (which gives wrong positions due to grid aliasing).
- **Screen data is real-space, not k-space**: LEED screens record density
  at a plane. FFT the screen to get the diffraction pattern.
- **Spectral resolution**: short TDDFT runs (~10 a.u.) give frequency
  resolution of ~1 eV — too coarse for spectral analysis. Use dedicated
  long-propagation runs (>1000 a.u.) for loss functions and plasmon spectra.
  State this explicitly in the report.
- **Asymmetric norms**: when positive and negative features differ by 3×+,
  use `TwoSlopeNorm(vmin=-vmax_positive, vcenter=0, vmax=vmax_positive)`.
  The stronger signal saturates; the weaker signal gets full colour range.
  This is honest — the colourbar labels make it clear.
- **Data loader consistency**: use `load_leed_pattern` (or equivalent
  library functions) rather than manual `np.loadtxt` + `reshape` to avoid
  axis transposition bugs.
- **v1 vs v2 runs**: filled markers = v2 (preferred), hollow = v1 (older).
  Document this convention for the report.

## Phase 3: Panel Composition in LaTeX

After all individual plots are approved, compose them into panels.

### Layout patterns:

```latex
% 1×2 side-by-side
\begin{figure*}[!t]
\centering
\begin{minipage}[t]{0.48\textwidth}\centering
  \includegraphics[width=\linewidth]{...}
\end{minipage}\hfill
\begin{minipage}[t]{0.48\textwidth}\centering
  \includegraphics[width=\linewidth]{...}
\end{minipage}
\caption{...}
\end{figure*}

% 2×2 grid
% Use 0.48\textwidth for each minipage, \vspace{4pt} between rows.

% 2-row with different widths (e.g. top: 2 plots, bottom: 1 full-width)
% Top row: two 0.48\textwidth minipages
% Bottom row: one 0.95\linewidth includegraphics

% n-column × m-row (e.g. 2×4 for density panels)
% Use 0.48\textwidth per column, \vspace{2pt} between rows.
```

### Rules:

- Use `figure*` for panels wider than a single column.
- Use `figure` with `[H]` for standalone plots.
- Use `\linewidth` fractions, not absolute widths.
- Tsubonoya or external reference images may need `width=0.75\linewidth`
  to match the aspect ratio of generated plots.
- Always fix orphaned `\ref{}` when removing or merging figures.
- Grep for all references to removed labels before deleting.

## Phase 4: Integration and Compilation

1. **Duplicate the .tex file** — never modify the original. Work in a `_v2.tex` copy.
2. **Replace figures systematically** — go section by section.
3. **Fix cross-references** — grep for every `\ref{fig:...}` that pointed
   to a removed or merged figure. Update to the new label or panel letter.
4. **Add appendix figures** — supplementary plots go in a new
   `\section{Supplementary figures}` before the bibliography.
5. **Compile twice** — first pass for layout, second for cross-references.
6. **Check for errors** — grep the `.log` for `Error` lines specific to
   remake figures (ignore pre-existing missing files).
7. **Verify all remake files exist** — grep `includegraphics.*remake/`
   and check each path.

## Phase 5: Feedback and Iteration

After compilation, the user reviews the PDF. Common feedback:

- **Text too small** → increase font sizes in `STYLE_CONFIG`, regenerate all plots.
- **Panel too cramped** → change layout (4×2 → 2×4, split into separate figures).
- **Features not visible** → adjust colour scale, norm, or floor.
- **Move to appendix** → cut from main text, paste into appendix section,
  fix all `\ref{}` cross-references.

Each round: apply changes → regenerate affected plots → recompile → show user.

## Production Log

Every session must produce a `plots_*_log.md` documenting:
- Each plot: script name, data source, key parameters, learnings.
- Style decisions: font sizes, colour scales, norms.
- What was tried and rejected (e.g. "FFT peak detection gave wrong positions").

## Integration with Other Skills

- **tufte-viz**: invoke for critique at each plot (TufteCritic class).
- **grill-with-docs**: invoke at Phase 0 for spec resolution.
- **journal-entries**: if plots relate to a specific run, cross-reference
  the journal entry for that run.
