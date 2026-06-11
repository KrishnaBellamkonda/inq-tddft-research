---
name: report-figures
description: |
  End-to-end workflow for producing publication-quality figures for the report.
  Use when the user has a plot specification (what figures to make) and wants to
  go from spec → individual plots → LaTeX panels → compiled PDF. Covers:
  (1) Grilling the user to resolve all ambiguities in the spec
  (2) Making plots one at a time with iterative feedback
  (3) Composing panels in LaTeX with minipage layouts
  (4) Integrating into the report .tex file
  Invoke tufte-viz skill for critique at each plot.
---

# TODO: Last I remember of these report figures usage, I remember it working very well. However, I wonder if this skill, tufte, and others are just overweight skills. Does each serve a specific and unique purpose. Or are there redundant instructions?  

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

Before any plots, configure `_shared_style.py`:

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
