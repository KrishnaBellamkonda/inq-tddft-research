---
name: tufte-viz
description: |
  Ideate and critique data visualizations using Edward Tufte's principles from "The Visual Display of Quantitative Information." Use this skill when:
  (1) Designing new data visualizations or charts
  (2) Critiquing or improving existing visualizations
  (3) Reviewing dashboards or reports for graphical integrity
  (4) Deciding between visualization approaches
  (5) Reducing chartjunk or improving data-ink ratio
  (6) Planning small multiples or high-density displays
  Applies principles: data-ink ratio, chartjunk elimination, graphical integrity, lie factor, small multiples, and data density.
---

# Tufte Visualization Ideation

Apply Edward Tufte's principles to design clear, honest, high-density data visualizations.

## Workflow

### For new visualizations:

1. **Clarify the data story**
   - What comparisons matter?
   - What's the key insight to communicate?
   - Who's the audience?

2. **Select approach** using Tufte principles:
   - High comparison need → Small multiples
   - Dense data → Consider data tables, sparklines
   - Time-series → Line charts with minimal grid
   - Part-to-whole → Avoid pie charts; prefer bar/table

3. **Design with data-ink in mind**
   - Start minimal, add only what's necessary
   - Every element must earn its ink
   - Default to grayscale; use color purposefully

4. **Apply the Tufte test** (see references/tufte-principles.md)

### For critiquing visualizations:

1. **Check graphical integrity**
   - Calculate lie factor if proportions seem off
   - Verify baselines and scales
   - Look for 3D distortion

2. **Identify chartjunk**
   - Decorative elements
   - Heavy grids
   - Unnecessary 3D effects
   - Moiré patterns

3. **Evaluate data-ink ratio**
   - What can be erased?
   - What's redundant?

4. **Suggest improvements** with specific before/after recommendations

## Key Principles Reference

- `references/tufte-principles.md` — core principles from *Visual Display of Quantitative Information*: lie factor, data-ink, chartjunk, small multiples, integrity.
- `references/analytical-design.md` — extensions from *Envisioning Information*, *Visual Explanations*, and *Beautiful Evidence*: the 6 principles of analytical design, sparklines, layering & separation, micro/macro, range-frames, causality, confections. Load when designing dashboards, dense displays, sparklines, or explanatory graphics.

## Layout change (2026-05-25)

The report template was switched from **two-column** (`elsarticle 3p,twocolumn`)
to **one-column** (`elsarticle 3p`) on 2026-05-25. Figure widths changed:

| Width class | Old (two-col) | New (one-col) |
|---|---|---|
| `single` | 3.50 in (88 mm) | 5.00 in (127 mm) |
| `1.5col` | 4.69 in (119 mm) | 5.50 in (140 mm) |
| `full`   | 7.09 in (180 mm) | 6.30 in (160 mm) |

All figures must be re-rendered after any width change. The executable
source of truth is `_shared_style.py:column_widths_in`.

## Project annotation rules (codified 2026-05-25)

These rules apply to every figure produced for this project. They override defaults when in conflict.

1. **No analytical or interpretive text on figures.** Explanatory paragraphs, physical arguments, and derivation results (e.g. "For electron projectiles: m_e/M_nucleus ≪ 1 ⇒ nuclear stopping negligible") belong in the LaTeX caption, never on the figure canvas.
2. **Short quantitative annotations are permitted** (e.g. E_cross ≈ 5 keV/u, ω_p = 3.53 eV) but must be placed in whitespace that does not obstruct any data ink. If no such whitespace exists, the annotation goes in the caption instead.
3. **No long leader lines.** Arrows or lines pointing from an annotation to a data feature are forbidden. If the annotation cannot sit adjacent to its feature, move it to the caption.
4. **Curve labels must contrast their background.** A label placed atop a filled region or a dark curve must be offset into a clear area, or given a minimal white bbox with no visible border.
5. **Standard units on all axes.** Stopping power: eV/Bohr. Energy: eV. Length: Bohr. Time: **fs** (femtoseconds, not a.u.). Momentum: Bohr⁻¹. No mixed-unit axes. Unit labels must appear on every axis.

**Quick checklist:**
- [ ] Lie Factor ≈ 1.0 (no visual distortion)
- [ ] Maximum data-ink ratio
- [ ] Zero chartjunk
- [ ] Clear labeling
- [ ] Answers "compared to what?"
- [ ] Shows causality or mechanism where relevant
- [ ] Multivariate (not over-reduced)
- [ ] Words, numbers, images integrated — not segregated
- [ ] Reveals multiple levels of detail (micro + macro)
- [ ] Layering: primary data dominates, secondary recedes
- [ ] Appropriate data density
- [ ] No analytical/interpretive text on the figure (rule 1)
- [ ] Quantitative annotations in whitespace only (rule 2)
- [ ] No long leader lines (rule 3)
- [ ] Curve labels contrast their background (rule 4)
- [ ] Standard units on all axes (rule 5)
