# Handover: Report 1 draft5 figure remake

## Current status
All panel plots are **produced** and **wired into `draft5.tex`**, which compiles
clean (pdflatex exit 0, no undefined refs, no overfull boxes, ~4.7 MB PDF).
The authoritative per-panel record is
`docs/reports/report1/drafts/draft5/panels_plan.md` (see its
"LaTeX wiring status — draft5.tex (2026-05-29)" block).

## What changed (this session)
- Finished **A10** 1D momentum: legend → upper-left, time labels in a.u.
- Made **B2** (`make_fig_jellium_gs.py` → two square `fig_jellium_gs_n162/n138.png`).
- Made **B3** (`make_fig_leed_backscatter_ccbond.py` → `fig_leed_backscatter_ccbond.png`,
  exact clone of main-body backscatter style, `run_cc_bond` screen 14 step 330).
- Rewired `draft5.tex`: A4 (1×2 setup), A5 (2×2, 3D placeholder), A6 (reordered 2×2),
  A9 (split → `fig:loss-function` main + `fig:plasmon-fft` appendix), B2 (1×2).
- Copied legacy `fig12_pseudopotential.png` into `drafts/draft5/figures/` (A3 kept).
- Commented out (reversible) REMOVE panels: A13 gantt, B1 scenario-ab,
  B4 leed-validation, B5 momentum-2d-full, B6 plasmon-realspace; repaired all
  dangling `\ref`s.

## Files touched
- `docs/reports/report1/drafts/draft5/scripts/make_fig_jellium_gs.py` (new)
- `docs/reports/report1/drafts/draft5/scripts/make_fig_leed_backscatter_ccbond.py` (new)
- `docs/reports/report1/drafts/draft5/scripts/make_fig_momentum_1d.py` (edited earlier)
- `docs/reports/report1/drafts/draft5/draft5.tex` (figure envs A4/A5/A6/A9/B2 + removals)
- `docs/reports/report1/drafts/draft5/figures/remake/*.png` (B2/B3/A10 outputs)
- `docs/reports/report1/drafts/draft5/figures/fig12_pseudopotential.png` (copied)
- `docs/reports/report1/drafts/draft5/panels_plan.md` (PRODUCED markers + wiring block)

## Commands run
- `venv/bin/python3 .../make_fig_momentum_1d.py` → `fig_momentum_1d.png` (2100×1650)
- `venv/bin/python3 .../make_fig_leed_backscatter_ccbond.py` → 2100×2100
- `venv/bin/python3 .../make_fig_jellium_gs.py` → n162/n138, 2100×2100 each
- `pdflatex -interaction=nonstopmode -halt-on-error draft5.tex` ×2 → exit 0

## Tests and validation
- All active `\includegraphics` paths resolve; no dangling `\ref`.
- Dimension checks: B2/B3 = 2100×2100 (3.5×3.5 square); A10 1D/2D = 2100×1650.
- Physics sanity: N=162 GS modulation 0.51% (flat closed shell) vs N=138 2.51%
  (broken-symmetry partial shell) — matches expectation.

## Trusted sources used
- Existing project run data only; styling from `inqview/report1/_shared_style.py`.

## Attribution notes
- B3 styling cloned from `make_fig_leed_backscatter_centre.py` (main body).
- B2 data-loading adapted from `inqview/report1/fig_jellium_gs.py`.

## Known issues / blockers
- **A5 panel (a) 3D render — DONE.** `render_setup3d.py` works via the full
  ParaView install (`Misc/ParaView-6.1.0.../bin/pvbatch`; system pvbatch lacks
  `paraview.simple`). Fixed a PV-6.1 API break in `inqview/paraview.py`
  (`ResetCameraClippingRange` → try/except `ResetCamera`). Frame 4 chosen →
  `fig_free_3d.png`, wired into A5(a). Candidates remain in
  `figures/remake/_setup3d_frames/` if a different timestep is preferred.
  Now rendered at 2100×2100 (3.5×3.5 in) with the **cell box** drawn (reusable
  `VolumeRenderSpec.draw_outline` option in `inqview/paraview.py`).
- **Tsubonoya distortion — FIXED.** `make_fig_tsubonoya.py` now fits-and-pads the
  397×319 extract to a 2100×2100 white square (no stretch). Matches the other A6
  cells at 3.5×3.5 in.
- All A5 and A6 panel cells are now uniform 2100×2100 (3.5×3.5 in) squares.
- **Raster "render-properly" fix (2026-05-29):** bare ParaView/PIL rasters are
  full-bleed and looked larger than inset matplotlib neighbours. Both
  `fig_free_3d` and `fig_tsubonoya` are now **embedded in 3.5×3.5 matplotlib
  figures** with the neighbour's data-box margins (3D → A5 line-plot `SQ`;
  Tsubonoya → A6 heatmap layout + blanked colourbar slot + `aspect='equal'`).
  Verified in-PDF pages 20/22. Both now carry 600 dpi. See
  [[reference-fixed-dimension-plot-pitfalls]] items 3–4.
- **A8 energy-decomposition (2026-05-29):** time axis was in fs → corrected to
  **atomic units** in both `make_fig_energy_decomp_{system,wp}.py`; sizing
  confirmed canonical one-column 3.5×3.0 (2100×1800). Verified in-PDF page 26.

## Assumptions still in play
- REMOVE panels (A13/B1/B4–B6) stay removed; commented (not deleted) for provenance.
- A3 pseudopotential kept as-is (not remade this round; figure flagged "too verbose").

## Exact next steps
1. User inserts the 3D render PNG and uncomments the A5 (a) `\includegraphics`.
2. Optional: clean Tsubonoya re-extraction to avoid the stretch distortion.
3. Later phases (per the original 8-step plan): caption-writing skill, then a
   per-panel caption documentation markdown.
