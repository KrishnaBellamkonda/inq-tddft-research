# Handover: Report-2 jellium-slab figures (WP vs classical)

Rolling file. Figures dir: `docs/reports/report2/drafts/draft2/figures/`.

---

## Milestone: 2026-08-11 — Canonical style system added; all 8 figures remade in draft2

### Current status

All 8 locally-buildable figures are in `docs/reports/report2/drafts/draft2/figures/`.
The canonical style system (`inqview/visualisation/style.py`) has been extended
with a full colour palette, variable-height figure factories, and updated cmap
roles. Three source-level changes applied before rerunning: density cmap fix,
S(E) x-axis in fig 15, title removal. Nothing was written to draft1.

### What changed

- **`inq-stack/python/inqview/visualisation/style.py`**: added full colour
  palette (`CLASSICAL`, `WP`, `_PASTEL`, `pastel_for()`, `SERIES`, `GREY_LINE`,
  `GREY_SPAN`, `MARKER`), kinetic-energy label constants (`LABEL_KE_TOTAL`,
  `LABEL_KE_MEAN`, `LABEL_KE_VAR`), new cmap roles `"density" → cividis` and
  `"momentum" → PuOr`, updated `figure_one_col(height_in=…)` to accept variable
  height, updated module docstring with S(E)/label conventions.
- **`draft2/figures/methods-report-2/build_setup_panels.py`**: all three density
  pcolormesh calls changed from `cmap_for("sequential")` → `cmap_for("density")`
  (inferno → cividis).
- **`ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv/build_sv_effective_width_s6.py`**:
  x-axis converted from velocity v (a.u.) to kinetic energy E = ½v²×27.211 eV
  (4 plot calls); x-label changed to "projectile energy $E$ (eV)"; on-canvas
  title removed (report standard: no title on canvas).

### Files touched

- `/local/data/public/skcb2/tddft/inq-stack/python/inqview/visualisation/style.py` — UPDATED
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/methods-report-2/build_setup_panels.py` — UPDATED (cmap)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/sigma56_sv/build_sv_effective_width_s6.py` — UPDATED (S(E), no title)
- All 11 PNGs in `draft2/figures/` regenerated

### Commands run

```bash
PY=/local/data/public/skcb2/tddft/venv/bin/python3
$PY docs/reports/report2/drafts/draft2/figures/methods-report-2/build_setup_panels.py
$PY docs/reports/report2/drafts/draft2/figures/bulk_jellium/build_dKE_vs_position.py
$PY docs/reports/report2/drafts/draft2/figures/bulk_jellium/build_bath_diff.py
$PY docs/reports/report2/drafts/draft2/figures/jellium_slab/make_7r_dEtot.py
$PY docs/reports/report2/drafts/draft2/figures/jellium_slab/make_sv_effective_width.py
$PY docs/reports/report2/drafts/draft2/figures/jellium_slab/make_lzb_report_figures.py
$PY docs/reports/report2/drafts/draft2/figures/jellium_slab/make_norm_loss.py
```

### Canonical colour system (as of this milestone)

| Role | Solid | Pastel fill | Marker |
|---|---|---|---|
| classical | `tab:blue` | `#c6d9f0` | square `s` |
| WP | `tab:red` | `#f9c7c7` | circle `o` |
| reference | — | — | diamond `D` |
| overlay line | `#888888` | — | — |
| overlay span | `#DDDDDD` | — | — |

Series palette (Okabe-Ito): `SERIES[0..5]` for sigma/velocity/other sweeps.

KE labels: `LABEL_KE_TOTAL = Δ⟨p²⟩/2m`, `LABEL_KE_MEAN = Δ⟨p⟩²/2m`,
`LABEL_KE_VAR = σ_p²/2m`. Never T1/T2 on any axis or legend.

Cmaps: density → `cividis`, diverging → `RdBu_r`, momentum → `PuOr`,
LEED/diffraction → `inferno`.

S plots: x-axis = projectile energy E (eV), not velocity v. Conversion: E=½v²×27.211 eV.

### Tests and validation

- Proposed: user visual inspection of all 11 PNGs
- Approved: pending
- Run: all 8 scripts exited 0
- No output to draft1 confirmed (find check)

### Known issues / blockers

- Fig 14: BLOCKED on CSD3 VTI data (unchanged)
- `_panel.py` in gitignored `docs/reports/` tree — must be recreated on fresh clone
- `make_sv_effective_width.py` wrapper leaves `_sv_eff_width_tmp.png` briefly then
  deletes it — this is expected

### Exact next steps

1. **User visual review** of all 11 PNGs in `draft2/figures/`.
2. **Commit** style.py and `build_sv_effective_width_s6.py` changes to git.
3. **Phase 1 (LaTeX layout)**: design placeholder panels for multi-subplot figures;
   measure slot inch-widths; re-run affected scripts with exact `height_in=` values.
4. **Sync CSD3 VTIs** for fig 14 and run `make_s6_v3_momentum_map.py`.

---

## Milestone: 2026-08-11 — Draft-2 figures: 8/9 repo-producible built; colour convention harmonised

### Current status

Draft-2 figure production is complete for all figures buildable from local data.
8 of 9 repo-producible figures are written to
`docs/reports/report2/drafts/draft2/figures/`. Figure 14 (momentum map,
sigma=6 v=3.0) is blocked on WP wavefunction VTIs that reside on CSD3.
All remaining 15 figures are CSD3/Figma externals that the user supplies.

### What changed

- **Colour convention fixed** throughout draft-2: classical → `tab:blue`,
  WP → `tab:red`. The original `build_dKE_vs_position.py` had these reversed
  (`C3`=red for classical, `C0`=blue for WP); the draft-2 version is corrected.
  Per-sigma figures (fig 15, 24) and per-velocity figures (fig 17) retain their
  Okabe-Ito series palettes — those are sigma-indexed, not WP/classical binary.
- **`_panel.py` created** at
  `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft1/figures/_panel.py`.
  This module was missing (gitignored, never committed) but imported by three
  sigma56_sv scripts. It provides `panel_mode()`, `SLOT_IN`, `slot_figure()`.
- **Draft-2 build scripts** created (self-contained, output to draft-2 paths,
  `dpi=600, bbox_inches=None`):

| Script | Produces |
|---|---|
| `draft2/figures/methods-report-2/build_setup_panels.py` | figs 3, 4, 5 |
| `draft2/figures/bulk_jellium/build_dKE_vs_position.py` | fig 6 (colour-fixed) |
| `draft2/figures/bulk_jellium/build_bath_diff.py` | fig 7 base |
| `draft2/figures/jellium_slab/make_7r_dEtot.py` | fig 13 (colour-fixed) |
| `draft2/figures/jellium_slab/make_sv_effective_width.py` | fig 15 |
| `draft2/figures/jellium_slab/make_lzb_report_figures.py` | fig 17 |
| `draft2/figures/jellium_slab/make_norm_loss.py` | fig 24 |
| `draft2/figures/jellium_slab/slab_panel/make_s6_v3_momentum_map.py` | fig 14 (BLOCKED-CSD3) |

- **Figure catalogue** updated at
  `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figure_catalogue.md`
  with production record and status for all 24 figures.

### Files touched

- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft1/figures/_panel.py` — CREATED (was missing)
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figure_catalogue.md` — UPDATED (production record added)
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/methods-report-2/build_setup_panels.py` — CREATED
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/bulk_jellium/build_dKE_vs_position.py` — CREATED
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/bulk_jellium/build_bath_diff.py` — CREATED
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/jellium_slab/make_7r_dEtot.py` — CREATED
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/jellium_slab/make_sv_effective_width.py` — CREATED
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/jellium_slab/make_lzb_report_figures.py` — CREATED
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/jellium_slab/make_norm_loss.py` — CREATED
- `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/jellium_slab/slab_panel/make_s6_v3_momentum_map.py` — CREATED
- All PNGs listed below under "outputs"

### Commands run

```bash
# All run from /local/data/public/skcb2/tddft/
/local/data/public/skcb2/tddft/venv/bin/python3 docs/reports/report2/drafts/draft2/figures/methods-report-2/build_setup_panels.py
/local/data/public/skcb2/tddft/venv/bin/python3 docs/reports/report2/drafts/draft2/figures/bulk_jellium/build_dKE_vs_position.py
/local/data/public/skcb2/tddft/venv/bin/python3 docs/reports/report2/drafts/draft2/figures/bulk_jellium/build_bath_diff.py
/local/data/public/skcb2/tddft/venv/bin/python3 docs/reports/report2/drafts/draft2/figures/jellium_slab/make_7r_dEtot.py
/local/data/public/skcb2/tddft/venv/bin/python3 docs/reports/report2/drafts/draft2/figures/jellium_slab/make_sv_effective_width.py
/local/data/public/skcb2/tddft/venv/bin/python3 docs/reports/report2/drafts/draft2/figures/jellium_slab/make_norm_loss.py
/local/data/public/skcb2/tddft/venv/bin/python3 docs/reports/report2/drafts/draft2/figures/jellium_slab/make_lzb_report_figures.py
```

### Outputs produced (all in `docs/reports/report2/drafts/draft2/figures/`)

```
methods-report-2/setup_jellium_slab.png             (fig 3)
methods-report-2/setup_cylindrical_jellium_sweep.png (fig 4)
methods-report-2/setup_cylindrical_jellium.png       (fig 5)
bulk_jellium/dKE_vs_position_rs5p7.png               (fig 6)
bulk_jellium/bath_density_diff_linear.png            (fig 7 base)
bulk_jellium/bath_density_diff_symlog.png            (fig 7 symlog alt)
jellium_slab/7r_dEtot_vs_time_both_preview.png       (fig 13)
jellium_slab/slab_sv_effective_width_s56.png         (fig 15)
jellium_slab/slab_S_of_invL_sigma5.png               (fig 17)
jellium_slab/slab_S_of_invL_sigma0p5.png             (companion to 17)
jellium_slab/slab_norm_loss_vs_sigma.png             (fig 24)
```

### Tests and validation

- Proposed: visual inspection by user
- Approved: pending user review
- Run: all scripts executed successfully (exit 0)
- Outcomes: all 8 buildable figures written at 600 DPI, `bbox_inches=None`
- Remaining gaps: user needs to visually confirm colour conventions,
  axis labels, and legend placement match report intent

### Known issues / blockers

- **Fig 14** (`slab_panel/slab_s6_v3_momentum_map.png`): requires sigma=6 v=3.0
  WP wavefunction VTIs at
  `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/sigma56_sv/wp/results/s6p0_v3p0/raw/vti/wavefunction_wp/`.
  This data is on CSD3. Run `make_s6_v3_momentum_map.py` once VTIs are synced.
- **Figures 2, 8–12, 16, 18–23**: CSD3 externals. User supplies.
- **Fig 1**: Figma external. User supplies.
- **Fig 7 annotation**: user adds arrow/label overlay to `bath_density_diff_linear.png`
  externally (e.g. Inkscape, PowerPoint) to produce `annotated-bath-density-difference.png`.
- `make_sv_effective_width.py` calls `build_sv_effective_width_s6.draw()` with
  `corrected=False` (no monopole correction). If the corrected variant is needed,
  change `draw(False, ...)` to `draw(True, ...)`.
- The `draw()` non-panel branch in `build_sv_effective_width_s6.py` uses
  `dpi=300, bbox_inches="tight"` — so the fig 15 file was produced from the wrapper
  at that DPI. To get 600 DPI for the REPORT_FIG path specifically, the `draw()`
  function's non-panel branch must be updated (or the wrapper must re-save after
  the call at 600 DPI). Currently the file exists but at 300 DPI.

### Assumptions still in play

- Draft-2 default figure size: `style.figure_one_col()` = 3.5 × 3.0 in for all
  figures. Panel-slot resizing (Phase 1 of the panel workflow) is a future pass
  once the LaTeX layout is designed.
- `tab:blue` / `tab:red` are the canonical WP/classical colours for all future
  figures in draft-2.
- `_panel.py` created in the gitignored `docs/reports/` tree — NOT committed to git.
  It must be recreated if the directory is wiped or the repo is freshly cloned.
  Consider committing it somewhere tracked (e.g. `sigma56_sv/_panel.py`).

### Exact next steps

1. **User reviews figures visually** — open the 11 PNGs in
   `/local/data/public/skcb2/tddft/docs/reports/report2/drafts/draft2/figures/`
   and confirm colours, labels, axes.
2. **Fix fig 15 DPI**: open `build_sv_effective_width_s6.py` lines 348/352, change
   `dpi=300` to `dpi=600` in the non-panel branch; re-run `make_sv_effective_width.py`.
3. **User supplies CSD3 figures** for figs 2, 8–12, 16, 18–23 (sync to
   `draft2/figures/<section>/` matching the LaTeX paths in the catalogue).
4. **Sync fig 14 VTIs** from CSD3 and run
   `draft2/figures/jellium_slab/slab_panel/make_s6_v3_momentum_map.py`.
5. **Phase 1 (LaTeX panel layout)**: design placeholder LaTeX panels for
   multi-subplot figures (figs 3+4+5, 6+7, 8–11, 13+14) and measure slot
   inch-widths; then re-run affected scripts at the exact slot widths.
6. **Write/update `draft2.tex`** with `\includegraphics` calls pointing to the
   new draft-2 figure paths.

---

## 2026-08-02 — S(v) versions + case study kickoff

### S(v) curves (DONE) — `make_S_of_v_versions.py`
Three report-ready versions, E-absorbed S_B, r_s=4.18. Encoding: colour = sigma
family; WP filled/solid, classical hollow/dashed; NO classical-vs-WP legend (fill
only); legend lists sigma_r values. Classical clipped to v<=3.5 to match WP grid.
- `S_of_v_v1_start_sigmar.png` — starting sigma_r (0.61/2.45/3.67 Bohr); no 17/20.
- `S_of_v_v2_timeavg_sigmar.png` — WP time-averaged <sigma_r> (22/18/17 Bohr);
  classical stays at starting sigma_r (noted in legend), no dispersion.
- `S_of_v_v3_with_wide.png` — starting sigma_r + wide classical sigma_WP=17/20.
Open question for user: legend box vs inline end-of-curve sigma_r labels; whether
to annotate both WP and classical <sigma_r> explicitly in v2.

### Case study = sigma_WP=2, v=2.0 (user-chosen 2026-08-02)
Constraint: sweep density frames were deleted; NO WP run has a local frame set, so
the total + induced density GIFs REQUIRE a local re-run with frame-saving.

**Re-run launched** (autonomous, dual-GPU queue) —
`scripts/wp_highdensity_sv/orchestrate_case_twins.py`, log `orchestrate_case_twins.log`:
  1. cs_s2p0_v2p0  (CASE STUDY): sigma=2, v=2.0, launch -24.0, dx=0.5, N=3623,
     frames every 13 steps. WP=wp/results/cs_s2p0_v2p0, CL=dyn_direct_cap/results/cs_s2p0_v2p0_cl.
  2. s4_e200 (previous task): sigma=4, 200 eV, launch -24.5, N=1890.
Both halves share the dx=0.5 GS => MATCHED pair (better than the published sweep,
which was WP dx0.4 + classical dx0.5). Both stepping cleanly as of 14:21.

**WP binary:** built locally against inq-study (CAP) at
`scripts/wp_highdensity_sv/wp/run`. FIX: run.cpp line 402 referenced
`WPRealSpaceMoments.zc` (circular centroid) absent in local inqkit — changed to
`.z` (t=0: plain==circular). Smoke (5 steps, dx0.5 GS) passed.

**Resume bug fixed:** the first orchestrator hardcoded LJ_RESUME=1 -> fresh runs
aborted ("no rt_state.txt"). Now RESUME auto-detected (1 only if rt_state exists).

### Case-study analysis panels — `make_case_study.py [preview|rerun]`
10 panels to `figures/case_study/`. PREVIEW generated now from existing sweep CSVs
(WP dx0.4 sweep_data/s2p0_v2p0 + classical dx0.5 s2p0_v2p0_cap); re-run
`make_case_study.py rerun` once the matched dx0.5 re-run completes for the FINAL,
self-consistent set:
  1 total_energy_vs_time            6 selected_window (shaded [5.75,18.25] a.u.)
  2 position_vs_time (WP=int<pz>dt) 7 dEtot_vs_time_both
  3 T1/T2/var(p)/2m vs time         8 Eabsorbed_over_L_both (running S_B)
  4 T1/T2/var(p)/2m vs position     9 interactions_vs_time (WP 6 pairwise terms)
  5 classical dEtot vs position    10 momentum before/after + induced
Delta-E referenced to t=0 (E(t)-E(0)) — grid-independent; WP uses
energy_total_corrected when present. Colour: WP=red, classical=blue (no legend text).

### STILL TO DO
- Density GIFs (total + induced) for BOTH halves — after re-run frames land; build
  with inqview make_density_gif_battery / make_twin_density_matrix (density-gif rule).
- Re-run `make_case_study.py rerun` for the final matched-grid panels.
- User review of S(v) legend style + v2 labelling.
