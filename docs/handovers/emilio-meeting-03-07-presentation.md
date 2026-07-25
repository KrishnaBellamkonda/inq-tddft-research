# Handover — 03-07-2026 Emilio meeting deck (draft 1)

Rolling handover for the supervisor-meeting presentation built from
`docs/reports/03-07-2026-meting-emilio/presentation_plan.md`.

## 2026-07-03 — draft 1 built

### State: done / partial / not done
- **DONE — AE skill edits** (`docs/presentations/method/assertive-evidence-presentation.md`):
  (1) inline numeric `[n]` citations, grey footer removed, numbered References
  slide is sole collation; (2) all on-slide text black (grey `666666` killed;
  connector lines may stay grey); (3) new **Caption convention** (every plot: what
  is visualised + quantity/units, ≤2 lines, black); (4) relaxed to **one assertion
  per slide, one _or more_ panels** when they build one argument; (5) reconciled
  build spec to **python-pptx** (not pptxgenjs).
- **DONE — new `scientific-tables` skill** (`.claude/skills/scientific-tables/`):
  `SKILL.md` + `make_table.py`. `add_native_table` (editable PPTX, coloured header
  only, black text, no zebra, thin borders) and `table_to_png` (engine `mpl`
  default / `latex` via pdflatex+booktabs). All three paths smoke-tested OK.
- **DONE — draft 1 deck** `docs/reports/03-07-2026-meting-emilio/draft1.pptx`
  (19 slides, 32 images incl. 8 animated GIFs, 5 native tables). Build script:
  `docs/reports/03-07-2026-meting-emilio/build/build_draft1.py`. Figure extractor:
  `.../build/extract_nb_figs.py` (pulls embedded PNGs from H0 notebook →
  `assets/extracted/`, with `manifest.txt`).
- **PARTIAL — two intentional placeholders** (data genuinely absent):
  (a) Section 1 slide "analytical model" — model-vs-classical-vs-WP comparison
  figure not built (H1 edge-model figure IS shown); (b) Section 2 loss-function
  slide — `L(q,ω)` figure to render from the 23 eV run notebook (S(E) localised +
  bulk context ARE shown).
- **NOT DONE (draft-2 refinements noted on-slide):** density-convergence plot
  needs y-axis reformatted to 2 s.f. + on-plot y-label (plan); 23 eV GS ledger
  wants E_GS/E_total(0)/E_H(0) extracted from the run notebook (table currently
  has run params + S from `se_state.csv`); ML section (4) deferred by user.

### Locked design decisions (grill)
1. Tooling first, then deck. 2. Dense working-meeting layout + assertion
headlines. 3. PPTX output, **animated GIFs** for movies. 4. Sections 1–3;
placeholders for gaps; no fabrication; ML section skipped.
5. **"ps6" resolved** → it meant `p5_wp_v6` (490 eV), which **aliases** (dx=0.5,
k0=6 ≈ k_max=π/dx=6.28). Replaced by the corrected top-velocity run
`p5_wp_v5p0_h035` (340 eV re-run at dx=0.35, k_max=8.98). See se_state.csv:
coarse v5/v6 show S=9.8/18.9 eV/Bohr (artefact); trusted point is the h035 rerun.

### Verified vs unverified
- **Verified:** every asset path referenced by the build script exists on disk
  (checked); deck opens, images+GIFs+tables embed (zip media = 21 png + 8 gif);
  both table engines + native table run; run→energy map from `run_summary.txt`;
  S(E) numbers from `hypotheses/qsp_phase5/se_state.csv`.
- **Unverified (needs user eye — do NOT self-preview figures):** whether each
  auto-extracted H0 panel (energy-decomposition c009, periodicity-2 c011) is the
  intended figure — captions carry "auto-extracted — verify".

### Key asset locations
- Section 1: `localised_jellium/scripts/campaign_autorun/runs/h{0..5}/H*.png`,
  `runs/extend_r160/extend_r160_excess_vs_r.png`; extracted H0 panels in
  `.../assets/extracted/`.
- Section 2: `localised_jellium/hypotheses/qsp_phase5/` (`figs/`,
  `p5_wp_v1p3_run_notebook_figs/`, `p5_wp_v5p0_h035_run_notebook_figs/`) and
  `qsp_phase4/{figs,p4wp_run_notebook_figs}/`. Bulk S(E):
  `systems/jellium/hypotheses/stopping_power_vs_energy_all.png`.
- Section 3: `cylindrical_jellium/hypotheses/annular_sv/{gs_validation,per_run_figs}/`,
  `Sv_beta.png`, `Sv_results.csv`.

### Next actions (post draft 1)
1. Reformat density-convergence axes; extract 23 eV GS energy ledger.
2. Rebuild via `build_draft*.py` (idempotent; robust to missing files).

## 2026-07-03 (later) — draft 2 built (user feedback folded in)

`draft2.pptx` (24 slides). Builder `build/build_draft2.py`; equations
`build/render_equations.py` → `assets/equations/*.png`; title crops in
`assets/notitle/`.

### Feedback applied
- **Titles inverted (anti-AE, user request):** every slide now has a SIMPLE
  descriptive title; the AE assertion is demoted to an on-slide `Takeaway:` line.
- **Section-divider slides** added (S03, S09, S18).
- **Section 1 reworked** around the new `hypotheses/extend_r160/` study:
  - S04 injection energy vs r (`excess_vs_r.png` + `wp_minus_cl_gap.png`, r→60).
  - S05 = the reworked old slide 3 ("fair comparison"): `component_decomposition.png`
    + `right_formula_compare.png` + rendered formulae (`decomp`, `fair` from
    `assets/equations/`) + 2-line method note. **The fair-comparison formula is
    classical ΔU_ext vs WP ΔU_H** (they cross near r≈22) — from the extend_r160
    notebook, grounded on KS decomposition [Parr & Yang, ref 3].
- **Wide-WP (p0b) new slides** in Section 2: S15 case study (WP `lead_density.gif`
  + stopping) + S16 gate review (6-panel c1..c6). NOTE: the WP run only rendered
  `lead_density.gif` (no `*_total_density.gif`; that exists only for the classical
  twin `p0b_classical`).
- **On-plot titles stripped** from the flagged stopping figures via
  `strip_title` (crops top 8% → `assets/notitle/`). Applied to 23 eV, 54 eV,
  h035, wide-WP stopping plots. **Unverified: the 8% crop fraction — user must
  eyeball it didn't clip axes.**
- **"ps6"** confirmed earlier = aliased v6; the corrected-top-velocity slide (S14)
  uses `v5_h035`.

### Draft-2 state
- Real data on every slide except **two intended placeholders**: S08 (analytical
  S vs classical/WP comparison — not built) and S17 (loss function L(q,ω) — to
  render from the 23 eV run). Media: 30 png + 8 gif.
- New data that arrived 2026-07-03 and is NOT yet used: `p5_wp_v6p0` full
  post-processing (momentum/KS-energy/occupation gifs, spectra); `p0b_classical`
  full observables. Available if wanted.

### Ongoing-work slide (S24, added on request)
Final content slide "Current directions I am exploring" (after Summary, before
the References backup). Five active campaigns with **audience-friendly titles**
(campaign names deliberately hidden): no-CAP quantum-vs-classical twin → "Do a
real electron and a point charge stop the same way?"; loss-function feasibility →
"Reading the stopping power off the energy-loss spectrum"; ml-patterns →
"Teaching a machine to spot the quantum fingerprint"; QKE Li kick →
"Cross-checking our plasmon response against real lithium"; TD-HF orbital →
"Is the simulated electron a faithful quantum stand-in?". Deck now 25 slides.
User to trim/confirm which threads to keep.

### Next actions (post draft 2)
1. User previews `draft2.pptx`; check the title crops, the fair-comparison slide
   density, and confirm the takeaway wording.
2. Render loss-function L(q,ω) (S17) + analytical S comparison (S08).
3. Reformat H2 density-convergence y-axis to 2 s.f. + on-plot label.
4. Extract 23 eV GS ledger (E_GS, E_total(0), E_H(0)) for the S11 table.
