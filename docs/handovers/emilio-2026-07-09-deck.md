# Handover — 09 July 2026 Emilio meeting deck (first draft)

Rolling handover for the assertion-evidence deck built from
`docs/reports/09-07-2026-meetng-emilio/presentation-plan.md`.

## 2026-07-09 — first draft assembled

### Done
- **Skill edits** (`docs/presentations/method/assertive-evidence-presentation.md`):
  1. New **§3b Section-divider slide** type — a single centred title, nothing else.
  2. New **§4b House body-slide layout** — topic title + captioned plots + a
     bottom **learning-point** line (blank placeholder where the plan gave none).
     Recorded as a deliberate, user-chosen deviation from strict AE (which puts the
     message in the headline).
- **Framing decisions (locked):** F1 sections-only (no title/mapping/conclusion),
  F2 empty Section 4 divider placeholder, F3 embed GIFs (deliverable is a live PPTX).
- **Roadmap:** `docs/reports/09-07-2026-meetng-emilio/figure-roadmap.md` — slide →
  figures → source → exists/generate → caption → learning blank. Authoritative.
- **New figures generated** (venv python):
  - `assets/make_section1_assets.py` → `assets/section1/` : infinite-plate &
    stacked-plates schematics + infinite-plate/slab equation PNGs.
  - `assets/make_extra_assets.py` → `assets/section2/cutoff_radial_potential.png`
    (V(r) from the 4 UPFs) and `assets/section3/muon_workflow.png` (per-orbital
    inverse-mass fork box-flow).
- **Deck built:** `docs/reports/09-07-2026-meetng-emilio/build_emilio_deck.py` →
  **`emilio_deck_draft1.pptx`** (20 slides, 16:9). Native params table (mass bold
  + red), GIF density triptych on the wide-WP slide, learning-point blanks
  throughout, section dividers, empty Section-4 divider.

### Verified vs unverified
- Verified: all 17 referenced asset paths exist; deck writes 20 slides; XC=LDA
  from `effmass_12h/gs/results/run_summary.txt`; muon key change read from
  `inq-study/src/{systems/electrons.hpp,real_time/propagate.hpp}` +
  `inqkit/wavepacket/wavepacket.hpp`.
- NOT previewed (user owns figure/deck preview per house rule). User to open the
  PPTX and review figures.

### Known refinements for the next pass (NOT blockers)
1. **S2.1 per-component waterfall** bar chart + bookkeeping table + equations. The
   waterfall is coded in `campaign_autorun_study/build_h0_p2_interpretation.py`
   but renders inside the notebook, not a standalone PNG. Currently the slide uses
   `extend_r160/component_decomposition.png` (where-energy-lives) +
   `wp_minus_cl_gap.png`. Export the waterfall + build the Δ-bookkeeping table.
2. **S2.3 ΔE_total(r) vs cutoff** (Fig B) from the 24 completed runs in
   `scripts/campaign_autorun/runs/cutoff_test/rc{10,20,30,40}/cl_r*_p2/` — coded in
   `build_theoretical_model.py`. Slide currently shows only the radial-potential
   V(r) figure (Fig A).
3. **Learning points** are blank on most slides by design — user fills.
4. **Triptych on other run-specific slides** (S2.1, high-mass densities) — only the
   wide-WP slide has the 3-GIF triptych so far; high-mass uses `p2_xz_density.png`.
5. Consider re-rendering the `extend_r160` PNGs at presentation scale (larger fonts).

## 2026-07-09 — consolidation + slide-4 revisions

- **Single downloadable figures folder:** ALL deck assets now live flat in
  `docs/reports/09-07-2026-meetng-emilio/figures/` with slide-keyed names
  (`s1_1_*`, `s2_3_*`, …). The three generators (`assets/make_section1_assets.py`,
  `assets/make_extra_assets.py`, `assets/make_s1_3_field_potential.py`) all WRITE
  into `figures/`; `build_emilio_deck.py` READS only from there. Stale
  `assets/section*/` subfolders removed. **Convention going forward: every new /
  remade plot lands in `figures/` under its slide-keyed name.**
- **Slide 4 (S1.3) revised per user:**
  - Projectile is now a **Gaussian charge, not a point charge** — the field/potential
    is the Gaussian convolution (s = sigma_WP/√2, sigma_WP=0.5). Point-charge
    `U(r)` panel dropped.
  - New plot `s1_3_semiempirical_field_potential.png`: field E(z) and potential
    phi_ext(z) vs z from the semi-empirical net density (same validated L_z=160 p3
    source as plate_model).
  - New equation PNG `s1_3_eq_phi_ext.png` rendered and placed on the slide.
  - **Open:** sigma_WP=0.5 assumed for slide 4 — confirm if a different sigma is wanted.

## 2026-07-09 — Section 2 real figures + bookkeeping + analytical calcs

- **S2.1 decomposition:** extracted the executed **"Full energy deconstruction —
  classical ghost vs WP electron"** figure from
  `campaign_autorun_study/theoretical_slab_model.ipynb` →
  `figures/s2_1_energy_decomposition.png` (now the S2.1 plot). Candidates also
  extracted: `s2_1_exc_total_difference.png`, `s2_1_screening_wp_potential.png`.
- **Bookkeeping table (S2.1b):** NEW native-pptx slide with the real per-r
  decomposition (r, dE_WP, dE_CL, WP-CL, dKin, dXC, d(H+E)) + the identity
  equation `s2_1_eq_bookkeeping.png`. Numbers are REAL (extracted from the
  notebook), not invented.
- **Analytical calculations** (`assets/make_analytical_calcs.py` →
  `analytical-calculations.md`): (1) S2.1 zero-point `3/(4 sigma^2)=81.6 eV`
  matches measured dKin; residual after subtracting analytic zero-point + Coulomb
  = dXC = -16 eV (self-XC) — the identity closes. (2) S1.3 expected Coulomb
  repulsion <~1 eV (from plate_model VALIDATION).
- **S2.3 empirical:** extracted **"Empirical cutoff sweep — 4 projectile UPFs
  (KS runs)"** → `figures/s2_3_cutoff_empirical_sweep.png`; S2.3 now shows the UPF
  V(r) AND the empirical dE_total(r) sweep (what actually happens in the sim).
- Deck now **21 slides**.
- **Notebook-figure extraction pattern:** executed PNGs pulled from `.ipynb`
  outputs by heading match (base64 image/png), saved into `figures/` with
  slide-keyed names. Reusable for future notebook figures.

## 2026-07-09 — refreshed stale U_ext / E_hartree plots

- User flagged the S2.2/S2.5 U_ext & E_hartree plots as **stale** (they came from
  `extend_r160/` figs dated 2026-07-03). The current analysis is
  `campaign_autorun_study/theoretical_slab_model.ipynb` (rebuilt 2026-07-08).
- **Regenerated from the CURRENT h0_p2 KS runs** (`scripts/campaign_autorun/runs/
  h0_p2/{wp,cl}_r{4..40}_p2/.../observables.csv`, E_GS from gs_p2_lz120) via
  `assets/make_s2_energy_vs_r.py`:
  - `figures/s2_2_energy_vs_distance.png` — dE_total(r), WP vs classical.
  - `figures/s2_energy_components_vs_r.png` — U_ext(r) & U_H(r) referenced to r=40
    (removes the p2 G=0 shift), WP & classical.
  - `figures/s2_5_hartree_external_sum.png` — (U_H+U_ext) vs r (G=0-robust), WP,
    classical, and their difference.
- **Deleted the stale files** from `figures/`: `s2_2_excess_vs_r.png`,
  `s2_component_decomposition.png`, `s2_5_right_formula_compare.png`,
  `s2_1_wp_minus_cl_gap.png` (all 07-03).
- Deck rebuilt, still 21 slides. **Lesson: prefer the newest run/notebook figures;
  check mtimes — extend_r160 (07-03) is superseded by theoretical_slab_model (07-08).**

## 2026-07-09 — Section 2 remade with actual standalone plots (ghost→projectile)

User asks: (1) turn the `theoretical_slab_model.ipynb` section-2 plots into ACTUAL
standalone figures, relabelling every "ghost" → "classical projectile"; (2) remake
ONLY Section 2 of the deck with the latest plots; (3) add a screening-baseline t=0
comparison (WP vs classical projectile + their difference).

### Done
- **New consolidated generator** `assets/make_s2_figures.py` (venv) ports the
  notebook's Section-2 analysis into presentation-scale PNGs in `figures/`, all
  labelled **classical projectile** (never "ghost"), **no baked-in annotation text**
  (plan rule — overlays/dashed lines only, numbers moved to captions):
  - `s2_1_energy_decomposition.png` — 3-panel at r=28: classical-projectile
    component decomposition | WP component decomposition (bars sum to dashed
    E_total line) | WP−CL ledger (ΔKinetic +82 eV zero-point, ΔXC −16 eV self,
    Δ(Hartree+External) −40 eV electrostatic; net +26 eV). Sum(parts)−total = 4.6e-13 Ha.
  - `s2_3_cutoff_empirical_sweep.png` — ΔE_total(r) from the 24 KS cutoff runs
    (`runs/cutoff_test/rc{10,20,30,40}/cl_r*_p2`), dotted = each cutoff radius.
  - `s2_4_potential_comparison.png` — poisson(n_WP) overlaid on the analytic
    classical-projectile Gaussian potential erf(r/(√2 σ_ρ))/r (σ_ρ=0.354).
  - `s2_6_screening_baseline.png` — **NEW (user request):** planar-mean densities;
    slab response n(t=0)−n_GS for WP AND classical projectile (both ≡0) + their
    difference; and the genuine WP−classical electron-density difference = n_WP.
    *Classical density is NOT saved on disk; used the exact fact that the classical
    projectile (z_valence=0, pure external potential, no evolution at t=0) has total
    density bit-identical to the GS — verified max|n_slab(WP,t=0)−n_GS| = 0.0, ∫n_WP=1.0000.*
- **Ledger validated:** `make_s2_figures.py` stdout reproduces the deck bookkeeping
  table row-for-row (dKin=81.7, dXC=−16.5, d(H+E) all match) → `bookkeeping_table()`
  is confirmed against real h0_p2 data, unchanged.
- **Relabelled + regenerated:** `assets/make_s2_energy_vs_r.py` ("classical ghost"→
  "classical projectile") → s2_2 / s2_energy_components / s2_5. `make_extra_assets.py`
  cutoff radial V(r) already said "classical projectile", refreshed.
- **Notebook itself:** `campaign_autorun_study/build_theoretical_model.py` had all 17
  "ghost" → "projectile"; **rebuilt AND re-executed** (`nbconvert --execute`, exit 0)
  → `theoretical_slab_model.ipynb` now consistent.
- **Deck Section 2 rebuilt** (`build_emilio_deck.py`): captions ghost-free, S2.4 now
  uses `s2_4_potential_comparison.png`, NEW screening-baseline slide added. Deck is
  now **22 slides** (was 21). Slides 1–4 (S1), 13–22 (S3+S4) **untouched**.

### Verified vs unverified
- Verified: all 9 Section-2 figures exist and resolve; deck writes 22 slides; no
  residual "ghost" in deck/asset scripts; notebook executes to exit 0; energy ledger
  matches the bookkeeping table; screening baselines exactly 0; ∫n_WP=1.
- NOT previewed (house rule — user owns figure/deck preview). User to open the PPTX.

### Pending / deliberate omissions
- **Density triptychs on S2.1/S2.6** (plan's run-specific rule) NOT added — Section 2
  is the energetics/electrostatics story; triptychs live better on the dynamics
  (wide-WP) slide. Flag for user if they want them.
- **Learning-point lines left blank** on the new Section-2 body slides (house default;
  user owns the scientific verdict).
- `s2_4_wp_weighting.png` (old plate-model figure) now unused but left in `figures/`.
- Optional empirical confirm: a 1-step classical density-saving run to replace the
  by-construction classical baseline in `s2_6` with a loaded density — offer to user.

### Key paths
- Section-2 generator: `docs/reports/09-07-2026-meetng-emilio/assets/make_s2_figures.py`
- Plan: `docs/reports/09-07-2026-meetng-emilio/presentation-plan.md`
- Roadmap: `.../figure-roadmap.md`
- Deck builder: `.../build_emilio_deck.py` → `.../emilio_deck_draft1.pptx`
- Driving notes: `docs/notes/localised-jellium-parameter-study.md`
