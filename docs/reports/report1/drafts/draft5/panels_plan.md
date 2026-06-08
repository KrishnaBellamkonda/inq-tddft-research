# Report 1 — Panels Plan (draft5 replot)

> **Living document.** Phase 1 output of the `replot` session (2026-05-28).
> I extract every panel as it currently exists in `draft5.tex`; **you** edit
> this file (change specs, add per-plot requests, plan new panels). I then
> **re-read it fresh and `git diff` it** to absorb your changes before Phase 2.
>
> **Authority split:** `draft5.tex` is authoritative for *which panels exist
> and their composition* (the figure environments). `structure.md` is
> authoritative for *section/subsection placement* — and it has already
> diverged from `draft5.tex` (Theory dissolved, Results subsections renamed).
> Phase 2 maps each panel onto the `structure.md` TOC, flagging any panel
> whose home section is mid-rename.
>
> **Status legend:** `REWORK` = PNG exists in `figures/remake/`, restyle to
> global rules + per-plot changes · `RELOOK` = non-remake figure, decide if
> restyle needed · `NEW` = designed/redesigned this round · `EXTERNAL` = user
> supplies, no script.

---

## 0. Decisions locked in the Phase-1 grill (2026-05-28)

1. **Canonical low-density jellium case study = `run_wp_n162_L50_E25_sigma1_v2` (25 eV).**
   The "20 eV" in `draft3/plots_remake.md` is stale and superseded.
2. **Mechanism panel (plasmon FFT + loss function) uses a *different* run** —
   the dedicated long-propagation `run_plasmon_n162_L50_E15` (E=15 eV,
   T=2000 a.u., same N=162/L=50/r_s=5.69 system). The E25 run (T≈11 a.u.,
   Δω≈1.2 eV) cannot resolve any spectral feature. **The report must state
   this explicitly** in the caption/text.
3. **Stopping power = TWO definitions**, with definition 2 decomposed into two
   component curves:
   - **Definition 1 — S₁**: total WP KS kinetic-energy change, `−ΔE_kin^WP / Δz`.
   - **Definition 2 — S₃**: total momentum-space KE change `Δ⟨p²⟩/2m`, split into
     its **components**: directed `⟨p⟩²/2m` (**S₂**) + spreading `σ_p²/2m` (**S₄**).
   - Master plot headline metric = **S₂** (directed KE). Definition-comparison
     panel shows S₂ vs S₃ vs S₄.
4. **`fig:bohr-schematic` (classical Coulomb scattering) is DROPPED.**
5. **Schematics remaining:** joint setup (`fig:leed-geometry`) + regime diagram
   (`fig:regime-diagram`). Regime diagram gets full per-plot grilling in Phase 3+.
6. **Sizing ground truth** (rule-of-thumb, tweak per plot, confirm in LaTeX):
   - One-column figure region: **3.5 in wide**; default one-column plot **3.5 w × 3.0 h** (square 3.5×3.5 for heatmaps).
   - Two-column figure region: **7 in wide**; two-column plot height ≈ **7 in**.
   - Axis labels / annotations / legends in the **8–12 pt** range on the page.
7. **Sizing workflow = save-at-final-width** (may pivot if it fights us): each
   plot is saved at the width it occupies on the page; `\includegraphics` uses
   that width with no rescaling, so on-page fonts are uniform across panels.

---

## 1. Global rules (apply to every panel unless a panel declares a deviation)

**Single source of truth = `inq-stack/python/inqview/report1/_shared_style.py`.**
`docs/reports/report1/figures/global_style.md` is its human-readable companion
and must be regenerated to match. They are currently **out of sync** and must
be reconciled (see issue below).

### 1.1 Sizing & layout
- Save each plot at its **final on-page width** (see §0.6/0.7). No LaTeX rescaling.
- Default aspect: one-column **3.5×3.0** (landscape) for line plots; **square** for 2D heatmaps/density. Two-column **~7×7**.
- Panel composition (rows/cols, minipage widths) handled in LaTeX per the draft5 layout below.

### 1.2 Fonts (TARGET: 8–12 pt on page — confirm in LaTeX in an implementation phase)
- ⚠️ **OPEN ISSUE — font drift to fix.** `global_style.md` specifies 10 pt label / 9 pt tick / 9 pt legend, but `_shared_style.py` `STYLE_CONFIG` was bumped to **14 / 13 / 13 pt** ("scaled up for panel compositions"). Combined with uniform 6.5 in saves + per-panel LaTeX scaling, this is the likely cause of inconsistent/wrong font sizes. **Decide the canonical on-page sizes (8–12 pt band) and set them once in `STYLE_CONFIG`, then regenerate `global_style.md`.**
- `text.usetex: True`, Computer Modern (`lmodern` + `mlmodern`), `amsmath`, `siunitx`.
- No figure titles (captions live in LaTeX). Panel labels `(a)`,`(b)` plain, top-left inside axes at `(0.02, 0.95)`, same size as axis label.

### 1.3 Annotations
- ⚠️ **User concern: remove unnecessary annotations.** Per-panel "Changes I want" blocks will list which in-axes annotations to cut. Default policy: keep only annotations that carry data meaning (reference lines, central-momentum markers, regime labels); drop decorative text.

### 1.4 Colour & lines
- Palettes from `global_style.md` §4: `palette_sweep5`, `palette_sweep3`, `palette_regime3`, regime tints, reference styles. Do not use matplotlib defaults. One palette per panel unless a stated semantic reason.
- Reference lines: asymptote `#808080` `--` lw 0.9; theory `#000000` `--` lw 1.0; fit `#881818` `-` lw 1.0.

### 1.5 Output
- PNG only, 600 DPI, white background, `bbox_inches="tight"`, `pad_inches=0.02`.
- Output dir: `docs/reports/report1/drafts/draft5/figures/remake/` (currently sourced from `draft3/figures/remake/` — confirm the draft5 figure dir to write to).

### 1.6 Tufte
- `TufteCritic` enabled by default; per-figure exemptions justified in the panel's "Changes I want" block.

---

## 2. Canonical runs

| Role | Run directory | Key params |
|------|--------------|------------|
| Jellium low-density case study | `run_wp_n162_L50_E25_sigma1_v2` | L=50, N=162, σ=1, E=25 eV, k₀=1.356, dt=0.01 |
| Matched free propagation | `run_free_wp_L50_E25_sigma1_v2` | same WP, empty box |
| Classical (low density) | `run_classical_n162_L50_E25` | Ehrenfest point charge |
| Mechanism / spectral (E15) | `run_plasmon_n162_L50_E15` | T=2000 a.u., 501 pts, Δω≈0.085 eV |
| High-density family | `run_wp_n162_L30_E{50,100,200,300}_highdens_sigma1_v2` | L=30, r_s=3.41 |
| Master plot aggregation | via `inqview/report1/stopping_power_data.py` | all energies/σ, v1+v2 |

**Timesteps for density panels:** t = 0.001, 0.089, 0.177, 0.266 fs (full run is interference-free). Subject to tweaking.

---

# PART A — Main-body panels

## A1. `fig:nuc-vs-elec` — Nuclear vs electronic stopping  ·  RELOOK
- **Section:** Introduction (§1, line 143). **Env:** `figure` (1-col), `\figwidthmed` (0.65 tw).
- **Plots:** `figures/fig01_nuclear_vs_electronic.png` · script `inqview/report1/fig01_nuclear_vs_electronic.py` · not in remake track.
- **Purpose:** Sₙ vs Sₑ vs energy; crossover ≈5 keV/u; motivates Sₑ-only focus.
- **Changes I want:**
  -

## A2. `fig:regime-diagram` — Stopping regimes + arbitrary stopping curve  ·  REWORK + NEW sub-plot
- **Section:** Literature Review §5.1 "Regime analysis" (line 260). **Env:** `figure*` (2-col), two-up: `[0.52 tw] (a) regime map` + `[0.46 tw] (b) stopping curve`.
- **Plots:**
  - (a) `figures/remake/fig_schematic_regime_map.png` · script `draft3/scripts/make_fig_schematic_regime.py` · REWORK.
  - (b) `figures/remake/fig_schematic_stopping_curve.png` · script (verify: likely same `make_fig_schematic_regime.py`) · NEW (added per spec to show arbitrary S(v) shape).
- **Spec notes:** make regime diagram ~1.2× larger; (b) shows linear rise → Bragg peak → Bethe v⁻²lnv² fall. Full per-plot grilling deferred to Phase 3+.
- **Changes I want:**
  -

## A3. `fig:pseudopotential` — Classical projectile pseudopotential  ·  RELOOK
- **Section:** Theory §"Classical Projectiles in TDDFT" (line 513) — *home section being dissolved per structure.md; reroute TBD*. **Env:** `figure` (1-col), `\figwidthmed`.
- **Plots:** `figures/fig12_pseudopotential.png` · script `inqview/report1/fig12_pseudopotential.py` · not in remake track.
- **Caption flagged in tex:** "TODO: CHANGE, too verbose, boring."
- **Changes I want:**
  -

## A4. `fig:leed-geometry` — Joint setup schematic (3 setups)  ·  NEW (redesigned 2-row)
- **Section:** Methods §"Simulation Setups" (line 688). **Env:** `figure*` (2-col). Row 1: coronene full width (`0.85\linewidth`, centred). Row 2: free-prop `[0.48 tw]` + jellium `[0.48 tw]`.
- **Plots:**
  - `figures/remake/fig_schematic_coronene.png` · `make_fig_schematic_coronene.py`
  - `figures/remake/fig_schematic_free.png` · `make_fig_schematic_free.py`
  - `figures/remake/fig_schematic_jellium.png` · `make_fig_schematic_jellium.py`
- **Spec notes:** annotate Lₓ,L_y,L_z and −L/2…L/2 axis; jellium WP at 4σ clearance; coronene WP start per Tsubonoya coords. Master-parameters table supplies dimensions (coronene 34.8×34.8×59.9, others 50³). t_IFW defined where centroid+4σ_r = L/2.
- **Changes I want:**
  -

## A5. `fig:free-wp` — Free wave-packet spreading  ·  RELOOK
- **Section:** Results §"Free propagation" (line 1080). **Env:** `figure` (1-col), `\figwidthmed`. Two sub-panels (a) σ_r(t) vs analytic, (b) residual.
- **Plots:** `figures/fig_free_wp_spreading.png` · script `inqview/report1/fig_free_wp_spreading.py` · not in remake track.
- **Changes I want:**
  -

## A6. `fig:leed-backscatter` — Coronene LEED panel (4 plots)  ·  REWORK + EXTERNAL
- **Section:** Results §"Coronene LEED study" (line 1109). **Env:** `figure*` (2-col), 2×2 minipages (`0.48 tw` each; (c) inner `0.75\linewidth`).
- **Plots:**
  - (a) `figures/remake/fig_coronene_target.png` · `make_fig_coronene_target.py` · GS density `run_save_gs_paper_replica`. REWORK (C–C crosshair cyan).
  - (b) `figures/remake/fig_leed_backscatter_centre.png` · `make_fig_leed_backscatter_centre.py` · `run_propagate_paper_replica` screen14 step330, log. REWORK.
  - (c) `figures/tsubonoya-comparison.png` · EXTERNAL (user adds, own linear colourbar, no script).
  - (d) `figures/remake/fig_leed_transmission_fft.png` · `make_fig_leed_transmission.py` · screen07 step330, FFT, log, with analytic |F(q)|² overlay (cyan circles). REWORK (FFT dots distinct colour + legend; log scale; optional point filter).
- **Note:** (b) and (d) share style; user wanted backscatter + transmission to share a colourbar — confirm in per-plot phase.
- **Changes I want:**
  -

## A7. `fig:delta-density-comparison` — Jellium induced density (4 rows × 2)  ·  REWORK
- **Section:** Results §"Jellium … representative case" (line 1185). **Env:** `figure*` (2-col), 4 rows; each row `[0.48 tw] 2D slice` + `[0.48 tw] z-profile`.
- **Run:** E25 case study minus matched free (`run_wp_n162_L50_E25_sigma1_v2` − `run_free_wp_L50_E25_sigma1_v2`).
- **Plots:**
  - Left col `fig_density_diff_2d_t{1,2,3,4}.png` · `make_fig_density_diff_2d.py` — y-midplane Δn_induced(x,z,t).
  - Right col `fig_density_profile_t{1,2,3,4}.png` · `make_fig_density_profile.py` — z-profile; solid dark-red = (jel−free), dashed blue = jel Δn, dash-dot grey = free Δn; shared y-range; explicit scale factor; legend in first only.
- **Caveat (from spec):** the script that made the old `fig_DD1_density_diff_grid.png` may be wrong/changed — verify formula & timesteps when reworking.
- **Changes I want:**
  -

## A8. `fig:energy-decomposition` — Jellium energetics (3 plots)  ·  REWORK
- **Section:** Results §jellium (after density, line 1252). **Env:** `figure*` (2-col). Row 1: `[0.48] system` + `[0.48] WP`. Row 2: GS decomposition full-width (`0.95\linewidth`).
- **Run:** E25 case study.
- **Plots:**
  - (a) `fig_energy_decomp_system.png` · `make_fig_energy_decomp_system.py` — ΔE_total (conserved), ΔE_kin, ΔE_H, ΔE_xc vs t. Source `observables.csv`.
  - (b) `fig_energy_decomp_wp.png` · `make_fig_energy_decomp_wp.py` — ΔE_kin^WP, directed Δ⟨p⟩²/2m (stopping signal), spreading Δσ_p²/2m. Source `wp_momentum_stats.csv`.
  - (c) `fig_gs_decomposition.png` · `make_fig_gs_decomposition.py` — δn_i^GS(t_end) bars + charge-balance sidebar. Source `gs_projected_occupations/…t_end.csv`.
- **Changes I want:**
  -

## A9. `fig:plasmon-fft` — Mechanism: plasmon + loss function (3 plots)  ·  REWORK
- **Section:** Results §jellium mechanism (line 1279). **Env:** `figure*` (2-col). Row 1: plasmon FFT (`0.85 tw`). Row 2: `[0.48] loss 2D` + `[0.48] loss 1D`.
- **Run: `run_plasmon_n162_L50_E15` (E15 long run) — NOT E25.** State in caption.
- **Plots:**
  - (a) `fig_plasmon_fft.png` · `make_fig_plasmon_fft.py` — |FFT[n_q_m(t)]| m=1,2,3, log-y, Bohm–Gross dotted lines. Source `…/n_q_spectrum.csv`.
  - (b) `fig_loss_function_2d.png` · `make_fig_loss_function.py` — L(q_z,ω) heatmap, Bohm–Gross (white dashed), e–h continuum (cyan dotted), Landau cutoff q_c (red dash-dot). Source `n_q_vs_time.csv` (Hann window). e–h transitions from E25 eigenvalues.
  - (c) `fig_loss_function_1d.png` · `make_fig_loss_function.py` — 1D cuts m=1,2,3, log-y, ω_p marked.
- **Changes I want:**
  -

## A10. `fig:momentum-before-after` — WP momentum redistribution (2 plots)  ·  REWORK
- **Section:** Results §jellium momentum (line 1315). **Env:** `figure*` (2-col), `[0.48] 1D` + `[0.48] 2D diff`. (Full before/after 2D → appendix B4.)
- **Run:** E25 case study.
- **Plots:**
  - (a) `fig_momentum_1d.png` · `make_fig_momentum_1d.py` — n_WP(|k|) at t=0 (dark red) & t_end (blue); dashed line at k₀=1.356; **central momenta of before & after both marked** (spec). Source `momentum_distribution.csv`.
  - (b) `fig_momentum_2d_diff.png` · `make_fig_momentum_2d.py` — Δ|ψ̃|² in (k_z,k_⊥), diverging cmap (after − before).
- **Changes I want:**
  -

## A11. `fig:master-stopping` — Master stopping power (2 plots)  ·  REWORK
- **Section:** Results §"Parameter sweep and the classical limit" (line 1340). **Env:** `figure*` (2-col), `[0.48] low-dens` + `[0.48] high-dens`.
- **Headline metric: S₂ (directed KE).** Filled = v2 (preferred), hollow = v1 — **report text must explain this convention** (legend doesn't).
- **Plots:**
  - (a) `fig_master_stopping_lowdens.png` · `make_fig_master_stopping.py` — r_s=5.69 (L=50): classical (squares), WP σ=1 (circles), WP σ=5 (triangles), supplementary σ=0.5/3/8 at E=100, Bethe v⁻² (grey dashed).
  - (b) `fig_master_stopping_highdens.png` · same script — r_s=3.41 (L=30): same layout; high-density family runs.
- **Data:** `stopping_power_data.py` pipeline. v2 preferred over v1 where both exist.
- **Changes I want:**
  -

## A12. `fig:definition-comparison` — Stopping-definition comparison  ·  REWORK
- **Section:** Results §stopping definitions (line 1359). **Env:** `figure` (1-col), `0.85\linewidth`. Single combined plot.
- **Run:** low-density only. Classical transient skip (first 10%); interference-free window enforced for both classical & WP.
- **Plot:** `fig_stopping_defs_combined.png` · `make_fig_stopping_defs_combined.py` — S₂ (directed, red filled, S>0), S₃ (total p-space KE, purple hollow |S|), S₄ (spreading, blue hollow |S|), classical (black squares), Bethe (grey dashed).
- **Changes I want:**
  -

## A13. `fig:gantt` — Report-2 plan Gantt  ·  RELOOK
- **Section:** §"Outlook and plan for Report 2" (line 1417). **Env:** `figure*`, `\figwidthwide`.
- **Plots:** `figures/fig_gantt.png` · script `inqview/report1/fig_gantt.py` · not in remake track.
- **Changes I want:**
  -

---

# PART B — Appendix / supplementary panels

## B1. `fig:scenario-ab` — KE decomposition Scenario A vs B  ·  RELOOK
- **Section:** Appendix "Kinetic-energy decomposition" (line 1489). **Env:** `figure*`, `\figwidthwide`.
- **Plots:** `figures/fig13_scenario_ab.png` · `inqview/report1/fig13_scenario_ab.py` · not in remake track.
- **Changes I want:**
  -

## B2. `fig:jellium-gs` — Ground-state density slices  ·  RELOOK
- **Section:** Appendix "Supplementary figures" (line 1584). **Env:** `figure*`, `\figwidthmed`. (a) N=162 closed shell, (b) N=138 partial shell.
- **Plots:** `figures/fig_jellium_gs.png` · `inqview/report1/fig_jellium_gs.py` · not in remake track.
- **Changes I want:**
  -

## B3. `fig:leed-ccbond` — Back-scattering LEED, C–C bond impact  ·  REWORK
- **Section:** Appendix (line 1591). **Env:** `figure` (1-col), `0.75\linewidth`.
- **Plot:** `figures/remake/fig_leed_backscatter_ccbond.png` · `make_fig_leed_backscatter_ccbond.py` · `run_cc_bond` screen14 step330, log floor 5e-6.
- **Changes I want:**
  -

## B4. `fig:leed-validation` — LEED validation (4-panel)  ·  REWORK
- **Section:** Appendix (line 1598). **Env:** `figure` (1-col), `0.85\linewidth`.
- **Plot:** `figures/remake/fig_leed_validation.png` · `make_fig_leed_validation.py` — (a) sim backscatter FFT, (b) analytic |F(q)|², (c) GS density FFT, (d) azimuthal I(θ). sim↔GS-FFT correlation r=0.80.
- **Changes I want:**
  -

## B5. `fig:momentum-2d-full` — Full 2D momentum before/after  ·  REWORK
- **Section:** Appendix (line 1605). **Env:** `figure*` (2-col), `[0.48] before` + `[0.48] after`, shared colourbar.
- **Run:** E25 case study.
- **Plots:** `fig_momentum_2d_before.png`, `fig_momentum_2d_after.png` · `make_fig_momentum_2d.py` — |ψ̃|² in (k_z,k_⊥) at t=0 and t_end.
- **Changes I want:**
  -

## B6. `fig:plasmon-realspace` — Real-space density wake  ·  REWORK
- **Section:** Appendix (line 1617). **Env:** `figure` (1-col), `0.85\linewidth`.
- **Run: E15 long run.**
- **Plot:** `fig_plasmon_realspace.png` · `make_fig_plasmon_fft.py` (or dedicated) — Σ_y Δn(x,z) at t=1000 a.u., SymLogNorm linthresh=1e-4·vmax.
- **Changes I want:**
  -

---

# PART C — High-density mirror track (do after low-density locked)

The jellium data panels (A7–A10) have **L=30 high-density (r_s=3.41) twins**.
Master plot A11(b) already consumes the high-density family. Decide per panel
whether a full high-density mirror figure is needed in the report or only the
master plot. Candidate runs: `run_wp_n162_L30_E{50,100,200,300}_highdens_sigma1_v2`.
Legacy draft3 scripts exist (`*_highdens.py` in `inqview/report1/`).

- **Changes I want / which panels to mirror:**
  -

---

# Open items for you to resolve while editing

1. **Font sizes** — pick the canonical on-page band (8–12 pt) and the exact
   label/tick/legend values; I'll set them once in `STYLE_CONFIG`.
2. **draft5 figure output dir** — confirm plots are written to
   `drafts/draft5/figures/remake/` (vs continuing to source from `draft3/`).
3. **Per-panel annotation cuts** — list unnecessary annotations to remove.
4. **High-density scope** (Part C) — full mirror panels or master-plot only.
5. **Tsubonoya panel** — confirm you place it manually with its own linear colourbar (no script).
6. **Backscatter/transmission shared colourbar** (A6) — confirm whether (b) & (d) must share one.
