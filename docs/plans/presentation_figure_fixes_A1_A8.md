# Plan: presentation figure fixes A1–A8 (interim-presentation deck)

Disambiguated 2026-06-02 via `/grill-with-docs`. Target deck:
`docs/presentations/drafts/KrishnaBellamkonda_InterimPresentation_Draft2 (2).pptx`
(10.0×5.625 in, 24 slides). **The user's slide numbers are unreliable** (a
~+3/+4 offset vs the deck's 1-indexed order); tasks are matched by *content*,
not slide number. The deck only *consumes* the regenerated PNGs.

## Glossary of resolved terms (read before editing)

- **Stopping-power definition clash.** `stopping_power_data.py` numbers them
  `S1_eV_per_bohr = momentum/expectation` (−Δ(⟨p⟩²+σ_p²)/2m/Δz, from
  `wp_momentum_stats.csv`) and `S2_eV_per_bohr = KS-orbital energy`
  (−Δ⟨E_KS⟩/Δz). `make_fig_stopping_defs_combined.py` uses the **opposite**
  numbering. **Always edit by physical meaning, never by S-number.**
  - "KS orbital energy-loss def" (A7 #1) = `S2_eV_per_bohr` ← vaxis fig uses this now.
  - "expectation-value def ⟨p⟩²/2m+σ_p²/2m" (A7 #2) = `S1_eV_per_bohr`.
  - A7's literal "⟨p²⟩/2m + σ_p²/2m" is redundant (⟨p²⟩=⟨p⟩²+σ_p²); intent is
    the Yao–Schleife drift+spread split = total kinetic = `S1`.
- **"Target plot"** = `fig_coronene_target.py` (coronene GS xy-density heatmap,
  "Target centre"/"Target C--C bond" markers). Confirmed by user.
- **r_s choice** = 5.69 Bohr (low-density L=50, Cs-like). Bold Cs in the metals
  table, footnote r_s=5.69. High-density L=30 (r_s=3.41, Li-like) is secondary.

## Task → producing script → action

| Task | Script | Action |
|---|---|---|
| A1 | `report1/make_fig_free_wp_panel.py` | Right panel: absolute E → ΔE (E−E₀). Combine σ(t)+ΔE into ONE figure, two vertically-stacked subplots, shared time axis. Match interference-free shading across both. |
| A2 | `report1/make_fig_coronene_setup.py` | Move k₀ label+arrow to WP centre at t₀ pointing −z (entry point, not z=+15). Remove floating WP-formula box + overlapping width/energy annotation boxes; use clean text labels *outside* the density-plot boundary. |
| A3 | `report1/make_fig_gs_decomposition.py` | (a) Fix clipped +496.8 bar label (raise y-upper or log). (b) Restructure `GridSpec(1,2)` (charge-balance sidebar *beside*) → charge balance as a subplot *below* the occupation bars, shared width (`GridSpec(2,1)`). |
| A4 | `report1/fig_coronene_target.py` | Report-quality polish: fonts, line weights, axis labels (Bohr units), clean ticks, well-placed legend, no clutter. Reconcile with `_shared_style` (note: currently uses `bbox_inches="tight"` — fixed-dim pitfall). |
| A5 | NEW (write into remake dir) | (a) Jellium schematic: uniform positive background slab, `n⁺ = const`, sample electron-density curve. (b) Metals table (Li 3.25, Na 3.93, K 4.86, Rb 5.20, Cs 5.62, Al 2.07): cols Metal | r_s (Bohr) | n (e/Bohr³); **bold Cs**; footnote "simulation r_s = 5.69". |
| A6 | `report1/make_fig_momentum_2d.py` + `make_fig_momentum_1d.py` | Jellium E25 run. Extend k_z to full range incl. negative (un-clip `set_xlim(K0±2)`). Compute F/B ratio = ∫_{k_z>0}Δ \|ψ̃\|² vs \|∫_{k_z<0}Δ\|ψ̃\|²\| ; print ratio on figure. |
| A7 | `draft5/scripts/make_fig_master_stopping_vaxis.py` | **TWO plots, one per WP definition.** Each: WP σ=1 trace (that def) + simulated loss-function (red) + classical datapoints. Keep analytical Lindhard box-q (grey) anchor. Drop σ=5/0.5/3/8. r_s=5.69 (L=50). Confirmed WP currently uses KS-energy def (`S2`). |
| A8 | — resolved | Analytical reference = Lindhard box-q (already in vaxis fig; docstring deliberately omits Bethe). No SRIM/Echenique addition. |

## Data availability (verified)
- A6: `wavefunction_wp` VTIs are full-box FFTs → full k-space available; only the
  plot xlim clips it. 1D source = `momentum_distribution.csv` (E25 run).
- A7: `collect_L50_data()` populates BOTH `S1_eV_per_bohr` and `S2_eV_per_bohr`;
  `wp_momentum_stats.csv` carries `sigma_pz2`, `e_kin_ha`, `pz_mean`.

## Validation (per project rules — known-case before integrate)
- A1: ΔE(t=0)=0 by construction; drift magnitude printed.
- A3: bar-sum sanity (Σ_occ + Σ_virt + untracked) unchanged after restructure.
- A5: table n = 3/(4π r_s³) per row (known-case: Cs r_s=5.62 → n≈0.00134 e/Bohr³).
- A6: ∫all Δ ≈ 0 (norm conserved); F/B ratio finite; t=0 peak at |k|=k₀.
- A7: S1,S2 curves both monotone-ish; loss-fn peak v & classical points unchanged.

## Notes
- venv python only: `/local/data/public/skcb2/tddft/venv/bin/python3`.
- Do NOT preview/Read generated PNGs (user previews).
- Deviation from grill-with-docs default: glossary captured here in `docs/plans/`
  (project file-placement rule) rather than a root `CONTEXT.md`.
