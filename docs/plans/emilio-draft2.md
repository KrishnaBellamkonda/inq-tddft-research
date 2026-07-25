# Plan — Emilio meeting deck, DRAFT 2

Resolved via grill-with-docs (2026-06-26). Deliverable: **new file**
`docs/reports/26-06-2026-meeting-emilio/emilio_meeting_draft2.pptx`, evolved from
the draft-1 builder. Draft 1 stays untouched as a fallback. **User edits draft 2
manually after this build.**

Scope = three additive tasks on top of draft 1 (whose FB-corrected plots are
KEPT as-is). Nothing in draft 1 is redone.

## Task 1 — add new plots (additive; draft-1 plots unchanged)
All four are new analysis plots; all land in **Section 3 (slab)**.

| New plot | Source | Spec |
|---|---|---|
| Electric field | `inqview.analysis.efield.electric_field` on a WP density frame | `|E(r)|` on an **xz slice**, mid-transit frame (qsp_phase2 `p2_wp`); canonical loader + colorbar idiom; dashed slab/CAP extents |
| KL divergence | `inqview.analysis.kl_divergence` over WP density | time series (KL of WP density vs reference, over the run) |
| KS-eigenstate energies | `state_energies.csv` | static bar chart, ΔE_i and ±WP variants |
| Momentum \|k\|-vs-time | `momentum_distribution.csv` | radial \|k\|-vs-time carpet, shown **alongside** the FB-017 k_z–k⊥ map (complementary: time vs direction) |

## Task 2 — gifs
- **BUILD (the §3.3 centerpiece):** 4 xz density gifs — total-density xz (WP,
  classical) + Δn xz (WP, classical) — from the qsp_phase2 VTIs, **frame-strided
  (~every 4th, ≈50 frames)** to dodge the draft-1 render timeout. Laid out
  **2×2** on the comparison slide. Keep the draft-1 `wake_compare` still as an
  appendix fallback.
- **BUILD:** animate the norm-vs-time plot into a gif (user override of the
  keep-static recommendation).
- **REUSE (no build):** KS-eigenstate-energy gifs (`ks_energies_delta(_no_wp).gif`,
  both runs); momentum `momentum_distribution.gif` (\|k\|-vs-time).
- **KEEP as-is:** the CAP/mask 1D z-profile gifs (NOT swapped to xz — user chose
  to keep the 1D versions).
- pptx note: python-pptx `add_picture` embeds `.gif`; it animates in slideshow.

## Task 3 — section divider slides
Keep the current natural flow / one deck. Between sections, a **plain
section-name slide** (title only), retitling the existing FB-009 dividers to
EXACTLY:
1. **Gaussian Potentials**
2. **Absorbing boundary conditions**
3. **Localised Jellium slab**

## Placement summary
- §1, §2: no new plots; dividers retitled; CAP/mask gifs unchanged.
- §3 absorbs ALL new content: comparison xz gifs (replace the still in-place,
  still→appendix), norm gif, KS-energy gifs, momentum \|k\|-vs-time (gif + carpet
  plot), E-field plot, KL-divergence plot — placed near their related draft-1
  slides.

## Open / carried from draft 1 (not draft-2 scope)
- Two physics flags (ΔE_WP negative; WP norm-loss ~24× classical) — still need
  user resolution.
- −Im[1/ε] absolute normalisation flagged uncertain.
- Wide-electron + quantum-S(v)-sweep slides remain placeholders (no runs).
