# Plan: Project website — Quantum Description of Electronic Stopping Power

Status: scaffold built (2026-10-08), prose placeholders pending user's own copy.

## Purpose
CV-linked project site to capture a PhD supervisor's attention and give a clear,
credible overview of the research. Headline = **method novelty (quantum
wavepacket projectile in rt-TDDFT) + the quantum correction to stopping power**.

## Decisions (from grill-with-docs brainstorm)
1. **Intent:** method-and-result site; the quantum correction is the whole story.
2. **Headline physics (grounded in `sE_results.csv`):** the correction is
   *energy-dependent*, largest in the low-energy (biological ~10 eV) regime and
   shrinking toward the classical limit. Bulk-extrapolated WP/classical ratio:
   E=31→0.24, 54→0.33, 100→0.56, 167→0.64, 276→0.65 eV. Frame as "suppressed by
   up to ~4× near 30 eV", NOT a flat "half". Exact wording: user writes after
   viewing plots.
3. **Architecture (one page):** Hero → §2 Why it matters (DNA ~10 eV) → §3 The
   idea (point charge vs wavepacket) → §4 Worked example (E(t)) → §5 Is it
   trustworthy? (S vs 1/Lz extrapolation) → §6 Result S(E) → §7 Secondary
   emission → §8 collapsed appendix (per-energy explorer).
4. **Hero:** option B now = animated corrected E(t) for E=100/Lz45 (classical vs
   WP). C (density split-screen) = later upgrade, needs a frame-saving re-run.
5. **Stack:** hand-built static HTML/CSS/JS, GitHub Pages.
6. **Location:** `website/` subfolder in this repo (compact).
7. **Publish:** push to `github.com/KrishnaBellamkonda/quantum-stopping-power`
   → `krishnabellamkonda.github.io/quantum-stopping-power`. User signs in + pushes.

## Caveats to keep visible (credibility)
- Each bulk S point is a quadratic extrapolation from only 3 slab lengths
  (Lz=15/30/45; 2 pts at E=31). Show the S-vs-1/Lz fits openly in §5.
- Classical S(E) is non-monotonic (peaks ~54 eV) — do not hide.

## Asset status
- **Have (real, embedded):** hero animation (built), setup schematic, E100/Lz45
  worked E(t), 5× S-vs-1/Lz panels, S(E) bulk + per-Lz, secondary emission,
  per-energy appendix E(t).
- **Missing / placeholder:** DNA-damage ~10 eV annotated plot (user has it, not in
  repo); density split-screen animation (C upgrade, needs re-run with frames).

## Hero provenance
`website/build_hero_animation.py` reuses `slab_stopping` skill loaders
(extensive_kinetic_total + norm-correction, +E_PB classical) and E_GS[45]=386.54
Ha — identical physics to the vetted `energy_vs_time.png`.
