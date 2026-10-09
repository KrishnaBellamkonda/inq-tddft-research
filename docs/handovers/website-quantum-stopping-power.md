# Handover: Project website — Quantum Description of Electronic Stopping Power

Rolling file. Task: build a CV-linked single-page site showcasing the quantum
(wavepacket) correction to electronic stopping power. Plan:
`/local/data/public/skcb2/tddft/docs/plans/website-quantum-stopping-power.md`.

---

## Milestone: 2026-10-09 — scaffold complete, real plots + hero animation wired

### DONE (verified)
- **Site scaffold** at `/local/data/public/skcb2/tddft/website/`:
  `index.html` (8 sections), `css/style.css`, `js/main.js`, `.nojekyll`.
  Total size 1.9 MB.
- **All 24 real plots copied + embedded**, each in its section
  (`website/assets/figures/{setup,worked,extrapolation,result,secondary,appendix}`).
  Every `src=` verified to resolve (checked programmatically).
- **Hero animation BUILT**: `website/assets/hero/hero_energy_E100_Lz45.{mp4,gif}`
  + `_poster.png`. Script: `website/build_hero_animation.py`. It reuses the
  `slab_stopping` skill loaders (`.claude/skills/slab-stopping-power/slab_stopping.py`:
  `extensive_kinetic_total`, `concat_obs`) and `E_GS[45]=386.538854026 Ha` — so the
  animated corrected E(t) curves are IDENTICAL physics to the vetted
  `results_gallery/.../E100_v2p71/Lz45/energy_vs_time.png`. MP4 is 32 KB.
- **Headline numbers** (section 05 table) grounded in
  `ResearchProject/systems/localised_jellium/scripts/slab_sE_sweep/sE_results.csv`:
  bulk WP/classical ratio E=31→0.24, 54→0.33, 100→0.56, 167→0.64, 276→0.65.
- **Publish-ready meta**: OpenGraph/twitter card (poster as og:image), inline SVG
  favicon, `.nojekyll` (lets GitHub Pages serve `assets/`).

### DESIGN DECISIONS (from grill-with-docs brainstorm — see plan)
- Headline = method + energy-dependent quantum correction (largest at low E,
  ~¼ of classical near 30 eV; → classical limit at high E). NOT a flat "half".
- Section order: Hero → Motivation → Idea → Worked example → Trustworthy
  (S vs 1/Lz) → Result S(E) → Secondary emission → collapsed appendix.
- Hero = option B (animated E(t)); density split-screen "C" = deferred upgrade.
- Stack = hand-built static HTML/CSS/JS on GitHub Pages.
- Location = `website/` subfolder in this repo.
- Publish target repo: `github.com/KrishnaBellamkonda/quantum-stopping-power`
  → URL `krishnabellamkonda.github.io/quantum-stopping-power`. USER pushes (their account).

### NOT DONE / waiting on user
- **All prose is placeholder** (violet dashed `.ph` boxes). User writes copy after
  reviewing the plots themselves. Hero headline sentence intentionally deferred.
- **DNA-damage ~10 eV plot** (section 01): not in repo; gold placeholder. User has
  the file → drop into `website/assets/figures/motivation/`, then wire the `<img>`.
- **§3 "The idea" density pair — DONE via option (b), 2026-10-09.** Side-by-side
  web MP4s (`website/assets/idea/idea_{classical,wp}_sigma1_density.mp4`, 65–70 KB
  each, re-encoded from 27 MB GIFs via ffmpeg). Source: matched
  `muon_mass_fork/effmass_sigma1` classical+WP runs — VERIFIED matched (both
  σ_WP=1, box 40×40×80, dx 0.333, v=2.711, E≈210 eV, same GS
  `shared_gs/slab_n52_L40x40x80_dx0p333`). **CAVEAT to tell user:** projectile
  effective mass is **2.1 m_e** (mass-fork study), not a bare electron, and
  E≈210 eV ≠ hero's 100 eV. Labelled "illustrative" in the HTML caption note.
  A heavier WP disperses *less*, so this UNDERSTATES the quantum visual, not
  overstates it. Swap for a true-electron E=100 density pair if/when (a) is run.
- **Density split-screen animation** ("C" hero upgrade; per-energy appendix): placeholder.
  Clean path = re-run ONE matched E=100/Lz45 classical+WP pair with frame saving
  (`LJ_SAVE_EVERY>0`), then composite. NOT launched — expensive GPU + disk-fill
  risk (density frames ~2.6 GB/run fill sda1; see memory
  `reference_sim_frames_fill_disk_crash_claude`). NEEDS user go-ahead on scope.
  - Candidate EXISTING slab density GIFs (localised_jellium) if an *illustrative*
    (non-E100) pair is acceptable instead of a matched re-run:
    `hypotheses/muon_mass_fork/effmass_sigma1_classical_run_figs/classical_total_density.gif`
    + `effmass_sigma1_wp_run_figs/sigma1_total_density.gif` — NOT verified to share
    box/energy; do NOT present as the E=100 pair without checking run_summary.txt.

### HOW TO VIEW
```
cd /local/data/public/skcb2/tddft/website
/local/data/public/skcb2/tddft/venv/bin/python3 -m http.server 8000   # → http://localhost:8000
```

### HOW TO PUBLISH (user does this)
1. Sign into GitHub as KrishnaBellamkonda.
2. Push `website/` contents to repo `quantum-stopping-power` root (or copy folder).
3. Settings → Pages → deploy from branch (root). `.nojekyll` already present.

### VALIDATION STATUS
- Hero physics provenance: reuses single-source-of-truth loaders (not re-derived). ✓
- Asset links: all resolve. ✓
- Numbers in section-05 table: copied from `sE_results.csv`. ✓
- No scientific prose asserted by assistant — all author-written placeholders.
