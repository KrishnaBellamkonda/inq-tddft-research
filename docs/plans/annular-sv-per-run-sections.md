# Plan: per-run deep-dive sections in the annular-tube report notebook

**Goal (user request, 2026-06-30):** In `hypotheses/annular_sv/annular_sv_report.ipynb`,
add **one section per projectile run** (9 classical + 1 WP) showing (a) the full
**matrix of density visualisations** and (b) the **high-value per-run observables**
from `docs/notes/plots_examples.md`, so each run can be understood in full.

## Runs (10)
`rs{6,4,2}_v{0p15,0p30,0p45}` (classical) + `wp_rs6_v0p30` (WP). Each results tree at
`annular_sv/<label>/results/<tag>/` (raw/ + run_summary.txt). τ=30–45 a.u. (light
projectile, by design); 301 density frames per stream.

## Geometry (NOT a slab — must stay tube-correct)
Periodic annular tube: wall at radius |x|∈[R_in=5, R_out=13], hollow bore |x|<5,
axis ∥ z, **no slab faces, no CAP**. xz-slice overlays = **vertical wall lines x=±5,±13**
(not the slab library's horizontal z-faces). S = **initial drag** −d(KE_ion)/ds over the
early vz≥0.85·v0 window (light-projectile rule), never the slab ΔE/L_z method.

## What each per-run section contains
1. **Density matrix** {density, Δn=n(t)−n(0), Δn=n(t+dt)−n(t)} × {total[, wp, bath]}
   — classical: 3 GIFs (total); WP: 9 GIFs (bath=total−wp). Wall-radius markers.
2. **z–t carpets** (total n, Δ-vs-0, per-step Δ).
3. **Energy decomposition** (total + kinetic/Hartree/xc) — pipeline.
4. **Induced current_z(t) + FFT**; **dipole**; **FFT-pipeline panel** — pipeline.
5. **Momentum** 1D incl/excl WP + difference-over-time (WP only) — pipeline.
6. **KL metric** KL(P_t‖P₀) (WP only) — pipeline.
7. **Initial-drag stopping**: KE_ion vs path, early window highlighted, fitted S(v0).

## Honest data-limitation notes (do NOT fabricate)
- **KS eigen-energy GIFs** (plots_examples item 4): `state_energies.csv` not stored
  by these runs → unavailable; note it.
- **GS KS-excitation decomposition / eigenvalues**: `eigenvalues.csv` not retrofitted.
- **E-field panel**: no pipeline phase; omit/note.
- **2D (k_z,k_⊥) scattering**: only 1D |k| stored; the 1D view is what exists.

## Architecture
- New `hypotheses/annular_sv/per_run.py` — tube-aware generator (self-contained,
  uses canonical `load_vti`; reuses campaign `extract_S` logic for the stopping plot).
  `generate(label, rs, v0, run_dir, results_dir, figroot) -> {group: [(cap, path)]}`.
- Pre-generate ALL figures to `hypotheses/annular_sv/per_run_figs/<label>/` (heavy
  GIF rendering once), then `build_report.py` embeds them **path-referenced** so the
  executed `.ipynb` stays fast + small (run-notebook img() convention).
- Extend `build_report.py`: after the synthesis sections, append a
  "## Per-run deep dives" section + one subsection per run.

## Validation
- Pipeline already verified to run on a WP run (figures land in results/analysis).
- Done = notebook executes 0 errors; every run subsection shows real figures or a
  labelled "not stored this run" note; density matrix has correct wall markers
  (verify via the GS tube geometry, not by eye on a flipped axis).
- Update `docs/handovers/cylindrical-jellium-projectile.md` + test-catalogue row.
