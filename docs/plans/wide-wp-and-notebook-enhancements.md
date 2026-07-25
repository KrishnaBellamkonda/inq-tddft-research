# Plan: wide-WP twin run + twin-analysis notebook upgrades

**Status:** planned 2026-07-15. Depends on the twin-run skills (built) + the
pairwise decomposition (verified, `reference_twin_pairwise_decomposition`).
Parent: `docs/plans/twin-run-rung2-dynamic-spec.md`.

## Motivation

The σ=0.5, k0=1 pair disperses 10× before reaching the slab, so classical↔WP
divergence is dominated by dispersion. A **wider, faster** WP stays coherent long
enough to actually **enter the slab**, making the classical-vs-WP difference
meaningful *during traversal*, not just free spreading.

## Run set — wide-WP twin pair

| param | value | rationale |
|---|---|---|
| σ_WP | 2 Bohr | ~4× slower dispersion than σ=0.5 |
| mass | 1 (m_e) | broadening over a fixed distance is **mass-independent** (D/(2k0σ_ρ²), mass cancels); mass fixed at 1 |
| k0 | 4 Bohr⁻¹ | most coherent within aliasing: k0+3σ_p = 4+1.06 = 5.06 < Nyquist 6.28 |
| dt | 0.025 au | k0=4 packet is fast; dt·E_cut ≈ 0.5 (stable) |
| N_STEPS | 300 | v0=4 → slab centre (24.5 Bohr) in ~245 steps; 300 covers "fully entered" |
| launch_z | −24.5 | density extent 3σ_ρ≈4.2 Bohr → tail at −20.3, clear of face −12.5 |
| GS | reuse `h2/gs_p2_lz120` | σ-independent slab GS; NO new GS |
| broadening | ~25% at face, ~80% at centre | acceptable (σ=0.5 was 10×); NOT the 4% target (infeasible at σ=2 — see below) |

- **Classical twin:** `proj_dyn` with `LJ_SIGMA=2 LJ_K0=4 LJ_MASS=1 LJ_DT=0.025
  LJ_N_STEPS=300` (Gaussian charge σ_pot=2/√2=1.414, moving Projectile Ehrenfest).
- **WP twin:** `phase5_wp` with the same, `LJ_SAVE_EVERY=25` (frames).
- **Pre-launch gate:** cutoff/aliasing guard (`reference_cutoff_aliasing_guard`) —
  MANDATORY. One GPU (GPU0), sequential, ~30–40 min each after the current 200-step run.
- **Why not 4% broadening (user's original target):** dispersion over a fixed
  traversal is D/(2k0σ_ρ²) — mass cancels; ≤4% at σ=2 needs k0≥10.5, far above the
  aliasing limit. ≤4% would require σ_WP≈4. σ=2/k0=4 is the practical compromise.

## Notebook upgrades (apply to EVERY twin pair, incl. the σ=0.5 pair)

Build into `twin_notebook_builder.py` (skill-local) + one run.cpp change:

1. **Density maps n(r,t)** — 2D (z–x) slices / z-lineouts of each run's density at a
   few times, from the saved VTI frames (via `inqview.load_vti`, physical order — NO
   fftshift, `vti-coordinate-mapping`). **Requires `proj_dyn` to save density frames**
   (add `LJ_SAVE_EVERY` + frame writer to `proj_dyn/run.cpp`, mirroring `phase5_wp`).
2. **Δn(r,t) = n_wp − n_classical** — the classical-vs-WP density difference map
   (shows the WP itself + the slab-polarisation difference). Shared colorbar +
   linear&log per `feedback_shared_colorbar_rule`.
3. **Pairwise-energy GIF** — a slow animation over steps: per frame, the pairwise
   Δ(WP−cl) (or absolute terms) as a bar plot, so the viewer watches E_PP, E_PS,
   E_PB, … diverge as the WP enters the slab. Saved as `.gif` (matplotlib animation /
   imageio). One frame per saved step.
4. **WP−classical bar plot** — a static bar chart of ΔE per component (ΔKE_total,
   ΔE_SS, ΔE_PP, ΔE_PS, ΔE_SB, ΔE_PB, ΔE_xc) at a chosen step (default final), the
   at-a-glance "where does the quantum energy difference live" figure.
5. Keep the existing pairwise ledger + gauge test + conservation + dispersion σ(t).

## Run notebooks (per run, deep-dive)

For each run (classical + WP of each pair), build a per-run deep-dive notebook via
the `run-notebook` skill (all inqview pipeline phases the run supports + the house
narrative). Output in `hypotheses/twin_dynamics/`.

## Build order

1. (No GPU) `proj_dyn/run.cpp` — add density-frame saving.
2. (No GPU) `twin_notebook_builder.py` — density maps + Δn, pairwise GIF, bar plot.
3. (No GPU) rebuild the σ=0.5 200-step notebook to validate the new cells.
4. (GPU, after 200-step) cutoff guard → wide-WP twin pair (σ=2, k0=4, 300 steps).
5. (No GPU) wide-WP analysis notebook + per-run notebooks.

## Files
- `ResearchProject/systems/localised_jellium/scripts/localised_jellium_dynamics/proj_dyn/run.cpp` (frames)
- `.claude/skills/twin-run-analysis/twin_notebook_builder.py` (new cells)
- runs → `…/proj_dyn/results/`, `…/phase5_wp/results/`; notebooks → `hypotheses/twin_dynamics/`
