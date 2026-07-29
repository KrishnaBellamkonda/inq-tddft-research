# 01 — Wave-packet energy spread (low / medium / high)

## Hypothesis

LEED pattern features (peak intensity, fringe spacing, angular spread)
depend systematically on the projectile electron's kinetic energy.
Specifically:

- At low energy (E ≪ 200 eV): the pattern should be diffuse / lacking
  sharp Bragg-like features; the scattering regime is closer to
  one-electron-bound-state hybridisation than diffraction.
- At high energy (E ≫ 200 eV): higher-frequency fringes appear because
  the de Broglie wavelength shortens, and the diffraction-dominated
  regime sharpens.

## Runs collated

- `run_E30`   — E = 30 eV (low).
- `run_base`  — E = 200 eV (medium / paper reference).
- `run_E800`  — E = 800 eV (high).

All other parameters (b, σ, cell) are at the Tsubonoya base.

## Produce comparison artefacts

```
python scripts/coronene_postprocess.py hypothesis \
  --hypothesis-dir hypotheses/01_wp_energy_spread \
  --runs run_E30=$(realpath run_E30) \
         run_base=$(realpath run_base) \
         run_E800=$(realpath run_E800)
```

Outputs:
- `leed_total_grid.png` — 3 × N_screens grid of total LEED.
- `peak_intensity_vs_label.png` — peak intensity at the central screen.
- `energy_drift_overlay.png` — energy conservation comparison.
- `current_z_overlay.png` — J_z(t) curves overlaid.
