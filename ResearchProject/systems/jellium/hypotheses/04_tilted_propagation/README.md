# 04_tilted_propagation

Pure +z vs 45° xz-tilted WP propagation.

## Runs in this comparison

- **run_01_base** — E=200.0 eV, σ=1.00 Bohr, direction=+z, N_e=38 (closed)
- **run_04_tilted_45** — E=200.0 eV, σ=1.00 Bohr, direction=45° xz, N_e=38 (closed)

## Artefacts

- `metadata.{md,csv}` — config snapshot for the runs above
- `observables_<col>.png` — overlay of current_z, dipole_z, energy_total vs time
- `residual_norm_bar.png` — WP overlap norm at t_final (capture indicator)
- `observables_*.png  (current_z, dipole_z, energy_total overlays)`