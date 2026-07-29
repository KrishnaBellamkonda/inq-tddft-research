# 05_electron_capture

Electron-capture diagnostic: slow WP (50 eV, σ=1) most likely capture candidate vs 200 eV baseline.

## Runs in this comparison

- **run_01_base** — E=200.0 eV, σ=1.00 Bohr, direction=+z, N_e=38 (closed)
- **run_02_low_energy** — E=50.0 eV, σ=1.00 Bohr, direction=+z, N_e=38 (closed)

## Artefacts

- `metadata.{md,csv}` — config snapshot for the runs above
- `observables_<col>.png` — overlay of current_z, dipole_z, energy_total vs time
- `residual_norm_bar.png` — WP overlap norm at t_final (capture indicator)
- `observables_*.png  (current_z, dipole_z, energy_total overlays)`
- `capture_diagnostic.png  (3-panel f_trap, J_z, GS overlap)`
- `capture_summary.txt   (numerical readout)`