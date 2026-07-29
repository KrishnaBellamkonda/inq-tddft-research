# 03_open_vs_closed_shell

Closed (N=38) vs open (N=40) shell jellium at otherwise identical WP.

## Runs in this comparison

- **run_01_base** — E=200.0 eV, σ=1.00 Bohr, direction=+z, N_e=38 (closed)
- **run_07_open_shell** — E=200.0 eV, σ=1.00 Bohr, direction=+z, N_e=40 (open)

## Artefacts

- `metadata.{md,csv}` — config snapshot for the runs above
- `observables_<col>.png` — overlay of current_z, dipole_z, energy_total vs time
- `residual_norm_bar.png` — WP overlap norm at t_final (capture indicator)
- `observables_*.png  (current_z, dipole_z, energy_total overlays)`