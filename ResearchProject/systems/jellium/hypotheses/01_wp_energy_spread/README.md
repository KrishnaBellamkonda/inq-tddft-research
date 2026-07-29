# 01_wp_energy_spread

WP kinetic-energy sweep at fixed σ=1 Bohr (50 / 200 / 400 eV). Includes Gaussian-broadening overlay.

## Runs in this comparison

- **run_02_low_energy** — E=50.0 eV, σ=1.00 Bohr, direction=+z, N_e=38 (closed)
- **run_01_base** — E=200.0 eV, σ=1.00 Bohr, direction=+z, N_e=38 (closed)
- **run_03_high_energy** — E=400.0 eV, σ=1.00 Bohr, direction=+z, N_e=38 (closed)

## Artefacts

- `metadata.{md,csv}` — config snapshot for the runs above
- `observables_<col>.png` — overlay of current_z, dipole_z, energy_total vs time
- `residual_norm_bar.png` — WP overlap norm at t_final (capture indicator)
- `observables_*.png  (current_z, dipole_z, energy_total overlays)`
- `gaussian_broadening_overlay.png  (measured σ_z(t) vs free-particle analytic)`