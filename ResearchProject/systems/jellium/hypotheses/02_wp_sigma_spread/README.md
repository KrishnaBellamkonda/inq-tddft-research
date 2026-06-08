# 02_wp_sigma_spread

WP-width sweep at fixed E=200 eV (σ=0.5 / 1.0 / 3.78 Bohr). Includes Gaussian-broadening overlay (narrow-σ shows fastest dispersion).

## Runs in this comparison

- **run_06_narrow_sigma** — E=200.0 eV, σ=0.50 Bohr, direction=+z, N_e=38 (closed)
- **run_01_base** — E=200.0 eV, σ=1.00 Bohr, direction=+z, N_e=38 (closed)
- **run_05_wide_sigma** — E=200.0 eV, σ=3.78 Bohr, direction=+z, N_e=38 (closed)

## Artefacts

- `metadata.{md,csv}` — config snapshot for the runs above
- `observables_<col>.png` — overlay of current_z, dipole_z, energy_total vs time
- `residual_norm_bar.png` — WP overlap norm at t_final (capture indicator)
- `observables_*.png  (current_z, dipole_z, energy_total overlays)`
- `gaussian_broadening_overlay.png  (measured σ_z(t) vs free-particle analytic)`