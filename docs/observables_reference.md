# Observables Reference — Jellium WP-RT Runs

Standard observable suite for all jellium wave-packet real-time TDDFT runs
(`jellium-wp-rt/run_01_base` through `run_07_open_shell`).

---

## Primary observables (written during simulation)

| Observable | C++ source | Output path | Frequency |
|---|---|---|---|
| GS total density | `density::total(electrons)` before WP injection | `results/density_gs/` | once |
| GS orbital densities | `density::orbital(electrons, i)` for i in [0, wp_idx) | `results/density_gs_orbitals/orbital_XXXX/` | once each |
| Total electronic density(t) | `density::total() + density::orbital(wp_idx)` | `results/density_rt_total/` | every `WRITE_EVERY` steps |
| Jellium-only density(t) | `density::total(electrons)` | `results/density_rt_jellium/` | every `WRITE_EVERY` steps |
| WP orbital density(t) | `density::orbital(electrons, wp_idx)` | `results/density_rt_wp/` | every `WRITE_EVERY` steps |
| KS overlap matrix O_ij(t) | `OrbitalOverlapMatrix::snapshot()` | `results/overlap/overlap_XXXXXX.csv` | every step |
| Overlap index | — | `results/overlap/index.csv` | appended each snapshot |
| Total energy, KE, Hartree, XC(t) | `ObservablesWriter` | `results/observables.csv` | every step |
| Current J_x, J_y, J_z(t) | `ObservablesWriter` | `results/observables.csv` | every step |
| Dipole μ_x, μ_y, μ_z(t) | `ObservablesWriter` | `results/observables.csv` | every step |
| Time-averaged LEED screens | `LeedPatternAccumulator::save()` | `results/screens/screen_NN.dat` | end of run (20 screens) |
| Instantaneous LEED snapshots | `PlaneScreen::extract()` + `save()` | `results/screens_snapshots/step_XXXXXX/screen_NN.dat` | every 3 steps |

### Notes

- `density::total()` returns the sum of all occupied KS orbital densities **excluding** the WP
  extra state. The WP is injected into an extra state beyond the occupied manifold.
- Total observable density = `density::total() + density::orbital(wp_idx)` at every frame.
- `WRITE_EVERY = 2` for runs 01–02 (slow WP); `WRITE_EVERY = 10` for runs 03–07.

---

## Overlap matrix O_ij(t)

**Definition:** `O_ij(t) = |dV × Σ_r conj(ψ_i^GS(r)) × ψ_j(r,t)|²`

- i ∈ [0, n_ref): index over GS reference orbitals (n_ref = wp_idx = number of occupied states)
- j ∈ [0, n_evolved): index over evolved orbitals (n_evolved = wp_idx + 1, includes WP)
- At t=0: O is approximately identity on the occupied block; O[i, wp_idx] ≈ 0

**CSV format per snapshot:**
```
row 0:  O[0,0]  O[0,1]  ...  O[0, n_evolved-1]
row 1:  O[1,0]  ...
...
row n_ref-1: O[n_ref-1, 0]  ...
```

**index.csv columns:** `step, time_au, file`

---

## Post-processing derivations (analysis.py)

| Derived output | Source data | Method |
|---|---|---|
| Observables summary PNG | `observables.csv` | `plot_observables_summary()` |
| Energy spectrum | `observables.csv`, `energy_total` column | `FourierTransform().transform_energy()` |
| Current spectrum (x,y,z) | `observables.csv`, `current_*` columns | `FourierTransform().transform_column()` |
| VTI series (3 × 50 frames) | `density_rt_*` directories | `convert_real_series_to_vti()` |
| Density slice GIFs (6) | `density_rt_*` directories | `plot_density_slice()` + imageio |
| LEED 4×5 grid PNG | `results/screens/*.dat` | `load_leed_pattern()` + matplotlib |
| Individual LEED PNGs (20) | `results/screens/*.dat` | `plot_leed_pattern()` |
| LEED time-evolution GIF | `screens_snapshots/*/screen_10.dat` | `load_leed_pattern()` + imageio |
| Overlap matrix bar-chart GIFs | `overlap/index.csv` + CSVs | matplotlib bar chart per evolved orbital j |
| ParaView 3D renders | `density_rt_total` VTIs | `ParaViewPipeline.render_density_from_meta_series()` |

---

## Screen positions (z, bohr)

20 screens, two at the cell boundaries plus 18 interior screens with small deterministic
offsets from evenly-spaced positions:

```
0.5, 2.53, 4.66, 6.78, 8.87, 10.95, 12.97, 15.03, 17.06, 19.09,
21.07, 23.11, 25.08, 27.12, 29.04, 31.03, 33.01, 34.97, 36.95, 39.5
```

Screen 10 (z ≈ 21.07 bohr, near the cell midpoint) is used for the time-evolution GIF.
