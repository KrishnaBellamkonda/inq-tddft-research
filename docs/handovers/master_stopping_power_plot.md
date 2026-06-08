# Handover: Master Stopping Power Plot

## Current status

**Complete — all 8 figures from the plan generated.**

Full pipeline implemented: computation module + 5 plotting scripts → 8 PNG figures at 600 DPI.

## What changed

Created the master stopping power S(E) plot pipeline from the plan at `docs/plans/master_stopping_power_plot.md`.

Six files:
- `inq-stack/python/inqview/report1/stopping_power_data.py` — computation module (new)
- `inq-stack/python/inqview/report1/fig_master_stopping.py` — master plot (rewritten from prototype)
- `inq-stack/python/inqview/report1/fig_matched_pair.py` — matched-pair figure (rewritten per plan)
- `inq-stack/python/inqview/report1/fig_definition_comparison.py` — definition bar chart (already existed, ran from pipeline)
- `inq-stack/python/inqview/report1/fig_sigma_rs.py` — σ/r_s scaling collapse (already existed, ran from pipeline)
- `inq-stack/python/inqview/report1/fig_sigma_sweep.py` — σ-dependence at E=100 (new)

## Files touched

| File | Action |
|---|---|
| `inq-stack/python/inqview/report1/stopping_power_data.py` | Created — reads CSVs, computes windows, extracts S₁/S₂/classical |
| `inq-stack/python/inqview/report1/fig_master_stopping.py` | Rewritten — produces four figures |
| `docs/reports/report1/figures/fig_master_stopping.png` | Generated — L=50, S₁ (momentum) |
| `docs/reports/report1/figures/fig_master_stopping_ks.png` | Generated — L=50, S₂ (KS orbital) |
| `docs/reports/report1/figures/fig_master_stopping_highdens.png` | Generated — L=30, S₁ |
| `docs/reports/report1/figures/fig_master_stopping_highdens_ks.png` | Generated — L=30, S₂ |
| `inq-stack/python/inqview/report1/fig_matched_pair.py` | Rewritten — 2-panel: S(E) + ΔE_H^bath(t) |
| `docs/reports/report1/figures/fig_matched_pair.png` | Generated — WP–classical comparison |
| `docs/reports/report1/figures/fig_definition_comparison.png` | Generated — L=50 S₁ vs S₂ bar chart |
| `docs/reports/report1/figures/fig_definition_comparison_highdens.png` | Generated — L=30 S₁ vs S₂ bar chart |
| `inq-stack/python/inqview/report1/fig_sigma_sweep.py` | Created — σ-dependence 2-panel |
| `docs/reports/report1/figures/fig_sigma_sweep.png` | Generated — |ΔE_KS| + Δσ_pz² vs σ |
| `docs/reports/report1/figures/fig_sigma_rs.png` | Generated — S_WP/S_cl vs σ/r_s scaling |

## Commands run

```bash
/local/data/public/skcb2/tddft/venv/bin/python3 -m inqview.report1.stopping_power_data   # standalone test
/local/data/public/skcb2/tddft/venv/bin/python3 -m inqview.report1.fig_master_stopping    # generate all 4 figures
```

## Tests and validation

### Data extraction verified

- L=50: 6 σ=1 pts, 6 σ=5 pts, 3 supplementary pts, 5 classical pts
- L=30: 4 σ=1 pts, 1 supplementary pt, 4 classical pts
- Full provenance printed for every data point (run dir, v1/v2, window, energies, S values)

### Bug found and fixed

- Classical E=50 CSV had truncated last row (NaN in velocity columns). Fixed by adding `dropna(subset=["vx","vy","vz"])` before regression.
- fig_matched_pair panel (b): raw ΔE_H was dominated by WP self-Hartree (-9.55 eV from spreading). Corrected by subtracting E_H_self = 1/(2σ_r√π) using actual σ_r(t) from `wp_real_space_stats.csv`. Bath-only response: +0.09 eV (WP) vs +0.87 eV (classical).

### Physics validation

1. **S₁ negative for σ ≤ 1**: Expected — free-particle spreading adds uncertainty KE faster than bath removes directed KE. Confirms S₂ is the reliable definition.
2. **σ=1 S₂ approaches classical**: At E=100, S₂=0.360 vs S_cl=0.325 (11% higher). Within expected range.
3. **σ=5 S ≪ classical**: ~10× lower at most energies — extended WP screens the Coulomb interaction.
4. **Classical S decreases with E**: Monotonic from E=20 (0.72) to E=600 (0.05), following Bethe-like falloff.
5. **L=30 (higher density) gives larger S**: S_cl(100eV, L=30) = 0.96 vs S_cl(100eV, L=50) = 0.33 — scales with electron density as expected.

## Known issues / blockers

1. **S₁ plots are sparse**: All σ=1 and σ=0.5 S₁ values are negative → filtered out. The L=30 S₁ plot shows only classical data. This is physically correct but makes S₁ less useful as a figure.
2. **Missing classical runs**: E=200 (L=50) not available. Plan lists 6 missing/needed classical runs.
3. **σ=8 point barely visible**: S₂=0.004 eV/Bohr is tiny on log scale, barely visible at bottom of KS plot.

## Assumptions still in play

- Interference-free window uses analytical free-particle spreading formula σ_r(t) = sqrt(σ²/2 + t²/(2σ²)), not the actual density from `wp_real_space_stats.csv`. The analytical approximation could underestimate spreading if the bath interaction accelerates it.
- Classical stopping uses the same time window as the paired σ=1 WP run. For energies with no σ=1 run, falls back to σ=5 window.
- The `E_expect_ha` column in `state_energies.csv` is the KS eigenvalue expectation value, not some other energy measure.

## Exact next steps

1. Review the generated figures and decide if the S₁ (momentum) figures are worth keeping in the report or should be dropped.
2. Dispatch missing classical runs (E=200 especially) to fill gaps in the classical line.
3. All dependent plots from the plan are now complete.
