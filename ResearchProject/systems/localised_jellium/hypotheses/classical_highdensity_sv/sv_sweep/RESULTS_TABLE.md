# Results — high-density classical S(v) benchmark (r_s=4.18, CAP-free)

Campaign `classical-highdensity-sv`. Mass-1 Gaussian-charge electron (sigma_WP=0.5),
Ehrenfest, 25-Bohr slab, z-open periodicity(2), no CAP. S = E_absorbed/L_slab
(slab excitation vs the projectile-free GS), cross-checked by projectile KE loss.

## Stopping power
| v_launch | v_final | v_mean | S (eV/Bohr) | S_keloss | E_abs (eV) | plateau flat (eV) |
| --- | --- | --- | --- | --- | --- | --- |
| 2.0 | 1.42 | 1.82 | 1.087 | 1.086 | 27.2 | 2.0e-05 |
| 2.5 | 2.11 | 2.41 | 0.97 | 0.97 | 24.3 | 2.0e-04 |
| 3.0 | 2.77 | 2.95 | 0.709 | 0.708 | 17.7 | 5.5e-04 |
| 3.5 | 3.36 | 3.47 | 0.509 | 0.509 | 12.7 | 7.8e-04 |
| 4.0 | 3.91 | 3.98 | 0.374 | 0.374 | 9.3 | 4.9e-04 |
| 4.5 | 4.44 | 4.49 | 0.283 | 0.283 | 7.1 | 3.0e-04 |

- Deposit vs KE-loss agree to ~0.1%; plateau flatness ~1e-4 eV (energy conserved).
- Bethe-tail power law: **S ∝ v^-1.72**.
- Electron decelerates (v_final < v_launch): S is at the mean in-slab velocity (v_mean).

## Component-ledger deltas (eV, t=0 → plateau) — Definition-1 staging
| v | ΔE_PP | ΔE_PS | ΔE_SS | ΔE_SB | ΔE_PB |
| --- | --- | --- | --- | --- | --- |
| 2.0 | -20.27 | 335.02 | -183.81 | 187.55 | -335.02 |
| 2.5 | -20.27 | 335.02 | -150.93 | 156.54 | -335.02 |
| 3.0 | -20.27 | 335.02 | -40.47 | 42.44 | -335.02 |
| 3.5 | -20.27 | 335.02 | -47.27 | 52.41 | -335.02 |
| 4.0 | -20.27 | 335.02 | -57.08 | 59.04 | -335.02 |
| 4.5 | -20.27 | 335.02 | -36.0 | 36.88 | -335.02 |

P=projectile, S=slab electrons, B=background; E_BB (background self) constant, omitted.
