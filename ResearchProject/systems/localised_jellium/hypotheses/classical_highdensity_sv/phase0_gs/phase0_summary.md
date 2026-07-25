# Phase 0 ground state — classical_highdensity_sv (slab_n100)

Denser localised jellium slab: 35x35x85 Bohr box, 25-Bohr slab (half-width 12.5),
N=100, dx=0.5, periodicity(2) (z-open), LDA, T=100 K.

| Quantity | Value |
|---|---|
| GS energy | 207.183 Ha |
| r_s | 4.1815 (target 4.18) |
| num_states | 74 (≈50 occupied + 24 extra) |
| ∫ n dV | 100.000 (target N = 100) |
| interior mean density | 3.2831e-03 a0^-3 |
| n0 (target) | 3.2653e-03 a0^-3 |
| interior / n0 | 1.005 |
| z-symmetry residual (alignment-corrected) | 4.32e-05 (reflection centre c=+0.243 Bohr = +dz/2) |
| z-symmetry residual (naive about z=0) | 1.15e-01 — spurious, from half-cell grid offset at the steep slab face |
| spill-out decay length λ | 1.06 Bohr |
| grid | 70x70x175, dx=0.50 Bohr |

**Verdict:** Yes — bulk interior sits at n0 (interior/n0 = 1.005), the erfc faces are symmetric to 4e-05 once the half-cell VTI grid offset is accounted for, and n(z) spills out exponentially into the vacuum; occupations are metallic (smeared top shell, no pathology).

Dashboard: `phase0_gs_dashboard.png`.
Source VTI: `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/scripts/classical_highdensity_sv/gs/results/density_gs/density_gs.vti` (loaded via `inqview.load_vti`, physical order, centered-z check passed).
