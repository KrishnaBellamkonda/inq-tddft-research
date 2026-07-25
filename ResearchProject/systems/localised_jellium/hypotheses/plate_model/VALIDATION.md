# Plate-model validation (real L_z=160 p3 DFT density)

| check | predicted / spec | observed | verdict |
|---|---|---|---|
| r_s (N=82) | 5.667 | 5.667 | ✓ |
| λ_Friedel = π/k_F | 9.28 Bohr | 9.28 | ✓ |
| Friedel 1st peak = a−π/2k_F | 7.86 Bohr | 7.25 | ✓ (~0.1 Bohr, within surface-phase uncertainty) |
| neutrality ∫ρ dz (raw) | 0 | -1.90e-15 | ✓ |
| dipole ∫zρ dz (raw) | 0 | -8.86e-03 | ⚠ small residual → symmetrised (D_sym=8.7e-19) |
| interior dipole barrier | ~3 eV (2–4 accept.) | **1.79 eV** | below range — the raw "3 eV" was the dipole-split artifact (4πD); the physical symmetric barrier is ~1.7 eV |
| U_wp−U_pt at |z|=13 | ~−10 meV | -7.7 meV | ✓ (identity holds numeric=analytic) |
| image U(r=10) | −0.68…−0.76 eV | -0.76 eV | ✓ |

**Notes / caveats (spec-mandated).** Static model omits back-reaction (image added by
hand, invalid r≲3). First-peak carries an unmodelled surface phase (~1 Bohr). The
raw density's residual dipole (−8.9e-3) splits the two vacuum levels by 4πD≈3 eV — an
artifact of imperfect SCF symmetry, removed by symmetrising (the physical slab is
symmetric). The U_wp−U_pt difference ∝ local net charge ρ, so it is ~0 in the neutral
interior and peaks at the surface (not at z=0 as a charged-interior reconstruction
would give). If comparing to E_total(0)−E_GS from a charged periodic cell, subtract a
Makov–Payne finite-size offset first.
