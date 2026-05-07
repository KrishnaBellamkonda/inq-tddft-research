# Monkhorst–Pack k-point sampling (PRB 13, 5188, 1976)

**Citation.** Monkhorst, H. J.; Pack, J. D. *Special points for Brillouin-zone
integrations.* Phys. Rev. B 13, 5188 (1976). [DOI:10.1103/PhysRevB.13.5188]

## What it is

A prescription for choosing a uniform mesh of k-points in the Brillouin zone (BZ)
that, after symmetry reduction, gives the most efficient sampling for integrating
periodic functions. For a `n1 × n2 × n3` mesh the points are

    k_{i,j,k} = (u_i, u_j, u_k),
    u_p = (2 p − n_p − 1) / (2 n_p),  p = 1 … n_p

with an optional half-cell shift to remove the Γ-point and centre the mesh on the
BZ. INQ exposes both unshifted and shifted forms via
`input::kpoints::grid({n1,n2,n3}, shifted=true)` (`inq/src/input/kpoints.hpp:21`).

## Why it matters here

The legacy QBall Li sweep
(`QuantumKickExtension/qball-codebase/Li/td_kicks/Li.54_td.inp:22-23`) uses
**Γ-only** sampling on a 3×3×3 BCC supercell. For a metallic system at the
edge of linear response, Γ-only on a supercell is equivalent to a `nx × ny × nz`
MP grid on the primitive cell only if the supercell sides are commensurate
multiples of the primitive cell sides — which they are here (3:1 in each
direction). In that limit Γ-only on a 3×3×3 supercell ≡ a 3×3×3 MP grid on the
primitive cell. The new INQ runs use a **2×2×2 shifted MP** on the supercell,
which is equivalent to a **6×6×6 shifted MP** on the primitive cell — denser
sampling that exposes Fermi-surface structure missed by QBall.

## Recommended use

For metals at low Fermi smearing, MP grids should be tested for convergence by
doubling/halving the density and checking that total energy per atom converges
to within a few meV. For Li in a 54-atom 3×3×3 supercell, 2×2×2 shifted MP
(equivalent to 6×6×6 primitive) is the standard starting density; 4×4×4 (≡
12×12×12 primitive) is the converged-reference choice if cost permits.

## Inferences (label as own)

- **Inference:** Γ-only QBall plasmon-like FFT peaks may be artefacts of Fermi-
  surface mis-sampling. The INQ 2×2×2 shifted comparison will test this.
