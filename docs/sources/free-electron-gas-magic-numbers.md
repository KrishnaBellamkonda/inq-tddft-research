# Free-electron gas: closed-shell magic numbers (Γ-only cubic)

## Why this note exists

The jellium TDDFT runs at `ResearchProject/systems/jellium/` use a Γ-only,
cubic, periodic LDA calculation with N electrons in volume V = L³. For such a
system the Kohn-Sham eigenstates of a perfectly uniform background are pure
plane waves ψ_G(r) = (1/√V) exp(i G·r) with G = (2π/L)(n_x, n_y, n_z),
n_i ∈ ℤ. Each spatial orbital holds two electrons (paramagnetic, doubly
occupied).

The total density is exactly uniform iff every degenerate shell is *fully*
occupied. With a partially-filled shell, the SCF can return any orthonormal
basis spanning the partial subspace, and the density is generically
non-uniform — see Ashcroft & Mermin Ch. 2 and the Phase-1 sign-off in
`docs/handovers/jellium_l60_observables.md`.

When choosing N for a jellium reference run, the practical rule is therefore:

1. Pick a target density n.
2. Compute N* = n · L³.
3. Snap N to the closed-shell magic number nearest N*.
4. Verify the corresponding r_s is close enough to the physical target.

## Cumulative magic-N table (verified 2026-05-03)

| |G|² | shell deg. (spatial) | cum. spatial states | cum. electrons |
|---:|---:|---:|---:|
|  0 |  1 |   1 |   2 |
|  1 |  6 |   7 |  14 |
|  2 | 12 |  19 |  38 |
|  3 |  8 |  27 |  54 |
|  4 |  6 |  33 |  66 |
|  5 | 24 |  57 | 114 |
|  6 | 24 |  81 | 162 |
|  8 | 12 |  93 | 186 |
|  9 | 30 | 123 | 246 |
| 10 | 24 | 147 | 294 |
| 11 | 24 | 171 | 342 |
| 12 |  8 | 179 | 358 |
| 13 | 24 | 203 | 406 |
| 14 | 48 | 251 | 502 |
| 16 |  6 | 257 | **514** |
| 17 | 48 | 305 | 610 |
| 18 | 36 | 341 | 682 |
| 19 | 24 | 365 | 730 |
| 20 | 24 | 389 | 778 |
| 21 | 48 | 437 | 874 |
| 22 | 24 | 461 | 922 |
| 24 | 24 | 485 | 970 |
| 25 | 30 | 515 | 1030 |
| 26 | 72 | 587 | 1174 |
| 27 | 32 | 619 | 1238 |
| 29 | 72 | 691 | 1382 |
| 30 | 48 | 739 | 1478 |

(no |G|² = 7, 15, 23, 28, … — Legendre's three-square theorem skips integers
of the form 4^a (8b + 7).)

## Choice for the L=60 / 4×-density rerun (2026-05-03)

- L = 60 bohr, target density = 4 × 5.926e-4 = **2.370e-3 bohr⁻³**
- Target N* = 2.370e-3 · 216000 = 511.9
- Closest magic-N: **N = 514** (|G|² ≤ 16)
- Achieved n = 514 / 216000 = 2.380e-3 bohr⁻³ (4.015× current — 0.4 % above target)
- r_s = (3 / (4π · n))^(1/3) = (100.34)^(1/3) ≈ **4.64 bohr** — between sodium
  (r_s ≈ 3.93) and potassium (r_s ≈ 4.86). `Source: Ashcroft & Mermin Table 1.1.`

## Enumeration script (reproducibility)

```python
# Cumulative closed-shell magic numbers, 3D Γ-only cubic.
# Each (n_x, n_y, n_z) lattice point holds 2 electrons.
from collections import Counter
nmax = 30
shells = Counter()
for nx in range(-nmax, nmax + 1):
    for ny in range(-nmax, nmax + 1):
        for nz in range(-nmax, nmax + 1):
            shells[nx*nx + ny*ny + nz*nz] += 1
cum_states = 0
print(f"{'|G|^2':>6} {'shell_deg':>11} {'cum_states':>11} {'cum_e':>8}")
for k in sorted(shells):
    if k > 30: break
    cum_states += shells[k]
    print(f"{k:>6d} {shells[k]:>11d} {cum_states:>11d} {2*cum_states:>8d}")
```

## Cross-check vs Sec. A of `docs/plans/in-this-task-we-lively-meerkat.md`

The plan's Phase-1.B table extension above |G|²=16 was wrong (e.g. it listed
N=518 at |G|²≤22 as the closest magic to 512, but the verified count is
N=922 at |G|²≤22). The verified closest magic is **N=514 at |G|²≤16**, the
same as the plan's original Section B candidate.
