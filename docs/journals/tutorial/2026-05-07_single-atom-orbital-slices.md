# 2026-05-07_single-atom-orbital-slices

**Title:** Single-atom orbital 2D slices: H, Li, Al
**Status:** complete

**Linked runs:**
- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/h/`   — H atom GS, LDA/60 Ry, L = 30 bohr cubic finite, atom at origin.
- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/li/`  — Li atom GS (all-electron pseudo: 1s, 2s valence).
- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/al/`  — Al atom GS ([Ne]-core pseudo: 3s, 3p valence).

**Linked spec / plan / handover:**
- `/local/data/public/skcb2/tddft/docs/superpowers/specs/2026-05-07-single-atom-orbitals-design.md`
- `/local/data/public/skcb2/tddft/docs/superpowers/plans/2026-05-07-single-atom-orbitals.md`
- `/local/data/public/skcb2/tddft/docs/handovers/single-atom-orbitals.md`

## Slice convention

Each panel is a 2D slice of the orbital density `|ψ_i(r)|²` through the cell centre at `(0, 0, 0)` bohr:

- **xy slice** at `z = 0`,
- **yz slice** at `x = 0`,
- **xz slice** at `y = 0`.

Cell extent in every panel is −15 .. +15 bohr (75 grid points, dx ≈ 0.4 bohr). Colour bar shared across the three panels of a single figure (one figure = one orbital). For non-negative densities the colour scale is `viridis` from 0 to the peak.

## Eigenvalues and per-figure peak amplitudes

Pure data (eigenvalue from SCF; peak from the slice generator). No interpretation in this section.

| Atom | Index | Suggested label | ε (Ha) | Occ | Peak |ψ|² (a.u.⁻³) |
|---|---|---|---|---|---|
| H  | 0 | 1s    | −0.2338 | 1.000 | 1.824 × 10⁻¹ |
| H  | 1 | 2s    | −0.0026 | 0.000 | 3.001 × 10⁻³ |
| H  | 2 | 2p_a  | +0.0151 | 0.000 | 4.790 × 10⁻⁴ |
| Li | 0 | 1s    | −1.8817 | 2.000 | 1.380 × 10⁰  |
| Li | 1 | 2s    | −0.1057 | 1.000 | 5.015 × 10⁻² |
| Li | 2 | 2p_a  | −0.0412 | 0.000 | 6.072 × 10⁻³ |
| Al | 0 | 3s    | −0.2853 | 2.000 | 1.108 × 10⁻² |
| Al | 1 | 3p_a  | −0.1020 | 0.333 | 1.363 × 10⁻² |
| Al | 2 | 3p_b  | −0.1020 | 0.333 | 1.298 × 10⁻² |
| Al | 3 | 3p_c  | −0.1020 | 0.333 | 1.469 × 10⁻² |
| Al | 4 | (4s candidate, ε=−0.0148 Ha) | −0.0148 | 0.000 | 7.306 × 10⁻⁴ |
| Al | 5 | (idx 5, ε=+0.0101 Ha)        | +0.0101 | 0.000 | 5.277 × 10⁻⁴ |
| Al | 6 | (idx 6, ε=+0.0130 Ha)        | +0.0130 | 0.000 | 2.897 × 10⁻⁴ |

For Al the labels at indices 4–6 are deliberately tentative — they sit close to ε = 0 and may mix with cubic-box modes. The point of this entry is to inspect the slices and decide.

Pseudopotential note for Al: the [Ne] core is frozen, so 1s/2s/2p shells of Al are *not* present in the calculation. The first computed valence state is 3s (index 0).

## Figures

### Hydrogen — first 3 orbitals

H (index 0) — 1s, ε = −0.234 Ha:

![H orbital 0 (1s)](attachments/2026-05-07_single-atom-orbital-slices/h_orbital_00_1s.png)

H (index 1) — 2s, ε = −0.003 Ha:

![H orbital 1 (2s)](attachments/2026-05-07_single-atom-orbital-slices/h_orbital_01_2s.png)

H (index 2) — 2p_a (one of the 2p triplet), ε = +0.015 Ha:

![H orbital 2 (2p_a)](attachments/2026-05-07_single-atom-orbital-slices/h_orbital_02_2p_a.png)

### Lithium — first 3 orbitals

Li (index 0) — 1s, ε = −1.882 Ha:

![Li orbital 0 (1s)](attachments/2026-05-07_single-atom-orbital-slices/li_orbital_00_1s.png)

Li (index 1) — 2s, ε = −0.106 Ha:

![Li orbital 1 (2s)](attachments/2026-05-07_single-atom-orbital-slices/li_orbital_01_2s.png)

Li (index 2) — 2p_a (one of the 2p triplet), ε = −0.041 Ha:

![Li orbital 2 (2p_a)](attachments/2026-05-07_single-atom-orbital-slices/li_orbital_02_2p_a.png)

### Aluminium — first 7 orbitals

Al (index 0) — 3s, ε = −0.285 Ha:

![Al orbital 0 (3s)](attachments/2026-05-07_single-atom-orbital-slices/al_orbital_00_3s.png)

Al (index 1) — 3p_a, ε = −0.102 Ha:

![Al orbital 1 (3p_a)](attachments/2026-05-07_single-atom-orbital-slices/al_orbital_01_3p_a.png)

Al (index 2) — 3p_b, ε = −0.102 Ha:

![Al orbital 2 (3p_b)](attachments/2026-05-07_single-atom-orbital-slices/al_orbital_02_3p_b.png)

Al (index 3) — 3p_c, ε = −0.102 Ha:

![Al orbital 3 (3p_c)](attachments/2026-05-07_single-atom-orbital-slices/al_orbital_03_3p_c.png)

Al (index 4) — 4s candidate, ε = −0.015 Ha:

![Al orbital 4](attachments/2026-05-07_single-atom-orbital-slices/al_orbital_04_n4_idx_4.png)

Al (index 5) — ε = +0.010 Ha:

![Al orbital 5](attachments/2026-05-07_single-atom-orbital-slices/al_orbital_05_n4_idx_5.png)

Al (index 6) — ε = +0.013 Ha:

![Al orbital 6](attachments/2026-05-07_single-atom-orbital-slices/al_orbital_06_n4_idx_6.png)

## Observations

By examining this, I want to understand the amplitude and the spread of these orbitals in space.
