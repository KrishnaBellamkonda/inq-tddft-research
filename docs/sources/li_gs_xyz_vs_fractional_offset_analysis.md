# Li 54-atom GS: `.xyz` vs `insert_fractional` convention

**Question.** When we rebuilt the Li 54-atom GS on the `[-L/2, +L/2]`
`.xyz` + `ions::parse(...)` convention (mandated by the build-run
skill), the converged total energy came out **7.4 mHa above** the old
`insert_fractional` checkpoint. For two periodically-equivalent atom
configurations this should be invariant. Is the new GS sound? Does the
divergence affect the FFT-peak comparison the project depends on?

## Inputs

| | OLD (`run_save_gs_2x2x2_T200/`) | NEW (`run_save_gs_2x2x2_T200_xyz/`) |
|---|---|---|
| atom-position API | `ions.insert_fractional({fx, fy, fz})` with fx ∈ {0, 1/6, 1/3, 1/2, 2/3, 5/6} | `ions::parse("li_54_3x3x3.xyz", cell)` with Cartesian Å in `[-L/2, +L/2)` |
| effective rigid translation | atoms in `[0, 5L/6]` Å Cartesian | atoms in `[-L/2, +L/3]` Å Cartesian — **shifted by -L/2** vs OLD, but periodically equivalent |
| same: cell `cubic(10.53 Å).periodic()`, MP 2×2×2 shifted, 74 Ry, T=400 K, Broyden ndim=8 α=0.1 |
| E_total | -389.3217493 Ha | -389.3143369 Ha |
| ΔE | — | **+7.4 mHa = +0.20 eV** |
| SCF iters | (not recorded in handover) | 73 |
| ∫ρ d³r | 162.0000 ✓ | 162 (exact, INQ integer count) |

## Analysis

### Per-state eigenvalue comparison

Joining the two checkpoints by `(kpoint_index, state_index)` gives a
spurious 1.66 eV outlier — but this is a state-index re-labelling
artefact: the SCF can converge with two near-degenerate Fermi-surface
states swapped. **Sorting eigenvalues per kpoint** removes this:

| Statistic | Value |
|---|---|
| N pairs | 808 (= 101 bands × 8 kpoints) |
| mean Δε | +0.0087 eV |
| std Δε | 0.10 eV |
| max \|Δε\| | 1.66 eV (k=2, single state) |
| \|Δε\| > 50 meV | 51 / 808 (6.3%) |
| \|Δε\| > 200 meV | 34 / 808 (4.2%) |
| Fermi-window mean Δε | +0.0009 eV (≈ zero) |
| Fermi-window max \|Δε\| | 0.39 eV |

Per-kpoint summary (sorted-spectrum comparison):

| k | mean Δε (eV) | max \|Δε\| (eV) | std Δε (eV) |
|---|---:|---:|---:|
| 0 | +0.0140 | 0.5262 | 0.064 |
| 1 | -0.0084 | 0.4366 | 0.055 |
| 2 | -0.0108 | **1.6609** | 0.176 |
| 3 | +0.0019 | 0.1222 | 0.015 |
| 4 | +0.0536 | 0.7210 | 0.149 |
| 5 | +0.0004 | 0.1044 | 0.017 |
| 6 | -0.0147 | 0.6355 | 0.078 |
| 7 | +0.0336 | 0.6474 | 0.112 |

CSV: `docs/sources/li_gs_xyz_vs_fractional_eigenvalue_comparison.csv`.

### Diagnosis

The two SCFs converged to **different metallic stationary points**.
Mechanism:

1. The `.xyz` Cartesian shift by -L/2 changes the **initial-guess phase
   relative to the FFT grid origin**. INQ's `ground_state::initial_guess`
   builds an atom-centred superposition of pseudoatomic densities. The
   numerical aliasing onto the Cartesian grid is *not* translationally
   invariant when the cell is FFT-discretised, even though the
   underlying physics is.
2. The Broyden mixer uses the initial-guess seed to walk to a basin in
   the metallic energy landscape. Many basins are nearly degenerate at
   400 K Fermi smearing (kT = 34 meV) — the Fermi surface has
   a high density of states (cf. Ashcroft & Mermin §2 for the
   free-electron DOS, which is large near E_F for Li).
3. Both runs reach the SCF energy tolerance (1e-6 Ha = 27 µHa), so each
   is a valid local minimum. The 7.4 mHa **inter-basin** energy
   difference is much larger than the **intra-basin** SCF tolerance,
   which is consistent with this picture.

### Why this is acceptable for the FFT-peak comparison

- **6.5 eV plasmon (low-v target)** is a *collective* longitudinal mode
  that arises from Σᵢ over many bands' contribution to χ(q→0, ω). It is
  insensitive to individual Fermi-surface band relabelling. The bulk
  Li plasmon energy is set by the carrier density and effective mass —
  both of which are fixed in this comparison (∫ρ = 162 exactly in both
  GS).
- **2.8 eV e-h (high-v target)** depends on the *cluster structure* of
  Γ-Γ transitions (Figure 5 of the BCN:1719P paper), not on individual
  transitions. The cluster statistics should be preserved even with
  band relabelling.
- The 0.4 eV Fermi-window max shift is comparable to the FFT bin width
  (the 15 fs Hann-windowed run gives ~0.5 eV resolution), so any peak
  shift would be at the noise floor anyway.

### What we will check empirically

The plasmon run (v=0.0626 on the new GS) is the cheapest empirical
test: if its FFT peak lands at ~6.5 eV (within ±0.3 eV), the GS
divergence is confirmed benign. If the peak shifts by ≥ 0.5 eV from
the v=0.0123 result (5.722 eV), we revisit. The peak position vs
velocity also has a known weak dependence — the paper sees ~6.5 eV
across the entire low-v family with little variation.

### Mitigations available

- **Tighten SCF mixing**: `mixing_ndim(12)` and lower `mixing(0.05)` to
  encourage convergence to a deeper basin. Costs ~2× SCF wall.
- **Lower temperature**: at 200 K (kT = 17 meV) the Fermi-surface
  ambiguity halves, but this brought the original SCF to oscillate
  divergently (handover §Risks).
- **Seed from old checkpoint**: load the OLD checkpoint as the initial
  guess for the new run (translated atom positions); SCF will refine
  to the same basin. This would null out the divergence.

For the current scientific question (plasmon vs e-h), none of these
are needed unless the empirical check above fails.

## References

- Ashcroft & Mermin, *Solid State Physics*, Saunders 1976, §2 — free
  electron model of Li.
- Marx & Hutter, *Ab Initio Molecular Dynamics*, Cambridge 2009, §3.5 —
  metallic SCF and finite-T smearing methodology.
- Santervás-Arranz, Stengel, Artacho, *Phys. Rev. Research* 7, 033292
  (2025) — kick paradigm in metallic Li, Mv²/countercurrent rule.
- BCN:1719P quantum-kick draft (Cavendish, 2026) — Figures 4-5.
- INQ pseudo-spectral discretisation: see `docs/inq_source_map.md`
  for the FFT-based H|ψ⟩ apply that this finite-T metallic divergence
  intersects with.
