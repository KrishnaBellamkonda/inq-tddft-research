# 2026-05-04_orbitals_at_each_band_at_each_k_point

**Title:** Orbitals at each band at each k point
**Run path:** `/local/data/public/skcb2/tddft/Tutorial/_inqkit_tests/orbital_per_kpoint_S2_li_2x2x2`
**Linked results:** `/local/data/public/skcb2/tddft/Tutorial/_inqkit_tests/orbital_per_kpoint_S2_li_2x2x2/results`
**Status:** complete

## Config snapshot

| Field | Value |
|---|---|
| run | orbital_per_kpoint_S2_li_2x2x2 |
| system | li_16_atom_bcc_supercell_2x2x2 |
| cell_angstrom | 7.02 7.02 7.02 |
| boundary | periodic |
| n_atoms | 16 |
| k_grid | 2 2 2 shifted |
| smearing | fermi_dirac |
| smearing_temperature_kelvin | 400 |
| xc | pbe |
| cutoff_ha | 30 |
| extra_states | 8 |
| scf_iters | 118 |
| ground_state_energy_ha | -115.4607060895939 |
| num_states | 32 |
| num_electrons | 48 |
| bands_dumped | 1,12,24,30 |
| smoke_pass | yes |

## Observations

We use the Li 2×2×2 supercell simulation as a case study for understanding
how Bloch orbitals look across the Brillouin zone. The configuration above
has been read off the `run_summary.txt`; the key numbers we will refer to
repeatedly below are **8 k-points** (2×2×2 shifted Monkhorst–Pack mesh),
**32 bands** per k-point (24 doubly-occupied + 8 extra states), and a
**semi-core ONCV PBE pseudopotential**. The pseudopotential keeps the 1s²
core inside the simulation rather than freezing it; effectively, only the
nuclear charge has been screened, and each Li atom contributes its full
**3 electrons** to the calculation. In a 2×2×2 supercell of BCC Li we have
2 atoms per primitive cell × 8 cells = 16 atoms, so the total electron
count is

> 16 atoms × 3 electrons/atom = **48 electrons**

(Note: an earlier draft said "24 electrons / 12 orbitals". That was a
typo — corrected to 48 electrons / 24 doubly-occupied orbitals here, with
the rest, `extra_states = 8`, sitting empty above the Fermi level.)

The eigenvalues for all 32 bands at all 8 k-points are tabulated in
`results/raw/observables/eigenvalues/eigenvalues.csv` (256 rows). The
band structure (E vs k-point index) plotted across the full eigenvalue
range and zoomed into the valence/conduction window is shown below:

![E vs k-point index for all 32 bands; right panel zooms on bands 16-31](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/e_vs_k_all_bands.png)

Bands 0–15 form the deep core manifold around −47.3 eV (16 narrow Bloch
combinations of the atomic 1s orbitals on the 16 atoms). Bands 16–22 are
the **fully-occupied lower valence** at −1.5…−2.7 eV. Band **23** is the
**partially-filled valence top** (occupation ≈ 0.13/0.125 ≈ half a doubly-
filled state). Band 24 onward is the conduction manifold, mostly empty
with a small Fermi-Dirac tail at 400 K.

### K points

K points represent a possibility that the system can stay in.
Specifically, each k point specifies the **phase twist** the Bloch
wavefunctions for each orbital would get. Consider the Γ point. The k
value here is 0. Hence, the phase twist at the boundaries of cells would
be none at this k point. However, if the k point considered is something
along the lines of (0, 0, 1), then there would be a phase twist of a
certain magnitude along the z axis, going from one cell to the other.
This phase twist is important as it might characterise important
properties in topological materials. Even in materials such as Li, we
find that the orbitals of the same band but different k point have
different density distributions.

(Note: the 2×2×2 shifted Monkhorst–Pack mesh used here actually has **no
exact Γ point** — the eight k-points sit at the corners ±(0.25, 0.25, 0.25)
of the BZ in 2π/a units. The argument above translates to "the k closest
to Γ has the smallest phase twist".)

The tight binding explanation will be given below. I wonder what specific
observables are affected by the k points, if the densities are not
affected.

### Orbitals in Li

Use python post-processing to make xz / yz / xy slices of the cell to
show the sample densities of different orbitals in different bands. Four
band-summary plots follow, each showing the eight k-points side by side
on a fixed colour scale; for every band we show three components of the
Bloch orbital ψ_{n,k}(r) — Re ψ, Im ψ, and the density |ψ|².

#### Band 1 — deep core 1s

![band 1 Re psi grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_001_re_psi_grid.png)
![band 1 Im psi grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_001_im_psi_grid.png)
![band 1 density grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_001_density_grid.png)

Consider the first band at index 001. Here, it can be seen that,
regardless of what the k point is, **all the orbitals look exactly the
same** in their density. Now, this might be surprising at first. However,
consider this fact. The first orbital is essentially a Gaussian — a 1s
orbital tightly bound to the Li nucleus. Here, we find that the density
does not reach the boundary, and the phase twist is of no use.

This can also be understood in a different manner. Consider a
tight-binding model for this system. We can say that the energy of this
band is not k-dependent (or is not substantially k-dependent), meaning
e(k) would have an expression

> e(k) = e₀ + 2 t · (overlap term)

where the overlap term encodes the hopping between neighbouring atomic
orbitals. Usually the overlap term would be k-point dependent. But for
this band, it is not. That means the overlap between the neighbouring
orbitals is almost zero — i.e., one Li atom's 1s orbital is so tightly
bound to its own nucleus that it barely overlaps with its neighbours'.
So we can confidently say that the first band is a **core 1s orbital**.

#### Band 12 — also 1s manifold

![band 12 Re psi grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_012_re_psi_grid.png)
![band 12 Im psi grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_012_im_psi_grid.png)
![band 12 density grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_012_density_grid.png)

The next band visualised in the results, band 12, also has similar energy
to that of the previous orbitals (E ≈ −47.27 eV vs −47.32 eV for band 1).
These must be **1s orbitals too** — bands 0–15 together form the
16-state core 1s manifold of the 16-atom supercell, with very small
splitting set by inter-atomic 1s–1s overlap.

#### Band 24 — bottom of the conduction band

![band 24 Re psi grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_024_re_psi_grid.png)
![band 24 Im psi grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_024_im_psi_grid.png)
![band 24 density grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_024_density_grid.png)

The next band would represent the **valence / conduction frontier**.
(Note (correction): the originally-anticipated "partially-filled valence
2s" is band **23**, which is what carries the Fermi-tail occupation
≈ 0.5 of a doubly-filled state. Band **24**, which we actually dumped, is
the *first empty* state — the bottom of the conduction band — with a
small Fermi-Dirac tail of ≈ 0.04 / 0.125 ≈ 0.34 of a state at 400 K.
So band 24 is essentially conduction; the tight-binding picture below
applies to both bands 23 and 24 because they are both built from the
same 2s atomic orbitals, half-filled in total.)

We must tabulate the occupation numbers too. Some numbers in the CSV
will be below 1 in absolute terms. **This is my new learning.** The
occupation in each orbital depends on the number of electrons in the
orbital — 1 or 2 — times 1/(number of k-points). The occupation is
effectively the overall weight given to a specific orbital. This means
that, in the calculation of observables, the occupation is used as a
weight.

Concretely, for any one-body observable Ô the expectation value is

```math
\langle \hat O \rangle \;=\; \sum_{n,\mathbf{k}}\, f_{n,\mathbf{k}}\;
\langle\psi_{n,\mathbf{k}}|\,\hat O\,|\psi_{n,\mathbf{k}}\rangle .
```

Here `f_{n,k}` is the occupation as INQ stores it: it already absorbs the
k-point weight `w_k` and the spin multiplicity (factor 2 for spin-paired
runs and Fermi–Dirac population). For a fully spin-paired, equally-
weighted shifted MP-2×2×2 mesh of 8 points, a doubly-occupied orbital
gets `f = 2 × 1/8 = 0.25`; that's the value we see for bands 0–22 at
every k-point in the eigenvalues CSV.

#### Band 30 — well inside the conduction band

![band 30 Re psi grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_030_re_psi_grid.png)
![band 30 Im psi grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_030_im_psi_grid.png)
![band 30 density grid](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/band_030_density_grid.png)

The final band visualised here is one of the conduction bands. The
visualisation shows that it is **quite widely spread**. Hence, an
electron, when in this state, is highly delocalised. This is why these
electrons contribute to current.

### Evidence of phase twist (bands 12 and 24)

User question: *what graph would I plot to highlight the twist?*

The cleanest single-figure answer is a **1D line cut along the cell
diagonal [111] of `Re ψ_{n,k}(s)` and the unwrapped phase
`arg ψ_{n,k}(s)`, overlaid for all 8 k-points** at fixed band. The
real-part panel shows the carrier wave that distinguishes one k from
another; the unwrapped-phase panel literally counts the radians of
twist accumulated along the diagonal — at Γ this is flat (zero twist),
at the BZ corners it accumulates ≈ ±π over the lattice scale.

#### Band 12 (deep 1s manifold)

![Phase twist evidence — band 12: Re ψ and unwrapped phase along [111]](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/phase_twist_band_012.png)

For band 12 the Re ψ trace looks like **localised spikes at each Li
site**, and the unwrapped phase is essentially flat between sites,
jumping in steps as the line cut passes through one tightly-localised
1s lobe to the next. The k-dependence here is mostly the *sign pattern
between sites*, not a smooth phase modulation — exactly the
tight-binding picture for a band where inter-site overlap is tiny.

#### Band 24 (conduction)

![Phase twist evidence — band 24: Re ψ and unwrapped phase along [111]](attachments/2026-05-04_orbitals_at_each_band_at_each_k_point/phase_twist_band_024.png)

For band 24 the Re ψ trace is **smooth and sinusoid-like across the
diagonal**, and the unwrapped-phase panel shows a clean ramp whose
slope changes with k-point — different k-points produce different
phase-twist rates along the diagonal. This is the free-electron-like
Bloch picture: the orbital wavefunction is well approximated by
`u_{n,k}(r) e^{i k·r}` with a smooth periodic envelope and the carrier
phase doing all of the k-dependent work.

These two figures together bracket the spectrum: deep cores (band 12)
have *sign-pattern* k-dependence, conduction bands (band 24) have
*phase-ramp* k-dependence.

## Open questions / next steps

- Why does the density `|ψ_{n,k}|²` look so similar across k for a
  given band, even when Re ψ and Im ψ differ wildly? The answer is that
  the Bloch phase is squared away, but the **observables that *do*
  depend on k** are momentum-like or current-like operators
  (`p̂`, `ĵ`, ...), where the gradient picks up the phase twist. This is
  the next sanity check to add to inqview.
- Apply the same `dump_orbitals_per_kpoint` driver to the 54-atom Li
  production checkpoint (`li_54_2x2x2_T200`) so the same intuition
  carries to the production-scale system.
