# understanding_pseudopotentials

**Title:** Understanding pseudopotentials (semi-core vs valence-only)
**Status:** complete
**Linked entries:**
- [`2026-05-04_orbitals_at_each_band_at_each_k_point`](2026-05-04_orbitals_at_each_band_at_each_k_point.md)
- [`2026-05-03_run_propagate_v0p0123_extensive`](../quantumkickextension/2026-05-03_run_propagate_v0p0123_extensive.md)

## Observations

Working with the Li 2×2×2 supercell simulation
([`2026-05-04_orbitals_at_each_band_at_each_k_point`](2026-05-04_orbitals_at_each_band_at_each_k_point.md))
gave concrete evidence that the ONCV PBE-1.2 pseudopotential we are
using is a **semi-core** pseudopotential rather than a smooth
valence-only one. A "smooth norm-conserving valence-only" Li
pseudopotential would freeze the 1s² core into the ionic potential and
expose only the 2s¹ valence electron to the simulation; a semi-core
pseudopotential keeps the 1s² electrons inside the simulation and only
screens the bare nuclear charge.

### Evidence

1. **Electron count.** The 16-atom Li supercell run reports
   `num_electrons = 48`. A valence-only pseudopotential would give
   `16 × 1 = 16` electrons. 48 is `16 × 3`, i.e. the full atomic
   electron count. Source: `run_summary.txt` at
   `Tutorial/_inqkit_tests/orbital_per_kpoint_S2_li_2x2x2/results/`,
   identical conclusion holds for the 54-atom production run with
   `num_electrons = 162 = 54 × 3`
   ([`2026-05-03_run_propagate_v0p0123_extensive`](../quantumkickextension/2026-05-03_run_propagate_v0p0123_extensive.md)).

2. **Band-structure shape.** The eigenvalue table shows two clearly
   separated manifolds: a tight cluster of 16 bands at E ≈ −47.3 eV
   (bands 0–15) and a broad cluster of 8 bands centred near the Fermi
   level at E ≈ −1 to −3 eV (bands 16–23). The deep cluster has nearly
   k-independent eigenvalues (≪ 0.01 eV spread across the BZ), and the
   per-(n,k) orbital plots in the linked entry show core-like 1s lobes
   confined to each Li site. A valence-only pseudopotential would have
   no deep manifold — only the dispersive 2s/2p valence band would
   appear. The **deep manifold *is* the 1s² core, included explicitly**.

3. **Per-orbital morphology.** The Re ψ / Im ψ / |ψ|² grids for band 1
   in the linked entry show tight Gaussian-like lobes pinned to atomic
   sites with negligible inter-site overlap — exactly what 1s-on-each-
   atom would look like. A pseudopotential that *removed* the 1s²
   would have band 0 instead start at 2s level (E ≈ −3 eV), and the
   morphology would already be delocalised at the lowest band.

### Why we settled on this conclusion

The combination of (a) the right total electron count for an "all
3-electron Li" picture, (b) the observed deep+shallow manifold split
in eigenvalues, and (c) the visibly atomic morphology of the deep
manifold orbitals, all point at the SG15 ONCV PBE-1.2 file shipped at
`inq/install/share/pseudopod/pseudopotentials/quantum-simulation.org/sg15/Li_ONCV_PBE-1.2.upf.gz`
being a **3-valence-electron semi-core pseudopotential**. For Li this
distinction matters: the 1s² electrons participate in the explicit KS
problem rather than living silently inside the ionic potential, so the
plane-wave / real-space cutoff must be high enough to resolve the
sharp 1s wave function (which is exactly why the QBall reference and
this work both run at 74 Ry plane-wave / ≈ 37 Ha real-space cutoff).

### Implication for downstream interpretation

Whenever we state "occupied bands" for any Li run in this project, we
mean the union of (1s manifold) ∪ (filled portion of the 2s manifold).
For a 16-atom supercell that's 16 + 8 = 24 doubly-occupied bands; for
the 54-atom production run that's 54 + 27 = 81 doubly-occupied bands
(matching `num_electrons / 2 = 162 / 2 = 81`).

## Open questions / next steps

- Cross-check by reading the SG15 ONCV PBE-1.2 UPF header for the
  electron-configuration field and recording it in
  `docs/sources/li-pseudopotential.md`.
- Confirm whether INQ's pseudopod loader exposes a way to drop the 1s²
  core to a frozen valence-only Li pseudopotential, in case we want to
  test a cheaper run for k-grid convergence studies.
