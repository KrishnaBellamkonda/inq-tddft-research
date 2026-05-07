# Plan: rerun L=50 with the professor-PDF orthonormalisation method

**Status:** drafted, awaiting user approval before execution.
**Target:** test whether the persistent density depression at the WP
launch site (z ≈ 0) seen in
`run_base_n138_L50_E1p5` and `run_base_n162_L50_E1p5` is an artefact
of the current modified-Gram-Schmidt (real-space) WP injection in
`inq-stack/include/inqkit/wavepacket/wavepacket.hpp`.

**Authoritative source:** the professor-PDF method is documented in
`docs/sources/orthonormalisation-professor.md` (with citations to
`ResearchProject/literature/misc/orthonormalization-by-professor.pdf`
and `viewables.hpp`).

---

## 1. The two algorithms in one table

| Aspect | inq-stack current (real-space MGS) | Professor PDF (momentum-space CGS) |
|---|---|---|
| Domain of inner product | real-space integral on grid | sum over reciprocal lattice $\mathbf G$ |
| Gaussian Fourier amplitude | implicit (lives in gridded $\psi_\text{wp}$) | analytic $F(\mathbf G - \mathbf k_0)$ |
| GS variant | **modified** (running residual) | **classical** (against fixed $\phi$) |
| Re-orthogonalisation pass | none (single sweep) | none (single sweep, by construction) |
| Final renormalisation | yes | yes |

For an exactly-orthonormal $\{\psi_n^\text{KS}\}$ the two are
algebraically identical (per source note §"Are the two algorithms
equivalent?"). They differ in floating-point rounding when the basis
is *not* exactly orthonormal — the SCF-converged INQ basis has
typical orthogonality defects $\sim 10^{-7}$, comfortably below the
0.13 max_overlap we observe, so the source-note conclusion is that
the two methods should give the *same* injection density to graphical
precision.

**This rerun therefore tests two things at once:**
1. Whether the hole is sourced by orthonormalisation at all (if it
   disappears under either re-implementation, yes).
2. Whether the real-space-vs-momentum-space difference matters
   numerically at our grid spacing.

---

## 2. Implementation plan (incremental, three options)

### Option A — Twice-is-enough (CGS2): minimal-change first test

The cheapest test. Add a second projection sweep after the first
renormalisation, in `inq-stack/include/inqkit/wavepacket/wavepacket.hpp`
at the end of `inject_into_last_extra_state`:

```cpp
// CGS2 / twice-is-enough re-orthogonalisation
for (int i = 0; i < n_occ; ++i) {
    auto overlap_i = real_space_overlap(psi_occ[i], psi_wp);  // existing helper
    psi_wp -= overlap_i * psi_occ[i];
}
renormalise(psi_wp);  // existing helper
```

If the hole disappears under CGS2, **the cause is residual non-orthogonality
of the inq-stack MGS** (the running residual accumulated a small piece
of each $\psi_i$ that the renormalisation alone could not fix). No
algorithmic change beyond the second sweep is needed.

**Expected outcome:** the `injection_report.txt` should report
`max_overlap` after the second sweep at machine epsilon (~10⁻¹⁵). If
the run-time density still shows the hole, Option A is **not** the
culprit and we move to Option B.

### Option B — Momentum-space classical Gram-Schmidt

The professor-PDF algorithm verbatim. Add a new method
`inject_into_last_extra_state_momentum_space` that:

1. FFT each $\psi_i^\text{KS}$ to plane-wave coefficients $C_i(\mathbf G)$.
2. Compute the analytic Gaussian Fourier amplitude
   $F(\mathbf G - \mathbf k_0)$ in closed form (it is a Gaussian of
   width $1/\sigma$ in $\mathbf G$-space, peaked at
   $\mathbf G = \mathbf k_0$).
3. For each $i$: $\langle \psi_i | \phi\rangle = \sum_\mathbf{G} C_i^*(\mathbf G)\, F(\mathbf G - \mathbf k_0)$.
4. Subtract the projection in plane-wave space:
   $C_\chi(\mathbf G) = F(\mathbf G - \mathbf k_0) - \sum_i \langle \psi_i | \phi\rangle\, C_i(\mathbf G - \mathbf k_0)$
5. iFFT $C_\chi$ → real-space $\psi_\text{wp}$, renormalise.

This is more invasive but tracks the reference algorithm exactly.

### Option C — Both

Run Options A and B, compare the resulting injected $\psi_\text{wp}$
densities point-wise. If they agree to floating-point precision (as
the source note predicts they should), the answer to "is the
real-space vs momentum-space distinction the culprit?" is **no** —
both are equivalent on our grid and both leave the hole, meaning the
hole is **not** an orthogonalisation artefact at all. It is a real
physical feature.

---

## 3. Run protocol

1. **Re-build inq-stack** with the new injection method (Option A
   first; promote to B if A is insufficient).
2. **Reuse the existing GS checkpoints** at
   `save_gs/gs_L50_cubic_N138_dx1p0/` and
   `save_gs/gs_L50_cubic_N162_dx1p0/` — the new orthonormalisation
   only affects WP injection, not the GS solve.
3. **Rerun two propagations**:
   - `run_base_n138_L50_E1p5_orthoX/` (X = "A" / "B" / "AB")
   - `run_base_n162_L50_E1p5_orthoX/`
   Total wall: ~25 min on a single A30 GPU per run.
4. **Postprocess** with the standard `pipeline.py`. Key artefacts:
   - `analysis/density/system_yz.gif` and `system_xz.gif` —
     visual hole-or-not.
   - `analysis/observables/density_fluctuation_l2.png` — should be
     order-of-magnitude smaller at t=0 if the hole is gone.
   - `analysis/observables/all_energies_vs_time.png` — energy-component
     swings should also shrink if the hole was sourcing the bath
     deformation artefact.

---

## 4. Companion test: per-orbital center-of-density vs time (resolves
"is cod_z slope really the WP velocity?")

Independently of the orthonormalisation test, add a postprocess phase
`inqview.postprocess.orbital_cod` that reads the per-step orbital
densities and computes the center-of-density of every KS orbital at
every snapshot. Outputs:

- `analysis/observables/orbital_cod_z_vs_time.png` — line plot, one
  trace per orbital, with the WP slot highlighted.
- `analysis/observables/orbital_cod_summary.csv` — long-format CSV.

**Predictions (pin down the slowdown mechanism):**

- **Bath orbitals stationary**, only WP moves → the cod_z slope
  interpretation is direct, the slowdown is real WP retardation.
- **Bath orbitals near the WP track the WP** → the bath is *being
  dragged* by the WP (consistent with the polarisation-cloud picture
  from the journal entries' §4 "Energy-component bookkeeping"). In
  this case `cod_z(total)` slope under-estimates the *individual* WP
  velocity, and the apparent slowdown is partly the bath catching up.
- **Bath orbitals drift opposite to the WP** → the "countercurrent"
  picture (consistent with the QBall Li 54-atom result, where
  $\langle P_{e,x}\rangle$ has the opposite sign to the kick).

This phase requires the per-orbital density VTIs that are *already*
emitted by the GS save (`raw/vti/density_gs_orbitals/orbital_NNNN.vti`)
plus a new RT extension that emits `density_rt_orbitals/orbital_NNNN_t<step>.vti`.
The latter is **not currently written**; it requires extending
`run_template.hpp` to loop over orbitals at every `WRITE_EVERY` step.
Cost: 81 (or 101) orbitals × 1500 / 2 steps × O(grid) writes — doable
but disk-heavy. Alternative: compute the cod *on the fly* in C++ and
write only a CSV (one row per orbital per step), no VTIs. Recommended.

---

## 5. Verdict criteria (what makes this rerun a "yes" vs "no")

A **report** under `docs/reports/orthonormalisation-rerun-verdict.md`
will be written using the report-writing skill, with the following
sections:

- **Methods**: Option A and/or B, GS reused, propagation parameters.
- **Results**: side-by-side `system_yz.gif` snapshots before/after;
  `density_fluctuation_l2(t=0)` numerical comparison;
  `max_overlap` final value.
- **Discussion**:
  - **Verdict YES (orthonormalisation was the culprit)** if:
    `density_fluctuation_l2(t=0)` drops by ≥ 1 order of magnitude AND
    the visual hole at z ≈ 0 in `system_yz.gif` disappears in the new
    method.
  - **Verdict NO (the hole is physical)** if both above measurements
    are unchanged within 10 %. In that case we promote the hole to a
    physical observable to interpret (likely the bath polarisation
    cloud's signature) and update the journal entries' §1 hypotheses.
  - **Mixed verdict** otherwise: document which features change and
    which don't.
- **Conclusion**: state whether the L=50 / N=162 baseline should be
  re-issued with the new method as the new canonical, or kept as-is.

---

## 6. Open questions before execution

- **CGS2 alone is cheap (~30 lines C++ + smoke test, ~20 min wall),
  Option B is invasive (~200 lines C++ + the FFT bookkeeping).**
  Recommend running A first; only escalate to B if A is insufficient.
  *User to confirm.*
- **Do we need to rerun the L=30 case for completeness?** The L=30
  hole was clearly visible in the partial-shell L=30 entry and it is
  arguably the more important diagnostic since the box is small enough
  for revival contamination to be relevant. *User to confirm whether
  the L=30 rerun is in scope.*
- **What does "professor-given orthonormalisation method" mean
  operationally?** I'm treating it as the momentum-space CGS in the
  professor's PDF (Option B). If the user means the existing real-space
  inqkit code with a CGS2 second sweep (Option A), the implementation
  effort is much smaller. *User to confirm.*
