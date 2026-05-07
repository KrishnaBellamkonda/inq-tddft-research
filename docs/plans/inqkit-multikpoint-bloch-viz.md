# Plan: Multi-k-point inqkit hardening + Bloch-orbital visualiser

Branch: `features/inqkit-multikpoint` (worktree of main repo, **not**
QuantumKickExtension). Independent of the `features/li-extensive` Phase 10b
work in QKE.

## Context

Phase 10b in `docs/handovers/li_extensive_kick.md` requires per-(k, n)
orbital diagnostics on the 54-atom Li supercell. Before producing those on
the production run, the user wants the underlying inqkit primitives
hardened for multi-kpoint use and **smoke-tested on simple systems** so
the visualisation can be trusted. The aim of the visualisation is
pedagogical — see what a Bloch orbital ψ_{n,k}(r) = e^{i k·r} u_{n,k}(r)
looks like in real space, how the phase twist depends on k, and how
|ψ|² compares with Re ψ / Im ψ.

Existing pieces (verified):
- `inq-stack/include/inqkit/fields/complex_field_3d.hpp` — struct.
- `inq-stack/include/inqkit/io/complex_field_3d_writer.hpp` — writer.
- `inq-stack/include/inqkit/observables/eigenvalue_dump.hpp` — multi-kpoint
  band CSV (already used by `run_save_gs_2x2x2_T200`).
- `inq-stack/include/inqkit/fields/density.hpp::orbital(electrons, n, k=0)`
  — single-orbital |ψ|², takes a kpoint_index but **defaults to 0**, so
  most call sites silently read only the first k-point.
- `inq-stack/python/inqview/postprocess/orbitals.py` — exists; need to
  audit whether it's k-aware.

## Goal

Stage A — produce a small, focused C++ helper and Python visualiser that
write **complex** ψ_{n,k}(r) for any (n, k), plus a bundle of
band-major figures comparing across k-points.

Stage B — run the bundle on two smoke-test systems before touching the
54-atom Li run:
- **Tier S1 — atomic N** in a cubic box, 2×2×2 shifted MP. Single atom,
  open-shell, gives a *molecular-like* (n,k)-orbital that is highly
  localised; the k-dependence shows mostly as a global phase twist on top
  of an essentially k-invariant envelope. Validates that the Re/Im VTIs
  capture the phase correctly.
- **Tier S2 — 2×2×2 Li supercell, 2×2×2 shifted MP**. Metallic 16-atom
  BCC chunk (one-quarter of the production system in each direction).
  Cheap (≤ 5 min on GPU 1, much smaller than the 54-atom run). Exposes
  the actual Bloch character: spatially extended u_{n,k}(r) with non-
  trivial k-dependent phase modulation across the unit cells.

After the user is satisfied with the S2 visualisation, the same code is
applied to the existing `li_54_2x2x2_T200` checkpoint without re-running
SCF.

## Recommended approach

### Phase A — inqkit multi-kpoint primitives

#### A1. Audit (read-only)

Walk every header under `inq-stack/include/inqkit/{fields,observables,wavepacket,real_time,screens,io}/`. For each public function, tag with one of:
`kpoint-summed` (already correct), `requires-explicit-kpoint`,
`single-kpoint-assumed-needs-fix`, `n/a`. Write the audit to
`docs/handovers/inqkit_multikpoint_audit.md`. Single-rank-basis
assumptions are documented but not changed in this plan.

Expected findings (from earlier session):
- `density::total` → kpoint-summed ✓
- `density::orbital(electrons, n, k=0)` → requires-explicit-kpoint;
  drop the default to force callers to be explicit.
- `wavepacket::inject_into_last_extra_state` → single-kpoint-assumed.
  Add an explicit `kpoint_index` parameter (no default for multi-k cases).
- `observables::orbital_overlap` → single-kpoint; out of scope here, will
  be addressed when needed for propagation diagnostics.
- `eigenvalue_dump` → already kpoint-aware ✓
- `state_energy_writer` — needs an audit pass; out of scope unless
  trivially fixable.

#### A2. New helper: `inqkit::fields::orbital_complex(electrons, n, k)`

Signature:

```cpp
namespace inqkit::fields {
  // Returns ψ_{n,k}(r) as a complex 3D field. No default for k —
  // multi-k callers must be explicit.
  ComplexField3D orbital_complex(
      inq::systems::electrons const& electrons,
      int orbital_index,
      int kpoint_index);
}
```

Implementation reuses the FFT-shift pattern from
`density.hpp::total` (lines 86–101), but instead of squaring stores
`Re ψ + i Im ψ` in `ComplexField3D::values`. Throws on multi-rank basis,
on out-of-range orbital/kpoint indices, and when the requested orbital is
not on the local rank.

Update `density::orbital` to take `kpoint_index` with no default (call
sites already in the codebase pass 0 or a real index — verify with grep
and patch).

#### A3. Smoke test the complex orbital primitive

`Tutorial/_inqkit_tests/orbital_complex_smoketest/run.cpp` — single H atom in
a 6 Å cubic box, Γ-only. Run SCF, then call `orbital_complex(electrons, 0, 0)`,
write Re/Im as separate binary VTIs via `RealField3DWriter`, and write the
combined complex field via `ComplexField3DWriter`. Pass criteria:
- ∫|ψ|² d³r ≈ 1 (norm).
- Im L² ≪ Re L² (for Γ, ψ is real up to global phase).
- Visual: 1s-like Gaussian peak.

This is the minimal known-case test from `development-feedback-loop.md`.

### Phase B — Bloch orbital visualiser

#### B1. C++ driver: `dump_orbitals_per_kpoint.cpp`

Generic driver that:
1. Constructs ions + electrons identically to the GS run that produced the
   checkpoint (cell, kpoint grid, cutoff, smearing, extra_states).
2. Loads the checkpoint via `electrons.load`.
3. Reads `eigenvalues.csv` (already present in checkpoint dir per
   `eigenvalue_dump.hpp`) to discover the kpoint count and band ordering.
4. Iterates a user-supplied list of band indices; for each (band, k) writes
   - `band_NNN/re_psi_kKKK.vti` (Re ψ, binary VTI)
   - `band_NNN/im_psi_kKKK.vti` (Im ψ, binary VTI)
   - `band_NNN/density_kKKK.vti` (|ψ|², binary VTI)
5. Appends a row to `orbital_index.csv`:
   `band, kpoint_index, kx, ky, kz, weight, occ, evalue_ha, re_l2, im_l2, density_l2`.

Output rooted at `<run_dir>/results/analysis/ground_state/orbitals_per_kpoint/`.
Folder layout is **band-major** so opening any `band_NNN/` directory in
ParaView and animating the "step" index sweeps across the BZ at fixed band
— the headline visualisation the user asked for.

A copy of this driver lives in each smoke-test tier (S1, S2) and one in the
production run's directory; they differ only in the cell + kpoint grid +
band list, so use a small helper header
`Tutorial/_inqkit_tests/_orbital_dump_helpers.hpp` to keep the shared logic
DRY without forcing a new inqkit module.

#### B2. Python visualiser: `inqview.postprocess.orbitals_per_kpoint`

Reads `orbital_index.csv` and the per-band VTI folders, produces:
- `band_NNN/re_psi_grid.png` — for each band, an N_k-panel grid of mid-cell
  slices of Re ψ. Fixed colour scale across the panels of one band.
- `band_NNN/im_psi_grid.png` — same for Im ψ.
- `band_NNN/density_grid.png` — same for |ψ|², on a separate (positive)
  scale.
- `bands_summary.png` — eigenvalue vs k-point index for each chosen band
  (toy band-structure-fragment; supplements the
  `gamma_transitions.py` histogram already in inqview).
- `paraview_recipe.md` — short note on how to open `band_NNN/re_psi_*.vti`
  in ParaView and animate.

Wired into `inqview.postprocess.pipeline.PHASES` as a new optional phase
`orbitals_per_kpoint`, gated on the presence of `orbitals_per_kpoint/`
under `analysis/ground_state/`.

### Phase C — Smoke tests

#### Tier S1 — atomic N, Γ-only + 2×2×2 shifted MP

`Tutorial/_inqkit_tests/orbital_per_kpoint_S1_nitrogen/`:

- `run.cpp` — single N atom in 8 Å cubic box, periodic, MP 2×2×2 shifted,
  PBE, cutoff 30 Ha (Tutorial-grade, not production), Fermi smearing
  300 K, extra_states 4.
- After SCF + `electrons.save`, calls
  `inqkit::observables::dump_eigenvalues` and the new
  `dump_orbitals_per_kpoint` driver for bands {0, 1, 2, 3}.
- Pass criteria:
  - ∫ρ d³r = 5 (N has 5 valence electrons with ONCV PBE; verify exact value
    from pseudopod).
  - For every (n, k): ∫|ψ_{n,k}|² d³r = 1.000 ± 1e-3.
  - At Γ-shifted (the central kpoint of the shifted 2×2×2 grid is not
    exactly Γ — record its position from `eigenvalues.csv`), the Im L² of
    bound orbitals is small but non-zero (atom is heavy enough that this
    won't be zero).
  - Visualisation shows a 1s-like core and three p-like lobes.

Cost: ≪ 1 GPU-minute (open-shell single atom). Run on either GPU.

#### Tier S2 — 2×2×2 Li BCC supercell, 2×2×2 shifted MP

`Tutorial/_inqkit_tests/orbital_per_kpoint_S2_li_2x2x2/`:

- 2×2×2 supercell of BCC Li (16 atoms, 48 valence electrons), cubic edge
  2 × 3.51 Å = 7.02 Å.
- Same fractional-coord generator as the production GS, just with the
  outer loop bound 2 instead of 3.
- MP 2×2×2 shifted (same as production), PBE, cutoff 30 Ha (lighter than
  74 Ry production — this is a teaching artefact, not a quantitative
  benchmark), Fermi 400 K, extra_states 8.
- Bands selected for visualisation: 1 (deep), 12 (mid), 24 (Fermi-region
  for 24 doubly-occupied bands), 30 (empty). Tunable.
- Pass criteria:
  - ∫ρ = 48 (electron count).
  - For every (n, k): ∫|ψ_{n,k}|² d³r = 1.000 ± 1e-3.
  - For k away from Γ: Re ψ shows clear sinusoidal phase modulation
    across the supercell with wavelength ≈ |k|^{-1}; Im ψ is shifted by π/2
    of the same wavelength. |ψ|² is much more uniform than either Re or
    Im — this is the headline pedagogical observation.
  - Eigenvalues align with a standard Li BCC band-structure literature
    reference (Splines, Mahan; see `docs/sources/`).

Cost target: SCF ≤ 5 minutes on GPU 1, single rank. Negligible compared to
the production 20 h propagation.

#### Tier production-apply — apply to existing 54-atom Li GS

After S1 and S2 pass and the user signs off:

- Add a new directory
  `QuantumKickExtension/inq-codebase/Li/run_dump_orbitals_2x2x2_T200/`
  with a `dump.cpp` that loads the existing checkpoint at
  `inq-codebase/Li/checkpoints/li_54_2x2x2_T200/` (no SCF; the
  checkpoint already exists from Phase 1).
- Same band selection logic as S2; user picks the band list (default:
  {1, 40, 81, 100} per earlier plan).
- Output under
  `run_dump_orbitals_2x2x2_T200/results/analysis/ground_state/orbitals_per_kpoint/`
  so it is distinct from the GS run's results tree.
- Run `inqview.postprocess.orbitals_per_kpoint` to produce the figure
  bundle.

### Phase D — wire-up

- Add `dump_orbitals_per_kpoint` and the new pipeline phase to the
  postprocess driver `inq-codebase/Li/scripts/postprocess_run.py`.
- Update `docs/journals/quantumkickextension.md` with a new entry once
  the production-apply step is done (uses the `journal-writing` skill).
- Update `docs/handovers/li_extensive_kick.md` Phase-10b section to
  reflect the dependency: production-apply blocked on S2 sign-off.

## Critical files

### To read

| File | Why |
|---|---|
| `inq-stack/include/inqkit/fields/density.hpp` (43–102) | FFT-shift pattern, k-sum reference. |
| `inq-stack/include/inqkit/fields/complex_field_3d.hpp` | Field struct. |
| `inq-stack/include/inqkit/io/complex_field_3d_writer.hpp` | Writer interface. |
| `inq-stack/include/inqkit/observables/eigenvalue_dump.hpp` | Existing kpoint-aware dump for band-list discovery. |
| `inq/src/systems/electrons.hpp` (600–720) | `kpin()`, `kpoint_index()`, save/load with k-points. |
| `Tutorial/li-bcc/li_bcc.cpp` | Reference primitive Li BCC; minimal multi-k example. |

### To create

| File | Phase |
|---|---|
| `inq-stack/include/inqkit/fields/orbital.hpp` (or new `orbital_complex.hpp`) — add `orbital_complex(...)` | A2 |
| `Tutorial/_inqkit_tests/orbital_complex_smoketest/run.cpp` | A3 |
| `Tutorial/_inqkit_tests/_orbital_dump_helpers.hpp` | B1 |
| `Tutorial/_inqkit_tests/orbital_per_kpoint_S1_nitrogen/run.cpp` | C-S1 |
| `Tutorial/_inqkit_tests/orbital_per_kpoint_S2_li_2x2x2/run.cpp` | C-S2 |
| `inq-stack/python/inqview/postprocess/orbitals_per_kpoint.py` | B2 |
| `QuantumKickExtension/inq-codebase/Li/run_dump_orbitals_2x2x2_T200/dump.cpp` | C-prod |
| `docs/handovers/inqkit_multikpoint_audit.md` (rolling) | A1 |
| `docs/handovers/inqkit_multikpoint_bloch_viz.md` (rolling) | all |

## Verification

Each step passes through the **write → known-case-test → fix → confirm**
loop required by `.claude/rules/development-feedback-loop.md`:

| Step | Test | Pass criterion |
|---|---|---|
| A1 | Audit completeness | Every public function in scope is tagged. |
| A2 | `orbital_complex` round-trip | For Γ on H atom: ∫|ψ|² = 1; Im L² < 1e-6. |
| A3 | Smoke test compiles + runs | Outputs Re/Im VTIs; Python reader can open them. |
| C-S1 | N atom, multi-k | ∫ρ = 5 (or pseudopod-reported value); ∀(n,k): ∫|ψ|² = 1. |
| C-S2 | Li 2×2×2 supercell | ∫ρ = 48; ∀(n,k): ∫|ψ|² = 1; visible Bloch phase modulation across cells in Re ψ for k ≠ Γ. |
| B2 | Visualiser on S2 | Grid figures readable; colour scale fixed within a band; eigenvalue summary plot matches `eigenvalues.csv`. |
| Prod-apply | 54-atom Li | Same checks, no SCF re-run. |

End-to-end:
```bash
cd Tutorial/_inqkit_tests/orbital_complex_smoketest && inq-run --cpu
cd ../orbital_per_kpoint_S1_nitrogen && CUDA_VISIBLE_DEVICES=1 inq-run
cd ../orbital_per_kpoint_S2_li_2x2x2 && CUDA_VISIBLE_DEVICES=1 inq-run
python -c "import inqview.postprocess.orbitals_per_kpoint as m; m.run('.')"
# user reviews S2 figures → on sign-off:
cd /local/data/public/skcb2/tddft/QuantumKickExtension/inq-codebase/Li/run_dump_orbitals_2x2x2_T200 && inq-run
```

## Risks and fallbacks

1. **`ComplexField3DWriter` not yet emitting binary VTI.** Fallback: write
   Re ψ and Im ψ as separate `RealField3DWriter` outputs (the path already
   tested with `density_rt_total/`) and skip the combined complex VTI.
2. **MP-shifted 2×2×2 has no exact Γ point.** The "k closest to Γ"
   convention from `gamma_transitions.py` carries over: the smoke-test
   pass criteria are stated in terms of "the k closest to Γ", not Γ
   itself.
3. **N atom open-shell.** If pseudopod ships only spin-restricted ONCV N,
   either accept restricted treatment for the smoke test or substitute He
   (closed shell, 2 electrons). Decide based on the first SCF attempt; a
   substitute is acceptable as long as the multi-k test still exercises
   non-trivial Bloch phases.
4. **State distribution across ranks.** `dump_orbitals_per_kpoint` runs
   single-rank for these systems by construction (k-parallelism with one
   GPU = one rank handles all k). For the production 54-atom dump where
   INQ launches multiple ranks, the dump driver is invoked with
   `mpirun.openmpi -np 1` to keep all kpoints on one rank.

## Out of scope

- Modifying INQ itself.
- Multi-rank gather/scatter for the orbital visualiser (future work).
- Spin-polarised electrons (no use case yet).
- Re-running the 54-atom Li GS at different k-grids.

## Approvals required before execution

1. Tier ladder S1 → S2 → production-apply, in that order, with user
   sign-off between S2 and production-apply.
2. Smoke-test home: `Tutorial/_inqkit_tests/` (proposed) — fits the
   existing pattern (`vti_writer_smoketest` already lives there).
3. Whether N atom (open-shell) or He (closed-shell) for S1. **Default
   proposal: N as requested**, fall back to He on pseudopod issue.
4. Bands to visualise on the production 54-atom run: default {1, 40, 81,
   100} per the earlier plan; tunable.
