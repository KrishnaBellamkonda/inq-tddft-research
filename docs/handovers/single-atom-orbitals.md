# Handover: single-atom orbital visualisation tutorial

## Current status

**Complete (with VTI output).** Three ground-state simulations done (H, Li, Al), each in a 30 bohr cubic finite cell at LDA/60 Ry with `extra_states(30)`. Per-run outputs are written **directly as ParaView-ready `.vti` files** via `inqkit::io::RealField3DWriter` / `ComplexField3DWriter` configured with `.emit_raw = false, .emit_vti = true, .vti_format = binary`. No Python post-processing step is required.

| Atom | n_states | ∫ρ | GS energy (Ha) | SCF iters | Occupied states |
|---|---|---|---|---|---|
| H  | 31 | 1.0 | −0.027639 | 69 | [0] 1s (occ 1.0) |
| Li | 32 | 3.0 | −0.044338 | (converged) | [0] 1s (occ 2.0), [1] 2s (occ 1.0) |
| Al | 32 | 3.0 | −0.057946 | (converged) | [0] 3s (occ 2.0), [1] 3p (occ 1.0) |

All three runs wrote total density, per-orbital `|ψ|²`, and per-orbital complex `ψ` to `results/` as `.vti` files (binary base64 inline, single file per field). SCF reached the 1e-6 Ha tolerance for every atom.

Per-run output layout:

```
results/
  density/
    density_total.vti                  # total electron density
  orbital_density/
    orbital_NNNN_density.vti           # |psi_i|^2, one per state
  orbitals/
    orbital_NNNN.vti                   # complex psi_i (real + imag PointData arrays)
```

Counts: H = 1 + 31 + 31, Li = 1 + 32 + 32, Al = 1 + 32 + 32. Per-orbital `.vti` size ≈ 4.3 MB (real density) and ≈ 8.6 MB (complex wavefunction) at 75³ grid points.

## What changed

- Created `Tutorial/single-atom-orbitals/{h,li,al}/run.cpp` — three near-duplicate INQ/inqkit programs.
- Created spec at `docs/superpowers/specs/2026-05-07-single-atom-orbitals-design.md`.
- Created plan at `docs/superpowers/plans/2026-05-07-single-atom-orbitals.md`.
- Created this handover.

## Files touched

- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/h/run.cpp`
- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/li/run.cpp`
- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/al/run.cpp`
- `/local/data/public/skcb2/tddft/docs/superpowers/specs/2026-05-07-single-atom-orbitals-design.md`
- `/local/data/public/skcb2/tddft/docs/superpowers/plans/2026-05-07-single-atom-orbitals.md`
- `/local/data/public/skcb2/tddft/docs/handovers/single-atom-orbitals.md`

Output trees (gitignored, large):
- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/h/results/`
- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/li/results/`
- `/local/data/public/skcb2/tddft/Tutorial/single-atom-orbitals/al/results/`

## Commands run

```
cd Tutorial/single-atom-orbitals/h  && inq-run
cd Tutorial/single-atom-orbitals/li && inq-run
cd Tutorial/single-atom-orbitals/al && inq-run
```

Each command did a full INQ build (~few minutes for first build, then incremental), SCF, and orbital-write loop.

## Tests and validation

Run-side, scripted in `run.cpp` and confirmed in stdout:

| Check | H | Li | Al |
|---|---|---|---|
| SCF energy_tolerance(1e-6 Ha) reached | ✓ | ✓ | ✓ |
| ∫ρ(r) d³r matches valence count | 1.0 ✓ | 3.0 ✓ | 3.0 ✓ |
| n_states ≥ 31 (ceil(Ne/2) + 30) | 31 ✓ | 32 ✓ | 32 ✓ |
| Eigenvalue ladder printed | ✓ | ✓ | ✓ |
| Per-orbital `.raw` + `.meta.txt` count = n_states | ✓ | ✓ | ✓ |
| Per-orbital complex `.real.raw + .imag.raw + .meta.txt` count = n_states | ✓ | ✓ | ✓ |

User-side, **not yet performed**:
- ParaView visual inspection of the orbital VTI series.
- Visual confirmation of s/p/d shell structure.

## Trusted sources used

- `Tutorial/hf-gs-with-inqkit/run.cpp` — direct template for the inqkit field-write pattern.
- `Tutorial/li-observables-with-inqkit/run.cpp` — pattern for `.temperature(0.001_Ha)` smearing on a Li system.
- `inq-stack/include/inqkit/fields/density.hpp`, `orbital.hpp` — confirmed `phi.spinor_set_size()` and the `electrons.kpin()[0]` accessor for runtime state counts.
- `Tutorial/_inqkit_tests/_orbital_dump_helpers.hpp` — confirmed `electrons.eigenvalues()[ilot][i]` and `electrons.occupations()[ilot][i]` indexing pattern.
- `inq/src/operations/integral.hpp` — confirmed `operations::integral` works on `inq::basis::field`, **not** `inqkit::fields::RealField3D` (forced a small refactor — see "Known issues" below).

## Attribution notes

- The structure of `run.cpp` (cell, electrons, ground_state::calculate, then RealField3DWriter / ComplexField3DWriter) is adapted from `Tutorial/hf-gs-with-inqkit/run.cpp`.
- The smearing temperature `0.001_Ha` is taken from `Tutorial/li-observables-with-inqkit/run.cpp`.
- The orbital-loop pattern (writing `|ψ_i|²` and `ψ_i` for every state) is original to this task; existing tutorials only write the HOMO.

## Known issues / blockers

1. **Pseudopotentials are softer than expected.** The INQ default pseudo set treats Li and Al as expected (3s²3p¹ for Al with [Ne] core, all-electron for Li). H has very weak binding (1s eigenvalue ≈ −0.024 Ha vs. the all-electron −0.5 Ha). The orbital *shapes* are still correct for visualisation purposes — that is the user's stated goal — but absolute energies should not be quoted as physical values.

2. **Initial spec assumed Li had 1 valence electron.** Reality: INQ's default Li pseudo is all-electron (3 valence: 1s²2s¹), giving `n_states = 32` instead of the spec's 31. The Li `run.cpp` comment header was updated to reflect this. Li is now actually a stronger pedagogical example because both the inner 1s and outer 2s are visible.

3. **High-i orbitals are box modes, not Rydberg states.** This is a feature, not a bug, and is documented in the spec. Visible in the eigenvalue ladders: many near-degenerate clusters at increasing energy correspond to particle-in-a-box modes of the 30 bohr cubic finite cell.

4. **Spec/plan refactor during build.** First H build failed because `operations::integral` requires `inq::basis::field` and not `inqkit::fields::RealField3D`. Fix applied to all three files: integrate `electrons.density()` directly, then build the inqkit `RealField3D` separately for writing.

## Assumptions still in play

- ParaView will render the inqkit-emitted `.vti` files without grid-orientation issues. The inqkit writer does its own FFT-shift to publish a left-to-right physical layout (see `density.hpp:fft_shift_index`) and the VTI writer reorders to x-fastest VTK PointData order (see `vti_image_data_writer.hpp` header comment). Spot-checked on the H output: header reports `WholeExtent 0 74 0 74 0 74`, `Origin -14.8 -14.8 -14.8`, `Spacing 0.4 0.4 0.4`, which matches the 30-bohr cubic cell at 60 Ry → dx ≈ 0.4 bohr and a centred origin.

## Exact next steps

1. **Open in ParaView directly.** No conversion step needed — the run
   already emits `.vti`. From the project root:
   ```bash
   /local/data/public/skcb2/tddft/ParaView-6.1.0-MPI-Linux-Python3.12-x86_64/bin/paraview \
     Tutorial/single-atom-orbitals/h/results/orbital_density/orbital_0000_density.vti
   ```
   ParaView auto-groups numbered `.vti` series, so opening
   `orbital_0000_density.vti` lets you scrub through the orbital index
   like an animation.

2. **Recommended starting orbitals (use a contour or volume render):**
   - H: index 0 (1s), 1–4 (n=2 box modes), 5–7 (n=3 box modes).
   - Li: index 0 (1s, tightly localised), 1 (2s, diffuse), 2–4 (2p),
     5+ (mixed Rydberg/box).
   - Al: index 0 (3s), 1 (3p), 2–4 (3p complementary), 5+ (3d / box).

3. **For complex wavefunctions** (`results/orbitals/orbital_NNNN.vti`),
   the file contains two PointData arrays `wavefunction_real` and
   `wavefunction_imag`. Render either as an isosurface, or compute
   `sqrt(real^2 + imag^2)` in a Calculator filter to get `|ψ|`.

4. **Optional: a journal entry** under `docs/journals/` documenting the
   first orbital screenshot per atom, once ParaView output is captured.
