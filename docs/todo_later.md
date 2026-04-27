# To-Do Later

Deferred tasks for the coronene WP scattering project. Items here are not blockers for the current simulation but should be addressed before paper-quality results are claimed.

---

## Legacy run-directory cleanup

After the new replication framework (see `docs/plans/coronene-replication.md`)
has reproduced the important runs, lump every legacy buggy coronene run under a
single `legacy/` subtree:

- `ResearchProject/systems/coronene/04_leed_simulation/` (run_001–005, all `z=L/2` shifted geometry)
- `ResearchProject/systems/coronene/coronene-wp-rt/` (run_01–06, all `z=L/2` shifted geometry)
- `Tutorial/coronene-leed/run_diagnoses/run_01_tight_scf` … `run_05_quarter_coronene` (buggy z=L/2 xyz files)
- `ResearchProject/systems/coronene/run_propagate_paper_replica/` and `run_save_gs_paper_replica/` once the new `save_gs/`/`run_*/` framework supersedes them

Move under `ResearchProject/systems/coronene/legacy/` (or per-tree `legacy/`)
rather than deleting; the on-disk artefacts are useful for diff-vs-corrected
comparisons.

---

## Integration tests

- **End-to-end smoke test:** Run `04_leed_simulation` on CPU with a reduced parameter set (E_cut=20 Ha, 10 steps, 1 snapshot) and verify that all output files are produced with sensible values.
- **Restart consistency:** Save the ground state density at step 2, reload it at the start of step 5, and confirm that the first few TDDFT steps reproduce identical results. Tests the SCF + TDDFT hand-off.
- **GPU/CPU parity:** Run the same short test on CPU and GPU, compare `leed_pattern.txt` element-wise to tolerance 1e-5. Validates the CUDA backend.

---

## Absorbing boundary conditions

- The current simulation has no absorbing boundaries. After scattering, reflected electrons accumulate and interfere with the incoming WP. For a clean LEED pattern the boundary should absorb outgoing probability flux.
- **Action:** Add a complex absorbing potential (CAP) in the z-margins. INQ supports `perturbations::absorbing_walls` — investigate API and add to `04_leed_simulation/config.hpp`.
- **Reference:** de Boor & Saalfrank, J. Chem. Phys. 147, 224107 (2017) for CAP parameterisation in real-space TDDFT.

---

## MPI-aware slice extraction

- `utils::extract_density_slice()` currently assumes a single MPI rank owns all grid points. Under MPI decomposition, `local_to_global` is required but the global slice array is only allocated on root.
- **Action:** Add an `MPI_Allreduce` (sum) over the slice array before returning, so all ranks contribute their local portion. Guard the write inside `if(electrons.root())`.
- **Risk:** Without this fix, MPI runs with >1 rank will silently produce incomplete density slices and a wrong LEED pattern.

---

## Transmission vs. reflection planes

- Paper Fig. 2 shows the LEED pattern on the observation plane at z = +D (same side as the incident WP). A transmission pattern at z = −D (far side) would also be physically meaningful.
- **Action:** Add a second accumulator `trans_accum` at `z = -cfg::WP_D_IMPACT_BOHR` and write `results/transmission_pattern.txt` alongside the LEED pattern.

---

## E_cut non-monotonicity investigation

- Our sweep shows E_total rises from 40→60 Ha (+90 meV), which is atypical. The leading hypothesis is pseudodojo_pbe projector completeness / XC aliasing above the natural cutoff.
- **Action:** Re-run the sweep with a different pseudopotential family (e.g. ONCV-PBE from Schlipf & Gygi) to check whether the non-monotonicity is pseudopotential-specific.
- **Status:** Non-blocking — 40 Ha (energy minimum) is justified for the current simulation. Must be resolved before claiming paper-quality accuracy.

---

## Geometry relaxation

- The coronene geometry is from an idealised XYZ file (bond lengths from a crystallographic template), not from a DFT-relaxed structure. The max force is 0.063 Ha/bohr (C outer ring), well above the typical 0.001 Ha/bohr convergence criterion.
- **Action:** Run a geometry relaxation (ground_state with force convergence) and compare forces before/after. Use the relaxed XYZ for the LEED simulation.
- **Note:** Forces in the plane will affect the molecular charge density and therefore the scattering potential. For a quantitative LEED comparison this matters.

---

## Analysis pipeline

- `04_leed_simulation/analysis.py` (Fig. 1 and Fig. 2 replication) is not yet written.
- **Action:** After a successful LEED run, write `analysis.py` to:
  - Animate `snapshot_t*.txt` frames (Fig. 1 comparison)
  - Plot `leed_pattern.txt` as a 2D intensity map and compare to paper Fig. 2

---

## VESTA geometry visualisation

- No visualisation of the coronene atomic structure has been made.
- **Action:** Convert the `.xyz` file to a VESTA-readable format and verify the molecular geometry matches Fig. 1 of Tsubonoya et al.

---

## Screen-z FFT-natural mapping bug (Phase 3)

**Bug:** `PlaneScreen::iz_nearest`
(`inq-stack/include/inqkit/screens/plane_screen.hpp:50-57`) clamped negative
physical z to grid index 0 instead of wrapping into INQ's FFT-natural
ordering (where index 0 = physical centre, upper-half indices = negative z).
Every transmission screen at `z < 0` was therefore sampling the molecule
plane (the static coronene electron cloud), not the requested z-plane.

**Symptom:** in every `all_screens_grid.png` prior to Phase 3, the
transmission half of the grid showed the static coronene electron density
rather than a diffraction pattern.

**Affected runs:** every run produced before Phase 3 — all 10 in the
framework. Backscattering screens (z ≥ 0) were unaffected because their
physical z mapped directly to a non-negative array index.

**Resolution:** FFT-natural wrap added in Phase 3:
`int iz_nat = ((iz % Nz) + Nz) % Nz;`. The branch
`coronene-fft-fixed-base` re-runs `run_base` for verification;
`coronene-fft-fixed-rerun` then re-runs the other 9. After both branches
merge, every total/, paper-window, and physics-window screen on the
transmission side is physically meaningful.

---

## Complex orbital wavefunction export was not fftshifted (Phase 3)

**Bug:** `inqkit::fields::orbital::wavefunction`
(`inq-stack/include/inqkit/fields/orbital.hpp:107-118`) iterated `(ix, iy, iz)`
and read `hc[ix][iy][iz][...]` directly **without** the `fft_shift_index`
mapping that the sibling `density::total` / `density::orbital` exporters
apply. The `ComplexField3D` it returned was therefore in FFT-natural order
while metadata claimed left-to-right physical layout (origin = −L/2).

**Symptom:** any VTI viewer or Python slice consumer of the complex
orbital wavefunction saw the WP at scrambled spatial positions relative
to the metadata origin. `density::orbital` (real density) was correctly
shifted, so total / system / wp **density** outputs were correct; only the
**complex wavefunction** export was scrambled.

**Affected runs:** every run that wrote
`results/raw/wavepacket/wavefunction_wp_initial.vti` (all 10).

**Resolution:** `fft_shift_index` mapping added in Phase 3, mirroring
`density.hpp:88-99`. The Phase-3 re-run regenerates the file correctly
on every run.
