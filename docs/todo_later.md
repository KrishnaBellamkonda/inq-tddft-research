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

---

## Jellium runs need to be moved to the right folder

The 8 jellium WP-RT runs currently live at
`ResearchProject/jellium/jellium-wp-rt/run_0[1-8]_*/`, predating the
unified `ResearchProject/systems/<material>/` layout used by the
coronene framework. They should be relocated under
`ResearchProject/systems/jellium/` (mirroring `systems/coronene/`)
so:

* The two material subtrees share the same conventions: `shared/`
  configs, `scripts/` for postprocess + dispatch, `run_*/` flat
  siblings, `hypotheses/` for cross-run comparisons.
* `inqview.postprocess.pipeline` can run unchanged on jellium results,
  not just coronene. Today the jellium runs use a flat
  `results/observables.csv` instead of the spec
  `results/raw/observables/observables.csv`, which is why a
  separate `jellium_spectra.py` script is needed.
* The cumulative handover at `docs/handovers/coronene-cumulative.md`
  can grow into a `wp-rt-cumulative.md` covering both materials.

**Action**:
1. Move `ResearchProject/jellium/jellium-wp-rt/run_0[1-8]_*` to
   `ResearchProject/systems/jellium/<run_name>/`.
2. Reshape each run's `results/` into the
   `docs/results_folder_structure_spec.md` layout
   (`results/raw/observables/observables.csv` etc.).
3. Drop `jellium_spectra.py` once the unified pipeline can be invoked
   directly via `coronene_postprocess.py run --results <jellium_run>/results`.

---

## Choice of drift-removal method for spectra: open question

The current postprocess builds three spectrum variants per quantity
(see `docs/observables_reference.md` §11):

| Variant | Subtraction |
|---|---|
| **raw_subtracted** | `s − s(0)` (initial value) |
| **mean_subtracted** | `s − ⟨s⟩` (mean) |
| **detrended** | `s − (linear fit)` |

**The open question is whether mean-subtraction or initial-value
subtraction is the more physically appropriate "DC removal" for a WP
scattering run.**

Arguments for **mean-subtraction (`s − ⟨s⟩`)**:

* Removes DC exactly. For a stationary or quasi-stationary signal the
  Fourier representation has no DC component, so removing the mean
  before windowing avoids a low-frequency artifact at ω = 0 that would
  leak into the first few non-zero bins.
* Standard in optical-response TDDFT post-processing.
* Symmetric in time: doesn't privilege t = 0.

Arguments for **initial-value subtraction (`s − s(0)`)**:

* Physically meaningful for an *induced* response: at t = 0 the system
  is in the ground state plus an idealised WP (or at the kick), so
  `s(t) − s(0)` is exactly the dynamical response to the perturbation.
* Equivalent to enforcing `(s − s(0))|t=0 = 0`, which is the natural
  initial condition for an induced dipole or current.
* For a current that should integrate to a small dipole change over
  the run, `J_z − J_z(0)` is more interpretable than `J_z − ⟨J_z⟩`
  (the mean is a windowed-time average, not a physically meaningful
  baseline).

For the **WP scattering runs specifically**:

* The WP injection at t = 0 is the stimulus; everything after is the
  response. So initial-value subtraction matches the linear-response
  framework most directly.
* The mean over the propagation window mixes the induced signal with
  the stimulus's own contribution, which biases the baseline.

**Recommendation (to be confirmed)**: prefer **initial-value
subtraction** (the `raw_subtracted` variant) as the default for
*current* and *dipole* spectra; keep mean_subtracted available for
diagnostic comparison; keep linearly-detrended as the variant most
robust to packet drift. For *energy*, neither is obviously better
because total-energy drift is a numerical-conservation artefact,
not a physical response — `mean_subtracted` is fine there.

**Decision pending** — once chosen, the spectrum block in
`inq-stack/python/inqview/postprocess/observables.py::_extended_spectra`
should mark the chosen variant as the canonical one (e.g. by emitting
it under a `spectrum_<col>.png` alias alongside the variant grid).

---

## Ground-state KS eigenvalues + occupations

The HOMO/LUMO gap, the WP slot's energy, and the seed for the overlap
analysis are all encoded in the GS KS eigenvalues `ε_i` and their
fractional occupations `f_i`. These are *not* currently emitted by any
run; the SCF log carries them but no post-run consumer can read them
back without parsing log text.

**Action — future runs**:

* Add `coronene::run_template::dump_eigenvalues(electrons, out_dir)`
  to `shared/cpp/run_template.hpp`, called once right after
  `electrons.load(<gs_checkpoint>)` and before the WP injection.
* Output CSVs at:
  ```
  results/raw/observables/eigenvalues/
    eigenvalues.csv      # state_index, eigenvalue_ha, eigenvalue_ev
    occupations.csv      # state_index, occupation
  ```
* Postprocess phase `observables` reads these and produces, in
  `results/analysis/observables/eigenvalues/`:
  - `eigenvalues_levels.png` (horizontal level diagram, HOMO/LUMO + WP slot marked)
  - `eigenvalues_dos.png` (Gaussian-broadened DOS, σ default 0.1 eV)
  - `eigenvalue_table.txt`

**Action — already-completed runs**: retrofit via
`scripts/extract_eigenvalues_from_log.py` which parses the SCF log
under each `save_gs/<sig>/run.log` (eigenvalues are printed per state
at SCF convergence) and copies the CSVs into every run that loaded the
matching checkpoint.

Documented in `docs/observables_reference.md` §12.

---

## Other useful additions to consider

Three things that would improve interpretability or reproducibility but
are not yet in the framework:

1. **HOMO/LUMO orbital character on screens** — for each transmission /
   backscattering screen, project the diffracted density onto the
   GS orbitals (HOMO, HOMO−1, LUMO) and plot which orbital character
   dominates the scattered signal. Tells us whether the WP excited a
   π-resonance, σ*, etc., rather than a generic continuum response.
   Cheap to add given the eigenvalues + the existing per-step
   `wp_overlap_with_gs_orbitals.csv`.

2. **Norm-conservation diagnostic** — write `norm_per_state.csv` at
   every WRITE_EVERY step (`Σ_r |ψ_i(r,t)|² · dV`) for every
   evolved orbital including the WP. Lets us track WP norm over the
   propagation directly (currently only `injection_report.txt` carries
   the t=0 norm). Also catches numerical leakage if any orbital's norm
   drifts beyond ~1e-3.

3. **Energy-component decomposition spectrum** — same Hann + zero-pad
   pipeline as §11 but applied to `energy_kinetic`, `energy_hartree`,
   `energy_xc` separately, not just `energy_total`. The KE spectrum in
   particular shows whether the post-collision excitation is dominated
   by single-particle motion or by collective modes. Trivial to add
   to the existing extended-spectra columns tuple.

All three are deferred until the GS-eigenvalue retrofit + post-Branch-3
verification cycle completes.
