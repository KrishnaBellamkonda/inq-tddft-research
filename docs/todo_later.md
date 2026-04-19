# To-Do Later

Deferred tasks for the coronene WP scattering project. Items here are not blockers for the current simulation but should be addressed before paper-quality results are claimed.

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
