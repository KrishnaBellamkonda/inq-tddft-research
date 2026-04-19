# INQ Source Code Map

A reference for understanding, extending, and modifying the INQ codebase.
All source files are under `inq/src/` and are **header-only** (`.hpp`).

---

## 1. Architecture Overview

INQ is a **header-only C++17 library** for GPU-accelerated DFT and real-time TDDFT.

Key design decisions:
- **Header-only**: all logic lives in `.hpp` files; the user's `.cpp` is the only
  translation unit that matters. This simplifies distribution but means compile times
  are proportional to what you include.
- **Template-heavy**: Hamiltonians, propagators, and solvers are templated on
  perturbation type, potential type, and basis type. Zero-overhead abstraction —
  the compiler specialises away unused paths at compile time.
- **GPU abstraction layer**: GPU kernels are written once using `GPU_LAMBDA` macros
  and `gpu::run()`. The same code compiles for CUDA, HIP (ROCm), or CPU depending on
  CMake flags.
- **MPI-parallel**: 3D Cartesian communicator decomposition across spatial domains,
  electronic states, and k-points. Communication is explicit and managed by the
  `parallel/` module.
- **Units system**: physical quantities carry unit information via `quantity<magnitude::X>`
  template wrappers. Internal storage is always atomic units.

### Directory inventory (27 subdirectories, ~178 header files)

```
inq/src/
├── inq/             ← master include + quantity wrapper
├── systems/         ← top-level data containers (ions, electrons, cell)
├── states/          ← Kohn-Sham state configuration
├── basis/           ← spatial representations (real-space, reciprocal, fields)
├── hamiltonian/     ← Kohn-Sham Hamiltonian and XC
├── ionic/           ← ion-ion interactions, Brillouin zone, ion dynamics
├── ground_state/    ← SCF driver and initial guess
├── real_time/       ← TDDFT propagation
├── eigensolvers/    ← orbital optimisation algorithms
├── mixers/          ← SCF density mixing
├── solvers/         ← linear algebra solvers (Poisson, CG, etc.)
├── operations/      ← field operators (overlap, gradient, FFT, I/O)
├── observables/     ← physical observables (density, forces, dipole)
├── perturbations/   ← external field definitions (kick, laser, etc.)
├── parallel/        ← MPI distribution and communication
├── interface/       ← CLI command parsing (shell interface)
├── options/         ← algorithm configuration (C++ API)
├── input/           ← MPI environment, parallelization config, k-points
├── math/            ← 3D vectors, complex, FFT size finding
├── matrix/          ← dense linear algebra (diagonalise, Cholesky, etc.)
├── magnitude/       ← physical units (energy, length, time, field)
├── physics/         ← fundamental constants
├── utils/           ← serialisation, string helpers, profiling
├── parse/           ← structure file parsers (CIF, POSCAR, XYZ)
├── config/          ← runtime path configuration
├── bomd/            ← Born-Oppenheimer MD
└── main/            ← CLI entry point + test runner main
```

---

## 2. Core System Containers

These three classes are the primary data containers. Almost every function takes
`ions` and/or `electrons` as arguments.

### `systems/cell.hpp` — Simulation Cell

Stores lattice vectors and derived quantities.

```cpp
// Construction
systems::cell::cubic(a)
systems::cell::orthorhombic(a, b, c)
systems::cell::lattice(v0, v1, v2)

// Periodicity
.periodic()        // 3D PBC (default for solids)
.finite()          // no PBC (molecules)
.periodicity(n)    // 0, 1, 2, or 3

// Key members
cell.lattice()      // 3×3 lattice matrix (rows = vectors), in bohr
cell.reciprocal()   // 3×3 reciprocal lattice (rows = G vectors), in 1/bohr
cell.volume()       // unit cell volume (bohr^3)
cell.periodicity()  // 0–3
```

### `systems/ions.hpp` — Ionic System (~1200 lines)

Manages atomic positions, velocities, species, and structure I/O.

```cpp
// Construction
systems::ions ions(cell);
systems::ions ions(cell, filename);  // read from CIF/POSCAR/XYZ

// Inserting atoms
ions.insert("Li", {x_b, y_b, z_b});           // Cartesian (bohr)
ions.insert_fractional("Li", {fx, fy, fz});    // fractional coords

// Access
ions.size()               // number of atoms
ions.positions()          // array of vector3<double> (bohr)
ions.velocities()         // array of vector3<double> (bohr/a.u.)
ions.species(i)           // chemical symbol string
ions.cell()               // reference to cell

// Supported file formats: CIF, POSCAR (VASP), XYZ
```

### `systems/electrons.hpp` — Electronic System (~1100 lines)

The central container for all electronic state. Created from `ions` and options.

```cpp
// Construction
systems::electrons electrons(
    ions,
    options::electrons{}.cutoff(30.0_Ry),
    kpoints   // optional; default = Gamma
);

// Key internal members (not usually accessed directly)
electrons.kpin_           // k-point orbital sets (the wavefunctions)
electrons.eigenvalues_    // KS eigenvalues per k-point per state
electrons.occupations_    // occupation numbers
electrons.atomic_pot_     // atomic (ionic) potential
electrons.states_         // ks_states config (spin, extra_states, etc.)

// Save/load (for TDDFT restart)
electrons.save("results/gs_save")
electrons.load(ions, "results/gs_save")
```

---

## 3. States and Basis

### `states/ks_states.hpp` — Kohn-Sham State Configuration

Stores spin config, number of states, temperature, number of electrons.

```cpp
// Spin configurations (states/spin_config.hpp)
states::spin_config::UNPOLARIZED   // 1 spin channel
states::spin_config::POLARIZED     // collinear (LSDA/GGA)
states::spin_config::NON_COLLINEAR // 4-component spinor
```

### `basis/real_space.hpp` — Real-Space Grid Basis

Discretises the simulation cell onto a uniform 3D grid.
Grid dimensions chosen to have only small prime factors (2, 3, 5, 7) for FFT
efficiency.

### `basis/field.hpp` — Single Scalar/Complex Field

```cpp
basis::field<Basis, Type>    // e.g. field<real_space, double>
```
Stores a linear array associated with a basis. Provides type-safe access and
skeleton wrappers.

### `basis/field_set.hpp` — Matrix of Fields

```cpp
basis::field_set<Basis, Type>  // e.g. field_set<real_space, complex>
```
Used for orbital sets (multiple states on the same grid). Distributed across
MPI ranks along the state dimension.

---

## 4. Hamiltonian and Self-Consistency

### `hamiltonian/ks_hamiltonian.hpp` — Kohn-Sham Operator (~800 lines)

The Hamiltonian operator `H = T + V_local + V_nl + V_xc + V_Hartree`.

```cpp
// Templated on potential type (double for GGA, complex for current-dependent)
hamiltonian::ks_hamiltonian<PotentialType> H(basis, states, atomic_pot, ions);

// Apply H to orbital set
H(phi)           // returns H|phi>

// Update projectors (call after ions move)
H.update_projectors(ions, atomic_pot)
```

Key members:
- `exchange_`: exact exchange operator (for hybrid functionals)
- `projectors_`: non-local pseudopotential projectors
- `relativistic_projectors_`: spin-orbit coupling

### `hamiltonian/self_consistency.hpp` — SCF Controller (~800 lines)

Manages XC potential, Hartree potential, ionic potential, and perturbation.
Updated each SCF iteration.

```cpp
hamiltonian::self_consistency sc(theory, electrons, perturbation);
sc.update_hamiltonian(H, energy, electrons, time);
```

### `hamiltonian/xc_term.hpp` — Exchange-Correlation (~1000 lines)

Calls libxc to compute XC energy and potential. Handles LDA, GGA, MetaGGA, and
hybrid (ACE approximation) functionals. Automatically computes density gradients,
kinetic energy density, or Laplacian as needed by the functional.

### `hamiltonian/atomic_potential.hpp` — Pseudopotentials (~800 lines)

Loads pseudopotentials via the `pseudopod` library. Manages non-linear core
corrections (NLCC) and double-grid interpolation for accuracy.

### `hamiltonian/energy.hpp` — Total Energy (~500 lines)

Container for all energy components: kinetic, Hartree, XC, ionic, exchange,
external. Computes the total DFT energy.

### `hamiltonian/paw.hpp` — PAW Method (~900 lines)

Projector Augmented Wave corrections for improved accuracy near nuclei.

---

## 5. Ionic Module

### `ionic/interaction.hpp` — Ion-Ion Interaction (~1200 lines)

Ewald summation for ion-ion Coulomb energy and forces. Handles both finite
and periodic systems. Also computes stress.

### `ionic/brillouin.hpp` — Brillouin Zone (~500 lines)

K-point management: Monkhorst-Pack grids, k-point weights, symmetry reduction
via spglib.

### `ionic/propagator.hpp` — Classical Ion Dynamics

Integrators for ionic motion: velocity Verlet (Born-Oppenheimer MD), fixed ions
(STATIC), or constant-velocity impulsive propagation (IMPULSIVE).

---

## 6. Solvers

### `ground_state/calculator.hpp` — SCF Driver (~700 lines)

Orchestrates the self-consistent field loop. Template parameter: perturbation type.

```
ground_state::initial_guess(ions, electrons)
ground_state::calculate(ions, electrons, theory, options)
  ↓
  for each SCF iteration:
    eigensolvers::steepest_descent(H, preconditioner, phi)
    observables::density::calculate(electrons)
    solvers::poisson(density) → V_Hartree
    hamiltonian::xc_term(density) → V_xc
    self_consistency::update_hamiltonian(...)
    mixers::broyden(density_in, density_out)
  ↓
  ground_state::results  (energies, eigenvalues, forces)
```

### `real_time/propagate.hpp` — TDDFT Propagation (~500 lines)

Drives the real-time loop. Template parameters: process function, ion sub-propagator,
perturbation type. Perturbation is passed as the 6th argument (optional):

```cpp
// Signature:
template<typename ProcessFunction,
         typename IonSubPropagator = ionic::propagator::fixed,
         typename Perturbation = perturbations::none>
void propagate(ions, electrons, func, theory, opts,
               Perturbation const & pert = {},
               int start_step = 0);
```

Time loop:
```
for each step:
  perturbation applied at time t
  self_consistency::update_hamiltonian(H, energy, electrons, t)
  crank_nicolson OR etrs (propagate phi → phi(t+dt))
  observables (dipole, current, density)
  ionic::propagator::step(ions, forces)   # if Ehrenfest
  process_function(viewables)
```

### `eigensolvers/steepest_descent.hpp` — Orbital Optimiser (~400 lines)

Preconditioned steepest descent for minimising the orbital energy functional.
Maintains orthogonality. Template on `(Hamiltonian, Metric, Preconditioner, Orbitals)`.

### `mixers/broyden.hpp` and `mixers/linear.hpp` — Density Mixers

Broyden (Pulay) mixing accelerates SCF convergence. Linear mixing as fallback.

### `solvers/poisson.hpp` — Hartree Potential (~600 lines)

Poisson equation solver for the Hartree (electrostatic) potential. Conjugate
gradient with various preconditioners. Supports periodic and finite systems.

---

## 7. Operations

`operations/` contains field-level operators — the building blocks for everything above.

| File | Operation |
|---|---|
| `overlap.hpp` | Overlap matrix ⟨φᵢ\|φⱼ⟩ |
| `overlap_diagonal.hpp` | Diagonal elements only |
| `orthogonalize.hpp` | Gram-Schmidt / Cholesky orthogonalisation |
| `laplacian.hpp` | −∇² operator (kinetic energy) in reciprocal space |
| `gradient.hpp` | ∇ operator |
| `divergence.hpp` | ∇· operator |
| `transform.hpp` | Real-space ↔ reciprocal-space (FFT wrapper) |
| `rotate.hpp` | Unitary rotation of orbital set |
| `exponential.hpp` | Matrix exponential (Crank-Nicolson step) |
| `integral.hpp` | Spatial integrals with MPI reduction |
| `transfer.hpp` | Transfer fields between bases or communicators |
| `io.hpp` | Save/load field sets to disk |
| `randomize.hpp` | Random orbital initialisation |
| `preconditioner.hpp` | Diagonal kinetic energy preconditioner |

---

## 8. Observables

```cpp
// Electron density from orbitals
observables::density::calculate(electrons)
// → field<real_space, double>  (or 2-component for spin)

// Forces and stress on ions
observables::forces_stress forces(ions, electrons, hamiltonian, energy)
forces.forces()    // std::vector<vector3<double>>  (Ha/bohr)
forces.stress()    // 3×3 tensor

// Dipole moment
observables::dipole(electrons)
// → vector3<double>  (a.u.)

// Electronic current density
observables::current(electrons, hamiltonian)

// Magnetisation (spin systems)
observables::magnetization(electrons)
```

---

## 9. Perturbations

All perturbations inherit from `perturbations::none` (base class — zero perturbation).
They are passed as the optional 6th argument to `real_time::propagate`.

| Class | Description |
|---|---|
| `perturbations::none` | No perturbation (default) |
| `perturbations::kick` | Instantaneous momentum kick at t=0 for optical spectrum |
| `perturbations::laser` | Monochromatic laser (sin field or cos vector potential) |
| `perturbations::simple_electric_field` | Static uniform electric field |
| `perturbations::magnetic` | Uniform magnetic field (Zeeman) |
| `perturbations::magnetic_pulse` | Pulsed magnetic field |
| `perturbations::absorbing` | Complex absorbing potential (CAP, open boundary) |
| `perturbations::sum` | Superpose multiple perturbations |
| `perturbations::blend` | Smooth on/off switching between perturbations |

```cpp
// Optical absorption spectrum kick (length gauge for finite, velocity for periodic):
perturbations::kick kick(ions.cell(), {0.01, 0.0, 0.0});

// Laser (monochromatic, velocity gauge by default):
perturbations::laser laser({0.0, 0.0, 1.0}, 0.1_Ha);

// Usage in propagation:
real_time::propagate(ions, electrons, output_fn, theory, opts, kick);
real_time::propagate(ions, electrons, output_fn, theory, opts, laser);
```

---

## 10. Parallel Infrastructure

INQ uses a **3D Cartesian MPI communicator** decomposition:

```
world communicator
    ↓ split by k-points
kpoint_communicator
    ↓ split by states
states_communicator
    ↓ split by spatial domains
domain_communicator
```

Key files:

| File | Role |
|---|---|
| `parallel/partition.hpp` | 1D partition of an index set across ranks |
| `parallel/communicator.hpp` | Thin MPI communicator wrapper |
| `parallel/transpose.hpp` | Tensor transpose between basis↔state distribution |
| `parallel/gather.hpp` | Gather distributed data to one rank |
| `parallel/alltoall.hpp` | All-to-all redistribution |
| `parallel/get_remote_points.hpp` | Fetch grid points from remote ranks |

Parallel distribution is managed automatically by `systems::electrons` during
construction based on `options::electrons` and `input::parallelization`.

---

## 11. GPU Abstraction Layer

GPU code is written once and compiles to CUDA, HIP, or CPU depending on CMake flags.

### Key macros and functions

```cpp
// Mark a lambda as a GPU kernel
GPU_LAMBDA   // expands to __device__ (CUDA), or empty (CPU)

// Launch a kernel over a 3D index space
gpu::run(Nz, Ny, Nx, GPU_LAMBDA (auto iz, auto iy, auto ix) {
    // per-element kernel code
});

// 1D version
gpu::run(N, GPU_LAMBDA (auto i) { ... });

// GPU-resident array (allocated in device memory)
gpu::array<Type, Rank>  // from boost-multi with GPU allocator

// Raw pointer for CUDA API calls
raw_pointer_cast(array.begin())
```

### Where GPU code lives

- `operations/transform.hpp` — FFT via cuFFT/FFTW
- `operations/laplacian.hpp`, `gradient.hpp` — reciprocal-space operators
- `operations/overlap*.hpp` — cuBLAS GEMM for overlap matrices
- `operations/rotate.hpp` — orbital rotations via BLAS
- `hamiltonian/projector.hpp` — non-local projector application
- `observables/density.hpp` — density construction from orbitals
- User `.cpp` files can also use `gpu::run()` directly (see `ResearchProject/`)

---

## 12. User Interface Layers

### C++ API (typical usage)

```cpp
#include <inq/inq.hpp>
using namespace inq;
using namespace inq::magnitude;

// options:: namespace — configure algorithms
options::electrons{}   .cutoff(30_Ry).temperature(0.1_eV)
options::theory{}      .pbe()
options::ground_state{}.energy_tolerance(1e-6_Ha).max_steps(100)
options::real_time{}   .num_steps(500).dt(0.04_atomictime).etrs()
```

### Shell CLI (`interface/`)

The compiled `inq` binary exposes a POSIX-style shell interface. Commands map to
the namespaces above:

```bash
inq cell cubic 10 bohr
inq electrons cutoff 30 Ry
inq ions insert Li 0 0 0
inq ground-state run
inq real-time run
inq results ground-state total-energy
```

Interface files in `interface/` parse these commands and call into the C++ API.
The CLI and C++ API are equivalent in capability.

---

## 13. Primary Data Flow

### Ground State DFT

```
systems::ions + systems::electrons
    ↓
ground_state::initial_guess(ions, electrons)
    ↓
ground_state::calculate(ions, electrons, theory, gs_options)
    ↓
  hamiltonian::self_consistency sc(theory, electrons, perturbation=none)
    ↓
  for each SCF iteration:
    eigensolvers::steepest_descent(H, preconditioner, phi)
      ← operations::laplacian (kinetic energy, in reciprocal space)
      ← hamiltonian::projector_all (non-local PP)
      ← operations::orthogonalize (maintain ⟨φᵢ|φⱼ⟩=δᵢⱼ)
    ↓
    observables::density::calculate(electrons)
      ← operations::transform (FFT φ → real space)
    ↓
    solvers::poisson(density) → V_Hartree
    hamiltonian::xc_term(density) → V_xc (via libxc)
    sc.update_hamiltonian(H, energy, electrons)
    ↓
    mixers::broyden(density_in, density_out)
    (check convergence: |E_new - E_old| < tolerance)
    ↓
  ground_state::results{energy, eigenvalues, forces}
```

### Real-Time TDDFT

```
systems::ions + systems::electrons  (from ground state)
    ↓
real_time::propagate(ions, electrons, func, theory, rt_options, perturbation)
    ↓
  apply perturbation at t=0 (kick: phase multiplied into φ; laser: vector potential)
    ↓
  for each time step t:
    sc.update_hamiltonian(H, energy, electrons, t)  // update V_xc[ρ(t)], V_laser(t)
    ↓
    if etrs:   real_time::etrs::step(H, phi, dt)    // U(t+dt,t)|φ⟩
    if cn:     real_time::crank_nicolson::step(...)  // (1 - iH·dt/2)/(1 + iH·dt/2)
    ↓
    observables::density::calculate(electrons)      // update ρ(t+dt)
    observables::dipole, current (if requested)
    ↓
    if Ehrenfest: observables::forces_stress → ionic::propagator::step
    ↓
    func(viewables)  // user callback: save data, print progress
```

---

## 14. Extension Points

### Adding a New XC Functional

1. Check if libxc already supports it (`libxc_funcs_list` or `xc-info` tool).
2. If yes: pass the libxc functional ID via `options::theory::set_libxc(id_x, id_c)`.
3. If adding a custom functional: implement in `hamiltonian/xc_term.hpp`, following
   the pattern of existing terms. Declare requirements (gradient, tau, laplacian) by
   setting the corresponding flags.

### Adding a New Perturbation

1. Create `perturbations/myperts.hpp`, inherit from `perturbations::none`.
2. Implement the relevant virtual methods:
   - `uniform_electric_field(double time)` → `vector3<double>`
   - `uniform_vector_potential(double time)` → `vector3<double>`
   - `zero_step(PhiType& phi)` → apply kick at t=0 to orbitals
3. Pass an instance as the 6th argument to `real_time::propagate`.
4. No template specialisation needed — the propagator is already templated.

### Adding a New Eigensolver

1. Create `eigensolvers/mysolver.hpp`.
2. Function signature: `void mysolver(H, metric, preconditioner, phi)` where
   `phi` is a `field_set` that gets updated in-place.
3. Must maintain orthogonality of `phi` on exit.
4. Wire it up in `ground_state/calculator.hpp` where eigensolvers are dispatched.

### Adding a New Observable

1. Create `observables/myobs.hpp`.
2. Use `operations::integral()` for spatial integrals (handles MPI reduction).
3. Use `operations::transform()` to go between real and reciprocal space.
4. Expose via the `real_time::viewables` struct if needed in the TDDFT callback.

### Modifying the Parallel Decomposition

1. Edit `parallel/partition.hpp` for custom 1D partitioning rules.
2. Edit `input/parallelization.hpp` to expose new options to the user.
3. Update `systems/electrons.hpp` constructor to use the new communicator layout.
4. Ensure `parallel/transpose.hpp` covers any new redistribution needed.

---

## 15. External Dependencies

| Library | Role | Notes |
|---|---|---|
| **boost-multi** | Multidimensional arrays | GPU-aware; basis for `field` and `field_set` |
| **boost-mpi3** | MPI C++ wrappers | Wraps MPI_Comm etc. |
| **pseudopod** | Pseudopotential parsing | Reads UPF format; provides `atomic_potential` data |
| **libxc** | XC functionals | Supports 600+ LDA/GGA/MetaGGA/hybrid functionals |
| **FFTW** | CPU FFT | Used when CUDA/HIP not enabled |
| **cuFFT / rocFFT** | GPU FFT | Linked when ENABLE_CUDA/HIP |
| **cuBLAS / rocBLAS** | GPU BLAS | Dense linear algebra on GPU |
| **BLAS / LAPACK** | CPU BLAS | Host-side linear algebra |
| **spglib** | Crystal symmetries | K-point reduction and Brillouin zone symmetry |
| **spdlog** | Logging | Structured JSON-compatible log output |
| **catch2** | Unit testing | Test framework for `src/` unit tests |
| **pybind11** | Python bindings | Optional Python interface |

All dependencies except BLAS/LAPACK/MPI/CUDA are fetched via CMake FetchContent.

---

## 16. File Size Reference

The largest files (by lines) give a sense of where complexity concentrates:

| File | ~Lines | Content |
|---|---|---|
| `systems/ions.hpp` | 1200 | Ion management, structure I/O, species |
| `ionic/interaction.hpp` | 1200 | Ewald summation, forces, stress |
| `systems/electrons.hpp` | 1100 | Electronic system container, MPI setup |
| `hamiltonian/xc_term.hpp` | 1000 | XC functional evaluation via libxc |
| `hamiltonian/projector_all.hpp` | 800 | Non-local pseudopotential projectors |
| `hamiltonian/ks_hamiltonian.hpp` | 800 | Kohn-Sham operator |
| `hamiltonian/atomic_potential.hpp` | 800 | Pseudopotential loading |
| `operations/transfer.hpp` | 800 | Basis/communicator transfer |
| `operations/transform.hpp` | 800 | FFT wrappers |
| `basis/field_set.hpp` | 600 | Distributed orbital storage |
| `solvers/poisson.hpp` | 600 | Hartree potential solver |
| `hamiltonian/paw.hpp` | 900 | PAW corrections |
