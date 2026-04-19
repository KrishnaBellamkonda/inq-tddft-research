# INQ Usage Guide

INQ is a GPU-accelerated DFT/TDDFT electronic structure engine written in C++17.
It is a header-only library (`#include <inq/inq.hpp>`); no separate API initialisation is required.

---

## 1. Minimal Structure of an INQ Program

```cpp
#include <inq/inq.hpp>

int main(){
    using namespace inq;
    using namespace inq::magnitude;   // enables _eV, _Ha, _Ry, _angstrom, _b, _atomictime

    // 1. Define the unit cell + ion positions
    systems::ions ions( systems::cell::cubic(5.0_angstrom).periodic() );
    ions.insert("Si", {0.0_b, 0.0_b, 0.0_b});

    // 2. Define electrons
    systems::electrons electrons(ions, options::electrons{}.cutoff(30.0_Ry));

    // 3. Initial guess
    ground_state::initial_guess(ions, electrons);

    // 4. Ground-state SCF
    auto gs = ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}.energy_tolerance(1e-6_Ha));

    // 5. (Optional) Real-time TDDFT
    real_time::propagate(ions, electrons, [](auto){},
        options::theory{}.lda(),
        options::real_time{}.num_steps(500).dt(0.05_atomictime));
}
```

---

## 2. Unit Cell — `systems::cell`

```cpp
// Cubic cell (single lattice parameter)
systems::cell::cubic(a)

// Orthorhombic cell (three distinct lengths)
systems::cell::orthorhombic(a, b, c)

// Fully general (three lattice vectors)
systems::cell::lattice(
    vector3<quantity<magnitude::length>>{a0x, a0y, a0z},
    vector3<quantity<magnitude::length>>{a1x, a1y, a1z},
    vector3<quantity<magnitude::length>>{a2x, a2y, a2z})

// Periodicity modifiers (chain after cell constructor)
.periodic()      // 3D periodic (default for solids)
.finite()        // isolated molecule (no PBC)
.periodicity(2)  // 2D slab
```

Units recognised by `quantity<magnitude::length>`: `_angstrom`, `_b` (bohr).

---

## 3. Ions — `systems::ions`

```cpp
systems::ions ions(cell);

// Insert atom at Cartesian position
ions.insert("Li", {x_b, y_b, z_b});          // positions in bohr

// Insert atom at fractional (crystal) coordinates
ions.insert_fractional("Li", {0.5, 0.5, 0.5});

// Access
ions.size();                  // number of atoms
ions.positions()[i];          // cartesian position (bohr)
ions.velocities()[i];         // velocity (bohr / a.u. time)

// Set velocity (ionic kick)
ions.velocities()[i] = vector3<double>{vx, vy, vz};
```

### Coordinate origin — the (0,0,0) corner rule

**The origin `(0,0,0)` is always the corner of the simulation box.** The grid runs
from `(0,0,0)` to `(Lx, Ly, Lz)`. All atom Cartesian coordinates must be positive
and within `[0, L]` in each direction.

The fractional API makes this explicit: `insert_fractional({0.0, 0.0, 0.0})` is the
corner, `insert_fractional({0.5, 0.5, 0.5})` is the centre.

**Centring a molecule correctly:**

```cpp
// WRONG — places H outside the box (z < 0); INQ wraps it to the opposite end
auto cell = systems::cell::cubic(8.0_bohr).finite();
ions.insert("H", {0.0_b, 0.0_b, -0.868_b});   // z < 0: outside [0, 8]
ions.insert("F", {0.0_b, 0.0_b, +0.868_b});   // near z=0 corner

// CORRECT — place both atoms at the geometric centre of the box
auto L = 16.0_bohr;
auto cell = systems::cell::cubic(L).finite();
auto cx = L / 2;   // 8 bohr
ions.insert("H", {cx, cx, cx - 0.868_b});
ions.insert("F", {cx, cx, cx + 0.868_b});
```

**What happens with negative or out-of-range coordinates:**
INQ wraps coordinates using the cell periodicity vectors, even for `.finite()` cells.
An atom at `z = -d` lands at `z = L - d` (opposite end of the box). This is
physically wrong for a finite (non-periodic) system — the molecule is split across
the box — and for LDA/GGA calculations the bond never forms, producing a density
concentrated at the box corners rather than at the molecular centre.

**The `Tutorial/HF/HF.cpp` anomaly:** that tutorial uses atoms at `(0, 0, ±0.459 Å)`
in a `cubic(8 bohr)` cell and appears to work, but only because it uses
`non_interacting` theory (no Coulomb bond) and the charge analysis happens to give
the right sign due to how `rvector_cartesian` returns centred coordinates. For any
interacting theory (LDA, PBE) those coordinates produce an unphysical molecule.

**Rule of thumb:** always place the molecule at `(Lx/2, Ly/2, Lz/2)` and verify
all atom coordinates are positive. The coronene geometry (`coronene_centered.xyz`)
is a correct reference example.

---

## 4. Electrons — `systems::electrons`

```cpp
systems::electrons electrons(
    ions,
    options::electrons{} /* see below */,
    kpoint_sampler          /* optional, default = Gamma */
);
```

### `options::electrons` methods

| Method | Description |
|---|---|
| `.cutoff(E)` | Plane-wave kinetic energy cutoff (`_Ry`, `_Ha`, `_eV`) |
| `.spacing(d)` | Real-space grid spacing (`_b`, `_angstrom`); alternative to cutoff |
| `.extra_states(n)` | Extra unoccupied KS states (needed for metals) |
| `.extra_electrons(n)` | Shift total electron count by n (charge) |
| `.temperature(kT)` | Fermi smearing width in energy units; required for metals |
| `.spin_unpolarized()` | Default; one spin channel |
| `.spin_polarized()` | Collinear spin (LSDA/GGA) |
| `.spin_non_collinear()` | Non-collinear spin |
| `.double_grid()` | Double grid for better XC accuracy |
| `.density_factor(f)` | Density grid factor (default 2.0) |

### k-point sampling — `input::kpoints`

```cpp
input::kpoints::gamma()                   // Gamma only
input::kpoints::grid({nx, ny, nz})        // Monkhorst-Pack, unshifted
input::kpoints::grid({nx, ny, nz}, true)  // Monkhorst-Pack, shifted (recommended)
input::kpoints::list()                    // manual list; use .add(k, weight)
```

---

## 5. Theory — `options::theory`

Controls the exchange-correlation functional:

| Method | Functional |
|---|---|
| `.non_interacting()` | No XC (free electrons) |
| `.lda()` | LDA (PZ81 parametrisation) |
| `.pbe()` | GGA-PBE (recommended for most calculations) |
| `.pbe0()` | Hybrid PBE0 |
| `.b3lyp()` | Hybrid B3LYP |
| `.scan()` | Meta-GGA SCAN |
| `.r2scan()` | Meta-GGA r²SCAN |
| `.hartree_fock()` | Pure HF exchange |
| `.hartree()` | Hartree only (no exchange or correlation) |
| `.induced_vector_potential(alpha)` | Add induced vector potential (TDDFT current response) |

---

## 6. Ground State — `ground_state::calculate`

```cpp
// Step 1: provide an initial guess for the density/orbitals
ground_state::initial_guess(ions, electrons);

// Step 2: run SCF
auto gs = ground_state::calculate(ions, electrons, theory, gs_options);

// Key results
gs.energy.total()     // Total DFT energy (Ha)
gs.energy.kinetic()   // Electronic kinetic energy
gs.energy.xc()        // Exchange-correlation energy
gs.energy.hartree()   // Hartree (Coulomb) energy
gs.energy.ion()       // Ion-ion interaction energy
gs.total_iter         // Number of SCF iterations
```

### `options::ground_state` methods

| Method | Description |
|---|---|
| `.energy_tolerance(dE)` | SCF convergence threshold on total energy (default 1e-5 Ha) |
| `.mixing(f)` | Linear mixing factor (default 0.3) |
| `.broyden_mixing()` | Use Broyden (Pulay) mixing — default |
| `.linear_mixing()` | Simple linear mixing |
| `.mixing_ndim(n)` | History depth for Broyden mixer |
| `.steepest_descent()` | Steepest-descent eigensolver |
| `.max_steps(n)` | Maximum SCF iterations |
| `.calculate_forces()` | Also compute atomic forces |
| `.silent()` | Suppress per-iteration output |

### Save / Load ground state

```cpp
electrons.save("path/to/dir");           // write orbitals + density to disk
electrons.load(ions, "path/to/dir");     // restore for a TDDFT restart
```

---

## 7. Real-Time TDDFT — `real_time::propagate`

```cpp
// Without perturbation (ions have been kicked via velocities, for example):
real_time::propagate(
    ions, electrons,
    output_function,          // callable: (auto data) -> void
    options::theory{}.pbe(),
    options::real_time{} /* see below */
);

// With an explicit perturbation (6th argument):
real_time::propagate(
    ions, electrons,
    output_function,
    options::theory{}.pbe(),
    options::real_time{}.num_steps(1000).dt(0.04_atomictime),
    perturbation_object          // e.g. perturbations::kick or perturbations::laser
);
```

### `options::real_time` methods

| Method | Description |
|---|---|
| `.dt(time)` | Time step (e.g., `0.04_atomictime`) |
| `.num_steps(n)` | Total number of time steps |
| `.propagation_time(t)` | Alternative to num_steps (total time) |
| `.etrs()` | ETRS propagator — default, efficient for metals |
| `.crank_nicolson()` | Crank-Nicolson propagator |
| `.parallel_transport()` | Parallel transport propagator |
| `.static_ions()` | Fix ions (no ionic dynamics) — default if no other set |
| `.impulsive()` | Ions move at constant velocity (no forces) |
| `.ehrenfest()` | Full Ehrenfest dynamics (forces computed) |
| `.observables_dipole()` | Compute dipole moment at each step |
| `.observables_current()` | Compute electronic current density |
| `.observables_clear()` | Remove all observables |

### Output callback

The second argument to `real_time::propagate` is a callable invoked at every time step:

```cpp
auto output = [&](auto data){
    if(data.every(10)){               // execute every 10th step
        data.iter();                  // current step number
        data.time();                  // current time (a.u.)
        data.energy().total();        // instantaneous total energy (Ha)
        data.energy().kinetic();      // electronic kinetic energy
        data.current();               // vector3<double>: J_x, J_y, J_z
        data.dipole();                // vector3<double>: dipole moment (a.u.)
        data.positions();             // std::vector<vector3<double>>: atom positions (bohr)
        data.forces();                // std::vector<vector3<double>>: forces (Ha/bohr)
    }
};
```

---

## 8. Perturbations

Perturbations are passed as the optional 6th argument to `real_time::propagate`.
All perturbation classes are in the `inq::perturbations` namespace.

```cpp
#include <inq/inq.hpp>
using namespace inq;
using namespace inq::magnitude;
```

### `perturbations::kick` — Momentum kick (optical spectrum)

Applies a uniform momentum kick at t=0. Automatically uses velocity gauge for
periodic systems and length gauge for finite systems.

```cpp
// Kick in the x direction with strength 0.01 a.u.
perturbations::kick kick(ions.cell(), {0.01, 0.0, 0.0});

// Explicit gauge
perturbations::kick kick(ions.cell(), {0.01, 0.0, 0.0},
                         perturbations::gauge::length);   // length gauge
perturbations::kick kick(ions.cell(), {0.01, 0.0, 0.0},
                         perturbations::gauge::velocity); // velocity gauge
```

### `perturbations::laser` — Monochromatic laser field

Monochromatic CW field. Default gauge: velocity (vector potential). Length gauge
is available but **only for finite systems**.

```cpp
// Linearly polarised along z, frequency 0.1 Ha
perturbations::laser laser({0.0, 0.0, 1.0}, 0.1_Ha);

// Explicit gauge
perturbations::laser laser({0.0, 0.0, 1.0}, 0.1_Ha,
                           perturbations::gauge::length);
```

### `perturbations::simple_electric_field` — Static electric field

```cpp
perturbations::simple_electric_field efield({0.0, 0.0, 0.01});
```

### Usage example: optical absorption spectrum

```cpp
// 1. Ground state
ground_state::initial_guess(ions, electrons);
auto gs = ground_state::calculate(ions, electrons,
    options::theory{}.lda(), options::ground_state{});

// 2. Kick + propagate
electrons.save("results/gs_save");

perturbations::kick kick(ions.cell(), {0.01, 0.0, 0.0});

real_time::propagate(ions, electrons,
    [&](auto data){
        if(data.every(1))
            output << data.time() << " " << data.dipole()[0] << "\n";
    },
    options::theory{}.lda(),
    options::real_time{}.num_steps(2000).dt(0.04_atomictime)
        .observables_dipole(),
    kick);

// 3. Fourier transform dipole signal → absorption spectrum
//    (post-process with Python)
```

---

## 9. Ionic Dynamics Modes

| INQ option | QBall equivalent | Description |
|---|---|---|
| `.static_ions()` | — | Ions frozen for entire TDDFT run |
| `.impulsive()` | `atoms_dyn IMPULSIVE` | Ions move at constant velocity; no forces |
| `.ehrenfest()` | `atoms_dyn EHRENFEST` | Ions follow Ehrenfest forces |

The ionic velocity kick (impulsive run) protocol:
```cpp
// After ground-state SCF, before real_time::propagate:
double v_kick = 0.04;   // bohr / a.u. of time
for(int i = 0; i < ions.size(); i++)
    ions.velocities()[i] = vector3<double>{v_kick, 0.0, 0.0};

real_time::propagate(ions, electrons, output,
    options::theory{}.pbe(),
    options::real_time{}.num_steps(2000).dt(0.04_atomictime)
        .impulsive().observables_current());
```

---

## 10. Magnitude Literals (units)

All available in `namespace inq::magnitude`:

| Literal | Quantity | SI equivalent |
|---|---|---|
| `_Ha` | Energy | Hartree |
| `_eV` | Energy | electron-volt |
| `_Ry` | Energy | Rydberg (= 0.5 Ha) |
| `_b` | Length | bohr |
| `_angstrom` | Length | Ångström |
| `_atomictime` | Time | atomic unit of time (≈ 24.19 as) |

---

## 11. Running on GPU

INQ automatically uses the GPU when compiled with CUDA/HIP. The log will print:

```
process 0 has gpu id <hash>
```

Environment variables required at runtime:
```bash
export INQ_SHARE_PATH=/path/to/inq/install/share
export PSEUDOPOD_SHARE_PATH=/path/to/inq/install/share/pseudopod
```

Use `inq-run` (see `docs/compilation.md`) to handle all of this automatically.

---

## 12. Common Patterns

### Metallic system (requires smearing)
```cpp
systems::electrons electrons(ions,
    options::electrons{}
        .cutoff(40.0_Ry)
        .extra_states(4)
        .temperature(0.086_eV),          // kB × 1000 K
    input::kpoints::grid({4, 4, 4}, true));
```

### Molecule in a box (no PBC)
```cpp
systems::ions ions(systems::cell::cubic(20.0_b).finite());
```

### Slab (2D periodic)
```cpp
systems::ions ions(systems::cell::orthorhombic(a, b, c).periodicity(2));
```

### Spin-polarised calculation
```cpp
systems::electrons electrons(ions,
    options::electrons{}.cutoff(30.0_Ry).spin_polarized());
```

### Restart TDDFT from saved ground state
```cpp
electrons.load(ions, "results/gs_save");
// then call real_time::propagate directly, skipping SCF
```
