# INQ C++ API Reference

This document covers every INQ API call found across your Tutorial and ResearchProject files.
It is written assuming you understand the physics, but not the C++ library's naming conventions or argument formats.

**Units convention used by INQ:**
All internal quantities are in Hartree atomic units (bohr, Hartree, a.u. time).
You pass values using *quantity literals* (e.g. `20.0_bohr`, `40.0_Ry`, `0.04_atomictime`)
and INQ converts internally. The suffix controls the unit — see the table at the end.

---

## Table of Contents

1. [System construction](#1-system-construction)
2. [Inserting atoms](#2-inserting-atoms)
3. [Electron options](#3-electron-options)
4. [Ground state](#4-ground-state)
5. [Ground state results](#5-ground-state-results)
6. [Real-time propagation options](#6-real-time-propagation-options)
7. [The propagation callback (data object)](#7-the-propagation-callback-data-object)
8. [Perturbations](#8-perturbations)
9. [Accessing orbitals and density](#9-accessing-orbitals-and-density)
10. [Unit literals quick reference](#10-unit-literals-quick-reference)

---

## 1. System Construction

### `systems::cell::cubic(L).finite()`

```cpp
// From n2.cpp
systems::ions ions(systems::cell::cubic(20.0_bohr).finite());
```

Creates a cubic box of side `L` with **no periodicity** (isolated molecule / cluster).
The electrons are contained in the box with zero-boundary conditions;
no interaction with periodic images.

- `L` — side length, a quantity with length units (e.g. `20.0_bohr`, `10.0_angstrom`)
- Returns a `systems::cell` object, which is immediately passed to the `systems::ions` constructor.

Use this for: molecules (H2, N2, coronene, HF).

---

### `systems::cell::cubic(L).periodic()`

```cpp
// From li_bcc.cpp and jellium runs
systems::ions ions(systems::cell::cubic(40.0_b).periodic());
```

Same cubic box but with **full 3D periodic boundary conditions**.
Plane waves are used as the basis set; the simulation repeats infinitely in all three directions.

Use this for: bulk crystals, jellium, anything with translational symmetry.

---

### `systems::cell::orthorhombic(Lx, Ly, Lz).finite()`

```cpp
// From coronene wp_scattering/run.cpp
auto ions = systems::ions(
    systems::cell::orthorhombic(Lx*1.0_b, Ly*1.0_b, Lz*1.0_b).finite()
);
```

Creates a rectangular box with **three independent side lengths**.
Use when your system is not cubic — e.g. a flat molecule where you want
a tall Z-dimension for a wavepacket.

- `Lx, Ly, Lz` — independent box dimensions, each a length quantity.

---

### `systems::cell::lattice(v1, v2, v3).periodic()`

```cpp
// From unit-cells.cpp
systems::cell lattice_cell = systems::cell::lattice(
    vector3<quantity<magnitude::length>>{lattice_param, zero, zero},
    vector3<quantity<magnitude::length>>{zero, lattice_param, zero},
    vector3<quantity<magnitude::length>>{zero, zero, lattice_param}
).periodic();
```

Creates an arbitrary Bravais lattice from three primitive vectors.
Use for non-cubic crystals or when you need explicit control over the lattice vectors.

- `v1, v2, v3` — three `vector3<quantity<magnitude::length>>` objects.
  Each is a 3-component vector like `{a, b, c}` where a/b/c are length quantities.

---

## 2. Inserting Atoms

### `ions.insert(element, position)`

```cpp
// From h2.cpp
double half_bond = 0.37_angstrom;
ions.insert("H", {0.0_angstrom, 0.0_angstrom, -half_bond});
ions.insert("H", {0.0_angstrom, 0.0_angstrom,  half_bond});
```

Adds an atom to the system.

- `element` — atomic symbol as a C++ string literal: `"H"`, `"C"`, `"N"`, `"Li"`, etc.
- `position` — a brace-initializer `{x, y, z}` where each component is a length quantity.

Coordinates are **Cartesian**, measured from the origin (corner of the cell for finite systems).

---

### `ions.insert_fractional(element, position)`

```cpp
// From li_bcc.cpp
ions.insert_fractional("Li", {0.0, 0.0, 0.0});
ions.insert_fractional("Li", {0.5, 0.5, 0.5});
```

Same as `insert`, but uses **fractional (crystal) coordinates**:
each component is a number between 0 and 1 that specifies position as a fraction
of the corresponding lattice vector.

- `position` — plain floating-point `{f1, f2, f3}`, no units needed.

Example: `{0.5, 0.5, 0.5}` puts the atom at the body-center of the cell.

---

### `systems::ions::parse(filename, cell)`

```cpp
// From coronene 02_ground_state_analysis/run.cpp
auto cell = systems::cell::orthorhombic(...).finite();
auto ions = systems::ions::parse(cfg::CORONENE_XYZ, cell);
```

Reads atomic geometry from an **XYZ file** instead of inserting atoms one by one.

- `filename` — path to the `.xyz` file (string).
- `cell` — a pre-constructed `systems::cell` object that defines the box.

The XYZ file format is:
```
36
coronene molecule
C    0.000000    1.228500    0.000000
C    1.063700    0.614250    0.000000
...
```
Line 1: number of atoms. Line 2: comment. Then one line per atom: element, x, y, z in **Angstrom**.

---

## 3. Electron Options

The `systems::electrons` object is constructed by passing `ions` plus a chain of option calls.
All the `.option()` calls return the options object, so they can be chained with dots.

### `systems::electrons(ions, options::electrons{}.option1().option2()..., kpoints)`

```cpp
// From h2.cpp — simplest possible call
systems::electrons electrons(ions, options::electrons{}.cutoff(80.0_Ry));

// From li_bcc.cpp — full example
systems::electrons electrons(ions,
    options::electrons{}
        .cutoff(40.0_Ry)
        .extra_states(4)
        .temperature(kT_1000K),
    input::kpoints::grid({4, 4, 4}, true)
);
```

The third argument (k-points) is optional. If omitted, only the Gamma point (k=0) is used.

---

### `.cutoff(E_cut)` — energy cutoff

```cpp
options::electrons{}.cutoff(40.0_Ry)
options::electrons{}.cutoff(54.0_Ha)
```

Sets the **plane-wave energy cutoff**. This controls the real-space grid:
a higher cutoff means a finer grid, more accuracy, and more memory/time.

The grid spacing `h` relates to the cutoff as: `E_cut = π² / (2 h²)`.
At 40 Ry, this gives roughly h ≈ 0.25 bohr.

- `E_cut` — energy quantity (`_Ry`, `_Ha`, or `_eV`).

---

### `.spacing(h)` — direct grid spacing

```cpp
// From jellium convergence study
options::electrons{}.spacing(0.5_b)
```

Alternative to `.cutoff()`. Sets the real-space grid spacing directly.

- `h` — length quantity, e.g. `0.5_b` for 0.5 bohr per grid point.

---

### `.extra_states(N)` — empty orbital slots

```cpp
// From jellium ground state
options::electrons{}.extra_states(8)

// From coronene LEED (1 for wavepacket + 2 for SCF buffer)
options::electrons{}.extra_states(3)
```

Allocates `N` extra orbital slots **above** the occupied states.
These are initially empty (occupation = 0).

This is needed for:
- Wavepacket injection (you inject the WP into the last extra state)
- SCF convergence buffer in metallic systems
- Visualizing unoccupied states

---

### `.extra_electrons(N)` — free electrons without ions

```cpp
// From jellium ground state (no atomic cores, only a uniform positive background)
options::electrons{}
    .extra_electrons(N_ELECTRONS)
    .extra_states(8)
```

For **jellium** (uniform positive background + free electrons), there are no ion insertions.
Instead, you tell INQ directly how many electrons to put in the box.
INQ automatically provides the compensating uniform positive background charge.

---

### `.temperature(T)` — Fermi-Dirac smearing

```cpp
// From li_bcc.cpp
double kT_1000K = 0.086 / 27.211;  // in Hartree
options::electrons{}.temperature(kT_1000K * 1.0_Ha)

// From jellium convergence
options::electrons{}.temperature(SMEAR_EV * 1.0_eV)
```

Enables **Fermi-Dirac thermal smearing** of orbital occupations.
Instead of sharp 0/2 occupations, states near the Fermi level get fractional occupations.
This helps convergence for metals and jellium.

- `T` — energy quantity equal to k_B × temperature. For 1000 K: 0.086 eV.

---

### `input::kpoints::gamma()` — single Gamma point

```cpp
// From jellium
input::kpoints::gamma()
```

Uses only the k=0 Gamma point for BZ integration.
Appropriate for large supercells, isolated systems, or when you are not doing
a periodic crystal calculation.

Pass as the third argument to `systems::electrons(...)`.

---

### `input::kpoints::grid({n1, n2, n3}, shifted)` — Monkhorst-Pack grid

```cpp
// From li_bcc.cpp
input::kpoints::grid({4, 4, 4}, true)
```

Creates a regular k-point mesh for periodic calculations.

- `{n1, n2, n3}` — number of k-points along each reciprocal lattice direction.
- `shifted` — `true` = Monkhorst-Pack (shifted off Gamma), `false` = Gamma-centered.

Pass as the third argument to `systems::electrons(...)`.

---

## 4. Ground State

### `ground_state::initial_guess(ions, electrons)`

```cpp
// From h2.cpp
ground_state::initial_guess(ions, electrons);
```

Initialises the electronic structure with a **superposition of atomic densities**.
This is a cheap starting point for the SCF loop — it gives a reasonable initial
density from known atomic data.

Must be called **before** `ground_state::calculate`.
Modifies `electrons` in-place.

---

### `ground_state::calculate(ions, electrons, theory, options)`

```cpp
// From h2.cpp
auto result = ground_state::calculate(ions, electrons,
    options::theory{}.pbe(),
    options::ground_state{}.calculate_forces()
);

// From jellium
auto gs = ground_state::calculate(ions, electrons,
    options::theory{}.lda(),
    options::ground_state{}.energy_tolerance(1e-8_Ha)
);
```

Runs the **self-consistent field (SCF)** loop to convergence.
This is the main ground-state calculation.

Arguments:
- `ions` — your ionic geometry.
- `electrons` — the electronic structure object (modified in-place; contains converged orbitals after).
- `options::theory{}.functional()` — which exchange-correlation functional to use.
- `options::ground_state{}.option1().option2()...` — convergence settings.

Returns a `ground_state::results` object (see Section 5).

---

### Theory options

```cpp
options::theory{}.pbe()         // PBE GGA — default for molecules and materials
options::theory{}.lda()         // LDA (Perdew-Zunger 1981) — default for jellium
options::theory{}.non_interacting()  // no Hartree, no XC — kinetic energy only
```

- **PBE** — Perdew-Burke-Ernzerhof (1996). Good for molecules and most materials.
- **LDA** — Local Density Approximation (exact for the uniform electron gas).
  Best choice for jellium.
- **non_interacting** — treats electrons as completely independent (no e-e interactions).
  Used for testing free-particle propagation.

---

### Ground state options (chainable)

```cpp
options::ground_state{}
    .energy_tolerance(1e-6_Ha)   // stop when ΔE < this between iterations
    .mixing(0.1)                  // fraction of new density to mix in each step
    .mixing_ndim(8)               // Broyden history depth
    .broyden_mixing()             // use Broyden mixer (better than linear for hard systems)
    .calculate_forces()           // also compute forces on ions
    .max_steps(300)               // hard limit on SCF iterations
```

**`.energy_tolerance(tol)`** — SCF convergence criterion.
The loop stops when the total energy changes by less than `tol` between iterations.
Typical values: `1e-6_Ha` (normal), `1e-8_Ha` (tight).

**`.mixing(alpha)`** — linear mixing parameter.
Each step, the new density is: `n_new = (1-α) n_old + α n_predicted`.
Lower α = more stable but slower. Typical: 0.05–0.3.

**`.mixing_ndim(N)` + `.broyden_mixing()`** — use Broyden's multisecant mixer
instead of linear mixing. Uses the last N density histories to extrapolate.
More expensive per step but far fewer iterations. Use for metals or molecules
that fail to converge with linear mixing.

**`.calculate_forces()`** — after SCF, compute the Hellmann-Feynman force on each nucleus.
Result goes into `result.forces`.

**`.max_steps(N)`** — abort if the SCF does not converge in N iterations.

---

## 5. Ground State Results

After `ground_state::calculate(...)` returns, the result object has these fields:

```cpp
auto result = ground_state::calculate(...);

// Energies
result.energy.total()      // total DFT energy: T_s + E_H + E_xc + E_ext + E_ion + E_nlpp  [Ha]
result.energy.kinetic()    // kinetic energy of Kohn-Sham electrons  [Ha]
result.energy.hartree()    // classical Coulomb self-energy of the electron density  [Ha]
result.energy.xc()         // exchange-correlation energy  [Ha]
result.energy.external()   // interaction with ionic pseudopotentials (local part)  [Ha]
result.energy.non_local()  // non-local pseudopotential energy (Kleinman-Bylander)  [Ha]
result.energy.ion()        // ion-ion (Ewald) repulsion  [Ha]
result.energy.nvxc()       // integral of n * V_xc (used in some decompositions)  [Ha]

// Other
result.forces    // vector<vector3<double>>, forces on each atom in Ha/bohr
result.dipole    // vector3<double>, electric dipole moment in e·bohr
result.total_iter  // int, number of SCF iterations to convergence
```

Example of printing total energy in eV:

```cpp
double eV = 27.211396;  // 1 Ha in eV
std::cout << "Total energy: " << result.energy.total() * eV << " eV\n";
```

The forces vector is indexed by atom, in the order the atoms were inserted:
```cpp
for (int i = 0; i < ions.size(); i++) {
    auto f = result.forces[i];
    std::cout << "Atom " << i << ": Fx=" << f[0] << " Fy=" << f[1] << " Fz=" << f[2] << "\n";
}
```

---

## 6. Real-Time Propagation Options

### `real_time::propagate(ions, electrons, callback, theory, options, ...perturbations)`

```cpp
// From li_bcc.cpp
real_time::propagate(
    ions, electrons, output,     // system + callback
    options::theory{}.pbe(),
    options::real_time{}
        .num_steps(2000)
        .dt(0.04_atomictime)
        .impulsive()
        .observables_current()
);
```

Runs the **TDDFT time propagation**.
Updates `electrons` at each step; calls `callback` at each step with observable data.

Arguments in order:
1. `ions` — atomic structure (positions updated if Ehrenfest dynamics is used)
2. `electrons` — electronic structure (orbitals are overwritten each step)
3. `callback` — lambda `[&](auto data){ ... }` called every step (see Section 7)
4. `theory` — same theory options as for ground state
5. `options::real_time{}....` — time step settings
6. (optional) perturbation objects — e.g. laser field

---

### Real-time options (chainable)

```cpp
options::real_time{}
    .num_steps(2000)          // total number of time steps to propagate
    .dt(0.02_atomictime)      // time step size
    .impulsive()              // ions move at constant velocity (no forces on ions)
    .ehrenfest()              // ions respond to electronic forces each step
    .observables_dipole()     // compute dipole moment at each step
    .observables_current()    // compute electronic current density J at each step
    .etrs()                   // use ETRS (time-reversal symmetry) propagator
```

**`.num_steps(N)`** — how many steps to take.
Total simulation time = N × dt.

**`.dt(time)`** — time step.
Units: `_atomictime` means atomic time units (ℏ/Ha ≈ 0.0242 femtoseconds).
Typical values: `0.02_atomictime` to `0.05_atomictime`.

**`.impulsive()`** — ions do not respond to the electronic forces.
They either stay fixed or move with a constant, preset velocity.
Use this when your ions are projectiles (LEED, ion scattering) or when you only care about
the electronic response.

**`.ehrenfest()`** — Ehrenfest dynamics: the ions move under the force `F = -∂E/∂R`
computed at each step. Use for coupled electron-ion dynamics.

**`.observables_dipole()`** — at each step, compute the total electronic dipole moment
`d(t) = -∫ r n(r,t) d³r` and make it available as `data.dipole()` in the callback.

**`.observables_current()`** — at each step, compute the total electronic current density.
Makes `data.current()` available. Used in the Li BCC absorption spectrum calculation.

**`.etrs()`** — selects the ETRS (Enforced Time-Reversal Symmetry) propagator.
This is a 4th-order symplectic integrator that exactly conserves a modified energy.
Used for non-interacting systems.

---

## 7. The Propagation Callback (data object)

The callback is a lambda you write yourself. It is called **once per time step**
during `real_time::propagate`. The `data` argument gives you access to current observables.

### Pattern

```cpp
auto output = [&](auto data) {
    // runs every step
    if (data.every(10)) {
        // runs every 10 steps
        double t = data.time();
        int step = data.iter();
        double E = data.energy().total();
        // write to file, print, etc.
    }
};
```

### Methods on `data`

**`data.every(N)`** → `bool`
Returns `true` on every Nth step (step 0, N, 2N, 3N, ...).
Use this to control how often you write output, since writing every step
creates enormous files.

**`data.time()`** → `double`
Current simulation time in atomic time units (a.u.).
To convert to femtoseconds: multiply by 0.02419.

**`data.iter()`** → `int`
Current step index, starting from 0.

**`data.energy().total()`** → `double` (Hartree)
Total energy at this step. Same sub-fields as the ground state result:
`.kinetic()`, `.hartree()`, `.xc()`, etc.

This is the main diagnostic for energy conservation —
in a correctly-running TDDFT simulation, `data.energy().total()` should stay
nearly constant throughout the propagation.

**`data.positions()`** → `vector<vector3<double>>`
Current atomic positions in bohr.
For fixed-ion runs, these are constant.
For Ehrenfest dynamics, these update each step.

**`data.current()`** → `vector3<double>` (only if `.observables_current()` was set)
Total electronic current density integrated over the cell: `J = ∫ j(r) d³r`.
Components `J[0]`, `J[1]`, `J[2]` in atomic units.

**`data.dipole()`** → `vector3<double>` (only if `.observables_dipole()` was set)
Electronic dipole moment in e·bohr.

### Full example from li_bcc.cpp

```cpp
std::ofstream dipole_file("results/dipole.csv");
dipole_file << "t_au,Dx,Dy,Dz\n";

auto output = [&](auto data) {
    if (data.every(1)) {
        auto d = data.dipole();
        dipole_file << data.time() << "," << d[0] << "," << d[1] << "," << d[2] << "\n";
    }
};

real_time::propagate(ions, electrons, output,
    options::theory{}.pbe(),
    options::real_time{}
        .num_steps(2000)
        .dt(0.04_atomictime)
        .impulsive()
        .observables_dipole()
);
```

---

## 8. Perturbations

### `perturbations::laser(E_field, frequency, gauge)`

```cpp
// From HF-laser-perturbation.cpp
auto laser = perturbations::laser(
    {0.0, 0.0, 1.0e-3},              // E-field vector (a.u.), pointing in z
    laser_frequency,                  // photon energy, e.g. 0.4911_eV
    perturbations::gauge::length      // dipole approximation, length gauge
);

real_time::propagate(ions, electrons, output,
    options::theory{}.pbe(),
    options::real_time{}.num_steps(2000).dt(0.04_atomictime).impulsive(),
    laser    // <-- add the perturbation here
);
```

Applies a time-dependent electric field `E(t) = E_field × cos(ω t)` during propagation.

Arguments:
- `E_field` — `{Ex, Ey, Ez}` plain doubles (not quantity literals), in atomic units.
  `1.0e-3` a.u. is a typical weak perturbation.
- `frequency` — photon energy as a quantity (e.g. `0.4911_eV`).
- `gauge` — `perturbations::gauge::length` (dipole, position coupling)
  or `perturbations::gauge::velocity` (momentum coupling).

The perturbation is passed as an additional argument to `real_time::propagate`
after the options object.

---

## 9. Accessing Orbitals and Density

These are used inside the TDDFT callback or after the ground state to read out or write
the wavefunction and charge density.

### `electrons.density()` — charge density field

```cpp
auto density = electrons.density();    // returns a field object
auto& basis  = density.basis();         // the real-space grid metadata
```

Returns the total electron density `n(r)` as a scalar field on the real-space grid.
To iterate over grid points:

```cpp
auto& basis = density.basis();
int Nx = basis.sizes()[0];
int Ny = basis.sizes()[1];
int Nz = basis.sizes()[2];
double dV = basis.volume_element();    // dx * dy * dz in bohr^3

auto hc = density.cubic();            // 3D array accessor
for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
        for (int iz = 0; iz < Nz; iz++) {
            double n = hc[ix][iy][iz];   // density at this grid point, in bohr^-3
        }
    }
}
```

---

### `electrons.kpin()` — Kohn-Sham orbitals

```cpp
auto& phi = electrons.kpin()[0];    // [0] = Gamma-point orbital set
```

`kpin()` returns a vector, one element per k-point.
For Gamma-point calculations (which is all your work), you only ever use `[0]`.

`phi` is the orbital set object. It holds all KS orbitals on the real-space grid.

**Accessing individual orbital values:**
```cpp
auto hc = phi.hypercubic();    // 4D array: [ix][iy][iz][state_index]

int Nst = phi.spinor_set_size();   // total number of states (occupied + extra)
int Nx = ...; int Ny = ...; int Nz = ...;   // from basis

for (int ist = 0; ist < Nst; ist++) {
    for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) {
            for (int iz = 0; iz < Nz; iz++) {
                auto psi = hc[ix][iy][iz][ist];  // complex<double>
                double re = psi.real();
                double im = psi.imag();
                double density = norm(psi);       // |psi|^2
            }
        }
    }
}
```

**State index convention:**
- States 0 through N_occupied-1 are the KS orbitals filled by electrons.
- States N_occupied through N_occupied + N_extra - 1 are the empty extra states.
- The **last** state (`ist = Nst - 1`) is where you inject wavepackets.

---

### `basis.point_op().rvector_cartesian(ix, iy, iz)` — grid coordinates

```cpp
auto& basis = phi.basis();
auto po = basis.point_op();

for (int ix = 0; ix < Nx; ix++) {
    for (int iy = 0; iy < Ny; iy++) {
        for (int iz = 0; iz < Nz; iz++) {
            auto r = po.rvector_cartesian(ix, iy, iz);
            double x = r[0];   // bohr
            double y = r[1];   // bohr
            double z = r[2];   // bohr
        }
    }
}
```

Converts grid index `(ix, iy, iz)` to Cartesian coordinates in bohr.

---

### `basis.sizes()`, `basis.volume_element()`, `basis.rspacing()`

```cpp
auto& basis = phi.basis();
int Nx = basis.sizes()[0];     // grid points along x
int Ny = basis.sizes()[1];
int Nz = basis.sizes()[2];
double dV = basis.volume_element();      // integration weight per grid point [bohr^3]
double dx = basis.rspacing()[0];         // grid spacing along x [bohr]
double dz = basis.rspacing()[2];         // grid spacing along z [bohr]
```

---

### `electrons.eigenvalues()` — KS eigenvalues

```cpp
auto evals = electrons.eigenvalues();   // 2D array [kpoint][state]

int n_states = ...;
for (int i = 0; i < n_states; i++) {
    double ev = evals[0][i];    // [0] = Gamma point; eigenvalue in Hartree
}
```

The eigenvalues are the Kohn-Sham orbital energies ε_i in Hartree.
For jellium, you can compare these to k²/2 + V_xc to validate the XC offset.

---

### `electrons.occupations()` — orbital occupations

```cpp
auto& occ = electrons.occupations();   // 2D array [kpoint][state]
occ[0][ist_wp] = 0.0;     // set the WP orbital to occupation 0 (exclude from density)
occ[0][ist_wp] = 1.0;     // set to occupation 1 (include in density as 1 electron)
```

Returns a reference to the occupation array (you can modify it directly).
Occupations are doubles: 2.0 for a spin-degenerate filled orbital, 1.0 for spin-polarized,
0.0 for empty.

Used in the LEED simulations to set the wavepacket orbital occupation.

---

### `electrons.save(path)` and `electrons.load(path)`

```cpp
// After ground state calculation:
electrons.save("results/gs_checkpoint");

// At start of a TDDFT run that restarts from saved ground state:
electrons.load("results/gs_checkpoint");
```

Saves/loads the complete electronic state (all orbitals, density, metadata)
to/from a directory on disk. The directory is created by `save`.
Use this to decouple the ground state run from the TDDFT run.

---

## 10. Unit Literals Quick Reference

| Literal | Quantity | Value |
|---|---|---|
| `1.0_b` or `1.0_bohr` | length | 1 bohr = 0.529177 Å |
| `1.0_angstrom` | length | 1 Å |
| `1.0_Ha` | energy | 1 Hartree = 27.211 eV |
| `1.0_Ry` | energy | 1 Rydberg = 0.5 Hartree |
| `1.0_eV` | energy | 1 eV = 0.036749 Ha |
| `1.0_atomictime` | time | 1 a.u. ≈ 0.02419 fs |

To use as plain doubles (for output):
```cpp
double eV_per_Ha = 27.211396;
double angstrom_per_bohr = 0.529177;
double fs_per_au = 0.024189;

double energy_eV = result.energy.total() * eV_per_Ha;
double spacing_angstrom = dx * angstrom_per_bohr;
double time_fs = data.time() * fs_per_au;
```
