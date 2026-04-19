# Custom Utility Functions Reference

This covers every function defined in your `utils.hpp`, `config.hpp`, and `jellium_utils.hpp` files.
These are functions *you* wrote (or that were written as part of the research project),
as distinct from the INQ library API.

---

## Table of Contents

1. [Jellium utilities (`jellium_utils.hpp`)](#1-jellium-utilities-jellium_utilshpp)
2. [Coronene configuration (`config.hpp`)](#2-coronene-configuration-confighpp)
3. [Coronene utilities (`utils.hpp`)](#3-coronene-utilities-utilshpp)

---

## 1. Jellium Utilities (`jellium_utils.hpp`)

File: `/local/data/public/skcb2/tddft/ResearchProject/jellium/01_ground_state/jellium_utils.hpp`

This file provides analytic formulas for the uniform electron gas (jellium).
These are used both to characterise the system (e.g. what is the Fermi energy for these parameters?)
and to validate the DFT output (e.g. do the eigenvalues shift by exactly V_xc?).

---

### `wigner_seitz_radius(N, L)`

```cpp
inline double wigner_seitz_radius(int N, double L);
```

**What it does:** Computes the Wigner-Seitz radius r_s for a jellium cell.
r_s is the radius of a sphere whose volume equals the volume per electron.
It characterizes how dense the electron gas is.

**Arguments:**
- `N` — number of electrons (integer)
- `L` — cell side length in bohr (double, not a quantity literal)

**Returns:** r_s in bohr

**Formula:** (4π/3) r_s³ = L³/N → r_s = (3/(4πn₀))^(1/3) where n₀ = N/L³

```cpp
double rs = wigner_seitz_radius(40, 40.0);   // → 7.256 bohr
```

---

### `mean_density(N, L)`

```cpp
inline double mean_density(int N, double L);
```

**What it does:** Returns the average electron number density.

**Returns:** n₀ = N/L³ in electrons per bohr³

---

### `fermi_wavevector(N, L)`

```cpp
inline double fermi_wavevector(int N, double L);
```

**What it does:** Computes the free-electron Fermi wavevector k_F.
This is the magnitude of the largest occupied plane-wave momentum at T=0.

**Returns:** k_F in bohr⁻¹

**Formula:** k_F = (3π²n₀)^(1/3)

---

### `fermi_energy(N, L)`

```cpp
inline double fermi_energy(int N, double L);
```

**What it does:** Free-electron Fermi energy (kinetic energy of the highest occupied state).

**Returns:** E_F = k_F²/2 in Hartree

---

### `plasmon_frequency(N, L)`

```cpp
inline double plasmon_frequency(int N, double L);
```

**What it does:** Drude model bulk plasmon frequency.

**Returns:** ω_p = √(4πn₀) in Hartree (atomic frequency units)

---

### `exchange_energy_pz81(rs)` and `exchange_potential_pz81(rs)`

```cpp
inline double exchange_energy_pz81(double rs);
inline double exchange_potential_pz81(double rs);
```

**What they do:** Analytic LDA exchange energy per electron and exchange potential,
from Perdew & Zunger (1981), which use the Dirac exchange formula.

**Argument:** `rs` — Wigner-Seitz radius in bohr

**Returns:** Energy/potential in Hartree

**Formulas:**
- ε_x = -0.4582 / r_s
- V_x = (4/3) ε_x = -0.6109 / r_s

---

### `correlation_energy_pz81(rs)` and `correlation_potential_pz81(rs)`

```cpp
inline double correlation_energy_pz81(double rs);
inline double correlation_potential_pz81(double rs);
```

**What they do:** LDA correlation energy per electron and correlation potential,
fitted to Ceperley-Alder Quantum Monte Carlo data by Perdew & Zunger (1981).

Uses two different fitting formulas depending on density:
- `rs < 1.0` (high density): logarithmic form
- `rs >= 1.0` (metallic density): Padé-approximant form

**Argument:** `rs` in bohr

**Returns:** Energy/potential in Hartree

These match what INQ's LDA functional (`options::theory{}.lda()`) uses internally.
You use them here to independently predict eigenvalue positions.

---

### `exc_pz81(rs)` and `vxc_pz81(rs)`

```cpp
inline double exc_pz81(double rs);
inline double vxc_pz81(double rs);
```

**What they do:** Total XC energy per electron and total XC potential.
These are just the sums of the exchange and correlation components above.

For a uniform electron gas, V_xc is a constant uniform shift applied to
all Kohn-Sham eigenvalues. Your eigenvalue validation checks that
ε_i - k_i²/2 = V_xc for all states.

---

### `free_electron_shells(L, n2_max)`

```cpp
struct Shell {
    int    n2;           // |n|² = nx² + ny² + nz² for integer quantum numbers
    int    degeneracy;   // how many distinct (nx,ny,nz) triples give this |n|²
    double energy_Ha;    // kinetic energy: (1/2) * (2π/L)² * n²
};

inline std::vector<Shell> free_electron_shells(double L, int n2_max = 6);
```

**What it does:** Enumerates the free-electron energy shells for a cubic periodic box.
In a periodic box of side L, the allowed momenta are k = (2π/L)(nx, ny, nz)
for integer (nx, ny, nz). States with the same |n|² = nx² + ny² + nz² are degenerate.

**Arguments:**
- `L` — cell side in bohr
- `n2_max` — maximum |n|² to enumerate (default 6 covers the first 4 shells)

**Returns:** Vector of Shell structs sorted by energy.

**Shell structure for a cubic box:**

| n² | Degeneracy | Example (nx,ny,nz) triplets |
|---|---|---|
| 0 | 1 | (0,0,0) |
| 1 | 6 | (±1,0,0), (0,±1,0), (0,0,±1) |
| 2 | 12 | (±1,±1,0) and permutations |
| 3 | 8 | (±1,±1,±1) |
| 4 | 6 | (±2,0,0) and permutations |

With spin degeneracy (factor 2): n²=0 holds 2, n²=1 holds 12, n²=2 holds 24, etc.
This is why magic numbers of electrons (closed shells) are 2, 14, 38, 54, 66...

---

### `kinetic_energy_shells(N_electrons, L)`

```cpp
inline double kinetic_energy_shells(int N_electrons, double L);
```

**What it does:** Computes the total kinetic energy at T=0 by filling shells in order,
up to the last occupied shell.

**Arguments:**
- `N_electrons` — must be a shell-closure number (2, 14, 38, 54, 66, ...)
- `L` — cell side in bohr

**Returns:** Total T_s in Hartree

Used to compare against INQ's computed `gs.energy.kinetic()`.

---

### `predicted_total_energy(N_electrons, L)`

```cpp
inline double predicted_total_energy(int N_electrons, double L);
```

**What it does:** Predicts the total DFT energy for a jellium system.
For jellium, the Hartree, external, and ion-ion terms all cancel exactly,
leaving only kinetic + XC.

**Formula:** E_tot ≈ T_s + N × ε_xc(r_s)

**Returns:** Predicted total energy in Hartree

This is compared against INQ's `gs.energy.total()` to validate the ground state.

---

## 2. Coronene Configuration (`config.hpp`)

Files: `ResearchProject/systems/coronene/04_leed_simulation/config.hpp` and per-run copies.

These files centralise all the numerical parameters for the coronene LEED simulation.
Everything is defined as `constexpr` constants or inline functions.
The idea is that you change only this one file to re-run with different parameters.

---

### Key constants

```cpp
// Cell dimensions
constexpr double LX_BOHR = ...;   // box width in x (bohr)
constexpr double LY_BOHR = ...;   // box width in y (bohr)
constexpr double LZ_BOHR = ...;   // box height in z (bohr)

// Energy cutoff
constexpr double ECUT_HA = 40.0;  // plane-wave cutoff in Hartree

// Wavepacket parameters
constexpr double WP_EKIN_EV    = 200.0;        // kinetic energy of incident electron [eV]
constexpr double WP_EKIN_HA    = WP_EKIN_EV / 27.211396;  // same in Hartree
constexpr double WP_D_A        = 1.4;          // Gaussian width σ in Angstrom
constexpr double WP_D_BOHR     = WP_D_A / 0.529177;  // same in bohr
constexpr double WP_D_IMPACT_A    = ...;       // initial z distance above molecule [Angstrom]
constexpr double WP_D_IMPACT_BOHR = ...;       // same in bohr
constexpr double WP_OCCUPATION    = 1.0;       // occupation assigned to WP orbital

// Time propagation
constexpr double DT_AU      = 0.02;            // time step in atomic time units
constexpr int    N_STEPS    = 516;             // total propagation steps
constexpr double T1_AU      = ...;             // time when LEED integration starts [a.u.]
constexpr double T2_AU      = ...;             // time when LEED integration ends [a.u.]

// Path to geometry file
constexpr const char* CORONENE_XYZ = "geometry/coronene_centered.xyz";
```

---

### `cfg::make_cell()`

```cpp
inline auto cfg::make_cell() {
    return inq::systems::cell::orthorhombic(LX_BOHR*1.0_b, LY_BOHR*1.0_b, LZ_BOHR*1.0_b).finite();
}
```

Factory function that creates the coronene simulation cell.
Encapsulates the cell construction so run.cpp just calls `cfg::make_cell()` without
having to repeat the dimensions.

---

### `cfg::wp_k0()`

```cpp
inline double cfg::wp_k0() { return std::sqrt(2.0 * WP_EKIN_HA); }
```

**What it does:** Computes the wavepacket momentum magnitude from kinetic energy.
This is the free-electron dispersion relation: E = k²/2 → k = √(2E).

**Returns:** k₀ in bohr⁻¹

For a 200 eV electron: k₀ ≈ 3.83 bohr⁻¹.

---

### `cfg::wp_norm()`

```cpp
inline double cfg::wp_norm() {
    return std::pow(M_PI * WP_D_BOHR * WP_D_BOHR, -0.75);
}
```

**What it does:** The normalization prefactor for a 3D Gaussian:
N = (π d²)^(-3/4) such that ∫ |N exp(-r²/(2d²))|² d³r = 1.

**Returns:** Normalization constant (units: bohr^(-3/2))

---

### `cfg::WP_BX()`, `cfg::WP_BY()`, `cfg::WP_BZ()`

```cpp
inline double cfg::WP_BX() { return LX_BOHR / 2.0; }
inline double cfg::WP_BY() { return LY_BOHR / 2.0; }
inline double cfg::WP_BZ() { return LZ_BOHR / 2.0 + WP_D_IMPACT_BOHR; }
```

**What they do:** Returns the initial centre position of the wavepacket.
X and Y are at the cell centre (above the centre of the coronene molecule).
Z is above the coronene plane by `WP_D_IMPACT_BOHR`.

**Returns:** Position in bohr.

---

### `cfg::Z_FLAKE_BOHR()` and `cfg::Z_OBS_BOHR()`

```cpp
inline double cfg::Z_FLAKE_BOHR() { return LZ_BOHR / 2.0; }
inline double cfg::Z_OBS_BOHR()   { return LZ_BOHR / 2.0 + WP_D_IMPACT_BOHR; }
```

`Z_FLAKE_BOHR()` — z-coordinate of the coronene molecule plane.
The molecule is centred in the z-direction of the box.

`Z_OBS_BOHR()` — z-coordinate of the observation plane
(where the LEED detector would be placed).
The WP starts here and the density is accumulated here to form the LEED pattern.

---

## 3. Coronene Utilities (`utils.hpp`)

Files: `04_leed_simulation/utils.hpp` and per-run copies (progressively more complex).

All functions live in the `leed_utils` namespace.

---

### `leed_utils::inject_wp(electrons, bx, by, bz, kx, ky, kz)`

```cpp
inline void inject_wp(
    inq::systems::electrons& electrons,
    double bx, double by, double bz,    // wavepacket centre [bohr]
    double kx, double ky, double kz     // momentum vector [bohr^-1]
);
```

**What it does:**
Writes a Gaussian wavepacket into the **last extra-state orbital**.
This modifies the orbital data directly in GPU memory.

The wavepacket has the form:
```
ψ(r) = N * exp(-|r - b|² / (2 d²)) * exp(i k·r)
```
where:
- `b = (bx, by, bz)` — packet centre
- `d` = `cfg::WP_D_BOHR` — Gaussian width parameter
- `N` = `cfg::wp_norm()` — normalization prefactor
- `k = (kx, ky, kz)` — momentum vector

**Implementation detail:**
It loops over all local grid points `(ix, iy, iz)`, computes the Cartesian coordinates
using `phi.basis().point_op().rvector_cartesian(ix, iy, iz)`, evaluates the Gaussian
amplitude and the phase factor `exp(ik·r)`, and stores the result as a `complex<double>`
in `phi.hypercubic()[ix][iy][iz][last_state_idx]`.

After calling this, you set `electrons.occupations()[0][last_state_idx] = WP_OCCUPATION`
to include the wavepacket electron in the density.

---

### `leed_utils::validate_wp(electrons)`

```cpp
inline std::pair<double, double> validate_wp(
    inq::systems::electrons const& electrons
);
```

**What it does:**
Checks the injected wavepacket by numerically integrating |ψ|² over the grid.
Returns the computed norm (should be ≈ 1.0).

**Returns:** `{norm, cfg::WP_EKIN_HA}` — the numerically computed norm plus
the expected kinetic energy (for reference).

```cpp
auto [norm, ke] = validate_wp(electrons);
std::cout << "WP norm: " << norm << " (should be ~1.0)\n";
```

---

### `leed_utils::iz_nearest(electrons, z_bohr)` *(run_002)*

```cpp
inline int iz_nearest(
    inq::systems::electrons const& electrons,
    double z_bohr
);
```

**What it does:**
Converts a physical z-coordinate (in bohr) to the nearest grid index `iz`.
Handles wrapping for periodic grids.

**Arguments:**
- `electrons` — to get the grid spacing from
- `z_bohr` — target z position in bohr

**Returns:** Grid index `iz` (integer, 0 ≤ iz < Nz)

Used to find which grid row corresponds to the coronene plane or the observation plane.

---

### `leed_utils::extract_density_slice(electrons, z_target)` *(run_002)*

```cpp
inline std::vector<std::vector<double>> extract_density_slice(
    inq::systems::electrons const& electrons,
    double z_target
);
```

**What it does:**
Extracts a 2D slice of the total electron density `n(x, y, z_target)` at a fixed z.

**Returns:** 2D vector indexed as `[iy][ix]` containing density values in bohr⁻³.

This is called repeatedly during propagation (every SNAPSHOT_INTERVAL steps)
to accumulate the LEED pattern: `I(x,y) = ∫ n(x, y, z_obs, t) dt`.

---

### `leed_utils::save_density_slice(slice, t, z, filename)` *(run_002)*

```cpp
inline void save_density_slice(
    std::vector<std::vector<double>> const& slice,
    double t,       // current time [a.u.]
    double z,       // z-plane position [bohr]
    const std::string& filename
);
```

**What it does:**
Writes a 2D density slice to a text file.
The output format is described in the output file formats document.

---

### `leed_utils::extract_z_profile(electrons)` *(run_002)*

```cpp
inline std::vector<double> extract_z_profile(
    inq::systems::electrons const& electrons
);
```

**What it does:**
Extracts a 1D density profile `n(z)` through the centre of the cell
(at `ix = Nx/2, iy = Ny/2`).

**Returns:** Vector of length Nz with density values in bohr⁻³.

Used to track the wavepacket's z-position over time.

---

### `leed_utils::compute_overlap_matrix(gs_electrons, electrons)` *(run_003)*

```cpp
inline std::vector<std::vector<std::complex<double>>> compute_overlap_matrix(
    inq::systems::electrons const& gs_electrons,
    inq::systems::electrons const& electrons
);
```

**What it does:**
Computes the N×N overlap matrix between the ground-state orbitals and the
current TDDFT orbitals:
```
S_ij = <ψ_i^GS | ψ_j(t)>  =  Σ_r  ψ_i^GS*(r) * ψ_j(t,r) * dV
```

This tells you how much each current orbital "looks like" each ground-state orbital.
The diagonal |S_ii|² gives the projected occupation of ground-state orbital i.

**Arguments:**
- `gs_electrons` — the saved ground-state electrons (reference)
- `electrons` — the current electrons at time t

**Returns:** 2D complex vector, indexed `[i][j]`, with N = number of occupied states.

---

### `leed_utils::save_orbital_3d(electrons, ist, filename)` *(run_003)*

```cpp
inline void save_orbital_3d(
    inq::systems::electrons const& electrons,
    int ist,                      // orbital state index
    const std::string& filename
);
```

**What it does:**
Saves the full 3D orbital density `|ψ_ist(r)|²` to a text file.
The file format is described in the output file formats document.

---

### `leed_utils::save_density_3d(electrons, filename)` *(run_003)*

```cpp
inline void save_density_3d(
    inq::systems::electrons const& electrons,
    const std::string& filename
);
```

**What it does:**
Saves the total 3D electron density `n(r)` to a text file.
Each line is one grid point value.

---

### `leed_utils::mkdir_p(dirname)` *(run_003)*

```cpp
inline void mkdir_p(const std::string& dirname);
```

**What it does:**
Creates a directory (and all intermediate parent directories),
like `mkdir -p` in the shell.
Does nothing if the directory already exists.

Used at the start of every run to ensure output directories exist
before attempting to write files.

---

### `leed_utils::save_grid_coords(electrons, dirname)` *(run_003)*

```cpp
inline void save_grid_coords(
    inq::systems::electrons const& electrons,
    const std::string& dirname
);
```

**What it does:**
Writes the grid coordinate arrays (x, y, z positions of all grid points)
to text files in `dirname`. This is written once at the start of the run.

The Python analysis scripts read these files to reconstruct real-space coordinates
when loading the 3D density or orbital snapshots.
