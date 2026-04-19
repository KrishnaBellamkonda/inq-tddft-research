# Output File Formats

This document describes every output file written by your C++ simulations.
Each section shows the exact format with real examples taken from actual output files.

---

## Table of Contents

1. [Jellium ground state: orbital slice files](#1-jellium-ground-state-orbital-slice-files)
2. [Jellium ground state: eigenvalue file](#2-jellium-ground-state-eigenvalue-file)
3. [Jellium convergence: `convergence_results.csv`](#3-jellium-convergence-convergence_resultscsv)
4. [Jellium WP propagation: 2D slice files (`slice_tNNN.txt`)](#4-jellium-wp-propagation-2d-slice-files-slice_tnnntxt)
5. [Jellium WP propagation: 3D density files (`density3d_tNNN.txt`)](#5-jellium-wp-propagation-3d-density-files-density3d_tnnntxt)
6. [Coronene: `sim_summary.txt`](#6-coronene-sim_summarytxt)
7. [Coronene: snapshot files (`snapshot_tNNNN.txt`)](#7-coronene-snapshot-files-snapshot_tnnnnntxt)
8. [Coronene: LEED pattern (`leed_pattern.txt`)](#8-coronene-leed-pattern-leed_patterntxt)

---

## 1. Jellium Ground State: Orbital Slice Files

**Location:** `ResearchProject/jellium/01_ground_state/results/orbitals/`

**Files:**
- `grid_slice.txt` — shared grid metadata (written once)
- `orbital_N_n2_M_real.txt` — real part of orbital N (shell quantum number n²=M)
- `orbital_N_n2_M_imag.txt` — imaginary part

### `grid_slice.txt`

The grid is an 80×80 slice at z = L/2 (the midplane of the box).
The file records the (x, y) coordinates of every point in that slice.

```
# 2D slice coordinates at z = 20 bohr
# Grid: 80x80 points,  h = 0.5 bohr
# Columns: ix  iy  x_bohr  y_bohr
0  0  0.000000  0.000000
0  1  0.000000  0.500000
0  2  0.000000  1.000000
0  3  0.000000  1.500000
...
0  79  0.000000  39.500000
1  0   0.500000  0.000000
1  1   0.500000  0.500000
...
```

**Columns:**
- `ix` — x-direction grid index (0 to Nx-1 = 79)
- `iy` — y-direction grid index (0 to Ny-1 = 79)
- `x_bohr` — x coordinate in bohr = ix × h
- `y_bohr` — y coordinate in bohr = iy × h

**Row ordering:** outer loop over ix (x is slow), inner loop over iy (y is fast).
Total rows: 80 × 80 = 6400.

---

### `orbital_N_n2_M_real.txt`

The real part Re[ψ_k(r)] of a single Kohn-Sham orbital, sampled on the 2D midplane slice.

```
# Jellium KS orbital — analytical plane wave  ψ_k(r) = exp(ik·r)/√Ω
# Part: real
# Shell |n|² = 0  k=(0,0,0)
# k = (0, 0, 0) bohr⁻¹
# Predicted KS eigenvalue  ε = k²/2 + V_xc = -0.111602 Ha  (-3.03684 eV)
# Cell Ω = 64000 bohr³,   norm = 1/√Ω = 0.00395285 bohr^{-3/2}
# Slice  z = 20 bohr  (ix=0..79, iy=0..79)
# Re[ψ_k(r)] = cos(k·r) / √Ω
# Columns: ix  iy  x_bohr  y_bohr  psi_value
0  0  0.000000  0.000000  3.95284708e-03
0  1  0.000000  0.500000  3.95284708e-03
0  2  0.000000  1.000000  3.95284708e-03
...
```

**Header lines** (all start with `#`):
- States which orbital this is (state index `N`, shell `n²=M`)
- The k-vector in bohr⁻¹
- The predicted KS eigenvalue (k²/2 + V_xc)
- The cell volume and normalization factor 1/√Ω

**Columns:**
- `ix`, `iy` — grid indices
- `x_bohr`, `y_bohr` — Cartesian coordinates in bohr
- `psi_value` — Re[ψ_k(x, y, z=L/2)] in bohr^(-3/2)

**`orbital_N_n2_M_imag.txt`** is identical in structure but the last column is Im[ψ_k].

**How to reconstruct |ψ|²:**
```python
import numpy as np
re_data = np.loadtxt("orbital_0_n2_0_real.txt", comments='#')
im_data = np.loadtxt("orbital_0_n2_0_imag.txt", comments='#')
psi_squared = re_data[:, 4]**2 + im_data[:, 4]**2   # |ψ|² at each grid point
# Reshape to 2D grid (ix outer, iy inner):
N = 80
psi_2d = psi_squared.reshape(N, N)
```

---

## 2. Jellium Ground State: Eigenvalue File

**Location:** `ResearchProject/jellium/01_ground_state/results/eigenvalues.txt`

Contains the Kohn-Sham eigenvalues for all states, plus the predicted values for validation.

```
# KS eigenvalues — jellium ground state
# V_xc (PZ81) = -0.111602 Ha  (-3.03684 eV)
# Slope-1 line:  eigenvalue_Ha = k2_over_2_Ha + V_xc_Ha
# Columns: state_idx  shell_n2  k2_over_2_Ha  eigenvalue_Ha  predicted_Ha  residual_Ha
0  0  0.00000000  -0.11159537  -0.11160177  0.00000640
1  1  0.01233701  -0.09924821  -0.09926476  0.00001655
2  1  0.01233701  -0.09924942  -0.09926476  0.00001534
3  1  0.01233701  -0.09925692  -0.09926476  0.00000785
...
```

**Header** (lines starting with `#`):
- The V_xc constant (exchange-correlation potential for this r_s value)
- The expected linear relationship: ε_i = k²/2 + V_xc

**Columns (space-separated):**

| Column | Name | What it is |
|---|---|---|
| 1 | `state_idx` | State number, 0-indexed |
| 2 | `shell_n2` | Shell quantum number n² = nx²+ny²+nz² |
| 3 | `k2_over_2_Ha` | Free-electron kinetic energy k²/2 in Hartree |
| 4 | `eigenvalue_Ha` | INQ-computed KS eigenvalue in Hartree |
| 5 | `predicted_Ha` | Analytic prediction k²/2 + V_xc in Hartree |
| 6 | `residual_Ha` | `eigenvalue - predicted` in Hartree (should be < 0.0001 Ha) |

**How to read in Python:**
```python
data = np.loadtxt("eigenvalues.txt", comments='#')
state_idx    = data[:, 0].astype(int)
shell_n2     = data[:, 1].astype(int)
k2_over_2    = data[:, 2]   # Hartree
eigenvalue   = data[:, 3]   # Hartree
predicted    = data[:, 4]   # Hartree
residual     = data[:, 5]   # Hartree (should be ~0)
```

---

## 3. Jellium Convergence: `convergence_results.csv`

**Location:** `ResearchProject/jellium/02_ground_state_convergence/results/convergence_results.csv`

This file contains the results of two convergence tests.
**Important:** The file also contains INQ runtime log output interspersed between data lines,
because INQ prints to stdout and the C++ code also writes to the same file.
Data lines start with `# TEST_A` or `# TEST_B` (yes, the data is in comment-style lines).

```
# Jellium convergence tests
# r_s = 7.25570000 bohr  (N=40, L=40 bohr reference)
# Smearing = 0.00862000 eV  (100 K Fermi-Dirac)
#
# TEST_A spacing_bohr,E_cut_Ha,E_total_Ha,T_s_Ha,E_xc_Ha,n_iter
[2026-04-13 15:04:30.640] [electrons:...] System information:
...   (INQ runtime output) ...
# TEST_A 0.66666700,5.55332826,1.08266017,-1.73484218,-1.02834063,35
[2026-04-13 15:05:12.482] [electrons:...] System information:
...   (INQ runtime output) ...
# TEST_A 0.50000000,9.86960440,0.45456183,-2.22736082,-1.02834063,40
...
# TEST_B N,L_bohr,k0_inv_bohr,Ts_Ha,Ts_per_N,T_TF_per_N,n_iter
# TEST_B 2,5.28835095,1.18937152,0.70730014,0.35365007,0.27011553,18
# TEST_B 14,10.96823788,0.57337617,1.98688219,0.14192016,0.14004099,24
```

**To extract data rows in Python:**
```python
test_a_rows = []
test_b_rows = []
with open("convergence_results.csv") as f:
    for line in f:
        line = line.strip()
        if line.startswith("# TEST_A "):
            vals = line.replace("# TEST_A ", "").split(",")
            test_a_rows.append([float(v) for v in vals])
        elif line.startswith("# TEST_B "):
            vals = line.replace("# TEST_B ", "").split(",")
            test_b_rows.append([float(v) for v in vals])

test_a = np.array(test_a_rows)
test_b = np.array(test_b_rows)
```

**TEST_A columns** (grid spacing convergence):

| Index | Name | Units | What it is |
|---|---|---|---|
| 0 | `spacing_bohr` | bohr | Real-space grid spacing h |
| 1 | `E_cut_Ha` | Hartree | Corresponding plane-wave cutoff (π²/2h²) |
| 2 | `E_total_Ha` | Hartree | Total DFT energy from INQ |
| 3 | `T_s_Ha` | Hartree | Kinetic energy from INQ |
| 4 | `E_xc_Ha` | Hartree | XC energy from INQ |
| 5 | `n_iter` | integer | SCF iterations to convergence |

**TEST_B columns** (shell-closure convergence):

| Index | Name | Units | What it is |
|---|---|---|---|
| 0 | `N` | integer | Number of electrons (closed-shell magic number) |
| 1 | `L_bohr` | bohr | Cell side length (scaled to keep r_s constant) |
| 2 | `k0_inv_bohr` | bohr⁻¹ | k₀ = 2π/L (fundamental BZ vector) |
| 3 | `Ts_Ha` | Hartree | Kinetic energy T_s from INQ |
| 4 | `Ts_per_N` | Hartree | T_s / N (kinetic energy per electron) |
| 5 | `T_TF_per_N` | Hartree | Thomas-Fermi reference value for comparison |
| 6 | `n_iter` | integer | SCF iterations |

---

## 4. Jellium WP Propagation: 2D Slice Files (`slice_tNNN.txt`)

**Location:** `ResearchProject/jellium/03_free_gaussian_wp_propagation/results/`

One file per time snapshot. The number `NNN` in the filename is the snapshot index (000, 001, 002, ...).

These files store a 2D cross-section of the wavepacket density `|ψ|²(x, y, z=L/2)`
at the midplane of the box.

```
# t=0 z_index=80 N=161 dx=0.248447
1.206664e-175 1.683325e-173 2.207710e-171 ... (161 values on one row)
1.683325e-173 2.347862e-171 3.079527e-169 ... (161 values on one row)
...
(161 rows total)
```

**Header** (single line starting with `#`):

| Field | Example | What it means |
|---|---|---|
| `t=` | `0` | Simulation time in atomic time units |
| `z_index=` | `80` | Grid index iz corresponding to z = L/2 |
| `N=` | `161` | Number of grid points per side |
| `dx=` | `0.248447` | Grid spacing in bohr |

**Data:**
- N rows, each with N space-separated floating-point values
- Row `iy` contains the density at y = iy × dx
- Column `ix` in that row contains the density at x = ix × dx
- Values are in bohr⁻³ (electron density)
- Very small values (near 0) appear as `1.2e-175` etc. — this is correct floating-point notation

**How to read in Python:**
```python
import numpy as np

def load_slice(filename):
    """Returns (t, dx, density_2d) where density_2d has shape (N, N)."""
    with open(filename) as f:
        header = f.readline()           # "# t=0 z_index=80 N=161 dx=0.248447"
    parts = header.split()
    t  = float(parts[1].split('=')[1])
    N  = int(parts[3].split('=')[1])
    dx = float(parts[4].split('=')[1])
    data = np.loadtxt(filename, comments='#')   # shape: (N, N)
    return t, dx, data
```

**Physical interpretation:**
At t=0, the packet starts at `z_init = WP_BZ()`, so the slice at z=L/2 shows
nearly zero density (the packet hasn't arrived at the midplane yet).
As time evolves and the packet propagates downward toward z=L/2, this slice
shows the packet cross-section growing.

---

## 5. Jellium WP Propagation: 3D Density Files (`density3d_tNNN.txt`)

**Location:** `ResearchProject/jellium/03_free_gaussian_wp_propagation/results/`

A subsampled 3D snapshot of the wavepacket density `|ψ|²(r)`.

```
# t=2.5 N=161 dx=0.248447 stride=4 NC=41
2.534996e-33
2.225555e-33
1.166894e-33
1.264347e-33
...
(one value per line, NC³ = 41³ = 68921 lines total)
```

**Header** (single `#` line):

| Field | Example | What it means |
|---|---|---|
| `t=` | `2.5` | Simulation time in a.u. |
| `N=` | `161` | Full grid size per side |
| `dx=` | `0.248447` | Grid spacing in bohr |
| `stride=` | `4` | Subsampling factor (every 4th grid point is saved) |
| `NC=` | `41` | Coarsened grid size = ceil(N / stride) |

**Data:**
- One floating-point value per line
- Values are stored in **Fortran order**: z-index varies fastest, then y, then x
  (i.e. iz = 0, 1, 2, ... NC-1 for each (ix, iy), then iy increments, then ix)
- Total values: NC³ = 41³ = 68,921 lines
- Values in bohr⁻³

**How to read in Python:**
```python
def load_density3d(filename):
    """Returns (t, dx_coarse, density_3d) where density_3d has shape (NC, NC, NC)."""
    with open(filename) as f:
        header = f.readline()
    parts = header.split()
    t      = float(parts[1].split('=')[1])
    stride = int(parts[4].split('=')[1])
    NC     = int(parts[5].split('=')[1])
    dx     = float(parts[3].split('=')[1])
    dx_coarse = dx * stride

    values = np.loadtxt(filename, comments='#')   # shape (NC^3,)
    density = values.reshape(NC, NC, NC)           # (ix, iy, iz)
    return t, dx_coarse, density
```

**Why subsampled?**
The full 161³ grid has ~4 million points. Saving it at every snapshot would use too much disk.
The stride=4 subsampling reduces to 41³ = 69k points, which is fast to write and sufficient
for visualization.

---

## 6. Coronene: `sim_summary.txt`

**Location:** `ResearchProject/systems/coronene/04_leed_simulation/results/sim_summary.txt`

A human-readable key-value summary of the simulation parameters.
Written at the start of the run so you know exactly what configuration produced these results.

```
# Coronene WP scattering — TDDFT LEED simulation
# Tsubonoya, Hu, Watanabe PRB 90, 035416 (2014)
GS_energy_Ha     -150.837
E_cut_Ha         40
WP_d_bohr        1.00155
WP_D_bohr        11.9998
WP_k0_bohr_inv   3.83402
WP_Ekin_eV       200
WP_occ           1
dt_au            0.0200092
t1_au            3.12981
t2_au            10.3353
n_steps          516
n_snapshots      11
```

**Format:** Each line is `key<whitespace>value`. No separator character. Lines starting with `#` are comments.

**Fields:**

| Key | Units | What it is |
|---|---|---|
| `GS_energy_Ha` | Hartree | Ground-state total energy before WP injection |
| `E_cut_Ha` | Hartree | Plane-wave energy cutoff used |
| `WP_d_bohr` | bohr | Gaussian width parameter σ |
| `WP_D_bohr` | bohr | Initial z-distance of WP above the molecule |
| `WP_k0_bohr_inv` | bohr⁻¹ | Wavepacket momentum magnitude k₀ = √(2E_kin) |
| `WP_Ekin_eV` | eV | Incident electron kinetic energy |
| `WP_occ` | — | Occupation assigned to WP orbital (usually 1) |
| `dt_au` | a.u. | Time step |
| `t1_au` | a.u. | Time at which LEED accumulation starts |
| `t2_au` | a.u. | Time at which LEED accumulation ends |
| `n_steps` | — | Total number of time steps propagated |
| `n_snapshots` | — | Number of 2D density snapshots saved |

**How to read in Python:**
```python
params = {}
with open("results/sim_summary.txt") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        key, val = line.split()
        try:
            params[key] = float(val)
        except ValueError:
            params[key] = val
```

---

## 7. Coronene: Snapshot Files (`snapshot_tNNNN.txt`)

**Location:** `ResearchProject/systems/coronene/04_leed_simulation/results/`
Also: `runs/run_004_.../results/density_obs_snapshots/`

One file per snapshot. The number `NNNN` is the step index (zero-padded to 4 digits).

These store a 2D cross-section of the total electron density `n(x, y, z_obs)`
at the observation plane z = z_obs, at the current time.

```
# t=0.000000 z=0.000000
2.336634e-02 2.810085e-02 4.423451e-02 7.727056e-02 ...  (Nx values, space-separated)
2.811090e-02 3.316593e-02 5.029041e-02 8.486844e-02 ...
4.426243e-02 5.023216e-02 7.013599e-02 1.087389e-01 ...
...
(Ny rows total)
```

**Header** (single `#` line):

| Field | Example | What it means |
|---|---|---|
| `t=` | `0.000000` | Simulation time in atomic time units |
| `z=` | `0.000000` | z-coordinate of the slice plane in bohr |

Note: For `run_004`, the z value in the header corresponds to the actual physical
z-coordinate of the observation plane (e.g. `z=84.848697` bohr).

**Data:**
- Ny rows, each containing Nx space-separated floating-point values
- Row index = iy (y-direction, 0 to Ny-1)
- Column index = ix (x-direction, 0 to Nx-1)
- Values are in bohr⁻³ (total electron number density)
- Includes both the ground-state coronene density AND the wavepacket density

**How to read in Python:**
```python
def load_snapshot(filename):
    """Returns (t, z, density_2d) where density_2d has shape (Ny, Nx)."""
    with open(filename) as f:
        header = f.readline()   # "# t=0.000000 z=0.000000"
    parts = header.split()
    t = float(parts[1].split('=')[1])
    z = float(parts[2].split('=')[1])
    data = np.loadtxt(filename, comments='#')   # shape (Ny, Nx)
    return t, z, data
```

**Difference from jellium slice files:**
- The coronene snapshot does not include N or dx in the header
  (those are inferred from the array shape and the sim_summary.txt file)
- The coronene density is much larger (because there are real molecular electrons),
  typically values of order 10⁻² to 10⁻¹ near the molecule

---

## 8. Coronene: LEED Pattern (`leed_pattern.txt`)

**Location:** `ResearchProject/systems/coronene/04_leed_simulation/results/leed_pattern.txt`

The time-integrated electron density at the observation plane —
this is the simulated LEED diffraction pattern.

```
# LEED pattern I(x,y) = integral_{t1}^{t2} n(x,y,z=D,t) dt
# z_obs=11.9998 bohr  t1=3.12981 a.u.  t2=10.3353 a.u.
# Rows: iy = 0..99  Cols: ix = 0..99
3.340575e-04 3.052934e-04 2.412094e-04 1.671907e-04 ...  (100 values)
3.059801e-04 2.813893e-04 2.244636e-04 1.552986e-04 ...
...
(100 rows total)
```

**Header** (3 lines, all starting with `#`):
- Line 1: what the quantity is: I(x,y) = ∫ n(x,y,z_obs,t) dt
- Line 2: `z_obs` = observation plane height (bohr), `t1` and `t2` = integration window (a.u.)
- Line 3: dimensions of the array

**Data:**
- Ny rows × Nx columns of space-separated floats
- Row iy, column ix → I(x_ix, y_iy)
- Units: electron density × time = bohr⁻³ × a.u.

**How the LEED pattern is accumulated in C++:**
```cpp
// In the callback, every step between t1 and t2:
if (t >= t1_au && t <= t2_au) {
    auto slice = leed_utils::extract_density_slice(electrons, cfg::Z_OBS_BOHR());
    for (int iy = 0; iy < Ny; iy++)
        for (int ix = 0; ix < Nx; ix++)
            leed_pattern[iy][ix] += slice[iy][ix] * dt_au;
}
// Written once at the end of propagation.
```

**How to read and plot in Python:**
```python
def load_leed_pattern(filename):
    """Returns (z_obs, t1, t2, pattern) where pattern has shape (Ny, Nx)."""
    with open(filename) as f:
        f.readline()   # skip first comment
        line2 = f.readline()
        f.readline()   # skip dimensions comment
    # Parse line2: "# z_obs=11.9998 bohr  t1=3.12981 a.u.  t2=10.3353 a.u."
    parts = line2.split()
    z_obs = float(parts[1].split('=')[1])
    t1    = float(parts[3].split('=')[1])
    t2    = float(parts[6].split('=')[1])
    pattern = np.loadtxt(filename, comments='#')   # shape (Ny, Nx)
    return z_obs, t1, t2, pattern

import matplotlib.pyplot as plt
z_obs, t1, t2, I = load_leed_pattern("results/leed_pattern.txt")
plt.imshow(I, origin='lower', cmap='hot')
plt.colorbar(label='I(x,y) [bohr⁻³ × a.u.]')
plt.title('Simulated LEED pattern')
plt.show()
```

**Physical interpretation:**
The LEED pattern shows where electrons arrive at the detector plane after scattering from the
coronene molecule. High intensity regions correspond to constructive interference (diffraction spots).
The integration window [t1, t2] is chosen so that the incident packet arrives during this window
and the scattered electrons reach the detector plane.

---

## Summary Table

| File | Written by | Dimensions | Format | Values |
|---|---|---|---|---|
| `grid_slice.txt` | jellium 01 | 6400 rows × 4 col | text, space-sep | ix, iy, x[bohr], y[bohr] |
| `orbital_N_n2_M_real/imag.txt` | jellium 01 | 6400 rows × 5 col | text, space-sep | ix, iy, x, y, ψ[bohr^-3/2] |
| `eigenvalues.txt` | jellium 01 | N_states rows × 6 col | text, space-sep | state, n², k²/2, ε, predicted, residual [Ha] |
| `convergence_results.csv` | jellium 02 | comment-rows | text, comma-sep in #-comments | spacing, E_cut, E_total, T_s, E_xc, n_iter |
| `slice_tNNN.txt` | jellium 03 | N×N array | text, space-sep rows | |ψ|²[bohr⁻³] |
| `density3d_tNNN.txt` | jellium 03 | NC³ values | text, one/line | |ψ|²[bohr⁻³], subsampled |
| `sim_summary.txt` | coronene 04 | ~12 rows | text, key-value | mixed (see table) |
| `snapshot_tNNNN.txt` | coronene 04 | Ny×Nx array | text, space-sep rows | n(r)[bohr⁻³] |
| `leed_pattern.txt` | coronene 04 | Ny×Nx array | text, space-sep rows | I(x,y)[bohr⁻³·a.u.] |
