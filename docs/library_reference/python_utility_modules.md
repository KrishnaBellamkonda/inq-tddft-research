# Python Utility Modules Reference

This covers the reusable Python modules in `Tutorial/angelo-jellium/jellium/`.
These are not stand-alone scripts — they are imported by the analysis scripts.

---

## Table of Contents

1. [`inq_io.py` — Reading INQ binary files and logs](#1-inq_iopy--reading-inq-binary-files-and-logs)
2. [`wavepacket.py` — Gaussian wavepacket math](#2-wavepacketpy--gaussian-wavepacket-math)
3. [`config.py` — Centralised system parameters](#3-configpy--centralised-system-parameters)

---

## 1. `inq_io.py` — Reading INQ Binary Files and Logs

**Location:** `Tutorial/angelo-jellium/jellium/inq_io.py`

**Purpose:**
INQ saves electronic states and densities as raw binary files with no header.
This module knows how to find those files in the INQ working directory,
read them into NumPy arrays, and parse INQ's text outputs (eigenvalue logs, energy tables, time series).

**Import:**
```python
from jellium.inq_io import read_orbital, read_density, read_time_series, ...
```

---

### Background: How INQ stores its state on disk

When INQ runs, it creates a directory structure under the current working directory:

```
default_ions/
    cell              ← cell side length (plain text)
default_electrons_options/
    spacing           ← real-space grid spacing (plain text)
default_orbitals/
    spin_density/
        0000000000    ← ground-state density (binary, N³ float64)
kpin0000000000/
    states/
        0000000001_0000000000   ← orbital 1 (binary, N³ complex128)
        0000000002_0000000000   ← orbital 2
        ...
default_checkpoint/
    real-time/
        orbitals/
            spin_density/
                0000000000      ← density at last checkpoint step
```

All binary files are raw C-contiguous arrays (no header, no shape information).
You must know N (the grid size) before reading. N = round(cell_side / spacing).

---

### `read_grid_params(inq_dir)`

```python
cell_side, spacing, N = read_grid_params("/path/to/inq/workdir")
```

Reads `cell_side` from `default_ions/cell` and `spacing` from `default_electrons_options/spacing`.
Computes N = round(cell_side / spacing).

**Returns:** `(cell_side: float, spacing: float, N: int)` all in bohr.

---

### `grid_N_from_config(cell_side, spacing)`

```python
N = grid_N_from_config(40.0, 0.25)   # → 160
```

Simple helper when you already know cell_side and spacing from your config.
Returns `int(round(cell_side / spacing))`.

---

### `read_orbital(path, N)`

```python
psi = read_orbital("kpin0000000000/states/0000000001_0000000000", N=161)
# psi.shape → (161, 161, 161), dtype complex128
```

Reads a single KS orbital from a binary file.

- `path` — path to the orbital binary file.
- `N` — grid points per dimension (same in all three directions).

**Returns:** complex128 array of shape `(N, N, N)`.

**Array layout:** `psi[ix, iy, iz]` = ψ(x_ix, y_iy, z_iz).
Axes: ix=x (outermost), iz=z (fastest-varying in memory).

Raises `ValueError` if the file size doesn't match N³ complex128 values.

---

### `read_density(path, N)`

```python
n = read_density("default_orbitals/spin_density/0000000000", N=161)
# n.shape → (161, 161, 161), dtype float64
```

Same as `read_orbital` but reads a real-valued density field (float64, not complex128).

- `path` — binary file path.
- `N` — grid size.

**Returns:** float64 array of shape `(N, N, N)`.

---

### `read_all_orbitals(inq_dir, N)`

```python
orbitals = read_all_orbitals("/path/to/inq/workdir", N=161)
# orbitals is a list of (161, 161, 161) complex128 arrays, one per orbital
```

Reads all KS orbitals from `kpin0000000000/states/`.
Files are sorted numerically, so `orbitals[0]` is state 0, `orbitals[1]` is state 1, etc.

**Returns:** List of complex128 arrays, each `(N, N, N)`.

---

### `read_spin_density(inq_dir, N)`

```python
density = read_spin_density("/path/to/inq/workdir", N=161)
# density.shape → (161, 161, 161), dtype float64
```

Convenience wrapper for reading the ground-state density from its standard location:
`default_orbitals/spin_density/0000000000`.

---

### `read_checkpoint_density(inq_dir, N)`

```python
density = read_checkpoint_density("/path/to/inq/workdir", N=161)
```

Reads the density stored at the last real-time checkpoint:
`default_checkpoint/real-time/orbitals/spin_density/0000000000`.

Use this to inspect the density mid-propagation if INQ was interrupted.

---

### `make_coords(cell_side, N)`

```python
x = make_coords(40.0, 161)   # → array of 161 values from 0 to ~39.75 bohr
```

Creates a 1D coordinate array: `np.linspace(0, cell_side, N, endpoint=False)`.

Used to build the physical axis when plotting slices:
```python
x = make_coords(cell_side, N)
y = make_coords(cell_side, N)
X, Y = np.meshgrid(x, y, indexing='ij')
```

---

### `read_time_series(path)`

```python
ts = read_time_series("run.log")
# ts is a dict: ts['time'] = array, ts['energy'] = array, etc.
```

Parses INQ time-series output files (tab-separated, `#`-commented header).

The header line looks like:
```
# time[atu]    energy[Ha]    Jx[au]    Jy[au]    Jz[au]
```

The function strips the unit strings in brackets and returns a dict mapping
column names to arrays.

**Returns:** Dict with column names as keys, 1D float64 arrays as values.
Also includes `'col0'`, `'col1'`, etc. as fallback integer-indexed access.

---

### `parse_eigenvalues_from_log(log_path)`

```python
evals = parse_eigenvalues_from_log("run.log")
# evals = [(state_idx, occupation, eigenvalue_Ha), ...]
```

Parses the final converged KS eigenvalues from an INQ log file.
INQ prints eigenvalue tables at each SCF iteration in the format:
```
st =    1  occ = 2.000  evalue =    -0.123456789012  res = 2e-09
```

The function finds the block after `"SCF ended after"` to get only the final converged values.
If that marker is absent (e.g. interrupted run), it returns the last complete block found.

**Returns:** List of `(state_index: int, occupation: float, eigenvalue_Ha: float)` tuples,
sorted by eigenvalue ascending. Empty list if no eigenvalues found.

---

### `write_eigenvalues(eigenvalues, path)`

```python
write_eigenvalues(evals, "results/eigenvalues.tsv")
```

Writes the parsed eigenvalues to a TSV file with this format:
```
# state	occupation	eigenvalue_Ha	eigenvalue_eV
1	2.0000	-0.123456789012	-3.35914
2	2.0000	-0.098765432109	-2.68754
```

- `ha_to_ev` argument: conversion factor (default 27.211386, NIST 2018 value).

---

### `read_eigenvalues(path)`

```python
evals = read_eigenvalues("results/eigenvalues.tsv")
# evals['eigenvalue_Ha'] → float64 array
# evals['state']         → int array
```

Reads a TSV written by `write_eigenvalues()`.

**Returns:** Dict with keys `'state'`, `'occupation'`, `'eigenvalue_Ha'`, `'eigenvalue_eV'`,
each a 1D NumPy array.

---

### `read_energy_components(path)`

```python
E = read_energy_components("energy_output.txt")
# E['total']   → float (Hartree)
# E['kinetic'] → float
# E['xc']      → float
```

Parses a text file containing INQ's energy component table:
```
Energy:
  total          =      -0.053015780605 Ha
  kinetic        =       0.000000000572 Ha
  ...
```

**Returns:** Dict mapping component name to value in Hartree.

---

## 2. `wavepacket.py` — Gaussian Wavepacket Math

**Location:** `Tutorial/angelo-jellium/jellium/wavepacket.py`

**Purpose:**
Provides a Python implementation of the 3D Gaussian wavepacket,
matching the C++ `inject_wp` function. Used to:
1. Test the wavepacket mathematically before injecting into INQ
2. Do analytical free-particle propagation (FFT-based) for comparison with INQ
3. Verify normalization and kinetic energy

**Import:**
```python
from jellium.wavepacket import GaussianWavePacket3D, propagate_free, evaluate_on_grid, ...
```

---

### `GaussianWavePacket3D` — dataclass

```python
from jellium.wavepacket import GaussianWavePacket3D
import numpy as np

wp = GaussianWavePacket3D(
    r0    = np.array([20.0, 20.0, 25.0]),  # centre in bohr
    sigma = 1.0,                            # width in bohr
    k0    = np.array([0.0, 0.0, -3.83]),   # momentum in bohr^-1 (pointing -z)
)
```

A dataclass that stores all wavepacket parameters.

**Fields:**
- `r0` — shape `(3,)` float64 array, centre position in bohr.
- `sigma` — float, Gaussian width (half-width parameter). Units: bohr.
- `k0` — shape `(3,)` float64 array, momentum vector in bohr⁻¹.

**Mathematical definition:**
```
ψ(r) = N * exp(-|r - r0|² / (4σ²)) * exp(i k0 · (r - r0))
```
where N = (2πσ²)^(-3/4) ensures ∫|ψ|² d³r = 1.

**Note on convention:** This uses `4σ²` in the exponent, not `2σ²`.
This is because σ here is the half-width parameter of the probability density
|ψ|² = N² exp(-|r-r0|²/(2σ²)), so σ is the standard deviation of the Gaussian density.

---

### `gaussian_norm_3d(sigma)`

```python
N_norm = gaussian_norm_3d(1.0)   # → (2π)^(-3/4) ≈ 0.4244
```

Returns the normalization constant N = (2πσ²)^(-3/4).

---

### `psi_gaussian_3d(r, wp)`

```python
# Evaluate at a single point:
r_point = np.array([20.0, 20.0, 20.0])
psi_val = psi_gaussian_3d(r_point[np.newaxis, :], wp)

# Evaluate on a 2D slice:
X, Y = np.meshgrid(x_coords, y_coords, indexing='ij')
Z_fixed = 20.0 * np.ones_like(X)
r_slice = np.stack([X, Y, Z_fixed], axis=-1)   # shape (Nx, Ny, 3)
psi_slice = psi_gaussian_3d(r_slice, wp)         # shape (Nx, Ny), complex
```

Evaluates ψ(r) at one or many points.

- `r` — array of shape `(..., 3)`, any number of leading batch dimensions.
- `wp` — `GaussianWavePacket3D` instance.

**Returns:** Complex array of shape `(...)` (leading dimensions only, no 3-component dimension).

---

### `make_grid(cell_side, spacing)`

```python
X, Y, Z, dx = make_grid(40.0, 0.248447)
# X, Y, Z each have shape (161, 161, 161)
# dx ≈ 0.248447
```

Creates a 3D uniform grid matching INQ's grid conventions.
Grid points run from 0 to cell_side (exclusive endpoint), with `N = round(cell_side/spacing)` points.

**Returns:** `(X, Y, Z, dx)` where:
- `X`, `Y`, `Z` — 3D meshgrid arrays of shape `(N, N, N)` containing the coordinate at each point.
- `dx` — actual grid spacing (slightly different from requested if rounding occurred).

**INQ convention:** Grid origin at 0, no endpoint, `np.linspace(0, L, N, endpoint=False)`.

---

### `evaluate_on_grid(wp, cell_side, spacing)`

```python
psi, dx = evaluate_on_grid(wp, 40.0, 0.248447)
# psi.shape → (161, 161, 161), complex128
```

Evaluates the full 3D wavepacket on the grid.
Calls `make_grid` + `psi_gaussian_3d` internally.

**Returns:** `(psi, dx)` where `psi` has shape `(N, N, N)`.

---

### `check_normalization(psi, dx)`

```python
norm = check_normalization(psi, dx)
print(f"Norm: {norm:.6f}")   # should print 1.000000 ± small error
```

Computes ∫|ψ|² dV ≈ Σ |ψ_ijk|² × dx³.

**Returns:** Float (should be close to 1.0 for a well-normalised packet that fits in the box).

---

### `propagate_free(wp, cell_side, spacing, t)`

```python
psi_t = propagate_free(wp, 40.0, 0.248447, t=2.5)
# psi_t.shape → (161, 161, 161), complex128
```

Free-particle time propagation using the split-operator FFT method.

**Algorithm:**
1. Evaluate ψ(r, t=0) on the grid
2. FFT to k-space: ψ̃(k)
3. Apply free-particle propagator: ψ̃(k) × exp(-i k² t / 2)
4. IFFT back to real space

This is **exact** for a free particle (no potential) with periodic boundaries.
Periodic boundaries mean the wavepacket wraps around at the box edge.

**Arguments:**
- `cell_side` — float, box size in bohr
- `spacing` — float, grid spacing in bohr
- `t` — float, time in atomic time units

**Returns:** Complex128 array of shape `(N, N, N)`.

**Use case:** Run this for a range of t values and compare σ(t) against the
analytical formula σ(t) = σ₀ √(1 + (t/2mσ₀²)²). This validates the INQ propagation.

---

### `write_wavepacket_params(wp, path)`

```python
write_wavepacket_params(wp, "results/wavepacket_params.txt")
```

Writes all wavepacket parameters to a text file:
```
# GaussianWavePacket3D parameters
r0_x    20.000000
r0_y    20.000000
r0_z    25.000000
sigma   1.000000
k0_x    0.000000
k0_y    0.000000
k0_z    -3.830000
norm    4.244131e-01
```

---

### `write_wavepacket_slice(psi, cell_side, spacing, path, axis=2)`

```python
write_wavepacket_slice(psi, 40.0, 0.248447, "results/slice_z_mid.tsv", axis=2)
```

Writes a 2D slice of |ψ|² through the centre of the box to a TSV file.
Slices through `axis=2` (z, default) at `iz = N//2`.
Slices through `axis=0` (x) or `axis=1` (y) if specified.

**Output format:**
```
# |psi|^2 slice at axis=2, index=80
# coord1	coord2	|psi|^2
0.000000	0.000000	1.234567e-08
0.000000	0.248447	3.456789e-07
...
```

Three columns: coord1, coord2, |ψ|² in bohr⁻³.

---

## 3. `config.py` — Centralised System Parameters

**Location:** `Tutorial/angelo-jellium/jellium/config.py`

**Purpose:**
Central configuration file for the angelo-jellium tutorial jellium simulations.
Import it to get all system parameters without hardcoding values in each script.

**Import:**
```python
from jellium.config import L, SPACING, N_ELECTRONS, ...
```

**Key parameters:**

```python
# Cell
L = 40.0          # box side length in bohr (a₀)
SPACING = 0.25    # grid spacing in bohr — gives N = 160 grid points per side

# Electrons
N_ELECTRONS = 40  # number of free electrons
SMEAR_EV = 0.0862 # Fermi-Dirac smearing = k_B × 1000 K, in eV

# Wavepacket (if used)
SIGMA_WP = 1.0    # Gaussian width in bohr
K0 = [0.0, 0.0, 1.5]   # initial momentum in bohr^-1
R0 = [L/2, L/2, L/2]   # initial centre position in bohr

# Time propagation
DT = 0.04         # time step in a.u.
T_MAX = 12.0      # total propagation time in a.u.

# Output
SNAPSHOT_EVERY = 25   # save a snapshot every N steps
```

These values match what was used in the C++ `run.cpp` files in the same tutorials.
