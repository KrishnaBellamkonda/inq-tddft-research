# Python Analysis Scripts Reference

This document covers every standalone Python analysis script.
For each script: what role it plays in the simulation workflow,
what files it reads, what it outputs, and what functions it defines.

---

## Table of Contents

**Jellium (ResearchProject/jellium/)**
1. [`01_ground_state/plot_results.py`](#1-01_ground_stateplot_resultspy)
2. [`01_ground_state/plot_shell_structure.py`](#2-01_ground_stateplot_shell_structurepy)
3. [`02_ground_state_convergence/plot_convergence.py`](#3-02_ground_state_convergenceplot_convergencepy)
4. [`03_free_gaussian_wp_propagation/plot_propagation.py`](#4-03_free_gaussian_wp_propagationplot_propagationpy)

**Coronene (ResearchProject/systems/coronene/)**
5. [`01_geometry/gen_geometry.py`](#5-01_geometrygen_geometrypy)
6. [`03_ecut_convergence/plot_convergence.py`](#6-03_ecut_convergenceplot_convergencepy)
7. [`04_leed_simulation/analysis.py` (base run)](#7-04_leed_simulationanalysispy-base-run)
8. [`04_leed_simulation/geometry_check.py`](#8-04_leed_simulationgeometry_checkpy)
9. [`runs/run_002/analysis.py`](#9-runsrun_002analysispy)

**Jellium system videos (ResearchProject/systems/jellium/)**
10. [`analysis/make_video.py`](#10-analysismake_videopy)
11. [`analysis/plot_spreading.py`](#11-analysisplot_spreadingpy)

**Angelo jellium tutorial (Tutorial/angelo-jellium/scripts/)**
12. [Script overview](#12-angelo-jellium-scripts-overview)

---

## 1. `01_ground_state/plot_results.py`

**Role in workflow:** Run after `inq-run` in the jellium ground state directory.
Validates the ground state by checking that INQ's eigenvalues match the
analytic free-electron prediction (k²/2 + V_xc).
Produces figures suitable for a report.

**Input files:**
- `results/eigenvalues.txt` — KS eigenvalue table (see output format doc)
- `results/orbitals/orbital_*_real.txt` and `*_imag.txt` — orbital slices

**Output files:**
- `results/shell_structure.png` — side-by-side energy level diagram
- `results/xc_offset.png` — eigenvalue scatter plot
- `results/orbitals.png` — 2D orbital visualizations

**Constants (hardcoded at top of file):**
```python
N  = 40       # number of electrons
L  = 40.0     # cell side in bohr
HA_TO_EV = 27.211386245988
```

All analytic quantities (n₀, r_s, k_F, E_F, k₀, V_xc) are derived from these.

---

### `vxc_pz81(rs)`

```python
def vxc_pz81(rs):
    """Perdew-Zunger 1981 LDA XC potential.
    
    Args:
        rs: Wigner-Seitz radius in bohr (float or array)
    
    Returns:
        V_xc in Hartree
    """
```

Computes exchange + correlation potential: V_x = (4/3)×(-0.4582/r_s) plus Broyden
correlation potential using the r_s ≥ 1 (metallic density) branch.

**Used for:** Computing the expected eigenvalue shift `Vxc_Ha` to compare against INQ.

---

### `compute_shells(L, N, n2_max=6)`

```python
def compute_shells(L, N, n2_max=6):
    """Enumerate free-electron shells and their occupancies.
    
    Args:
        L: cell side in bohr
        N: number of electrons to fill
        n2_max: maximum |n|² shell to enumerate
    
    Returns:
        List of dicts:
        {
            'n2':       int,    # |n|² quantum number
            'deg':      int,    # spatial degeneracy (k-states per spin)
            'Ek_Ha':    float,  # kinetic energy = k²/2 in Hartree
            'electrons':int,    # electrons actually placed in this shell
            'cap':      int,    # maximum electrons this shell can hold (2*deg)
            'frac':     float,  # fill fraction = electrons/cap
        }
    """
```

The shell structure enumerates all (nx, ny, nz) integer triplets with nx²+ny²+nz² ≤ n2_max
and groups them by |n|². States are filled from lowest energy up until all N electrons
are placed.

**Example output for N=40, n2_max=6:**

| n² | deg | cap | electrons | frac |
|---|---|---|---|---|
| 0  | 1  | 2  | 2  | 1.0 (full) |
| 1  | 6  | 12 | 12 | 1.0 (full) |
| 2  | 12 | 24 | 24 | 1.0 (full) |
| 3  | 8  | 16 | 2  | 0.125 (partial) |

---

### `shell_style(frac)`

```python
def shell_style(frac):
    """Return (color, alpha) for a shell based on its fill fraction."""
    if frac >= 0.999:   return '#2166ac', 0.90   # blue = fully occupied
    elif frac > 0.001:  return '#e6711f', 0.90   # orange = partially occupied
    else:               return '#aaaaaa', 0.55   # gray = empty
```

---

### `plot_shell_structure()`

Creates Figure 1: a side-by-side energy level diagram.
- **Left panel:** Free-electron levels ε_k = k²/2 (no XC shift)
- **Right panel:** KS levels ε_k = k²/2 + V_xc (shifted uniformly down)

Each energy level is drawn as a horizontal line, color-coded by fill fraction
(blue=full, orange=partial, gray=empty). Line width scales with degeneracy.

The right panel includes an annotation showing the V_xc shift with a double-headed arrow.

Saves to: `results/shell_structure.png`

---

### `load_orbital(orbital_idx)`

```python
def load_orbital(orbital_idx):
    """Load a single orbital's Re/Im slices.
    
    Globs for results/orbitals/orbital_{idx}_n2_*_real.txt
    and the matching _imag.txt.
    
    Returns:
        (Re_psi, Im_psi, density, x_vals, n2)
        where Re_psi, Im_psi, density are (N_g, N_g) arrays
        and x_vals is a (N_g,) array of x coordinates in bohr.
        Returns (None,)*5 if files not found.
    """
```

**How it parses the filename to get n²:**
```python
fname_re = "results/orbitals/orbital_0_n2_0_real.txt"
n2 = int(fname_re.split('_n2_')[1].split('_real')[0])
# → 0
```

**How it reshapes the data:**
The file has N_g² rows. Column 4 is the wavefunction value.
Column 0 gives ix, which goes from 0 to N_g-1 (the slow axis).
So: `Re_psi = raw_re[:, 4].reshape(N_g, N_g)`.

---

### `plot_xc_offset()`

Creates Figure 2: eigenvalue scatter plot.

Reads `results/eigenvalues.txt` (see output formats doc for column layout).
Plots ε_i (INQ eigenvalue) vs k²/2 (analytic kinetic energy),
coloured by shell n².

Overlays the reference line ε = k²/2 + V_xc.
Fits a linear model to get slope (should be 1.000) and intercept (should match V_xc).
The fit result is printed in a text box on the figure.

**Bottom panel:** Residuals (ε_i - (k²/2 + V_xc)) in milliHartree.

Saves to: `results/xc_offset.png`

---

### `plot_orbitals()`

Creates Figure 3: 2D orbital visualizations.

Loops over orbital files found in `results/orbitals/orbital_*_real.txt`.
Creates a 3-row × N_orbs-column figure:
- Row 0: Re[ψ_k] (colormap: RdBu_r, diverging, centred at 0)
- Row 1: Im[ψ_k] (colormap: PuOr_r, diverging)
- Row 2: |ψ_k|² = Re² + Im² (colormap: viridis)

For plane waves, |ψ_k|² is always exactly 1/Ω everywhere, so the density panel is nearly uniform.
The interesting panels are Re[ψ] and Im[ψ] which show the oscillation pattern of each plane wave.

Saves to: `results/orbitals.png`

---

## 2. `01_ground_state/plot_shell_structure.py`

**Role in workflow:** Standalone script for generating just the shell structure diagram,
independent of a simulation run. Uses only analytic formulas (no data files needed).

**Input:** None (everything is hardcoded or computed analytically)

**Output:** `results/shell_structure.png`

Essentially a simplified version of `plot_shell_structure()` from `plot_results.py`.
Includes the PZ81 V_xc formula directly.

---

## 3. `02_ground_state_convergence/plot_convergence.py`

**Role in workflow:** Run after `inq-run` in the convergence directory.
Reads the convergence CSV and generates two figures showing:
1. Total energy vs energy cutoff (E_cut convergence)
2. Kinetic energy per electron vs N (shell convergence towards TF limit)

**Input files:**
- `results/convergence_results.csv` — the convergence data

(See output formats doc for the unusual format of this file — data is in `# TEST_A` comment lines.)

**Output files:**
- `results/convergence_Ecut.png`
- `results/convergence_shells.png`

---

**How to parse the file:**
```python
test_a, test_b = [], []
with open("results/convergence_results.csv") as f:
    for line in f:
        if line.startswith("# TEST_A "):
            vals = line.replace("# TEST_A ", "").split(",")
            test_a.append([float(v) for v in vals])
        elif line.startswith("# TEST_B "):
            vals = line.replace("# TEST_B ", "").split(",")
            test_b.append([float(v) for v in vals])
test_a = np.array(test_a)   # columns: spacing, E_cut, E_total, T_s, E_xc, n_iter
test_b = np.array(test_b)   # columns: N, L, k0, Ts, Ts_per_N, T_TF_per_N, n_iter
```

---

**`convergence_Ecut.png`:** Two-panel figure.
- Top: E_total(Ha) and T_s(Ha) vs E_cut(Ry). Shows how energy converges.
- Bottom: ΔE_total (difference from the finest-grid value) on a log scale.
  Reference horizontal lines at 1 mHa and 10 mHa tolerance.

**`convergence_shells.png`:** Single panel.
- x-axis: N (number of electrons, discrete shell-closure values)
- y-axis: T_s/N (kinetic energy per electron)
- Data points: INQ values
- Reference line: Thomas-Fermi bulk limit T_TF/N = (3/10)(3π²n₀)^(2/3)

---

## 4. `03_free_gaussian_wp_propagation/plot_propagation.py`

**Role in workflow:** Run after the jellium free wavepacket propagation.
Reads the slice snapshots and the spreading data to visualize and validate
free-particle propagation.

**Input files:**
- `results/slice_tNNN.txt` — 2D density snapshots
- `results/density3d_tNNN.txt` — 3D density files
- `results/grid_info.txt` — simulation metadata (L, E_cut, dt, etc.)

**Output:** Animation frames or figures of wavepacket spreading.

**Typical functions:**
```python
def load_slice(filename):
    """Parse header (t, z_index, N, dx) and data array from slice_t*.txt."""
    # Returns: (t, dx, density_2d) where density_2d.shape = (N, N)

def compute_sigma(density_2d, dx):
    """RMS width σ = √(⟨x²⟩ - ⟨x⟩²) of the 2D density distribution."""
    # Returns: (sigma_x, sigma_y)

def sigma_analytic(t, sigma0):
    """Free-particle spreading: σ(t) = σ0 √(1 + (t/(2σ0²))²)."""
    # Valid for unit-mass particle (m=1 a.u.)
    return sigma0 * np.sqrt(1 + (t / (2 * sigma0**2))**2)

def plot_spreading(t_vals, sigma_num, sigma0):
    """Plot σ(t): numerical (squares) vs analytic (line)."""
    sigma_ana = sigma_analytic(np.array(t_vals), sigma0)
    plt.plot(t_vals, sigma_num, 'o', label='INQ')
    plt.plot(t_vals, sigma_ana, '-', label='Analytic')
```

---

## 5. `01_geometry/gen_geometry.py`

**Role in workflow:** Run once to generate the coronene XYZ geometry file.
Creates the atomic coordinates for coronene (C₂₄H₁₂) from the known bond lengths
and D6h symmetry, then writes a centered `.xyz` file.

**Input:** None (hardcoded geometry parameters)

**Output:** `coronene.xyz` or `coronene_centered.xyz`

**Typical operations:**
- Generates C positions using 6-fold rotational symmetry
- Appends H positions at each peripheral C
- Translates coordinates so the molecule is centred at the cell centre
- Writes in standard XYZ format (element, x, y, z in Angstrom)

---

## 6. `03_ecut_convergence/plot_convergence.py`

**Role in workflow:** Run after the coronene E_cut convergence study.
Reads the convergence CSV and plots E_total vs E_cut to identify the
optimal cutoff for the LEED simulations (typically 40 Ha).

**Input:** `results/ecut_convergence.csv`

**CSV format:**
```
# E_cut_Ha,E_total_Ha,E_total_eV,grid_points,scf_steps
30.0,-150.234,...
35.0,-150.789,...
40.0,-150.837,...
```

**Output:**
- `results/ecut_convergence.png` — E_total vs E_cut with convergence annotations

---

## 7. `04_leed_simulation/analysis.py` (base run)

**Role in workflow:** Post-processes the base LEED simulation (the simplest version).
Produces the two key figures: density evolution (Fig. 1 analogue) and LEED pattern (Fig. 2 analogue)
from Tsubonoya, Hu, Watanabe PRB 90, 035416 (2014).

**Input files:**
- `results/sim_summary.txt` — run metadata
- `results/snapshot_t*.txt` — 2D density slices at the coronene plane
- `results/leed_pattern.txt` — the accumulated LEED pattern

**Output:**
- `results/density_evolution.png` — grid of N_snapshot panels showing n(x,y,z_flake,t)
- `results/leed_pattern.png` — I(x,y) in linear + log scale

---

## 8. `04_leed_simulation/geometry_check.py`

**Role in workflow:** Quality check script. Run after generating or modifying
the coronene XYZ file, before running a simulation.

**What it checks:**
1. All 36 atoms are within the cell boundaries [0, Lx) × [0, Ly) × [0, Lz)
2. Bond lengths: C-C ≈ 1.421 Å (±0.05 Å), C-H ≈ 1.086 Å (±0.05 Å)
3. Atom count: exactly 24 C + 12 H
4. Molecule centred near (Lx/2, Ly/2, Lz/2)

**Input:** `geometry/coronene_centered.xyz` (or hardcoded path)

**Output:** Pass/fail messages to stdout. Non-zero exit code if checks fail.

---

## 9. `runs/run_002/analysis.py`

**Role in workflow:** The most complete analysis script. Processes the run_002 extended LEED
simulation and generates 9 figure files covering energy conservation, orbital dynamics,
wavepacket trajectory, density evolution, and the LEED pattern.

**Input files:**
- `results/energy/energy_vs_time.csv`
- `results/ks_overlaps/projected_occ_vs_time.csv`
- `results/wp_orbital/wp_slice_t*.txt` (WP orbital 2D slices)
- `results/wp_trajectory/density_z_profile_vs_time.csv`
- `results/density_snapshots/snapshot_t*.txt`
- `results/leed_pattern/leed_pattern.txt`

**Output files:**
- `results/energy/total_energy.png`
- `results/energy/kinetic_energy.png`
- `results/energy/all_energies.png`
- `results/energy/energy_fft.png`
- `results/ks_overlaps/projected_occ_heatmap.png`
- `results/wp_orbital/wp_orbital_snapshots.png`
- `results/wp_trajectory/density_trajectory.png`
- `results/density_snapshots/fig1_density_snapshots.png`
- `results/leed_pattern/fig2_leed_pattern.png`
- `results/leed_pattern/fig2_leed_reciprocal.png`

**Unit conversions used:**
```python
AU_TO_FS = 0.024188843   # 1 a.u. = 0.02419 fs
HA_TO_EV = 27.21138625
```

---

### `savefig(path, **kw)`

```python
def savefig(path, **kw):
    """Save current figure to path, creating intermediate directories."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    plt.savefig(path, dpi=150, bbox_inches="tight", **kw)
    plt.close()
```

Used throughout to save figures without needing `os.makedirs` calls everywhere.

---

### `load_energy(path)`

```python
data = load_energy("results/energy/energy_vs_time.csv")
# data is a structured NumPy array with named columns:
# data['step'], data['t_au'], data['E_total'], data['E_kinetic'],
# data['E_hartree'], data['E_xc'], data['E_external'],
# data['E_nonlocal'], data['E_ion']
```

Reads the energy CSV using `np.genfromtxt(path, delimiter=",", comments="#", names=[...])`.
The `names` list matches the CSV column order exactly.

Returns `None` if the file doesn't exist (all plot functions skip gracefully).

---

### `plot_total_energy(data)`

Plots E_total vs time in femtoseconds.
Computes total drift = `(E[-1] - E[0]) * 1e6` in μHa and includes it in the title.

---

### `plot_all_energies(data)`

3×2 subplot grid showing all 6 energy components (kinetic, Hartree, XC, external, non-local, ion-ion)
vs time, each in its own panel.

---

### `plot_energy_fft(data)`

Computes the FFT of E_total fluctuations (after subtracting a linear trend).
Plots |FFT| vs frequency in eV on a log scale.
The plasmon frequency or other resonances would appear as peaks.

**How it converts frequency:**
```python
freqs = np.fft.rfftfreq(N, d=dt)      # frequencies in 1/a.u.
omega_eV = freqs * HA_TO_EV * (2*np.pi)  # angular freq in a.u., converted to eV
```

---

### `load_slice(path)`

```python
t_au, z_bohr, arr = load_slice("results/density_snapshots/snapshot_t0003.txt")
# t_au: simulation time in a.u.
# z_bohr: z-coordinate of the slice plane in bohr
# arr: 2D NumPy array of shape (Ny, Nx), density in bohr^-3
```

Parses the header line `# t=0.000000 z=0.000000` and loads the data array with `np.loadtxt`.

Used by both `plot_wp_orbital_snapshots()` and `plot_density_snapshots()`.

---

### `plot_ks_overlap_heatmap(path)`

Reads the projected occupation CSV (`step, t_au, occ_0, occ_1, ..., occ_56`).
Creates a 2D heatmap: time on x-axis, KS state index on y-axis, color = occupation.

Uses `ax.imshow(..., cmap="inferno", vmin=0, vmax=1)`.

Marks the boundary between occupied and empty states with a cyan dashed line,
and the WP state (last state) with a yellow dashed line.

**What to look for:** The heatmap should show the WP state (state 56) having occupation 1
before t=t_scatter, then the occupation spreading to multiple states after scattering.
The ground-state orbitals (0–53) should have occupation ≈ 2 throughout.

---

### `plot_z_trajectory(path)`

Reads `density_z_profile_vs_time.csv` (shape: N_snapshots × (2+Nz) columns).
Creates a 2D image: z on y-axis, time on x-axis, color = density at cell centre.

This shows the wavepacket moving down in z, hitting the coronene plane, and then
the scattered component moving back up.

Uses hardcoded `LZ_BOHR = 59.904` to set the z-axis scale.

---

### `plot_leed_pattern(path)`

Reads `leed_pattern.txt` with `np.loadtxt(path, comments="#")`.

Creates a 2-panel figure:
- Left: `ax.pcolormesh(..., cmap="hot")` — linear colour scale
- Right: same data with `LogNorm` — log scale to see weak diffraction features

Also saves a Fourier transform of the LEED pattern:
```python
leed_fft = np.abs(np.fft.fftshift(np.fft.fft2(leed))) ** 2
```
This gives the autocorrelation / power spectrum of the diffraction pattern.

---

## 10. `analysis/make_video.py`

**Location:** `ResearchProject/systems/jellium/analysis/make_video.py`

**Role in workflow:** Creates an animation of the jellium wavepacket evolution.

**Input:** Sequence of `slice_tNNN.txt` or `density3d_tNNN.txt` files.

**Output:** Video file (`.mp4` or `.gif`) showing density evolving over time.

**Typical approach:**
```python
import matplotlib.animation as animation

fig, ax = plt.subplots()
im = ax.imshow(slices[0], cmap='hot', origin='lower')

def update(frame):
    im.set_data(slices[frame])
    ax.set_title(f"t = {times[frame]:.2f} a.u.")
    return [im]

ani = animation.FuncAnimation(fig, update, frames=len(slices), interval=100)
ani.save("wavepacket_evolution.mp4", writer='ffmpeg', dpi=150)
```

---

## 11. `analysis/plot_spreading.py`

**Location:** `ResearchProject/systems/jellium/analysis/plot_spreading.py`

**Role in workflow:** Validates the free-particle propagation by comparing
the numerically computed wavepacket width σ(t) against the analytical formula.

**Input:**
- Sequence of `slice_tNNN.txt` files from the propagation
- `grid_info.txt` for σ₀ (initial width) and time step

**Functions:**
```python
def compute_sigma_from_slice(density_2d, dx):
    """Compute RMS width of a 2D distribution.
    
    σ_x = √(⟨x²⟩ - ⟨x⟩²) = √(Σ_i x_i² ρ_i / Σ_i ρ_i - (Σ_i x_i ρ_i / Σ_i ρ_i)²)
    
    Returns: (sigma_x, sigma_y)
    """

def sigma_analytic(t, sigma0):
    """
    Free-particle spreading: σ(t) = σ0 √(1 + t² / (4σ0⁴))
    
    Derivation: for a Gaussian WP, the width in position space grows as
    σ(t)² = σ0² + (t/2σ0)² (in a.u., m=1)
    """
    return sigma0 * np.sqrt(1.0 + (t / (2.0 * sigma0**2))**2)
```

**Output:**
- `results/sigma_vs_time.png` — σ(t) comparison plot
- Points: INQ numerical σ computed from each slice
- Dashed line: analytical σ(t) formula

---

## 12. Angelo Jellium Scripts Overview

**Location:** `Tutorial/angelo-jellium/scripts/`

These scripts form a complete analysis pipeline for the angelo-jellium tutorial.
They rely on the utility modules in `Tutorial/angelo-jellium/jellium/`.

| Script | What it does |
|---|---|
| `apply_wavepacket.py` | Injects a Gaussian WP into the last extra-state orbital in INQ binary files; validates norm |
| `experiment.py` | Main workflow: load GS, inject WP, run propagation, collect snapshots |
| `plot_ground_state.py` | Plots density, eigenvalues, and shell structure from the INQ ground state |
| `plot_eigenvalues.py` | Histogram or level diagram of KS eigenvalue spectrum |
| `plot_wavepacket.py` | 3D slices of the initial Gaussian WP; validates against analytical formula |
| `plot_real_time.py` | Real-time observables: energy, dipole, current vs time |
| `plot_density_evolution.py` | Grid of density snapshots showing WP evolution |
| `plot_comparison.py` | Side-by-side: INQ-propagated density vs FFT free-particle propagation |
| `plot_eigenvalue_dynamics.py` | KS eigenvalue changes during propagation |
| `plot_floquet.py` | Floquet analysis: quasi-energy spectrum under periodic driving |

---

### Common patterns across all these scripts

**Reading the INQ binary state:**
```python
from jellium.inq_io import read_grid_params, read_spin_density, read_all_orbitals

cell_side, spacing, N = read_grid_params(".")
density = read_spin_density(".", N)           # (N, N, N) float64
orbitals = read_all_orbitals(".", N)          # list of (N, N, N) complex128
```

**Plotting a 2D slice of a 3D array:**
```python
from jellium.inq_io import make_coords

x = make_coords(cell_side, N)
mid = N // 2
slice_xy = density[:, :, mid]    # z = midplane

plt.figure()
plt.pcolormesh(x, x, slice_xy.T, cmap='viridis', shading='auto')
plt.colorbar(label='n (bohr⁻³)')
plt.xlabel('x (bohr)'); plt.ylabel('y (bohr)')
```

**Propagating a WP and comparing with INQ:**
```python
from jellium.wavepacket import GaussianWavePacket3D, propagate_free, check_normalization
import numpy as np

wp = GaussianWavePacket3D(
    r0=np.array([cell_side/2]*3),
    sigma=1.0,
    k0=np.array([0.0, 0.0, 1.5])
)

for step in range(N_steps):
    t = step * DT
    psi_t = propagate_free(wp, cell_side, spacing, t)
    density_t = np.abs(psi_t)**2
    # Compare density_t to INQ checkpoint density
```

**Reading the INQ log for eigenvalues:**
```python
from jellium.inq_io import parse_eigenvalues_from_log, write_eigenvalues

evals = parse_eigenvalues_from_log("run.log")
write_eigenvalues(evals, "results/eigenvalues.tsv")
```
