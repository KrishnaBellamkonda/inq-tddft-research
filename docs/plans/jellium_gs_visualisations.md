# Plan: Jellium Ground-State Visualisations

Related experiment: `ResearchProject/jellium/01_ground_state/`  
Status: plan — awaiting user confirmation before implementation

---

## 1. KS Orbital Visualisation   ← AWAITING CONFIRMATION

### Goal
Visualise individual KS orbitals and the total electron density in real space.
For jellium the KS orbitals are exact plane waves, so the orbital density
|ψ_k(r)|² is perfectly uniform and the real part Re[ψ_k] oscillates sinusoidally.
The visualisation should confirm both facts.

### Proposed approach — Option A (recommended): write a cube file from run.cpp

Add a post-SCF block to `run.cpp` that:
1. Reads the first few orbital grids from `electrons` (using INQ's field
   iteration API, to be confirmed by reading `inq/src/systems/electrons.hpp`
   and `inq/src/operations/io.hpp`).
2. Writes a Gaussian cube file (standard volumetric format readable by VESTA,
   VMD, Python with `ase.io`).

Then a Python script `plot_orbitals.py` reads the cube files and produces:

| Panel | What it shows |
|---|---|
| 2D heatmap (x-y plane, z=L/2) | `|ψ_k(r)|²` for the |n|²=0 and |n|²=1 shells — should be perfectly uniform and ~cos²-modulated respectively |
| 2D heatmap (x-y plane, z=L/2) | Total density `ρ(r) = Σ_i f_i |ψ_i(r)|²` — should be nearly constant (deviations show finite-grid noise) |
| 1D line plot (along cell diagonal) | `Re[ψ_k]` for each occupied shell — shows the plane-wave oscillation with wavelength λ=L/|n| |
| Colour-bar annotation | Maximum fractional deviation from uniform: (max ρ − min ρ)/n₀ |

### Alternative — Option B: use `electrons.save()` + Python reader

After SCF, call `electrons.save("jellium_gs")`. INQ writes binary files.
A Python script reads them (same format as Angelo's code — flat complex128
arrays of size N_x × N_y × N_z per orbital). This avoids adding C++ plotting
code but requires knowing the exact INQ save format (to be verified).

### Files to create
- `run.cpp` — add a post-SCF cube-file output block (Option A)
  OR `electrons.save("jellium_gs")` call (Option B)
- `plot_orbitals.py` — reads cube/binary, produces figure
- `orbitals/` subfolder — output cube files (large, not checked in)

### Decision needed from user
- Option A (cube from C++) or Option B (save + Python reader)?
- Which shells to visualise (all occupied, or just |n|²=0,1,2)?

---

## 2. Shell Structure Plot   ← DONE

Script: `plot_shell_structure.py`  
Output: `shell_structure.pdf`, `shell_structure.png`

Side-by-side level diagram:
- Left panel: free-electron energies ε_k = k²/2
- Right panel: KS eigenvalues ε_k = k²/2 + V_xc(n₀)

Shells colour-coded by occupancy (blue=full, amber=partial, grey=empty).
Line thickness proportional to shell degeneracy.
Fermi energy shown as a dashed red line.
A labelled double-headed arrow between the panels marks the V_xc offset.

Run with: `python3 plot_shell_structure.py`

---

## 3. XC Offset Verification Plot   ← AWAITING CONFIRMATION

### Goal
Numerically verify that every KS eigenvalue satisfies ε_i = k_i²/2 + V_xc(n₀)
by plotting the INQ eigenvalues against the predicted free-electron energies.

### Plan

**Data needed:**
After SCF, extract KS eigenvalues from `electrons.eigenvalues()`.
This returns a GPU array indexed [kpin][state] — shape (1, N_states) for
Gamma-only. Eigenvalues are ordered by energy ascending.

Map each eigenvalue to its shell |n|² by matching the energy ordering to the
analytical shell list (unambiguous because the shells are well separated for
L=13.89 bohr).

**The plot:**

```
x-axis:  k_i²/2  (Ha)  — free-electron kinetic energy of shell i
y-axis:  ε_i     (Ha)  — KS eigenvalue from INQ
```

Expected: all points lie on the line  y = x + V_xc,  with slope exactly 1.

Overplot:
- The line  y = x + V_xc(PZ81)  in dashed red (analytical prediction).
- Points colour-coded by shell |n|² (or by occupancy).
- A text annotation showing the fitted intercept vs PZ81 prediction.
- Possibly a residual panel beneath showing ε_i − (k_i²/2 + V_xc) vs state index.

**Code required:**

Add a block to `run.cpp` (or a separate `run_eigenvalues.cpp`):
```cpp
// Access eigenvalues from electrons object after SCF
auto evals = electrons.eigenvalues();   // gpu::array<double, 2>
// Copy to host and print or write to file for Python
```
The Python script `plot_xc_offset.py` then reads the eigenvalue file and
produces the scatter plot.

**Alternative:** Print eigenvalues directly in run.cpp and redirect stdout to a
file, then parse with Python.

### Decision needed from user
- Add eigenvalue output to existing `run.cpp`, or create a separate lightweight
  script that loads a saved ground state?
- Print to stdout (simple) or write to a structured file (cleaner)?

---

## 4. Finite-Size Error Analysis   ← AWAITING CONFIRMATION

### What finite-size errors are and how they arise

For a finite periodic cell of N electrons, the energy per electron E_total/N
differs from the infinite bulk limit for three reasons:

**4a. Kinetic energy shell oscillations (dominant for small N)**

In the bulk HEG, T_s/N = (3/5) E_F (Thomas-Fermi). In a finite box, the
discrete shell structure at Gamma gives a different value for each N because
some shells are partially filled. Only for "magic numbers" N where a shell is
exactly complete does T_s/N approach the bulk value cleanly.

For our system: |n|²=3 is partially filled for N=40, giving a kinetic energy
error of order (E_{|n|²=3} − E_F) × f_partial ~ 0.3 Ha.

Test: fix r_s (keep rs=2.52), vary N ∈ {2, 14, 38, 54, ...} (closed-shell
numbers), compute T_s/N and compare to (3/5)E_F. Shows convergence to bulk as
N increases.

**4b. Grid cutoff error (dominant for coarse spacing)**

INQ represents wavefunctions on a real-space grid of spacing h. Plane waves
with |k| > π/h are excluded. The kinetic energy cutoff is E_cut = π²/(2h²).

For a smooth free-electron system, the error in E_total scales as
exp(−c × E_cut / E_F) — exponentially fast convergence. Our current spacing
h=0.347 bohr gives E_cut ≈ 41 Ha >> E_F ≈ 0.29 Ha, so this error is
completely negligible.

Test: fix N=40, L=13.89, vary h ∈ {0.55, 0.45, 0.40, 0.347, 0.30, 0.25} bohr.
Corresponding E_cut ∈ {16, 24, 31, 41, 55, 79} Ha.
Plot E_total vs E_cut. Expect convergence by ~25 Ha (well below our default).

**4c. Finite-size XC correction (small, but present)**

The LDA XC energy E_xc = N × ε_xc(r_s) is exact for uniform jellium regardless
of N (it's a local functional applied to a uniform density). No finite-size
correction expected from XC itself.

The Hartree energy is exactly cancelled by the background for any N (the
neutralising background is exact). So E_H = 0 for any box size.

**Recommended tests to run (user to authorise):**

| Test | What to vary | What to measure | Expected result |
|---|---|---|---|
| Grid convergence | h ∈ {0.55, 0.45, 0.40, 0.347, 0.30, 0.25} bohr | E_total vs E_cut | Flat by E_cut=30 Ha |
| Shell closure (larger effort) | N ∈ {2, 14, 38, 54} (magic closed shells) | T_s/N vs N | Converges to (3/5)E_F |
| Spacing effect on T_s | same as grid convergence | T_s vs E_cut | Should be ~constant (T_s is insensitive to grid for smooth plane waves) |

**Implementation needed:**
- For grid convergence: add a spacing loop to a new `run_finite_size.cpp`
  in the same folder, or run `inq-run` multiple times with modified spacing.
- Plot with `plot_finite_size.py`.

### Decision needed from user
- Run grid convergence test (cheap: 6 SCF runs)?
- Run shell-closure test (requires new L values, more work)?
- Implement as a loop in C++ or as separate calls?

---

## Files produced so far

| File | Status |
|---|---|
| `plot_shell_structure.py` | Done — run `python3 plot_shell_structure.py` |
| `energy_summary.csv` | Done — fill Numerical INQ column after running `inq-run` |
| `plot_orbitals.py` | Pending confirmation (§1) |
| `plot_xc_offset.py` | Pending confirmation (§3) |
| `run_finite_size.cpp` or script | Pending confirmation (§4) |
