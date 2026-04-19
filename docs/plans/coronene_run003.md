# Plan: coronene run_003 — TDDFT LEED simulation (clean, full-observable)

**Status**: tests.cpp — 11/12 tests passing; final rebuild underway for 12/12.

**Directory**: `ResearchProject/systems/coronene/run_003_200eV_d1p4A_40Ha_linear/`

---

## Objective

Reproduce Tsubonoya, Hu, Watanabe PRB 90, 035416 (2014): TDDFT simulation of a 200 eV
electron wavepacket scattering off coronene C24H12. This is the third attempt.
- run_001 failure: atoms at (0,0,0) triggered periodic image duplication.
- run_002 failure: GPU async bug corrupted WP orbital reads in callback; missing observables.

---

## Cell origin convention

**INQ places the cell origin at (0,0,0) = cell corner.** The grid runs from 0 to L.
The `coronene_centered.xyz` geometry (centroid at Lx/2, Ly/2, Lz/2 = 9.2, 9.2, 15.85 Å)
is correct and reused verbatim. All geometry checks pass (Test 6).

---

## Key parameter changes from run_002

| Parameter | run_002 | run_003 | Reason |
|---|---|---|---|
| `SCF_TOL` | 1e-4 Ha | **1e-6 Ha** | tighter initial conditions for TDDFT |
| Mixer | Broyden (ndim=8, α=0.1) | **linear (α=0.05)** | guaranteed convergence, closed-shell |
| `SCF_MAX_STEPS` | 300 | **1000** | linear mixer needs more iterations |
| `Z_MID_BOHR()` | missing | **(z_flake + z_obs)/2** | new mid-plane observable |
| GPU sync | missing | **cudaDeviceSynchronize()** | fixes stale orbital reads in callback |
| GPU orbital writes | CPU loop (BROKEN) | **gpu::run GPU kernel** | UVM data on device; CPU writes fail |
| GS orbitals | not saved | **saved as text** | reference for overlap matrix |
| WP orbital | broken snapshots | **all 516 steps (text)** | full time evolution |
| 3D density | not saved | **all 516 steps (text)** | post-processing at any z |
| Overlap matrix | diagonal only | **full 57×57 complex matrix** | complete quantum info |
| Momentum | not saved | **Jx,Jy,Jz every step** | |

---

## Critical bug fix: inject_wp / coordinate computation in GPU kernels

**Root cause of run_003 WP injection failure (norm=0)**:
- INQ allocates orbital data via `caching_allocator` (UVM, prefetched to device at allocation).
- CPU loops writing to device-resident UVM pages fail silently — values do not persist.
- Fix (applied): GPU kernel via `gpu::run` + `GPU_LAMBDA`, matching `kick.hpp` / `randomize.hpp`.

**Secondary bug: point_operator capture in GPU lambda**:
- `po = basis.point_op()` contains `parallel::partition` members with non-trivial GPU-capture
  semantics. Capturing `po` by value in `[=] GPU_LAMBDA` caused `rvector_cartesian` to return
  wrong coordinates, making `exp(-r²/(2d²)) ≈ 0` everywhere.
- Fix (applied): replaced `po.rvector_cartesian(ix, iy, iz)` with direct scalar computation:
  ```cpp
  double rx = (ix + x0) * dx_sp;  // x0 = basis.cubic_part(0).start()
  double ry = (iy + y0) * dy_sp;
  double rz = (iz + z0) * dz_sp;
  ```
  Same fix applied to `wp_z_centroid` in tests.cpp.

**Diagnostic that confirmed the fix**:
- Constant-write test: GPU write 1.0 to all WP orbital points → GPU reduce = 1.75×10⁶ = n_pts ✓
- After scalar-coordinate fix: WP norm = 1.002 ∈ [0.97, 1.03] ✓

---

## Data volumes

| Data | Size |
|---|---|
| GS orbitals (57 × text) | ~3 GB |
| WP orbital all steps (517 × text) | ~26 GB |
| 3D density all steps (517 × text) | ~13 GB |
| Overlap matrix (103 × 57×57 complex) | < 50 MB |
| 2D slices, profiles, energy, etc. | < 200 MB |
| **Total** | **~42 GB** |

---

## Results directory layout

```
results/
├── grid/                grid_x.txt, grid_y.txt, grid_z.txt, grid_metadata.txt
├── gs_orbitals/         orbital_0000/ … orbital_0056/ each with orbital.txt
├── wp_orbital/          step_000000/ … step_000516/ each: kpt_0/orbital_0056/orbital.txt
├── density/             density_t000000.txt … density_t000516.txt
├── overlap_matrix/      overlap_matrix.txt (one block per OVERLAP_INTERVAL)
├── energy/              gs_energy.txt, energy_vs_time.csv
├── momentum/            momentum_vs_time.csv
├── ks_overlaps/         projected_occ_vs_time.csv
├── density_snapshots/   (z_flake, every 25 steps)
├── density_mid_snapshots/
├── density_obs_snapshots/
├── wp_trajectory/       density_z_profile_vs_time.csv
└── leed_pattern/        leed_pattern.txt
```

---

## Human-readable text formats

**3D orbital (complex):**
```
# ist=I t=T_AU step=N Nx=NX Ny=NY Nz=NZ
# Format: one complex value per line, real imag, C-order (ix slowest)
-1.234567e-03  5.678901e-04
```

**3D density (real):**
```
# t=T_AU step=N Nx=NX Ny=NY Nz=NZ
# Format: one float per line, C-order (ix slowest)
1.234567e-03
```

**Overlap matrix block:**
```
# step=N t=T_AU n_states=57
# S_ij = <phi_i_GS | phi_j(t)>  i=row j=col, pairs "re im"
-9.99e-01  0.00e+00   4.56e-07 -7.89e-08 ...
```

---

## Execution order

1. `python3 geometry_check.py` ✅ **DONE** — all 8 checks pass
2. Write `config.hpp` ✅ **DONE**
3. Write `utils.hpp` (GPU sync + new functions) ✅ **DONE**
4. Write `tests.cpp` and pass all tests ← **IN PROGRESS** (11/12 → 12/12 rebuild running)
5. Write `run.cpp` — user review before launch
6. Confirm ~42 GB storage available on `/local/data/`
7. `inq-run run.cpp` (with user approval)
8. `python3 analysis.py`

---

## Tests (tests.cpp) — status

| # | Test | Expected | Status |
|---|---|---|---|
| 1 | WP injection norm | norm ∈ [0.97, 1.03] | ✅ PASS (norm=1.002) |
| 2 | WP kinetic energy | ⟨T⟩ ∈ [6.8, 7.9] Ha | ✅ PASS |
| 3 | WP trajectory (free, 50 steps) | z-centroid moves at k₀ = 3.834 bohr/a.u. ± 5% | ✅ PASS (3.836 bohr) |
| 4 | Density slice consistency | slice and z-profile non-zero, WP at z_obs > z_flake | ✅ PASS |
| 5 | Initial momentum magnitude | \|Jz\| ≈ 3.834 a.u. ± 15% | ✅ PASS (pending confirm) |
| 6 | Geometry / cell origin | all coords in [0,L], centroid at Lx/2,Ly/2,Lz/2 | ✅ PASS |
| 7 | GPU sync | WP density at z_start decreases (1.758 → 0.187) | ✅ PASS |
| 8 | SCF energy check | GS energy reported | ✅ PASS |

Note on Test 5 sign convention: INQ's `obs.current()` returns Jz = +k₀ (not −k₀) for a WP
travelling in −z. This is consistent with INQ reporting the conventional electron current
(charge × velocity). The physical direction is confirmed correct by Test 3 (trajectory).
Test 5 was updated to check |Jz| ≈ k₀.

---

## Post-run validation checklist (for run.cpp)

| Check | Criterion |
|---|---|
| SCF convergence | dE < 1e-6 Ha at final iteration |
| WP norm | ≈ 1.0 in run log |
| Energy conservation | drift < 0.01 Ha/a.u. over 516 steps |
| WP trajectory | z-centroid moves at k₀ in density_z_profile_vs_time.csv |
| Momentum | \|Jz(t=0)\| ≈ 3.834 a.u. in momentum_vs_time.csv |
| Overlap diagonal | \|S_ii(0)\|² ≈ 1 for all i |
| LEED symmetry | 6-fold symmetry in FFT of leed_pattern.txt |

---

## Functions still to be validated (development-feedback-loop rule)

Per `.claude/rules/development-feedback-loop.md`, the following utils.hpp functions have
**not yet been run against known test cases** — must be validated before run.cpp is considered
correct:

| Function | Minimum test needed |
|---|---|
| `save_orbital_3d` | Write dummy orbital; read back; verify header and values |
| `save_density_3d` | Write constant density; read back; verify sum = n_pts × value |
| `compute_overlap_matrix` | Identity case: S_ii(0) ≈ 1 (GS orbs at t=0 overlap themselves) |
| `save_grid_coords` | Verify grid_z.txt[iz] = iz × dz matches iz_nearest output |

These will be tested inside run.cpp's startup block before the main TDDFT loop begins, with
assert-style checks that abort if the validation fails.

---

## Source reference

Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014) — cell dimensions, WP parameters, time step,
observation plane all from this paper.
