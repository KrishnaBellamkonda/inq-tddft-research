# Jellium Gaussian Wave Packet — Implementation Notes

**Status:** Task 1 complete (2026-04-01). Task 2 planned but not started.

**Location:** `ResearchProject/systems/jellium/`

---

## Overview

This document records the full plan, implementation decisions, bugs encountered, and fixes applied
for the Gaussian wave packet simulation in INQ. It covers both the completed Task 1 (Angelo
validation) and the planned Task 2 (Kononova-style interacting jellium).

### References

1. **Angelo** — "Quantum Wave-Packet Preparation and Electron Dynamics in Jellium," Cavendish Lab
   report, Candidate 3221L. Primary reference for all parameters, Eqs. (5)–(7).
2. **Kononova et al.** — "Electron dynamics in extended systems within real-time TDDFT," *MRS
   Communications* 12, 1002–1014 (2022). DOI: 10.1557/s43579-022-00273-7. Task 2 reference.
3. **Robinett** — "Quantum wave packet revivals," *Phys. Rep.* 392, 1–119 (2004). Revival physics
   (T_rev = L²/π).
4. **Castro, Marques, Rubio** — "Propagators for the TD-KS equations," *J. Chem. Phys.* 121, 3425
   (2004). ETRS method.

---

## Task 1: Angelo Validation Test — Single-Orbital Free-Particle Spreading

### Scientific Basis

Angelo §3.2: *"Before analysing the many-electron TDDFT data, it is useful to verify that a
genuinely localized orbital can be prepared in the same periodic geometry."*

Model: 1 electron, non-interacting theory, cubic periodic cell. The KS potential is identically
zero (Hartree=0, XC=0 for non-interacting). Dynamics reduce to the free-particle TDSE:

    i ∂ψ/∂t = −∇²/2 ψ

Analytical solution: the wave packet width evolves as

    σ(t) = σ₀ √(1 + t²/(4σ₀⁴))     [Angelo Eq. 6]

and the exact quantum revival in a periodic box occurs at

    T_rev = L²/π     [Angelo Eq. 7]

Initial wavefunction (Angelo Eq. 5):

    ψ_wp(r, 0) = N exp[−d²(r,r₀)/(4σ²_wp)] exp[ik₀(z−z₀)]

where d(r,r₀) is the minimum-image distance.

### Parameters

| Parameter | Symbol | Value | Units | Notes |
|-----------|--------|-------|-------|-------|
| Cell type | — | Cubic, 3D-periodic | — | — |
| Cell side | L | 40 | bohr | Enlarged from Angelo's 13.89; prevents re-entry |
| Electrons | N_e | 1 | — | Single-orbital validation |
| Theory | — | non-interacting | — | KS potential ≡ 0 |
| Spin | — | unpolarized | — | — |
| Grid cutoff | E_cut | 40 | Ha | → h ≈ 0.333 bohr, 120³ grid |
| Packet centre | r₀ | **(0, 0, 0)** | bohr | Cell centre in INQ's [-L/2,L/2) coords |
| Packet width | σ_wp | 1 | a₀ | Angelo §2.2 |
| Carrier momentum | k₀ | 1.5 | a₀⁻¹ | Angelo §2.2 (z-direction) |
| Propagator | — | ETRS | — | Castro et al. (2004) |
| Time step | Δt | 0.04 | atu | User request (Angelo used 0.02) |
| Total time | T | 50 | atu | ~1.9 cell traversals |
| Num steps | N_steps | 1250 | — | T/Δt |
| Snapshot interval | — | every 25 steps | — | Every 1 atu |
| Revival time | T_rev | L²/π ≈ 509 | atu | >> simulation window |
| Wrap time | T_wrap | L/k₀ ≈ 26.7 | atu | Centre crosses boundary |

**Important:** r₀ = (0,0,0), not (L/2, L/2, L/2). See error notes below.

### Results (Validated 2026-04-01)

| t (atu) | σ_num (a₀) | σ_ana (a₀) | Relative error |
|---------|-----------|-----------|----------------|
| 0 | 1.000000 | 1.000000 | 0% |
| 1 | 1.11803 | 1.11803 | 0% |
| 2 | 1.41421 | 1.41421 | 0% |
| 3 | 1.80277 | 1.80278 | < 0.001% |
| 4 | 2.23607 | 2.23607 | 0% |
| 5 | 2.69282 | 2.69258 | 0.009% |

Energy at t=0: 1.500000000000 Ha (= k₀²/2 + 3/(8σ₀²) = 1.125 + 0.375).
Energy drift at t=50: ~1×10⁻⁶ Ha (ETRS is essentially unitary for non-interacting).
Centre of mass ⟨z⟩ = k₀t for early t, wraps at t≈7 when leading edge crosses z=+L/2.

**Physical note on late-time divergence:**
After t≈7 atu, σ_num diverges from σ_ana. This is expected physics: the packet's leading edge
(at z ≈ k₀t + 3σ(t)) crosses the PBC boundary at z_max = +20 bohr when t≈7. Raw-z moments
⟨z⟩, ⟨z²⟩ then pick up wrapped density from z ≈ −20, inflating the computed σ. This is
correct periodic-cell behaviour, not a code error.

### Output Files

```
ResearchProject/systems/jellium/
├── gaussian_wave_packet.cpp          # C++ driver
├── analysis/
│   ├── plot_spreading.py             # σ(t) static figure
│   └── make_video.py                 # MP4 animation
├── results/
│   ├── sigma_t.txt                   # σ(t) data (51 rows, t=0–50 atu)
│   ├── grid_info.txt                 # grid dimensions & parameters
│   ├── sigma_comparison.pdf          # 4-panel σ(t) figure
│   ├── energy_conservation.pdf       # energy drift figure
│   ├── gaussian_spreading.mp4        # 5.1 s video (51 frames, 10 fps)
│   └── snapshots/
│       ├── slice_NNNN.npy            # 2D midplane density ρ(x,y) per atu
│       └── line_NNNN.npy             # 1D line ρ(x) per atu
└── run.log                           # full stdout from simulation
```

### Build & Run

```bash
cd ResearchProject/systems/jellium
/local/data/public/skcb2/tddft/shared/bin/inq-run          # GPU build (CUDA sm_80)
/local/data/public/skcb2/tddft/shared/bin/inq-run --cpu    # CPU build

# After simulation:
python3 analysis/plot_spreading.py
python3 analysis/make_video.py
```

GPU wall time: ~10 min (1250 steps × 0.47 s/step on A100/sm_80).
CPU wall time: ~29 min (1250 steps × 1.38 s/step).

---

## Errors Encountered and Fixes

### Error 1: `gpu::run` does not accept `std::array`

**Symptom:** Compile error —
```
error: no matching function for call to 'run(const std::array<int, 3>&, lambda)'
```

**Cause:** The code used `gpu::run(basis.local_sizes(), lambda)`, passing the whole array.
`gpu::run` requires individual `size_t` arguments.

**Fix:** Expand to individual indices in z-outermost order:
```cpp
// Wrong:
gpu::run(basis.local_sizes(), GPU_LAMBDA (auto iz, auto iy, auto ix) { ... });

// Correct:
gpu::run(basis.local_sizes()[2], basis.local_sizes()[1], basis.local_sizes()[0],
         GPU_LAMBDA (auto iz, auto iy, auto ix) { ... });
```

**Rule:** Always use `[2]` (z, outermost), `[1]` (y), `[0]` (x, innermost) to match the
`(iz, iy, ix)` lambda signature and the `cubic()[iz][iy][ix]` accessor convention.

---

### Error 2: `GPU_LAMBDA` inside a generic lambda (CUDA error)

**Symptom:** GPU build fails with —
```
error: An extended __device__ lambda cannot be defined inside a generic lambda expression("operator()")
```

**Cause:** The propagation callback `[&](auto data){ ... GPU_LAMBDA ... }` is a *generic lambda*
(has `auto` parameter). CUDA's extended lambda feature (`__device__` lambdas) cannot be used
inside generic lambdas.

**Fix:** Precompute all GPU-kernel results outside the callback. Specifically, the weight fields
`w_x`, `w_y`, `w_z`, `w_z2` are time-independent (just Cartesian coordinates) and were being
recomputed at every snapshot with a `GPU_LAMBDA`. Instead, compute them once before the propagation
call, then just call `operations::integral_product(density, w_z)` inside the callback — no GPU
kernel needed:

```cpp
// Outside the callback (not a generic lambda context):
basis::field<basis::real_space, double> w_z(dbasis), w_z2(dbasis);
gpu::run(dbasis.local_sizes()[2], dbasis.local_sizes()[1], dbasis.local_sizes()[0],
    [wz = begin(w_z.cubic()), wz2 = begin(w_z2.cubic()), dpoint_op]
    GPU_LAMBDA (auto iz, auto iy, auto ix) {
        auto r = dpoint_op.rvector_cartesian(ix, iy, iz);
        wz[iz][iy][ix] = r[2];
        wz2[iz][iy][ix] = r[2]*r[2];
    });

// Inside the callback (no GPU_LAMBDA):
real_time::propagate<>(ions, electrons, [&](auto data) {
    auto density = data.electrons().density();
    double sigma_z_sq = operations::integral_product(density, w_z2) / N
                      - pow(operations::integral_product(density, w_z) / N, 2);
    // ...
}, ...);
```

**Side benefit:** Also faster — weight fields are allocated and filled once, not 50 times.

---

### Error 3: Gaussian at wrong location (σ_num(0) ≈ 19 instead of 1)

**Symptom:** The t=0 snapshot gives σ_num = 19.034 instead of σ₀ = 1.0.

**Cause:** INQ uses a **symmetric coordinate system**: grid index `i` maps to coordinate
`(i - N/2) × h` via `to_symmetric_range`, so real-space coordinates span `[-L/2, L/2)`.
The code had `x0 = y0 = z0 = L/2 = 20 bohr`, which places the Gaussian centre at the *cell
boundary*, not the centre. The minimum-image convention splits the packet into two halves at
z = ±20 bohr. Raw-z moments of this split packet give ⟨z⟩ ≈ 0 and ⟨z²⟩ ≈ (L/2)² ≈ 400,
so σ_num ≈ √400 = 20 ≈ 19 (measured).

Simultaneously, the energy (4.72 Ha instead of 1.5 Ha) and its non-conservation revealed that
the split Gaussian was the wrong initial state.

**Fix:** Set `x0 = y0 = z0 = 0.0`. The cell centre in INQ is always at the origin.

**Rule:** In any INQ driver, `rvector_cartesian` returns coordinates in `[-L/2, L/2)`. The cell
centre is (0, 0, 0). Confirmed by `to_symmetric_range` in `src/basis/grid.hpp`:
```cpp
if(ii[idir] >= (sizes[idir] + 1)/2) ii[idir] -= sizes[idir];
```

---

### Error 4: `data.every()` skips iter=0

**Symptom:** No t=0 snapshot in the callback output.

**Cause:** From `src/real_time/viewables.hpp`:
```cpp
auto every(int every_iter) const {
    if(iter() == 0) return false;
    return (iter()%every_iter == 0) or last_iter();
}
```
The `iter()==0` guard is intentional — at iter=0 the propagator hasn't taken a step yet.

**Fix:** Add a standalone t=0 snapshot block *before* calling `real_time::propagate<>()`,
operating directly on `electrons.density()` (which is already set by
`observables::density::calculate(electrons)` right after orbital initialization).

---

### Error 5: Race condition — two processes writing to `results/`

**Symptom:** `sigma_t.txt` had only the t=0 row; snapshot files had gaps (0000, 0040–0050
but missing 0001–0039). Energy reported as 0% drift because all E_total values rounded to "1.5".

**Cause:** `inq-run` builds *and then runs* the binary (via `exec`). When `inq-run` was run as
a background task and `./gaussian_wave_packet` was also started manually, both processes wrote
to the same `results/` directory simultaneously. The `sigma_t.txt` was opened (truncated) by
the second process, discarding the first process's data.

**Fix:** `pkill -9 -f gaussian_wave_packet` to kill all instances, then `rm -rf results` and
restart with `nohup ./gaussian_wave_packet > run.log 2>&1 &` to get a single clean run. Use
`nohup` to prevent SIGHUP killing the process when the shell session closes.

---

### Error 6: `sigma_t.txt` E_total column shows "1.5" (no drift visible)

**Symptom:** The energy conservation plot shows ΔE = 0.000e+00 %.

**Cause:** The C++ code writes `sigma_log << E_tot`, which uses `std::ostream` default precision
(6 significant figures). At t=50 atu, E = 1.499999000610 Ha, which rounds to "1.5" in 6 digits.

**Workaround:** Energy conservation is confirmed from `run.log`, which shows full precision:
`e = 1.499999000610` at step 1250. Drift = (1.5 − 1.499999000610) / 1.5 = 6.7×10⁻⁵ %.

**Permanent fix (not yet applied):** Use `sigma_log << std::setprecision(12) << E_tot` in
the write_snapshot lambda in `gaussian_wave_packet.cpp`.

---

### Error 7: `\pmod` not supported in matplotlib mathtext

**Symptom:** `plot_spreading.py` crashes with `ValueError: Unknown symbol: \pmod`.

**Fix:** Replace `r"$z_0 + k_0 t \pmod{L}$"` with `r"$k_0 t$ mod $L$"`.

---

## INQ API Gotchas (Summary)

These apply to any INQ C++ driver, not just this project.

| Issue | Rule |
|-------|------|
| `gpu::run` array | Pass sizes as `[2], [1], [0]` individually; `std::array` not accepted |
| `GPU_LAMBDA` in generic lambda | Precompute outside `[&](auto data)` callbacks |
| Coordinate origin | INQ coords are `[-L/2, L/2)`, cell centre = (0,0,0) |
| `data.every(N)` | Returns `false` at iter=0; add explicit t=0 snapshot before propagation |
| `electrons.density()` | Returns `observables::density::total(spin_density_)` — valid after `electrons.spin_density() = observables::density::calculate(electrons)` |
| `kpin()` | Returns reference to `kpin_` (type `std::vector<states::orbital_set<basis::real_space, complex>>`); modifications persist |
| `data.root()` | Equivalent to `electrons.full_comm().root()` (see `viewables.hpp:59`) |
| `operations::integral_product(f1, f2)` | GPU-safe; no need for custom reduction loops |
| Stream precision | Use `std::setprecision(12)` when writing energy to text files |

---

## Task 2: Kononova-Style — Gaussian in 40-Electron Jellium (Planned)

### Scientific Basis

Kononova et al. Eq. (8)–(9): *"The time-dependent Kohn–Sham states include all the electrons
of the target material and the incident electron of the wave packet."*

The Gaussian is an **additional (41st) electron** added to a self-consistent 40-electron jellium
ground state and propagated under full TDDFT (LDA). This studies how a localised wave packet
evolves in a many-electron environment.

### Parameters (Differences from Task 1)

| Parameter | Task 1 | Task 2 |
|-----------|--------|--------|
| N_e | 1 | 41 (40 ground-state + 1 Gaussian) |
| Theory | non-interacting | LDA (Ceperley-Alder/Perdew-Zunger) |
| Temperature | none | 100 K Fermi smearing |
| Extra states | 0 | 5 (Gaussian + buffer for unoccupied) |
| r_s (L=40, N=40) | — | ≈ 7.45 a₀ (lower density than Angelo's r_s=2.52) |
| Δt | 0.04 atu | 0.02 atu (smaller for interacting stability) |

**Note on r_s:** Angelo used L=13.89 bohr with N=40 electrons, giving r_s=2.52 a₀ (metallic density).
At L=40 bohr with N=40, r_s ≈ 7.45 a₀ (much lower density). To exactly match Angelo's density
at L=40 would require N ≈ 955 electrons (prohibitively expensive). This is an accepted trade-off
from using the larger cell.

### Implementation Plan

```cpp
// Stage A: 40-electron jellium ground state
auto electrons = systems::electrons(env.par(), ions,
    options::electrons{}
        .cutoff(40.0_Ha)
        .extra_electrons(40.0)
        .temperature(0.00861_Ha)      // 100 K
        .extra_states(5),
    input::kpoints::gamma());
ground_state::calculate(ions, electrons,
    options::theory{}.lda(),
    options::ground_state{}.steepest_descent().mixing(0.3));

// Stage B: Append Gaussian as 41st orbital
// (Insert orbital at end of kpin()[0], Gram-Schmidt against existing 40)

// Stage C: Real-time TDDFT
real_time::propagate<>(ions, electrons, callback,
    options::theory{}.lda(),
    options::real_time{}.num_steps(2500).dt(0.02_atomictime).etrs());
```

### Key Observables to Add (vs Task 1)

- Induced density δn(r,t) = n(r,t) − n₀
- Electronic current J(t) (Kononova Eq. 4)
- Projected orbital energies ε̃_j(t) (Angelo Eq. 11)

### Status

Not started. Prerequisites: Task 1 validated (complete).

---

## Session History

| Date | Action |
|------|--------|
| 2026-03-28 (est.) | Initial plan written, `gaussian_wave_packet.cpp` created |
| 2026-04-01 | GPU build fixed (3 compile errors), coordinate origin fixed, t=0 snapshot added, simulation run to completion, analysis figures and video generated |
