# Plan: Jellium WP-RT Initial Exploration

## Scientific motivation

Uniform jellium (a positive background charge with no ionic structure) is the
simplest exactly-solvable DFT system. Using it as the scattering target before
returning to coronene or other real materials allows:

- Baseline observables with no structural artefacts (plane-wave symmetry expected)
- Controlled comparison between runs (energy, σ, angle, open vs closed shell)
- Full observable suite development without worrying about material-specific effects

A parallel concern: the coronene LEED run did not reproduce the expected cross
feature or diffraction pattern. The jellium runs provide a sanity check for the
observable pipeline before returning to that problem.

---

## Simulation parameters (common)

| Parameter | Value |
|---|---|
| Cell | L = 40.0 bohr cubic, periodic |
| Functional | LDA |
| Grid spacing | 0.50 bohr |
| k-points | Γ only |
| Time step dt | 0.02 a.u. |
| Smearing | kT = 0.00862 eV (Fermi-Dirac) |
| Extra states | 3 |

---

## Run matrix

| Run | E (eV) | σ (Å) | N_gs | direction | N_STEPS | T_sim (a.u.) | T_loop (a.u.) | % |
|---|---|---|---|---|---|---|---|---|
| 01_base       | 200 | 0.53  | 38 | +z     | 417 | 8.34  | 10.43 | 80% |
| 02_low_energy |  50 | 0.53  | 38 | +z     | 834 | 16.68 | 20.87 | 80% |
| 03_high_energy| 400 | 0.53  | 38 | +z     | 295 | 5.90  | 7.38  | 80% |
| 04_tilted_45  | 200 | 0.53  | 38 | 45° xz | 350 | 7.00  | 14.76 | 47% |
| 05_wide_sigma | 200 | 2.0   | 38 | +z     | 480 | 9.60  | 10.43 | 92% |
| 06_narrow_sigma| 200 | 0.265 | 38 | +z    | 480 | 9.60  | 10.43 | 92% |
| 07_open_shell | 200 | 0.53  | 40 | +z     | 417 | 8.34  | 10.43 | 80% |

Physics: k₀ = √(2 E_Ha), v = k₀ (a.u.), T_loop = L/v = 40/k₀.
All closed-shell runs use N_gs = 38 (r_s = 7.38 a₀, n = 5.94×10⁻⁴ e/bohr³).
run_07 uses N_gs = 40 (r_s = 7.26 a₀, n = 6.25×10⁻⁴ e/bohr³) to test
fractional orbital occupation effects on scattering.

---

## Observable suite (all 7 runs)

| Observable | Source | Output directory | Frequency |
|---|---|---|---|
| GS total density | `density::total()` pre-WP | `results/density_gs/` | once |
| GS orbital densities | `density::orbital(i)` pre-WP | `results/density_gs_orbitals/orbital_XXXX/` | once each |
| Total electronic density(t) | `total + orbital(wp_idx)` | `results/density_rt_total/` | every WRITE_EVERY |
| Jellium density(t) | `density::total()` | `results/density_rt_jellium/` | every WRITE_EVERY |
| WP density(t) | `density::orbital(wp_idx)` | `results/density_rt_wp/` | every WRITE_EVERY |
| KS overlap matrix O_ij(t) | `OrbitalOverlapMatrix::snapshot()` | `results/overlap/` | every step |
| Energy, current, dipole(t) | `ObservablesWriter` | `results/observables.csv` | every step |
| 20 LEED screens (time-avg) | `LeedPatternAccumulator::save()` | `results/screens/` | end of run |
| Instantaneous screens | `PlaneScreen::extract()` | `results/screens_snapshots/step_XXXXXX/` | every 3 steps |

WRITE_EVERY = 2 for runs 01–02 (slow WP, more frames needed); 10 for runs 03–07.

---

## Analysis outputs (analysis.py per run)

1. N-electron conservation (20 samples from density_rt_total)
2. Density consistency: max|total − jellium − wp|
3. Observables summary (3-panel PNG)
4. FFT spectra (energy + current_x/y/z)
5. VTI conversion (50 frames × 3 density series)
6. Density slice GIFs: xz + yz for total/jellium/wp (6 GIFs, fps=6)
7. LEED 4×5 grid PNG + individual screen PNGs
8. LEED time-evolution GIF (screen_10, midpoint screen)
9. Overlap matrix GIFs: one per evolved orbital j (fps=5)
10. ParaView 3D renders (low priority, skips if pvbatch absent)

---

## New inqkit header

`inq-stack/include/inqkit/observables/orbital_overlap.hpp`
- `OrbitalOverlapMatrix(electrons, n_ref, output_dir)` — stores n_ref GS wavefunctions
- `snapshot(electrons, time_au, step)` — computes n_ref × (n_ref+1) overlap matrix O_ij(t)
- Saves `results/overlap/overlap_XXXXXX.csv` + appends to `results/overlap/index.csv`

---

## Sequential launcher

`ResearchProject/jellium/jellium-wp-rt/run_all_wp_rt.sh`
Runs all 7 simulations in sequence using `inq-run`. Each run's stdout is tee'd
to `run.log` in its directory. Summary table printed at end.

Usage: `bash run_all_wp_rt.sh`

---

## Implementation status

| Component | Status |
|---|---|
| run_01_base/run.cpp | ✓ full observable suite |
| run_02_low_energy/run.cpp | ✓ full observable suite |
| run_03_high_energy/run.cpp | ✓ fixed N=38, full observable suite |
| run_04_tilted_45/run.cpp | ✓ fixed N=38, full observable suite |
| run_05_wide_sigma/run.cpp | ✓ fixed N=38, N_STEPS=480, full suite |
| run_06_narrow_sigma/run.cpp | ✓ fixed N=38, N_STEPS=480, full suite |
| run_07_open_shell/run.cpp | ✓ N=40, full observable suite |
| orbital_overlap.hpp | ✓ implemented |
| run_0[1-7]/analysis.py | ✓ all 7 written |
| run_all_wp_rt.sh | ✓ written |
| docs/observables_reference.md | ✓ written |
| docs/notes/future-todos.md | ✓ written |

---

## Validation pending (after runs complete)

- [ ] Build test: `inq-run` in run_01_base; check `norm_after ∈ [0.97, 1.03]`
- [ ] Overlap at t=0: diagonal ≈ 1, off-diagonal ≈ 0
- [ ] N-electron conservation: mean ≈ N_gs + 1 within 0.5%
- [ ] 20 .dat files in results/screens/
- [ ] No loop-back in run_05 (WP should not reappear from z=0 side at t=9.6 a.u.)
- [ ] Density consistency: max ratio < 1e-6 at mid-simulation frame
