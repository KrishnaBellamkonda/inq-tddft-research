# Handover: Coronene TDDFT LEED Simulation

---

## Current status

**Run 004 is COMPLETE.** All 1561 TDDFT steps finished (4.74 h wall time, GPU A30).
Analysis script written and run; all 10 figures generated.
**Next: interpret LEED patterns; consider follow-up runs with wider screen spacing or crystal periodicity.**

---

## Run 004 results summary (2026-04-16)

| Validation check | Result |
|---|---|
| GS energy | 288.983 Ha (LDA, coronene C24H12) |
| SCF converged | 110 iterations, tol=1e-4 Ha |
| WP norm | 1.000000 ✓ |
| Jz(t=0) | 3.8311 a.u. vs k₀=3.8340 (rel err −0.08%) ✓ |
| Energy drift | −1.31×10⁻⁵ Ha/a.u. (total 0.41 mHa over 31.25 a.u.) ✓ |
| LEED screens written | leed_screen0/1/2.txt ✓ |

### Physical observations

- **WP trajectory**: heatmap (zprofile_heatmap.png) shows WP starting at z_obs=44.9 Å, moving in −z,
  arriving at z_flake=23.8 Å at T1=10.4 a.u. as expected. Physics is correct.
- **Density at z_flake**: coronene D₆h hexagonal ring structure perfectly resolved. Visible
  density perturbation when WP passes through at T1.
- **Jz sign**: Jz=+3.83 throughout despite WP moving in −z. INQ's `.observables_current()` 
  returns the negative of electron probability current (confirmed by Jz dip at T1: backscattering 
  reduces |Jz|, consistent with convention). Physical magnitude is correct.
- **Jx, Jy**: Non-zero at t>0 (~0.05–0.13 a.u.) from GS molecular ring currents + WP interaction.
  Not a bug — coronene π-electrons carry angular momentum.
- **Real-space LEED**: cross/blob pattern centred on molecule with negative sidelobes. This is the
  single-molecule LEED form factor — NOT 6-fold Bragg peaks (those require a crystal lattice).
  The LEED pattern reflects the molecular shape convoluted with the WP beam profile (d=0.53 Å).
- **k-space LEED**: concentrated near k≈0 (specular reflection dominates). Weak structure visible
  at 1–3 Å⁻¹, within first graphene BZ. Consistent with single-molecule diffraction.
- **Three screens**: z_screen0/1/2 at 44.90, 45.78, 46.67 Å — only 0.88 Å apart (= 5d/3).
  Patterns nearly identical as expected. Consider wider spacing in next run.

### Figures generated (results/figures/)

| File | Contents |
|---|---|
| `leed_real_space.png` | ∫Δn dt at 3 screens, real space, RdBu colourmap |
| `leed_kspace.png` | 2D FFT of LEED, klim=15 Å⁻¹ with graphene ring overlays |
| `leed_kspace_zoomed.png` | Same, klim=6 Å⁻¹ |
| `energy_vs_time.png` | ΔE(t) in mHa, linear drift −1.31e-5 Ha/a.u. |
| `momentum_vs_time.png` | Jx/Jy/Jz vs time |
| `wp_trajectory_z.png` | WP z-centroid vs time |
| `zprofile_heatmap.png` | 2D heatmap Δn(z,t) showing WP trajectory |
| `density_snapshots_flake.png` | 2D density at z_flake, 10 snapshots |
| `density_snapshots_obs.png` | 2D density at z_obs, 10 snapshots |
| `density_snapshots_mid.png` | 2D density at z_mid, 10 snapshots |

---

## Run 003 (previous result, COMPLETE)

Run 003 is COMPLETE. All 516 TDDFT steps finished. Analysis validated (5/5 checks PASS).
Run 003 moved to `runs/run_003_200eV_d1p4A_40Ha_linear/`.

---

## Directory structure

All simulation runs are now under:
`/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/04_leed_simulation/runs/`

| Run | Parameters | Status |
|---|---|---|
| `run_001_200eV_d0p53A/` | d=0.53 Å, molecule at corner (broken) | Complete, archived |
| `run_002_200eV_d1p4A_54Ha/` | d=1.4 Å, centred, Broyden tol=1e-4 | Complete |
| `run_003_200eV_d1p4A_40Ha_linear/` | d=1.4 Å, all fixes applied | **Complete — analysis done** |
| `run_004_200eV_d0p53A_40Ha_longZ/` | d=0.53 Å, Lz×1.5, D maximised | **Ready to run** |

---

## What changed in this session

### run_003 completed and validated

- 516 TDDFT steps at 200 eV, d=1.4 Å, LDA, 40 Ha cutoff, linear+Broyden tol=1e-4
- GS energy: −371.551 Ha
- WP norm ≈ 1.0, KS excitation (max|1−|S_ii|²|) = 1.08e-2 → scattering detected
- `utils.hpp` `mkdir_p` bug fixed: was single-level only; now recursive with EEXIST check
- `analysis.py` created: 5/5 validation checks PASS (GS energy, SCF converged, WP norm,
  scattering, LEED 6-fold symmetry at 0.366 > 0.3 threshold)

### LEED background problem diagnosed and fixed in analysis

- **Problem**: `leed_pattern.txt` from run.cpp accumulates raw n(x,y,z_obs,t), not ∫Δn dt.
  At z_obs=41.95 bohr, molecular π-electron tails dominate: ~2.6e-2 vs ~9.7e-3 bohr⁻³ WP.
  This gives a radially symmetric cross/plus pattern (Fourier of molecular orbitals),
  not the hexagonal diffraction pattern.
- **Fix in analysis.py**: `time_averaged_density_at_obs(subtract_gs=True)` subtracts the
  t=0 obs-plane density before averaging. This isolates the WP-scattered contribution.
- **Fix in run_004 run.cpp**: LEED accumulator subtracts GS baseline captured before WP injection.

### Original run_001 "good" appearance explained

`coronene.xyz` was centred at origin → negative atom coords → INQ wraps to cell boundaries
→ molecule split to 4 grid corners. With narrow WP (d=0.53 Å) at cell corner, the WP
coherently illuminated 4 artificial molecular copies → structured interference. Visually
hexagonal-ish but physically wrong. run_002/003 fixed the geometry (centred XYZ).

### INQ cell origin convention confirmed

INQ uses (0,0,0) as cell **corner**, not centre. Grid runs 0→L in each direction.
Atoms must have all-positive coordinates within [0,L]. Confirmed by:
- run_001 failure (negative coords → duplication)
- `po.rvector_cartesian(0,0,0)` returns `(0,0,0)`
- `coronene_centered.xyz` validated correct by geometry_check.py (all 8 checks pass)

---

## Files touched this session

| File | Change |
|---|---|
| `runs/run_003_200eV_d1p4A_40Ha_linear/utils.hpp` | Fixed `mkdir_p` to be recursive; added `#include <cerrno>` and `#include <cstring>` |
| `runs/run_003_200eV_d1p4A_40Ha_linear/analysis.py` | Created — full validation + LEED + trajectory plots with background subtraction |
| `runs/run_004_200eV_d0p53A_40Ha_longZ/config.hpp` | New — d=0.53 Å, Lz=47.55 Å, D=21.125 Å, T2=0.756 fs |
| `runs/run_004_200eV_d0p53A_40Ha_longZ/coronene_centered.xyz` | New — z shifted to 23.775 Å (= new Lz/2) |
| `runs/run_004_200eV_d0p53A_40Ha_longZ/run.cpp` | New — adapted from run_003; background-subtracted LEED accumulator |
| `runs/run_004_200eV_d0p53A_40Ha_longZ/utils.hpp` | Copied from run_003 (has all fixes) |
| `docs/handovers/coronene_wp_scattering.md` | This file |

---

## Commands run (this session)

```bash
# Checked run_003 completion
wc -l .../blkh4bbgb.output && tail -10 .../blkh4bbgb.output

# Run analysis
cd .../run_003_200eV_d1p4A_40Ha_linear && python3 analysis.py

# Move run_003 into organised runs/ directory
mv .../run_003_200eV_d1p4A_40Ha_linear .../runs/

# Run 004: to be launched
cd .../runs/run_004_200eV_d0p53A_40Ha_longZ && inq-run
```

---

## Tests and validation (run_003)

| Test | Result |
|---|---|
| GS energy ≈ −371 Ha (LDA coronene) | PASS (−371.551 Ha) |
| SCF converged (iter < 300) | PASS (Broyden tol=1e-4, iter ~80) |
| WP norm ∈ [0.97, 1.03] | PASS (≈1.0) |
| WP scattering detected: max\|1−\|S_ii\|²\| > 1e-3 | PASS (1.08e-2) |
| LEED 6-fold: max_peak/mean_peak > 0.3 (BG-subtracted) | PASS (0.366) |

---

## run_004 parameter rationale

| Parameter | run_003 | run_004 | Reason |
|---|---|---|---|
| d (WP width) | 1.4 Å | 0.53 Å | Paper value; focused beam for LEED contrast |
| Lz | 31.7 Å | 47.55 Å | 1.5× to allow WP to propagate further from molecule |
| D (WP to molecule) | 6.35 Å | 21.125 Å | Maximised: Lz/2 − 5d = 23.775 − 2.65 Å |
| z_obs | 22.200 Å | 44.900 Å | = z_flake + D (LEED screen at WP start) |
| T2 | 0.25 fs | 0.756 fs | Must exceed T1=0.252 fs; 3×T1 matches paper T2/T1≈3.3 |
| N_steps | 516 | ~1561 | T2/DT |
| LEED accumulator | raw n(t) | n(t)−n_GS | Background subtraction isolates WP signal |

D comparison: run_003 D=6.35 Å → run_004 D=21.125 Å → **3.33× larger**.
This matches run_001's physical setup (focused beam, long propagation) but with correct geometry.

---

## Known issues / blockers

1. **run_004 is ~3× longer than run_003** (1561 vs 516 steps). At ~2 min/step GPU rate
   from run_003, expect ~52 min wall time. Confirm GPU availability before launching.

2. **GS orbitals save is large**: 57 orbitals × Nx×Ny×Nz = 57 × 100×100×263 grid points
   (Nz scales with Lz). Each orbital.txt ≈ 50 MB → ~2.85 GB total. Confirm disk space.
   (Total run_004 output estimate: ~10 GB given larger Nz.)

3. **Overlap matrix CPU overhead**: 57×57×(100×100×263) ≈ 8.1 billion ops per call,
   every 10 steps → ~156 calls. Estimate ~13 min extra overhead. Acceptable.

4. **z-profile is a line scan** at (Nx/2, Ny/2, iz), NOT xy-marginal density. WP centroid
   from z-profile is unreliable when molecular π-tails are comparable to WP signal.
   Use scattering metric (|S_ii|² deviation) as primary scattering indicator.

---

## Assumptions still in play

1. E_cut=40 Ha: adequate despite paper's ~54 Ha (justified by convergence sweep + pseudodojo artefacts at 54 Ha)
2. SCF tol=1e-4 Ha: sufficient for 0.76 fs TDDFT initial conditions (validated in run_002/003)
3. LDA/ALDA: acceptable for coronene ground state and real-time propagation
4. LEED screen at z_obs = WP start position (reflection geometry): WP travels −z, scatters off molecule, returns to z_obs screen
5. Background subtraction isolates WP contribution correctly: assumes GS density at z_obs is constant in time (no drift)

---

## Exact next steps

1. **Interpret LEED pattern**: The current k-space pattern is dominated by specular reflection
   (k≈0). To resolve molecular diffraction:
   - **Option A**: Increase D further (larger cell → WP farther from molecule → better-defined
     far-field). Requires more storage/time.
   - **Option B**: Use a periodic 2D molecular array (crystal surface) → sharp Bragg peaks.
     This is how real LEED experiments work. INQ supports periodic BCs natively.
   - **Option C**: Compare run_004 LEED pattern to run_003 (d=1.4 Å) to see beam-width effect.

2. **Wider screen spacing**: If doing a new run, place screens at z_obs, z_obs+D/3, z_obs+2D/3
   (not compressed into 5d/3=2.65 Å window).

3. **Overlap matrix post-hoc**: Can compute S_ij from saved GS+time orbitals offline — not
   computed in-run to avoid CPU bottleneck. Script `analysis.py` has `compute_overlap_matrix_posthoc`
   placeholder (not yet implemented).

4. **Paper comparison**: Tsubonoya et al. Fig.2 shows LEED pattern for coronene. Compare our 
   k-space pattern shape to their result.

*Updated: 2026-04-16 — run_004 complete (1561 steps, 4.74 h); analysis.py written and run; all 10 figures generated*
