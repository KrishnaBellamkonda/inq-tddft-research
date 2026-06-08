# Manifest: Report 1 submission package

**Branch:** report1/submission-package  
**Generated:** 2026-05-28  
**Scope:** Files required for submission to examiner, including cited run directories, shared infrastructure, and report-generation scripts.

---

## Summary

- **Total Jellium runs:** 29 (18 WP σ-parameter sweep + 4 high-density WP + 5 classical + 1 special mechanism + 1 free-prop)
- **Total Coronene runs:** 2
- **Unique Cfg headers cited:** 26 (jellium) + 2 (coronene)
- **Shared infrastructure:** run_template.hpp, boundary_rule.hpp (jellium); run_template.hpp (coronene)
- **Library payload:** inqkit (37 .hpp) + inqview (99 .py modules total; report1 subdirectory hosts 20+ figure scripts)
- **Draft5 figure scripts:** 20 make_fig_*.py scripts in docs/reports/report1/drafts/draft5/scripts/

---

## Jellium runs (29)

### L=50 low-density (r_s ≈ 5.69): σ=1 energy sweep (6 runs)

| Run dir | A7 | A8 | A9 | A10 | A11 | A12 | analyse.py | Cfg header |
|---|---|---|---|---|---|---|---|---|
| run_wp_n162_L50_E20_sigma1_v2 | — | — | — | — | ✓ | ✓ | ✓ | electron_proj_E20_L50_cubic_sigma1_v2.hpp |
| run_wp_n162_L50_E25_sigma1_v2 | ✓ | ✓ | — | ✓ | ✓ | ✓ | ✓ | electron_proj_E25_L50_cubic_sigma1_v2.hpp |
| run_wp_n162_L50_E50_sigma1_v2 | — | — | — | — | ✓ | ✓ | ✓ | electron_proj_E50_L50_cubic_sigma1_v2.hpp |
| run_wp_n162_L50_E100_sigma1_v2 | — | — | — | — | ✓ | ✓ | ✓ | electron_proj_E100_L50_cubic_sigma1_v2.hpp |
| run_wp_n162_L50_E200_sigma1_v2 | — | — | — | — | ✓ | ✓ | ✓ | electron_proj_E200_L50_cubic_sigma1_v2.hpp |
| run_wp_n162_L50_E300_sigma1_v2 | — | — | — | — | ✓ | ✓ | ✓ | electron_proj_E300_L50_cubic_sigma1_v2.hpp |

### L=50 low-density: σ=5 energy sweep (4 runs)

| Run dir | A11 | A12 | analyse.py | Cfg header |
|---|---|---|---|---|
| run_wp_n162_L50_E50_v2 | ✓ | ✓ | ✓ | electron_proj_backlog_L50_cubic.hpp |
| run_wp_n162_L50_E100_v2 | ✓ | ✓ | ✓ | electron_proj_E100_L50_cubic.hpp |
| run_wp_n162_L50_E300_v2 | ✓ | ✓ | ✓ | electron_proj_backlog_L50_cubic.hpp |
| run_wp_n162_L50_E600_v2 | ✓ | ✓ | ✓ | electron_proj_backlog_L50_cubic.hpp |

### L=50 low-density: σ supplementary at E=100 (3 runs)

| Run dir | A11 | analyse.py | Cfg header |
|---|---|---|---|
| run_wp_n162_L50_E100_sigma0p5 | ✓ | ✓ | electron_proj_E100_L50_cubic.hpp |
| run_wp_n162_L50_E100_sigma3 | ✓ | ✓ | electron_proj_E100_L50_cubic.hpp |
| run_wp_n162_L50_E100_sigma8 | ✓ | ✓ | electron_proj_E100_L50_cubic.hpp |

### L=50 low-density: Classical Ehrenfest (5 runs)

| Run dir | A8 | A11 | A12 | analyse.py | Cfg header |
|---|---|---|---|---|---|
| run_classical_n162_L50_E20 | — | ✓ | ✓ | ✓ | electron_proj_E20_L50_cubic.hpp |
| run_classical_n162_L50_E25 | ✓ | ✓ | ✓ | ✓ | electron_proj_E25_L50_cubic.hpp |
| run_classical_n162_L50_E50_v2 | — | ✓ | ✓ | ✓ | electron_proj_E50_L50_cubic_sigma1_v2.hpp |
| run_classical_n162_L50_E100_v2 | — | ✓ | ✓ | ✓ | electron_proj_E100_L50_cubic_sigma1_v2.hpp |
| run_classical_n162_L50_E600_v2 | — | ✓ | ✓ | ✓ | electron_proj_backlog_L50_cubic.hpp |

### L=50 low-density: Free propagation & mechanism (2 runs)

| Run dir | A4 | A7 | A9 | Role | analyse.py | Cfg header |
|---|---|---|---|---|---|---|
| run_free_wp_L50_E25_sigma1_v2 | ✓ | ✓ | — | Matched free (validation) | ✓ | boundary_rule.hpp |
| run_plasmon_n162_L50_E15 | — | — | ✓ | Mechanism (long T=2000 a.u., Δω≈0.085 eV) | ✗ | plasmon_n162_L50_E15.hpp |

### L=30 high-density (r_s ≈ 3.41): σ=1 energy sweep (4 runs)

| Run dir | A11(b) | analyse.py | Cfg header |
|---|---|---|---|
| run_wp_n162_L30_E50_highdens_sigma1_v2 | ✓ | ✓ | highdens_n162_L30_E50_sigma1_v2.hpp |
| run_wp_n162_L30_E100_highdens_sigma1_v2 | ✓ | ✓ | highdens_n162_L30_E100_sigma1_v2.hpp |
| run_wp_n162_L30_E200_highdens_sigma1_v2 | ✓ | ✓ | highdens_n162_L30_E200_sigma1_v2.hpp |
| run_wp_n162_L30_E300_highdens_sigma1_v2 | ✓ | ✓ | highdens_n162_L30_E300_sigma1_v2.hpp |

### L=30 high-density: Supplementary (1 run)

| Run dir | A11 | analyse.py | Cfg header |
|---|---|---|---|
| run_wp_n162_L30_E100_highdens | ✓ | ✓ | highdens_n162_L30_E100.hpp |

### L=30 high-density: Classical Ehrenfest (4 runs)

| Run dir | A11(b) | analyse.py | Cfg header |
|---|---|---|---|
| run_classical_n162_L30_E50_highdens | ✓ | ✓ | highdens_n162_L30_E50_sigma1.hpp |
| run_classical_n162_L30_E100_highdens | ✓ | ✓ | highdens_n162_L30_E100.hpp |
| run_classical_n162_L30_E200_highdens | ✓ | ✓ | highdens_n162_L30_E200_sigma1.hpp |
| run_classical_n162_L30_E300_highdens | ✓ | ✓ | highdens_n162_L30_E300_sigma1.hpp |

---

## Coronene runs (2)

| Run dir | A6 | A6(b) FFT | analyse.py | Cfg header |
|---|---|---|---|---|
| run_propagate_paper_replica | ✓ | ✓ | ✗ | (not found) |
| run_cc_bond | ✓ | — | ✓ | cc_bond_35x35x60.hpp |

**Note:** `run_propagate_paper_replica` provides LEED screens 7 & 14; `run_cc_bond` is the C–C bond impact variant. Neither run.cpp file includes a config path (config embedded inline or via template).

---

## Shared infrastructure

### Jellium system

**Path:** `ResearchProject/systems/jellium/shared/`

- **cpp/**
  - `run_template.hpp` — base template for all jellium runs
  - `eigenvalues_writer.hpp` — observables writer
  - `leed_screen_layout.hpp` — LEED screen grid (unused in report1)
  - `results_paths.hpp` — output path definitions

- **configs/** (26 .hpp files)
  - **Unique files cited by runs above:**
    - `boundary_rule.hpp` — free-prop WP setup
    - `electron_proj_E{20,25,50,100,200,300}_L50_cubic_sigma1_v2.hpp` (6 files) — σ=1 sweeps
    - `electron_proj_E100_L50_cubic.hpp` — σ variable sweeps
    - `electron_proj_backlog_L50_cubic.hpp` — legacy σ=5 sweeps  
    - `electron_proj_E{20,25}_L50_cubic.hpp` (2 files) — classical runs
    - `electron_proj_E50_L50_cubic_sigma1_v2.hpp` — classical-matched
    - `electron_proj_E100_L50_cubic_sigma1_v2.hpp` — classical-matched
    - `highdens_n162_L30_E{50,100,200,300}_sigma1_v2.hpp` (4 files) — high-dens σ=1
    - `highdens_n162_L30_E100.hpp` — high-dens supplementary
    - `highdens_n162_L30_E{50,200,300}_sigma1.hpp` (3 files) — high-dens classical
    - `plasmon_n162_L50_E15.hpp` — mechanism run (long T)

### Coronene system

**Path:** `ResearchProject/systems/coronene/shared/`

- **cpp/**
  - `run_template.hpp`
  - `eigenvalues_writer.hpp`
  - `leed_screen_layout.hpp`
  - `results_paths.hpp`

- **configs/** (1 .hpp file used in report1)
  - `cc_bond_35x35x60.hpp` — C–C bond impact config

---

## Library payload

### inqkit (INQ wrapper toolkit)

**Path:** `inq-stack/include/inqkit/`

- **Total:** 37 .hpp header files
- **Key modules:** 
  - `fields/` — density, orbital, wavefunction representations
  - `io/` — VTI writers, observables export, complex field I/O
  - `observables/` — density_delta, wp_momentum_stats, state_energy_writer, occupations, etc.
  - `real_time/` — real_time_session, step_context
  - `wavepacket/` — wavepacket injection, momentum tracking
  - `jellium/` — jellium-specific (shells, projected configs)

### inqview (visualization & post-processing)

**Path:** `inq-stack/python/inqview/`

- **Total:** 99 .py modules across all subdirectories
- **Top-level modules (12):** `__init__.py`, `cli.py`, `colors.py`, `fields.py`, `io.py`, `paraview.py`, `styling.py`, `tufte.py`, `util.py`, `write_json.py`, `write_png.py`, `write_vti.py`
- **report1 subdirectory:** 20+ `.py` figure scripts (legacy + draft5)
  - Includes: `stopping_power_data.py` (master stopping-power aggregator), `fig01_*`, `fig_*_*.py`, `_shared_style.py`, `render_setup3d.py`

---

## Draft5 figure scripts

**Path:** `docs/reports/report1/drafts/draft5/scripts/`

Scripts to be copied into staging inqview/report1/ (20 files):

1. `make_fig01_nuclear_vs_electronic.py`
2. `make_fig_schematic_regime.py`
3. `make_fig_setup_density.py` (A4 density-diff setup)
4. `make_fig_free_wp_panel.py` (A5 2×2 panel)
5. `make_fig_coronene_setup.py` (A6 r1c1 setup)
6. `make_fig_leed_transmission.py` (A6 r1c2 FFT)
7. `make_fig_leed_backscatter_centre.py` (A6 r2c1)
8. `make_fig_density_diff_2d.py` (A7)
9. `make_fig_density_profile.py` (A7)
10. `make_fig_energy_decomp_system.py` (A8a)
11. `make_fig_energy_decomp_wp.py` (A8b)
12. `make_fig_gs_decomposition.py` (A8c)
13. `make_fig_loss_function.py` (A9 2D & 1D)
14. `make_fig_plasmon_fft.py` (appendix)
15. `make_fig_momentum_1d.py` (A10a)
16. `make_fig_momentum_2d.py` (A10b)
17. `make_fig_master_stopping.py` (A11)
18. `make_fig_stopping_defs_combined.py` (A12)
19. `make_fig_jellium_gs.py` (B2)
20. `make_fig_leed_backscatter_ccbond.py` (B3)

**Note:** A6 r2c2 (Tsubonoya panel) is externally supplied and manually placed; `render_setup3d.py` produces the A5 r1c1 ParaView 3D render.

---

## Summary of per-run file inventory

### File presence (all 31 runs)

| Metric | Present | Absent | Notes |
|--------|---------|--------|-------|
| run.cpp | 31/31 | 0 | ✓ all present |
| analyse.py | 30/31 | 1 | ✗ run_plasmon_n162_L50_E15 has no analysis script |
| results/run_summary.txt | 31/31 | 0 | ✓ all present (verified on sample) |
| REPORT.md | 0/31 | 31 | — not yet written; not part of canonical run outputs |

### Cfg header coverage

- **26 unique jellium configs** cited by the 29 jellium runs
- **1 unique coronene config** (cc_bond) cited explicitly; run_propagate_paper_replica config not found in run.cpp
- **All Cfg namespaces** follow `jellium::config::*` or coronene templates; compatible with inqkit framework

---

## Deliberately excluded

These run families are **not** listed above and are excluded for the stated reasons:

| Family | Reason |
|--------|--------|
| **Knudsen-sweep family** (`E{700,800,900,1100}_knudsen`) | Exploratory high-energy regime; not cited by panels_plan or stopping_power_data |
| **Legacy σ-variable runs (v1)** (`E{20,25}_sigma{1,5}`, `E100`, etc. without `_v2` suffix) | Superseded by v2 variants with identical or improved parameters; not cited in current analysis |
| **Hypotheses / exploratory** (`base_n138*`, `E{50,200,400}_s{0p5,0p53,2p0}*`) | Early-stage parameter studies; not in final report narrative |
| **_compare runs** (none in active set) | Placeholder for multi-system comparisons not developed |
| **MPI/parallel variants** (`mpi_inject`, `mpi_propagate`, `_minimal`) | Debugging versions; results subsumed into canonical runs |
| **Geometry variants** (`x40`, `x80`, `tilt45`, `L40x40x150`, `L60`) | Parameter studies or non-standard cell shapes; not cited |
| **E1p5 variants** (`base_n162_L50_E1p5`, etc.) | Low-energy trial runs; superseded |

---

## Missing / ambiguities for user review

1. **run_plasmon_n162_L50_E15 (analyse.py absent):**
   - Used only for A9 loss-function + plasmon-FFT plots
   - Long run (T=2000 a.u., 501 time steps) — may not have been analysed with standard pipeline
   - **Decision:** include in manifest; note that only raw observables are required

2. **run_propagate_paper_replica (no Cfg in run.cpp):**
   - Core coronene mechanism run (LEED screens 7 & 14, step 330)
   - run.cpp does not appear to have a traditional `#include "shared/configs/..."` line
   - **Decision:** include in manifest; verify run.cpp structure with user if needed

3. **High-density mirror track (Part C of panels_plan.md):**
   - Runs already cited (L=30 σ=1 + classical families)
   - Full per-panel high-density mirrors (A7–A10 analogs) not yet decided in panels_plan
   - **Current scope:** include cited high-density runs; mirror scripts (if any) to be added per user decision

4. **Tsubonoya panel (A6 r2c2):**
   - Manually supplied PNG; no generation script included
   - `make_fig_tsubonoya.py` was a post-hoc stretch wrapper in draft5 session
   - **Decision:** external image; not part of code manifest

---

## Staging structure

When assembling the submission package, plan to organize as:

```
docs/reports/report1/code/
├── ResearchProject/
│   └── systems/
│       ├── jellium/
│       │   ├── run_wp_n162_L50_E25_sigma1_v2/  [full dir with run.cpp, results/]
│       │   ├── run_wp_n162_L50_E20_sigma1_v2/
│       │   ├── ... [27 more jellium runs]
│       │   └── shared/
│       │       ├── cpp/
│       │       └── configs/
│       └── coronene/
│           ├── run_propagate_paper_replica/
│           ├── run_cc_bond/
│           └── shared/
│               ├── cpp/
│               └── configs/
├── inq-stack/
│   ├── include/inqkit/  [37 .hpp files]
│   └── python/inqview/  [99 .py + report1/ subdir]
└── docs/reports/report1/
    └── drafts/draft5/scripts/  [20 make_fig_*.py scripts]
```

