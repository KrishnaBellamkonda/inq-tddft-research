# Handover: Jellium Ground State

## Current status

**Phase: Ground-state and convergence runs complete. All figures generated.**

Both INQ runs have finished and all Python visualisations have been produced.
The next substantive task is: interpret results, decide whether to refine anything,
and proceed to TDDFT / perturbation calculations.

---

## What changed (this session)

### Rules
- `.claude/rules/file-placement.md` — rule 7 added: always save figures as `.png`; never `.pdf`/`.svg` unless explicitly requested.
- `.claude/rules/testing.md` — rule 0 added: use GPU whenever available (`inq-run` not `inq-run --cpu`).

### 01_ground_state
- `run.cpp` — orbital section 7 rewritten: now writes **two separate files per shell**:
  `orbital_N_n2_M_real.txt` (Re[ψ_k]) and `orbital_N_n2_M_imag.txt` (Im[ψ_k]).
- `plot_results.py` — updated to read split files; added Im[ψ] as a third plot row; all PDF saves removed (PNG only).
- `plot_orbitals_3d.py` — **new file**: PyVista headless 3D isosurface renderer. Reads k-vectors from file headers, reconstructs Re[ψ_k] on 80³ grid, renders ±40% isosurfaces red/blue. Requires `quantum-wave-packet` pyenv.

### 02_ground_state_convergence
- `run_convergence.cpp` — output now written to `results/convergence_results.csv` via `std::ofstream` (with `out.flush()` after each row). No longer relies on stdout redirect. Also had comments added by user linter.
- `plot_convergence.py` — reads from `results/convergence_results.csv` by default (no CLI argument needed); parser hardened with try/except to skip column-header lines mixed in by INQ logger.

---

## Files touched

| File | Status |
|---|---|
| `/local/data/public/skcb2/tddft/ResearchProject/jellium/01_ground_state/run.cpp` | Complete; run succeeded |
| `/local/data/public/skcb2/tddft/ResearchProject/jellium/01_ground_state/plot_results.py` | Complete; figures generated |
| `/local/data/public/skcb2/tddft/ResearchProject/jellium/01_ground_state/plot_orbitals_3d.py` | Complete; figures generated |
| `/local/data/public/skcb2/tddft/ResearchProject/jellium/01_ground_state/jellium_utils.hpp` | Unchanged |
| `/local/data/public/skcb2/tddft/ResearchProject/jellium/02_ground_state_convergence/run_convergence.cpp` | Complete; run succeeded |
| `/local/data/public/skcb2/tddft/ResearchProject/jellium/02_ground_state_convergence/plot_convergence.py` | Complete; figures generated |
| `/local/data/public/skcb2/tddft/docs/handovers/jellium_ground_state.md` | This file |
| `/local/data/public/skcb2/tddft/.claude/rules/file-placement.md` | PNG rule added |
| `/local/data/public/skcb2/tddft/.claude/rules/testing.md` | GPU rule added |

---

## Commands run

```bash
# Ground-state run (~5 min, GPU)
cd ResearchProject/jellium/01_ground_state && inq-run > /tmp/gs_run.log 2>&1

# Convergence run (~18 min, GPU)
cd ResearchProject/jellium/02_ground_state_convergence
./run_convergence > results/convergence_results.csv 2>/tmp/conv_progress.log

# Visualisations (quantum-wave-packet pyenv)
cd ResearchProject/jellium/01_ground_state && python plot_results.py
cd ResearchProject/jellium/01_ground_state && python plot_orbitals_3d.py
cd ResearchProject/jellium/02_ground_state_convergence && python plot_convergence.py
```

---

## Tests and validation

**Proposed:** component tests in run.cpp (E_H≈0, E_ext=0, E_nl=0, E_ion=0, ΔT_s<0.5 Ha, ΔE_total<0.5 Ha).  
**Approved:** user said "Run the scripts."  
**Run:** all tests executed.  
**Outcomes:**

- `01_ground_state`: **All 6 tests PASSED.**
  - E_Hartree = 0 (tol 1e-4), E_ext = 0, E_non_local = 0, E_ion = 0
  - ΔT_s = small (within 0.5 Ha tol)
  - E_total (INQ) = −2.621561 Ha vs analytical −2.621752 Ha (Δ = 0.19 mHa) ✓

- `02_ground_state_convergence` (Test A — grid spacing):
  - E_total varies < 0.3 mHa across h = 0.70 → 0.30 bohr
  - **h = 0.50 bohr (E_cut ≈ 20 Ha) is confirmed converged** — no update to SPACING needed
  - All spacings gave E_total ≈ −2.6215 Ha

- `02_ground_state_convergence` (Test B — shell closure):

  | N | T_s/N (mHa) | T_TF/N (mHa) | Δ |
  |---|---|---|---|
  | 2 | ≈ 0 | 20.99 | N=2 has only k=0 state; T_s=0 |
  | 14 | 21.29 | 20.99 | +1.4% above T_TF |
  | 38 | 20.16 | 20.99 | −4.0% below |
  | 54 | 20.20 | 20.99 | −3.8% below |
  | 66 | 20.88 | 20.99 | −0.5% — nearly converged |

  Shell oscillations visible; converging to Thomas-Fermi bulk limit as N→∞.

---

## Output figures

All in `ResearchProject/jellium/`:

| Figure | Path |
|---|---|
| Shell structure diagram | `01_ground_state/results/shell_structure.png` |
| XC offset verification | `01_ground_state/results/xc_offset.png` |
| 2D orbital slices (Re, Im, density) | `01_ground_state/results/orbitals.png` |
| 3D isosurfaces per shell | `01_ground_state/results/orbital_3d_n2_{0..4}.png` |
| 3D summary (tiled) | `01_ground_state/results/orbitals_3d_summary.png` |
| E_cut convergence | `02_ground_state_convergence/convergence_Ecut.png` |
| Shell-closure convergence | `02_ground_state_convergence/convergence_shells.png` |

---

## Trusted sources used

- Perdew & Zunger (1981) PRB 23, 5048 — PZ81 LDA; V_xc and ε_xc used in jellium_utils.hpp
- Dreizler & Gross, *Density Functional Theory* — Thomas-Fermi T_TF = (3/5)E_F
- INQ tutorial (`docs/inq_tutorial.md`) — API usage

---

## Known issues / blockers

- `run_convergence.cpp` source has been updated to write to file directly, but the **binary has not been rebuilt** with this change. Next `inq-run` will rebuild automatically. The current CSV was produced by the old binary via stdout redirect — results are valid.
- INQ's internal logger writes to stdout, which contaminated the CSV when using stdout redirect. Fixed in both source (ofstream) and parser (try/except).
- Clang LSP false positives in all `.cpp` files (inq/inq.hpp not in IDE include path) — not real errors.

---

## Assumptions still in play

- SPACING = 0.50 bohr is confirmed converged (ΔE < 0.3 mHa). No update needed.
- LDA (Perdew-Zunger) is appropriate for the homogeneous electron gas.
- PyVista 3D renders use `quantum-wave-packet` pyenv; must NOT be active when running `inq-run`.
- N=2 T_s≈0 is physically correct: only the k=(0,0,0) state is occupied, which has zero kinetic energy.

---

## Exact next steps

1. Review figures (particularly `orbitals_3d_summary.png` and `convergence_shells.png`).
2. Decide whether to proceed to TDDFT / real-time propagation or refine any ground-state aspect.
3. When running `inq-run` again: ensure `quantum-wave-packet` pyenv is NOT active.
4. When running Python visualisation scripts: activate `quantum-wave-packet` pyenv first.
