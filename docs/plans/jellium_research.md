# Plan: Jellium Research Project

## Goal

Reproduce and extend Angelo's jellium TDDFT study using the INQ **C++ API**
(not the CLI) in `ResearchProject/jellium/`.  Each sub-directory is one
self-contained experiment compiled and run via `inq-run`.

## Physical setup (all experiments)

- N = 40 electrons, cubic cell L = 13.89 a₀ → r_s = 2.52 a₀ (Al-like)
- Spacing h = 0.347 a₀ → 40³ grid, E_cut ≈ 41 Ha
- LDA (Perdew–Zunger)
- Gamma-point only (jellium is a uniform system; PBC with Γ gives exact
  free-electron shell structure)
- Fermi smearing at 100 K (needed for open-shell |n|²=3 manifold)
- No ions inserted → INQ treats system as jellium with uniform background

## Directory layout

```
ResearchProject/jellium/
├── 01_ground_state/        ← this plan
│   ├── run.cpp             ← jellium SCF + rigorous validation tests
│   ├── jellium_utils.hpp   ← analytical HEG functions (PZ81, shell structure)
│   └── results/
├── 02_kick_tddft/          ← planned: uniform δ-kick, dipole, eigenvalue dynamics
├── 03_wavepacket/          ← planned: Gaussian localised perturbation
└── lib/                    ← future: shared library when utils stabilise
```

## Experiment 01 — Ground-state benchmark

### What we do
1. Set up jellium (no ions, extra_electrons = 40).
2. Run LDA SCF to convergence.
3. Apply four rigorous tests (see Tests section).

### Tests (quantitative)

| Test | Expected | Tolerance | Rationale |
|---|---|---|---|
| `E_hartree ≈ 0` | 0 Ha | < 1×10⁻⁴ Ha | Uniform ρ → all G≠0 density components ≈ 0 |
| `E_external ≈ 0` | 0 Ha | < 1×10⁻⁸ Ha | No ionic pseudopotentials |
| `E_non_local ≈ 0` | 0 Ha | < 1×10⁻⁸ Ha | No non-local PP |
| `E_ion ≈ 0` | 0 Ha | < 1×10⁻⁸ Ha | No nuclei |
| `T_s ~ HEG shell sum` | ~6.75 Ha | < 0.5 Ha | KS orbitals should be plane waves |
| `E_total ~ T_s + N·ε_xc` | ~ −2.1 Ha | < 0.5 Ha | HEG analytical estimate |

### Analytical reference (jellium_utils.hpp)

Source: Perdew & Zunger, PRB 23, 5048 (1981) for LDA exchange-correlation.

Functions to implement:
- `wigner_seitz_radius(N, L)` — r_s in bohr
- `fermi_wavevector(N, L)` — k_F in bohr⁻¹
- `fermi_energy(N, L)` — E_F in Ha
- `exc_pz81(rs)` — ε_xc per electron in Ha
- `vxc_pz81(rs)` — V_xc in Ha (LDA potential)
- `kinetic_energy_shells(N, L)` — T_s from discrete Gamma-point shell sum

## Planned experiments (not yet implemented)

- **02_kick_tddft**: uniform δ-kick → dipole response → absorption spectrum
- **03_wavepacket**: Gaussian localised kick, compare to uniform
- **04_rs_sweep**: vary r_s ∈ {2,3,4,5,6}, compare plasmon peak to ω_p

## Library migration trigger

When the same `jellium_utils.hpp` functions appear in ≥ 2 experiments,
consolidate into `ResearchProject/jellium/lib/jellium_utils.hpp` and add
Catch2 unit tests.

## Sources

- Perdew & Zunger (1981), PRB 23, 5048 — PZ81 LDA parametrisation
- Angelo's report: Cavendish candidate 3221L — system parameters and reference values
- INQ documentation: `docs/inq_tutorial.md`, `docs/inq_source_map.md`

## Validation status

| Test | Proposed | Approved | Run | Outcome |
|---|---|---|---|---|
| 01 ground-state benchmark | yes | pending | no | — |
