# Plan — CAP thin-absorber (L=5) reflectivity tuning

Date: 2026-06-15
Status: design locked (grilled 2026-06-15); awaiting launch approval + Task #7.
Owner task / handover: `docs/handovers/absorbing-boundary.md`

## Goal

Using **INQ's in-built CAP** (`perturbations::absorbing`, the team's region-restricted
sin² imaginary potential — NOT the MFA), find the CAP parameters that give a **low
reflection error ε across the 1–100 eV (10⁰–10²) domain with the minimum near 10 eV**,
under the constraint of a **thin** absorber (`L = 5` Bohr) and **shallow** depth
(`η ∈ [−0.01, −0.30]` Ha).

This is a re-parameterisation of the already-working `cap_sweep` pipeline
(`ResearchProject/systems/vacuum/scripts/cap_sweep/`), not new physics code.
ε = WP norm left in the inner region at τ (CONTEXT.md "reflection error ε(E,L)";
`run.cpp`: `eps = inner_tau / N0`). Lower ε = better absorption.

## Scientific design (the run menu)

Absorber = inq CAP over the last `L` of the box (fractional slab
`[0.5−L/Lcell, 0.5]`); propagator **ETRS** (CN renormalises → kills absorption);
geometry scales with energy (`σ=4√2/k₀`, `Lcell=6σ+L`, WP at `z₀=−L/2`,
`τ=2(3σ+L)/k₀`, `dx=clamp(0.75/k₀,0.18,0.30)`, `dt=0.01`).

| Knob | Values | Notes |
|---|---|---|
| `L` (width) | **5 Bohr**, fixed | the "thin" target |
| `η` (depth) | **−0.01, −0.05, −0.30** Ha | log-even over [−0.01, −0.30]; one ε(E) curve each |
| `E` (energy) | **1.01, 1.87, 3.48, 6.46, 7, 10, 15, 22.24, 41.28, 76.62, 100 eV** | MFA-comparable log ladder + densified around 10 eV |
| `k₀` | derived `√(2E/27.2114)` | 0.27 … 2.71 |
| Runs | **33** (3 η × 11 E) | — |

`dt = 0.01` fixed across all runs (ε is dt-sensitive → varying dt would break
curve comparability). Low-E runs are the costly corner (1.01 eV ≈ 50k steps,
box ≈ 131 Bohr); 10 eV ≈ 5.75k steps.

**Observables.** Every run emits the full free-WP **minimum observable set** +
`observables_manifest.json` (ADR 0006), exactly as the current `cap_sweep/run.cpp`
does. **Showcase/bulk split** (CONTEXT.md): only the **3 runs at E≈10 eV** (one per η)
write `density_wp` VTI frames (for the density GIF); the other 30 write CSV +
`epsilon.txt` only (33 × 60 3-D frames would be GBs).

**Expectation (calibration, not a target).** Prior autonomous data: L=5, η=−0.5 →
ε≈0.14. A thin 5-Bohr absorber has a *higher floor* than L=20–50 (which hit
1e-5…1e-7). So expect best ε here ~10⁻²–10⁻¹. The deliverable is the **curve shape
and where the minimum sits**, not a deep null. Physics: ε rises at low E
(wavelength ≫ L → edge reflection) and at high E (too little dwell time in a thin
absorber) → a minimum at intermediate E that η slides along the axis.

## Folder reorganisation (decided 2026-06-15 — amends ADR 0007)

Runs are **grouped by sweep**, not flat: `systems/vacuum/<sweep_name>/<run_name>/`.
Analysis stays in `hypotheses/<sweep_name>/`, renamed to **bare sweep names** (drop
the `NN_` prefix) so runs and analysis share an identical `<sweep_name>`.

Target tree (full migrate):

```
systems/vacuum/
├── cap_real/        run_cap_*        ← was flat run_cap_* (autonomous knobs study)
├── mfa_sweep/       run_mfa_*        ← was runs/ (73 MFA runs + epsilon_grid.csv)
├── cap_thin_L5/     run_cap_k*_L5_eta*  ← NEW (33 runs)
├── scripts/
│   ├── cap_sweep/   (machinery: run.cpp, dispatch.py, build/)   ← exists
│   └── mfa_sweep/   (machinery)      ← MOVE from top-level vacuum/mfa_sweep/
├── shared/  shared_gs?  tests/
└── hypotheses/
    ├── cap_real/        ← rename of 01_cap_real/ (ipynb, figs, build_cap_report.py, tests/)
    ├── mfa_sweep/        ← MFA combined CSV + analysis (new home for epsilon_grid.csv)
    └── cap_thin_L5/      ← NEW: combined CSV, build script, study ipynb, figs, tests/
```

Migration map:
- `vacuum/run_cap_*`  → `vacuum/cap_real/run_cap_*`
- `vacuum/runs/run_mfa_*` → `vacuum/mfa_sweep/run_mfa_*`; `runs/epsilon_grid.csv` → `hypotheses/mfa_sweep/`
- `vacuum/mfa_sweep/` (machinery: run.cpp, run, build/, *.log, profile.dat) → `vacuum/scripts/mfa_sweep/`
- `vacuum/hypotheses/01_cap_real/` → `vacuum/hypotheses/cap_real/`
- Fix paths in `hypotheses/cap_real/build_cap_report.py` (globs `run_cap_*/results/…`
  → `../../cap_real/run_cap_*/results/…`) and re-run to confirm the notebook still builds.

## Ecosystem reconciliation (keep the layout decision consistent)

The flat-`run_*` rule is encoded in 3 places; all must change together:
1. `docs/adr/0007-system-folder-structure.md` — add a dated **Amendment (2026-06-15)**:
   runs group under `<sweep_name>/`; `hypotheses/<sweep_name>/` uses bare names;
   jellium/coronene stay grandfathered-flat; vacuum migrated as reference.
2. `CONTEXT.md` — "System folder structure" section: `run_<type>_<params>/ (FLAT)`
   → `<sweep_name>/<run_name>/`; `hypotheses/<NN_purpose>/` → `hypotheses/<sweep_name>/`.
3. `.claude/rules/file-placement.md` — the ADR-0007 table + canonical-structure rows.

(Standard binds NEW systems going forward; jellium/coronene NOT migrated.)

## Validation gate

- **Provisional.** All ε remain provisional until **Task #7** (inq-study engine
  ctest validating the scalar-potential complexification). Same precedent as the
  autonomous `cap_real` study. Notebook + handover must carry the provisional flag.
- **Pre-launch (simulation-validation, Tier A):** the `cap_sweep` binary is already
  built + validated (autonomous run); only the parameter menu changes. Smoke = 1
  cheap run (e.g. 10 eV, η=−0.05) before the full 33.
- **Cost:** 33 runs on 2 GPUs; low-E corner dominates (~50k steps × 3). Rough ETA
  to size before launch.

## Execution checklist

1. [ ] Update CONTEXT.md + ADR 0007 amendment + file-placement.md (ecosystem).
2. [ ] Migrate existing runs/machinery per the map; fix `build_cap_report.py`; re-run to confirm.
3. [ ] New `cap_thin_L5` dispatch menu: clone `scripts/cap_sweep/dispatch.py`,
       fix `L=5`, set η={−0.01,−0.05,−0.30}, E-ladder above, run dir →
       `vacuum/cap_thin_L5/run_cap_k{..}_L5_eta{..}`, VTI only at E≈10 eV.
4. [ ] Smoke 1 run → validate manifest (`python -m inqview.validation <run>`).
5. [ ] **Get launch approval (cost gate)** → dispatch 33 runs on GPU 0,1.
6. [ ] Build `hypotheses/cap_thin_L5/`: combined ε CSV, ε(E) curves (3 η overlaid,
       log-y, min-marker near 10 eV), density GIF, study ipynb (provisional banner).
7. [ ] Update handover `docs/handovers/absorbing-boundary.md`.

## Open risks
- Thin+shallow CAP may not reach "low" ε anywhere (floor ~0.1); the study still
  answers *where the min sits* and *which η is best* — report honestly.
- Low-E runs slow; if wall-clock is unacceptable, trim the 1.01/1.87 eV points
  (keep ≥6.46 eV) — decide after the smoke timing.
