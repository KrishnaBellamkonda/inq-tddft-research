# Plan: phase-2 completion — autonomous overnight run

Opened 2026-06-10 (grilling session). Executes the **remaining phase-2**
revitalisation tasks unattended, one after another, behind a validation loop,
ending with bit-identical complex-run validation. Survives compaction — the
executing session reads this + `docs/handovers/inqview-restructure.md`.

## Preconditions (MET)

- 4-package structure complete (io/analysis/visualisation/pipeline); suite green.
- **inqkit behaviour-preservation confirmed** by 3 GPU sanity runs vs golden:
  `run_free_wp_L50_E100` Δ=0 · `run_free_wp_L50_E100_sigma1` Δ=0 · jellium
  `run_wp_n162_L30_E300_highdens_sigma1` ≈1e-13 (FP). (NVML mismatch is NOT a
  blocker — CUDA compute verified working; memory `reference_gpu_driver_mismatch`.)

## Operating rules (locked in grill)

- **Validation gate after every step.**
  - Python steps: `pytest python/tests python/inqview/pipeline/test_lindhard.py`
    + `compileall python/inqview` + deps-clean probe (`import inqview.analysis`
    pulls no matplotlib/vtk). Must stay **green** (currently 86 pass / 5 xfail /
    1 xpass — the xfail/xpass are pre-existing Lindhard, treat as baseline).
  - inqkit C++ steps: pure-tier `ctest -L pure` + a GPU build check; the
    feature's engine test for behaviour (observables byte-identical).
  - Complex-run steps: bit-identical comparison of `raw/observables/*.csv` to
    backed-up golden (np.allclose rtol1e-6; report max|Δ|).
- **Failure policy:** each step is a git checkpoint. On a red gate → bounded
  auto-fix (≤2 attempts) → if still red, **revert the step** (restore green),
  log a deferred-blocker with diagnostics, **continue to the next independent
  task**. The tree is NEVER left red.
- **Git:** commit each green step to `unit-tests/inq-stack` with a scoped
  `action(scope): …` subject + body per `.claude/rules/commit-messages.md`
  (no claude/ai/anthropic; no Co-Authored trailer). **No push.**
- **Autonomy:** no user questions mid-run (all forks resolved here). Long
  builds/GPU runs go to background tasks that re-invoke on completion.
- A running **morning report** is appended to this file: done / deferred / blockers.

## Task list (ordered; one after another)

### Step 0 — commit the validated baseline (do first)
106 uncommitted paths (the whole migration + glue + inqkit refactor + docs).
Commit as scoped, rule-compliant commits to give the revert-policy a clean base:
- `refactor(inqkit): behaviour-preserving compute()-split writers + grid_layout`
  (validated bit-identical by the 3 sanity runs).
- `refactor(inqview): 4-package split (io/analysis/visualisation/pipeline) + lazy init`
- `feature(inqview): energy_components + wp_integrity renderers/glue + tests`
- `docs: ADRs 0003-0006, observable spec, overnight plan, catalogue updates`
Then the overnight steps commit individually on top.

### Group A — finish new-feature glue (pure Python)
- **A1** `plasmon_spectrum` 3d_binned mode + a field-frame/VTI loader + test
  (undamped-plasmon phasor → δ-peak; 1/q² scaling already covered).
- **A2** `center_of_density` from-run: field-frame loader + `compare(total,wp)`
  on a real frame; test the python-COD vs inqkit-CSV **dx/2** offset (documents E04).
- **A3** `gs_projected_occupations` t=0 identity test (IV-M09): `n_i^GS(0)=f_i(0)`.

### Group B — Step-2 phase splitting (pure Python)
- **B1** split `pipeline/wake.py` → bath math to `analysis.wake`, `shared_clim`
  + movies to `visualisation.wake`, thin `pipeline` run. Retarget `test_wake`.
- **B2** split `pipeline/energy_balance.py` → `analysis` compute + `visualisation`
  render (the energy_components renderer already exists; wire the band-sum ledger).
- **B3** split `pipeline/kl_divergence.py` → pure helpers to `analysis`, render to
  `visualisation`. Retarget `test_kl_divergence`.

### Group C — cleanup
- **C1** `config.py`/`defaults.py`: confirm superseded by `visualisation.style`;
  fold needed bits in, cut the rest from the public surface.
- **C2** `email.py` → internal `_notify` (out of the public API / lazy map).
- **C3** `visualisation/vti.py` verify-then-cut: confirm C++ owns VTI writing and
  no live consumer needs the Python writer; cut if dead, else keep + document.
- **C4** relocate `report1/**` + `scripts/**` out of the importable package
  (applications, not library); fix any import fallout.
- **C5** purge the 28 stale `.understand-anything/.trash-*` KG JSON files from
  the package.
- **C6** update CLAUDE.md inqview section → 4-package structure (drop the stale
  flat-module map + dead `features/python-paraview` branch ref).
- **C7** rename `inqview/pipeline/pipeline.py` → `runner.py` (kill the
  `pipeline.pipeline` double-name); update `pipeline/__init__` + the `postprocess`
  shim re-export.
- **C8** minimal `inq-stack/README` refresh + a short inqview public-API section
  (4-package entry points + lazy top-level). Full polish = phase 3.

**Decision (2026-06-10):** KEEP the `inqview.postprocess` deprecated shim
(back-compat for ~30 run analyse.py); removing it = a phase-3 migration task.

### Group D — inqkit features (C++ / GPU)
- **D1** current+dipole as `Vec3` in `observables_writer.hpp` + `StepContext`
  (CSV columns UNCHANGED — back-compat). Engine test: observables byte-identical
  to a known vector; pure `Vec3` parts in `test_vec3.cpp`.
- **D2** N-dim `plane_screen.hpp` (axis x/y/z) + time-averaged variant. Engine
  test: extract a known slice along each axis; `Σ_t ρ·dt/T` == mean of constant frames.

### Group E — complex-run validation (GPU finale)
Back up each golden `raw/observables/*.csv` first; rebuild against current
inqkit; run; compare bit-identical. Use **both GPUs** for concurrency.
- **E1** coronene `run_cc_bond` (LEED screens + overlap; golden REPORT).
- **E2** L30-highdens **E50** matched pair: `run_classical_n162_L30_E50_highdens`
  + `run_wp_n162_L30_E50_highdens_sigma1` (fast, ~15 min each).
- **E3** one true low-E L50 matched pair (E20 or E25): classical + WP, concurrent.
- Each must reproduce golden to FP noise (non-interacting bit-identical). Any
  deviation beyond tolerance = blocker (the refactor changed physics) → STOP E,
  report (this is the one place a failure is not auto-deferred — it's a finding).

## Out of scope (phase 3, designed not built)
- Minimum-observable-set manifest + validator (ADR 0006, spec
  `docs/observables/minimum-set-spec.md`). Designed in this grill; implemented
  in phase 3.

## Morning report
_(appended during execution: completed steps, commits, deferred blockers, complex-run table.)_
