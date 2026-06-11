# Plan: deferred-changes cleanup (release-grade, validated)

Opened 2026-06-11 (grilling session). Completes the items DEFERRED by the
overnight run (`docs/plans/inqview-phase2-overnight.md` FINAL SUMMARY) to reach a
**clean, release-ready** codebase. Survives compaction.

## Workflow (locked in grill)
- **Branch:** `deferred-cleanup/inq-stack` off the current committed HEAD — so
  `/code-review` sees ONLY this effort's diff. The user reviews ONCE at the end.
- **Autonomous + gated**, commit each green step, defer-on-failure (bounded
  auto-fix → revert → log → continue). Tree never left red. NO push.
- Live morning-report log appended here.

## Scope (everything deferred; observable-set impl stays phase-3)
Migration removes the back-compat shims for a CLEAN end state (user: "I'm going
to release this package — it needs to be clean"). Validation = import-resolution
check over ALL run dirs + RUN 1–2 `analyse.py` on existing data (not all 90).

## Group P — low-risk pure-Python cleanup
- **B2** split `pipeline/energy_balance.py`: extract the pure energy-ledger compute
  → `analysis/energy_balance.py` + a test (ΔE_WP/ΔE_bath/unaccounted on a known
  observables frame); render stays in pipeline/visualisation. Gate: pytest.
- **C1** dispose `config.py`/`defaults.py`: fold the value `plots.py` still needs
  (`DEFAULT_THEME`) into `visualisation.style`; retarget `plots`; cut the rest
  from the public surface. Gate: pytest + import-check.
- **C3** `visualisation/vti.py` verify-then-cut: `paraview`/`defaults` use it, so
  keep the Python writer but DOCUMENT it as legacy (C++ owns VTI); only trim dead
  exports. Gate: pytest.

## Group M — the migration (clean end state; removes shims)
1. **Sweep 90 run `analyse.py` + 8 report1/scripts importers**: `inqview.postprocess.X`
   → `inqview.pipeline.X` (submodules: pipeline→runner, density_fourier, wake,
   lindhard, spectral_weight, compare, _common); `inqview.report1/scripts` →
   their new homes.
2. **C7** rename `pipeline/pipeline.py` → `pipeline/runner.py`; update
   `pipeline/__init__` + all `from . import pipeline as _pipeline` → `runner`.
3. **C4** relocate `inqview/report1/**` + `inqview/scripts/**` →
   `inq-stack/python/applications/**` (out of the importable package).
4. **Remove** the `inqview/postprocess/` shim.
- **Validation (the key concern):**
  - new `test_import_resolution`: a script asserting every run dir's `analyse.py`
    imports resolve under the new layout (AST/import probe — cheap, deterministic).
  - **RUN 1–2 `analyse.py` end-to-end on existing run data** (e.g. a completed
    jellium WP + a coronene run) → confirm post-processing still produces REPORT.md
    + figures. (User: "a few checks is good enough.")
  - pytest + compile + deps-clean green.

## Group F — inqkit C++ features (GPU)
- **D1** current/dipole as `inqkit::detail::Vec3` (consistency with
  `center_of_density`, which already returns Vec3). CLEAN: assignment is
  centralized in `real_time/real_time_session.hpp` (NOT in run.cpp).
  - add `operator[]`/`operator[] const` to the PURE `detail/vec3.hpp` (no INQ dep)
  - `StepContext.current/dipole` → `inqkit::detail::Vec3`
  - convert `inq::vector3 → Vec3` at the 1 callback site in real_time_session.hpp
  - `observables_writer` reads `ctx.current[0..2]` unchanged → CSV byte-identical
  - update `test_observables_writer_engine` + `test_vec3` (operator[] cases)
- **D2** N-dim `plane_screen.hpp`: add `axis` (0=x,1=y,2=z, default 2 = back-compat)
  + a time-averaged variant. Generalise the `iz_nearest` FFT-wrap per axis.
- **Validation bar (locked):** engine tests (new behavior + byte-identical
  round-trip assertions; existing tests still green) **+ one bit-identical
  production re-run each** vs golden — D1: any current/dipole run; D2: coronene
  `run_cc_bond` (uses plane_screen; z default ⇒ must stay bit-identical).

## Order & gates
P (fast Python) → M (migration, the risky sweep) → F (C++, GPU). Each step gated;
Python gate = pytest+compile+deps-clean; migration gate adds import-resolution +
analyse.py run; C++ gate = engine ctest + bit-identical re-run. Commit per green
step. END: user runs `/code-review` on the branch.

## Out of scope
Minimum-observable-set manifest+validator (ADR 0006) — phase-3.

## Morning report
_(appended during execution.)_
