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

## Morning report (branch deferred-cleanup/inq-stack)

### Group P — DONE (commits b4a5b2c, 7c363b7, c5824e6)
- **B2** split energy_balance ledger → `analysis/energy_balance.py` + 3-case test.
- **C1** dropped legacy config/defaults from the public API (kept internal).
- **C3** marked legacy `vti.py` as deprecated (still used by paraview/defaults).
- Gate green throughout (94 pass).

### Group M — DONE + VALIDATED (commits 471b96f, 7260986, 96ff45c, e3ec45e)
- **M1/C7** renamed dispatcher `pipeline.py`→`runner.py` + all internal refs.
- **M2** swept 90 run files `inqview.postprocess.*`→`inqview.pipeline.*`.
- **M4** removed the postprocess shim entirely (clean public surface).
- **C4** relocated report1/scripts → `applications/` (excluded from the wheel);
  rewrote 51 figure files + 8 ResearchProject importers (a wrong-path miss caught
  by the import-check, fixed in e3ec45e).
- **Validation (key concern):** import-resolution check over all run dirs (11
  distinct inqview/applications targets, ALL resolve) + **ran 2 real analyse.py
  end-to-end on existing data** — jellium WP (all phases incl. energy_balance =
  B2) and coronene (incl. screens) — both wrote REPORT.md, exit 0.

### Group F — in progress
- **D1 (Vec3) DONE + engine-tested:** `Vec3::operator[]`; `StepContext.current/
  dipole` → `inqkit::detail::Vec3`; convert at the 1 callback site
  (real_time_session.hpp); observables_writer unchanged. Engine test
  `test_observables_writer_engine` PASSED (9 assertions — CSV byte-identical).
  Pure `test_vec3` operator[] case added.
- **D2 (N-dim plane screen):** generalised `PlaneScreen` to axis 0/1/2 (default
  z = byte-identical; axis loop outermost) + `TimeAveragedScreen`. Engine test
  (x/y/z symmetry + time-average) building.
- **D1 + D2 committed** (1a8d1d8, a125e86), both engine-tested (D1 9 assertions
  CSV-byte-identical; D2 811 assertions). 
- **Coronene production re-run BIT-IDENTICAL ✓** (wall 2884s): D1 observables.csv
  max|Δ|=**0.00e+00**; D2 all 20 LEED screens worst|Δ|=**0.00e+00** (z default).
  Validation bar fully met (engine tests + bit-identical production re-run).

## FINAL SUMMARY — deferred-cleanup COMPLETE (13 commits)
Everything deferred is DONE on `deferred-cleanup/inq-stack` (off post-overnight
HEAD; `/code-review` scope = this branch's diff). Gated + committed per green step.
- **Group P:** energy_balance split + test; config/defaults & vti out of the
  public surface.
- **Group M (clean release):** dispatcher renamed `runner`; 90 run analyse.py +
  8 importers swept `postprocess`→`pipeline`; **shim removed**; report1/scripts →
  `applications/` (excluded from wheel). Validated by import-resolution (all
  resolve) + 2 real analyse.py runs (jellium + coronene, exit 0).
- **Group F:** D1 current/dipole→Vec3 (CSV byte-identical), D2 N-dim plane screen
  + TimeAveragedScreen. Engine tests pass; coronene production re-run Δ=0.
- **Nothing deferred this round.** The minimum-observable-set manifest+validator
  (ADR 0006) remains the only parked item, intentionally (phase-3).
- **READY for the user's single `/code-review`.**
