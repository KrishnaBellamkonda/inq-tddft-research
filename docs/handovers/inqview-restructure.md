# Handover: inqview restructure + test suite

## Current status
Design/architecture phase of the inqview (Python) rejuvenation. The inqkit
(C++) round is COMPLETE (see `docs/handovers/inqkit-unit-testing.md`, 26 test
files green). This task is its successor.

**Architecture locked (ADRs 0003/0004/0005 + CONTEXT.md). ALL method TODOs
resolved** (`docs/validation/inqview-findings.md`, IV-M01..M12 + IV-E01..E03).
Three validation dossiers ran (loss function, FFT normalization, FFT
drift-removal). **Planning phase COMPLETE. Test harness now STOOD UP.**

### Milestone — harness built (2026-06-10)
- `inq-stack/pyproject.toml`: added `[tool.pytest.ini_options]` (testpaths,
  markers analysis/io/integration/theme) + `[project.optional-dependencies]`
  (analysis=scipy/pandas, viz=matplotlib/vtk/imageio, test=pytest) — encodes
  the ADR-0003 deps layering.
- `inq-stack/python/tests/`: `conftest.py`, `_signals.py` (Tone/make_tone),
  `test_fourier.py`, `test_kl_divergence.py`, `test_screens_io.py`,
  `test_fields.py`, `test_wake.py`, `test_deps_clean.py`.
- **Run (venv pytest): `python/tests/` = 29 passed, 2 xfailed** (IV-E03
  coherent-gain + ADR-0003 deps-clean, both strict-xfail captured bugs).
- **`docs/validation/test-catalogue.md`** — full index of ALL tests (inqkit
  pure 26 + engine 19, inqview 6 modules + lindhard) with what/asserted per test.

### Milestone — free-WP integration test + features backlog (2026-06-10)
- `inq-stack/tests/cpp/engine/test_free_wp_engine.cpp` (IV-M11): free Gaussian
  WP, non-interacting (empty ions + ghost via extra_electrons(2.0), WP in the
  extra state), L=20/σ=2/k₀=0.8ẑ/T=3 a.u. Asserts analytic free-particle laws:
  ⟨z⟩=k₀T, Var=σ²/2+T²/(2σ²), ⟨p⟩=k₀, norm & E_kin conserved + qualitative
  moved/spread. Registered in engine `CMakeLists.txt`. **BUILT + PASSED**
  (ctest, 3.5 s, 0 failures, first run) — confirms the WP injector is
  minimum-uncertainty and non-interacting propagation reproduces exact free-
  particle spreading. inqkit engine tier now 18 files / 20 cases.
- `docs/plans/new-features-backlog.md`: consolidated NEW-FEATURE ideas from
  both libraries' TODOs (inqview 9 incl. energy_components/wp_integrity/
  plasmon_spectrum/theme/subtract; inqkit 8 incl. N-dim screen/any-k-point/
  projected-occupation snapshot/Gaussian potential/norm-per-state), each with
  its proving test. This is the next-phase build roadmap.

### Milestone — first 4 features built (2026-06-10)
User selected 10 features; **4 built + tested** (testing rule #6 added: every
feature ships with a test):
- canonical theme → `inqview/visualisation/style.py` (semantic cmap roles +
  fixed-dim factory; test_theme 6 cases). **Starts the new package tree.**
- fourier `subtract=` + coherent-gain `/win.sum()` → `inqview/fourier.py`
  (IV-M12 + IV-E03; the IV-E03 strict-xfail is now a real PASS; test_fourier 17).
- `energy_components` → `inqview/analysis/energy_components.py` (Σ==E_total
  invariant; test 6). Renderer (bars/lines/GIF) still to build in visualisation.
- `wp_integrity` → `inqview/analysis/wp_integrity.py` (ipr/momentum_kl/variance
  + dataclass; test 8). From-run time-series assembly follows (VTI fixtures).
- New packages created: `inqview/analysis/`, `inqview/visualisation/` (with
  minimal `__init__`). `inqview/__init__` NOT yet cleaned → deps-clean test
  stays xfail until the migration finishes that.
- **inqview suite: 56 passed, 1 xfailed** (was 2 xfail; IV-E03 flipped green).

### Milestone — remaining features (2026-06-10, batch 2)
**All 7 inqview features DONE + tested** (suite now **68 passed, 1 xfailed**):
+ `plasmon_spectrum` (`analysis/plasmon_spectrum.py`, complex-FFT + |n_q|²/q² +
  axial extraction; fixes IV-E01/E02; test 4)
+ KL drift-rate `kl_series` (in `wp_integrity`) + `(k,t)` carpet renderer
  (`visualisation/carpets.py`, untested per IV-M10); wp_integrity test now 11
+ Python `center_of_density` (`analysis/center_of_density.py`, node convention,
  WP/total/bath compare, E04 dx/2 cross-check; test 5)
- Two gotchas fixed mid-build: submodule-vs-function name shadowing in
  `analysis/__init__` (don't re-export bare `center_of_density`); plasmon axial
  layout is x-slowest/z-fastest (F[0,0,m] = z-mode).
**inqkit (1 of 3):** norm-per-state → `observables/state_norm_writer.hpp` +
`test_state_norm_engine.cpp`. **BUILT + PASSED** (ctest 2.67 s). Engine tier now
19 files / 21 cases.
**inqkit REMAINING (2) — deferred as a focused non-breaking task** (specced in
backlog): current+dipole as Vec3 modifies `io/observables_writer.hpp` +
`StepContext` (current/dipole are ALREADY vector-accessed `ctx.current[0..2]`,
so the change is mostly type-naming + `ObservableSelection` ergonomics — but it
is depended on by EVERY production run.cpp + `test_observables_writer_engine`,
so it needs backward-compat + a full engine-suite regression). N-dim plane
screen modifies `screens/plane_screen.hpp` (coronene LEED depends on it).
Not rushed at session end to avoid a production regression.

## FINAL STATUS (2026-06-10): 8 of 10 selected features BUILT + verified
- inqview (7/7): theme, fourier subtract=+gain, energy_components, wp_integrity,
  plasmon_spectrum, kl_series+carpet, center_of_density. Suite 68 pass / 1 xfail.
- inqkit (1/3): norm-per-state (engine, passed). 2 remaining specced + deferred.
- `test_fourier.py` is analytic (tone → exact spectrum); independently
  re-confirms agent findings (boxcar amplitude exact, zero-pad invariance,
  DC-hijack = todo #2 bug) and pins **IV-E03 as `strict xfail`** (flips
  red→green when `/win.sum()` lands). `test_kl_divergence.py` validates the KL
  math underpinning IV-M05.
- Tests import CURRENT module paths (`inqview.fourier`,
  `inqview.postprocess.kl_divergence`); assertions are migration-invariant —
  only import lines change when modules move to `inqview.analysis`.

### Milestone — GPU sanity gate DEFERRED (2026-06-10)
Before the migration the user wanted 2 WP-in-vacuum + 1 WP-in-jellium sanity
re-runs to confirm the inqkit characterization edits are behaviour-preserving on
real production runs. **Blocked by an environment fault** and DEFERRED by the
user. Details:
- The inqkit header edits are behaviour-preserving refactors (verified by diff):
  `compute()`-split observables writers emit identical CSV columns; `fft_shift_
  index` extracted unchanged into `detail/grid_layout`; `rho_bath`→`rho_system`
  rename only. Low numerical risk; the re-runs were to confirm end-to-end.
- **GPU blocker: NVML "Driver/library version mismatch"** — loaded kernel module
  535.288.01 vs userspace/DKMS-installed 535.309.01 (unattended update, no
  reboot). No userspace workaround (old libs deleted). Fix needs root (kill stale
  GPU kernel pid 1816189 = user's 33-day idle weather-climate ipykernel → stop
  lightdm → reload modules, or reboot). Assistant has no sudo. Recorded in
  memories `reference_gpu_driver_mismatch` + `feedback_gpu_default_expectation`.
- **Chosen run menu (parked until GPU back):** vacuum `run_free_wp_L50_E100`
  (σ=5, 462 steps, golden norm 0.99999999 + python_toy analytic) + `run_free_wp_
  L50_E100_sigma1`; jellium `run_wp_n162_L30_E300_highdens_sigma1` (E=300, 190
  steps, loads GS `checkpoints/gs_L30_cubic_N162_dx0p40`, golden observables.csv).
  Back up golden CSVs before each re-run; compare via `np.allclose`; finish all
  then report a consolidated table.

### Milestone — MIGRATION STARTED: keystone landed (2026-06-10)
Package migration (ADR-0003 step 3) is underway. First, behaviour-preserving move:
- `git mv inqview/fourier.py → inqview/analysis/fourier.py` (numpy/scipy/pandas
  only — clean for the analysis layer). `analysis/__init__` now exports
  `FourierResult/FourierTransform/WindowSpec`; `plots.py` type-hint import
  retargeted to `.analysis.fourier`.
- **`inqview/__init__.py` rewritten as a LAZY (PEP 562 `__getattr__`) module** —
  `_LAZY_EXPORTS` maps every public name to its submodule and imports it only on
  attribute access. So `import inqview` / `import inqview.analysis` pull in NO
  matplotlib/VTK, while `from inqview import X` and `inqview.plot_*` still work.
- **deps-clean invariant now ENFORCED** — `test_deps_clean` flipped from strict
  xfail to a live parametrized test over (`inqview`, `inqview.analysis`,
  `inqview.analysis.fourier`); all import matplotlib/VTK-free. `test_fourier`
  retargeted to `inqview.analysis.fourier`.
- **Suite: 71 passed, 0 xfailed** (was 68 pass/1 xfail). Lazy API smoke-tested:
  analysis clean, `from inqview import RealField3D/FourierTransform` clean,
  `plot_energy_vs_time` lazily loads matplotlib on access, bad name → AttributeError.

### Milestone — 3 of 4 layers migrated (2026-06-10)
Bottom-up package moves, suite GREEN (71) + byte-compile clean after each:
- **io** (numpy-only): `git mv fields.py→io/fields.py`, `data.py→io/data.py`,
  `screens.py→io/leed.py` (rename resolves the loader-vs-`postprocess.screens`
  phase collision). New `io/__init__`. `io/leed.py` ifft import depth fixed
  (`..postprocess._ifft`). Consumers retargeted: vti/paraview/defaults/plots +
  lazy map + TYPE_CHECKING + tests (`inqview.io.fields`, `inqview.io.leed`) +
  4 report1 figs + 2 scripts. Verified `import inqview.io` pulls no mpl/vtk.
- **visualisation** (the ONLY mpl/VTK layer): `git mv plots.py/paraview.py/vti.py
  → visualisation/`. Relative-import depth bumped (`..io`/`..config`/`..analysis`;
  `paraview→.vti` sibling kept). `visualisation/__init__` deliberately imports
  ONLY `style` (NOT vti/paraview) so `import inqview.visualisation` doesn't force
  VTK and `test_theme` stays green. Consumers retargeted: defaults + lazy map +
  TYPE_CHECKING + render_density_series script.
- Finding: **zero** postprocess phases import `inqview.plots`/`_common` — phases
  plot inline, so plots/vti were NOT entangled with the phase code.
- **Still flat at top level:** `config.py`, `defaults.py`, `email.py` (low
  priority: config/defaults superseded by `visualisation.style`; email→internal).

### Milestone — 4-PACKAGE STRUCTURE COMPLETE: postprocess → pipeline (2026-06-10)
User chose **relocate-first, split-incrementally** for the 34 phases.
- `git mv inqview/postprocess → inqview/pipeline` (34 modules + dispatcher
  `pipeline.py` + `_common`/`layout`/`test_lindhard`). Canonical 4th package name.
- **Backward-compat shim** `inqview/postprocess/__init__.py`: sets
  `__path__ = inqview.pipeline.__path__` and re-exports `run/PHASES/PipelineResult`
  + emits `DeprecationWarning`. So existing run `analyse.py`
  (`from inqview.postprocess import pipeline/density_fourier`, deep
  `from inqview.postprocess.kl_divergence import _kl`) keep working unchanged —
  VERIFIED they forward to the same pipeline/ source. ~30 ResearchProject run
  dirs NOT broken.
- Internal refs retargeted to `inqview.pipeline`: `io/leed.py` (`..pipeline._ifft`),
  `test_lindhard` (path + import), `tests/test_wake`, `tests/test_kl_divergence`,
  `pyproject` testpaths. report1 had 0 postprocess importers.
- **Suite GREEN: 78 passed, 5 xfailed, 1 xpassed** (the xfail/xpass are
  PRE-EXISTING lindhard dynamical-sign characterizations, deferred — not from
  this change; `python/tests` alone still 71 passed). Package byte-compiles.

**4 packages now exist & are clean:** `io` (numpy), `analysis` (numpy/scipy),
`visualisation` (mpl/VTK), `pipeline` (orchestration; holds the still-unsplit
phases). deps-clean enforced. Lazy top-level API + `postprocess` shim both intact.

### Milestone — energy_components RENDERER wired up (2026-06-10)
User prioritised "wire up new-feature renderers/glue". First (their detailed
energy-flow ask) DONE:
- `inqview/visualisation/energy_components.py`: `render_initial_vs_final_bars`
  (grouped bars, Ha), `render_flow_lines` (ΔE(t) per component, eV — the flow),
  `render_breakdown_gif` (animated ΔE bars → GIF via PillowWriter). Consumes the
  `EnergyComponents` dataclass; NEVER recomputes. Canonical theme (ADR-0004).
  Shared `COMPONENT_COLORS` so a component keeps its colour across all 3 views.
  Exported from `visualisation/__init__`.
- `tests/test_energy_components_render.py` (4 cases) — **data-contract** test
  (rule #6 + ADR-0005): bar heights == `breakdown(...)`, line y-data ==
  `dE_*·HA_TO_EV`, GIF writes a non-empty file. Catches recompute/mis-wire
  without pixel comparison. Catalogue row added; deps-clean row updated.
- **Suite: 82 passed, 5 xfailed, 1 xpassed.** analysis still mpl/VTK-clean.

### Milestone — wp_integrity from-run assembly (2026-06-10)
Verified real run output formats first (E300 run):
`momentum_distribution.csv` = long `step,time_au,k_bohr_inv,n_total,n_wp` (WP
spectrum = `n_wp` per step); `wp_real_space_stats.csv` has `sigma_x2/y2/z2`;
**no WP-only density VTI is saved** (only total/system/delta).
- `analysis/wp_integrity.py::assemble_from_run(run_dir)`: reads the two CSVs →
  `kl_mom` (kl_series of per-step n_wp vs initial/previous) + `sigma_r`
  (√Σσ²). `ipr`=NaN (WP density not saved; run-vintage system-vs-bath ambiguity,
  memory `reference_canonical_bath_density`). Deps-clean (numpy/pandas).
  Exported as `analysis.assemble_wp_integrity`.
- `tests/test_wp_integrity_from_run.py` (4 cases): synthetic run dir, analytic
  KL=½ln2 + σ_r=√1.5 known up front. Catalogue row added.
- **Suite: 86 passed, 5 xfailed, 1 xpassed.** analysis still mpl/VTK-clean.

### REMAINING — Step 2 (incremental split) + glue + cleanup
- **Split high-value phases** into `analysis.compute()` + `visualisation.render()`
  + thin `pipeline.run()`: wake (bath math→analysis, shared_clim/movies→viz),
  energy_balance, density_fourier (compute already mirrored in
  `analysis.plasmon_spectrum`), kl_divergence (pure helpers→analysis). When a
  phase splits, retarget its test import (`pipeline.X`→`analysis.X`/`viz`).
- **New-feature renderers/glue:** energy_components bars/lines/GIF; wp_integrity
  from-run assembly; plasmon 3d_binned + VTI loader; COD VTI loader + real-run
  inqkit-CSV cross-check; gs_projected_occupations t=0 identity test.
- **Top-level leftovers:** `config.py`/`defaults.py` (superseded by
  `visualisation.style`; decide keep/move/cut), `email.py` (→ internal `_notify`).
- **Legacy cuts:** `visualisation/vti.py` (verify-then-cut; C++ owns VTI),
  verify `visualisation/paraview.py` still used; relocate `report1/`+`scripts/`
  applications out of the importable package.
- **2 inqkit features** (Vec3 current/dipole, N-dim plane screen) — engine work,
  GPU-gated. **GPU sanity gate** (2 vacuum + 1 jellium) still parked on driver fix.
- **END:** user triggers `/code-review ultra` on the full suite (IV-M10).

## What changed (this session)
Ran `/grill-with-docs` on the inqview reviews. Locked, via AskUserQuestion:
1. **Package split** (ADR 0003): `inqview.io` / `analysis` / `visualisation`
   / `pipeline`; **deps-clean invariant** — `analysis`/`io` import no
   matplotlib/VTK (testable).
2. **Scope pruning**: `email.py`→internal/relocated; `vti.py` Python writer→
   legacy verify-then-cut; `paraview.py`→verify still used; `config.py`/
   `defaults.py` generic theme→superseded.
3. **Canonical theme** (ADR 0004): promote `report1/_shared_style.py` into
   `inqview.visualisation.style`; **semantic cmap roles**
   (sequential→inferno, diverging→RdBu_r, phase→twilight); **fixed-dimension
   factory** (one-col 3.5×3.0, two-col 7.0); individual plots only.
4. **Data contract**: `compute(...) -> frozen dataclass (numpy+units)`;
   `render_*(result) -> Figure`.
5. **Test strategy** (ADR 0005): whole suite **pure-tier, portable** (no GPU/
   INQ, fixtures <5 MB, `np.allclose` not bit-exact). Combination of
   analytic (reduced systems) + characterization (I/O golden) + a
   **free-space-WP** integration anchor (analytic σ(t)/⟨p⟩/centroid).
   **Renderers untested**; only a numeric theme-config test.
6. **Finding**: the "duplicate `screens.py`" is NOT a dup — `inqview/screens.py`
   is the `LeedPattern` loader (→io), `postprocess/screens.py` is the render
   phase (→pipeline). Same for `overlap.py`. Name collision to fix, not delete.

## Files touched
- `/local/data/public/skcb2/tddft/CONTEXT.md` — glossary: 4 packages,
  deps-clean, theme/roles/fixed-dim, suite portability.
- `/local/data/public/skcb2/tddft/docs/adr/0003-inqview-package-split.md` (new)
- `/local/data/public/skcb2/tddft/docs/adr/0004-inqview-canonical-visualisation-theme.md` (new)
- `/local/data/public/skcb2/tddft/docs/adr/0005-inqview-test-strategy.md` (new)
- `/local/data/public/skcb2/tddft/docs/plans/inqview-restructure-and-tests.md` (new — full plan, module mapping, test matrix, OPEN TODOs)
- `/local/data/public/skcb2/tddft/docs/code-revitalisation/inqview-todo-catalogue.md` (input, from prior session)

## Commands run
Read-only exploration only (find/sed/du over `inq-stack/python/inqview/` and
`ResearchProject/.../results/`). No code changed, no tests run, no git ops.

## Tests and validation
None run — design phase. The inqview suite does not exist yet. First tests to
lock (against current pre-move code) per the plan: `test_deps_clean`,
`test_theme`, `test_fourier`, io-parsing, free-space-WP integration.

## Trusted sources used
- `report1/_shared_style.py` + `.claude/skills/report-figures/SKILL.md` +
  `docs/reports/report1/figures/global_style.md` — the designed viz standard.
- Memories: `reference_fixed_dimension_plot_pitfalls`,
  `reference_loss_function_method`, `reference_canonical_bath_density`.

## Attribution notes
Theme values (inferno/RdBu_r, 3.5×3.0 fixed-column scheme, dpi/fonts) are the
user's own designed standard from the report1 figure work; credit
`report1/_shared_style.py` near the promoted `inqview.visualisation.style`.

## Known issues / blockers
- Fixture-build mechanism for the pipeline orchestration smoke test is OPEN
  (real runs are 145 MB–27 GB; must trim or lean on the free-space-WP combo).
- Apparent earlier contradiction (theme "must be tested" vs "visualisations
  don't need tests") RESOLVED: renderer untested, theme *numbers* tested.

## Assumptions still in play
- inqview test suite stays entirely pure-tier (no engine) — relies on
  committed fixtures generated once by us.
- `paraview.py` Python pipeline and `vti.py` writer may be dead; not yet
  verified — do not delete until consumers confirmed.

## Method decisions locked this session (findings IV-M01..M09, IV-E01..E03)
- IV-M01/M04 loss-fn → `PlasmonSpectrum` (axial+3d_binned, peak-locator NOT
  −Im[1/ε]); validated by agent; memory corrected. IV-E01 real-only FFT (fix),
  IV-E02 relabel.
- IV-M02 COD recomputed in Python (node convention); E04 → dx/2 cross-check.
- IV-M03 minimum set = global `(summary, observables, density)` + extras.
- IV-M05 `WPIntegrity(kl_mom, σ_r, ipr)`; free-WP analytic test.
- IV-E03 fourier window coherent-gain fix `/win.sum()`; signal-agent confirmed.
- IV-M07 `energy_components` (kinetic/H/xc/ext flow + bars + ΔE(t) + GIF)
  primary; band-sum ledger kept caveated; orbital×component deferred.
- IV-M08 pipeline strictly sequential (parallel = future).
- IV-M09 band-structure deferred to multi-k QKE runs; gs_projected kept + t=0 test.
- IV-M10 validation = code-review at END (user-triggered `/code-review ultra`);
  expected values still derived analytically up front (anti-circularity).
- IV-M11 free-WP tests SEPARATE per library: inqkit `test_free_wp_engine.cpp`
  (non-interacting, analytic free-particle laws — inqkit's first integration
  test); inqview independent free-space-WP fixture.
- IV-M12 FFT subtraction: `subtract={'initial','mean','detrend','none'}`,
  default 'detrend'; canonical per column = initial (dipole/current), detrend
  (energy); fixes todo.txt #2. Combine with IV-E03 coherent-gain fix.

3 validation dossiers in `docs/validation/`: loss-function, FFT-normalization,
FFT-drift-removal. TODO-list (todo.txt/todo_later) examined + cross-referenced
in the plan. **Planning phase 100% complete — all method TODOs closed.**

## New analysis kernels to build (with their portable tests)
`plasmon_spectrum` (peak@ω_p, 1/q²) · `center_of_density` (E04 dx/2 cross-check)
· `wp_integrity` (free-WP σ_r(t)) · `energy_components` (Σ==E_total) ·
`fourier` (coherent-gain, sinusoid==A) · `gs_projected_occupations` (t=0 identity).

## Exact next steps
1. Decide the free-space-WP pipeline-fixture build mechanism (ADR 0005).
2. Implement: stand up `inq-stack/python/tests/`; lock first tests green vs
   CURRENT code (test_deps_clean, test_theme, test_fourier coherent-gain,
   plasmon-peak, io-parsing, free-space-WP integration).
3. Execute the package migration (io/analysis/visualisation/pipeline) keeping
   the suite green; split each phase into compute+render+thin-run.
4. THEN apply deferred fixes as separate red→green: IV-E01 (FFT real-only),
   IV-E03 (coherent-gain), relabels, energy_components build.
