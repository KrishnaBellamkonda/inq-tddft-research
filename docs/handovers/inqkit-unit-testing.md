# Handover: inqkit / inqview unit-testing rejuvenation

## MILESTONE: CI/CD pipeline shipped — phase deliverable complete (2026-06-11)
GitHub Actions CI added (`.github/workflows/ci.yml`, commit `9f3280c` on
`test-reorg/inq-stack`). Runs the **two portable tiers** on every push/PR
touching `inq-stack/**` (+ workflow_dispatch, concurrency cancel-in-progress):
- **Python suite** — `pip install -e ./inq-stack[analysis,test] matplotlib
  pillow imageio`, `pytest -q` over `tests/python`, plus the deps-clean import
  check (`tests/python/test_deps_clean.py`). No VTK, no GPU, no INQ.
- **C++ pure tier** — apt cmake/g++/ninja, `cmake -S inq-stack/tests/include`,
  `ctest -L pure`. Catch2 v3.5.4 fetched via FetchContent (added to
  `tests/include/CMakeLists.txt`) since stock runners have no INQ build tree;
  local builds still reuse INQ's Catch2 source.
- `.gitignore`: negated the broad `.*` rule for `.github` so the workflow is
  tracked. README test paths refreshed.
**Engine/GPU tiers deliberately NOT in CI** — they need INQ + a CUDA GPU; doc'd
as local / self-hosted-runner (ADR 0001/0002). Validated locally pre-commit:
YAML parses (10 steps), pure tier builds + 10/10 ctest pass, portable pytest
needs no VTK (vtk-block test → 111 pass).
**Status of "this phase" (unit-testing/restructure):** test reorg (mirrored
`tests/include` + `tests/python`), minimum-observable-set (ADR 0006, C++
manifest + Python validator, real jellium-wp run PASSES), and CI/CD all DONE.
Remaining is user-gated (push 4 stacked branches to origin; merge) + minor
rollout (wire `write_manifest` into run templates; test-catalogue rows).
NEXT PHASE = claude-ecosystem (`task_calude_ecosystem.md`), interview-driven.

## NEXT PHASE: inqview Python tests (2026-06-10)
inqkit round done → pivot to inqview. **29 TODO comments aggregated** into
`docs/code-revitalisation/inqview-todo-catalogue.md` (themes: Φ-imports,
Φ-minimum-set, Φ-viz-rule [1a/1e/1f], Φ-redundancy [dup screens.py, vti.py],
Φ-cod-reuse [wake.py ↔ inqkit E04/bath]). Testable surface = numeric
post-processing (fourier FFT, density_fourier loss function, KL divergence,
overlap, energy/bath balance). Likely tiers: pure-numpy (CI) vs VTK/data-dependent
(needs run outputs). Full grilling-based test planning = next session; catalogue
is the input. Inqkit + inqview density/COD semantics (E02/E04) must be coordinated.

## Current status — 26 TEST FILES GREEN, inqkit ROUND COMPLETE (2026-06-10)
**26 green: pure 9 + engine 17** (incl. 1 mpirun -np 2 multi-rank test), 0
failures. ALL inqkit characterization tests + the fix round done & verified:
- Fixes: **E01** (plane_screen all_reduce), **E03** (iterated GS) — red→green.
- Refactors: **Vec3 swap** (COD→Vec3 + callers), **stats compute()→Moments**
  (both wp_momentum_stats + wp_real_space_stats), each with a direct unit test.
- Deferred (documented): **E02** (cached density t=0; proper fix = orbital-based
  `total_from_orbitals`), **E04** (COD half-cell — production trajectory bias).
- Coverage: detail (grid_layout, vec3), fields (density, orbital), io (real/
  complex/vti writers, observables), observables (COD, density_delta, eigenvalue,
  occupations, orbital_overlap, state_energy, momentum_distribution, both stats
  compute()), screens (plane_screen extract + multi-rank), wavepacket (inject/
  ortho/T02), jellium/shells. SKIPPED: 11 stubs, trivial PODs, density_delta
  coarse_grain (binary VTI), config data.
- ⚠ production .cpp Vec3 migration compile UNVERIFIED (next inq-run).
- Harness optimisation noted: test_engine_main recompiled per target (make it an
  OBJECT lib to ~halve build times).
**inqkit characterization + fix round substantially complete.** 23 test files
green (pure 9 + engine 14, incl. 1 multi-rank mpirun test). Harness: pure tier
(`tests/cpp/`) + engine tier (`tests/cpp/engine/`, links INQ via add_subdirectory,
GPU). 4 bugs found: **E01 + E03 FIXED & verified (red→green)**; **E02 + E04
documented/deferred** (production-affecting). Refactors done: **Vec3 swap**
(CenterOfDensityResult→Vec3 + callers); **stats compute()→Moments**
(wp_momentum_stats). Full design + matrix in `docs/validation/inqkit-tests.md`;
findings in `docs/validation/inqkit-errors.md`.
**Remaining (precisely scoped, see Fix-round section):** wp_real_space_stats
compute() refactor; momentum_distribution + state_energy tests; verify production
Vec3 migration compiles (next inq-run).

## (historical) TODO-driven grilling (2026-06-10). Flow pivoted: the user's
inqkit review lives as `TODO:` comments in the headers + a written review
(`docs/code-revitalisation/inqkit-review-and-next-steps.md`). I catalogued all
30 substantive TODOs as T01–T31 (`docs/code-revitalisation/inqkit-todo-
catalogue.md`) and am grilling them cluster-by-cluster to produce **specific
tests + a proposed architecture**, recorded in `docs/validation/inqkit-tests.md`.
inqkit and inqview are handled separately (parallel to the user's Python review).

**Clusters grilled & locked (2026-06-10):**
- **Θ-coord** (T01/T05/T06/#12) — DONE. Move `fft_shift_index` →
  `detail/grid_layout.hpp`; Test A (pure, gates refactor) + Test B (engine,
  off-centre He, both parities, VTI→Python). Baselines BL-coord-1a/1b.
- **Θ-density-semantics** (T02/T03/T04) — DONE. T02 direct WP-inclusion test
  (BL-dens-1, He+WP); T03 resolved by grep (function is load-bearing in jellium
  `_wf` runs → keep+document); T04 momentum-stats complex-ψ unit test (folds
  with T28).
- **Θ-parallel** (T21/T22/T25/T27) — DONE. **E01 confirmed bug:** plane_screen
  missing two `all_reduce`s (reference-correct pattern already in
  wp_real_space_stats:209-214). Fix via red→green parallel-invariance test
  (np=1 frozen baseline BL-parallel-1 vs np=2 elementwise). **E01 is LATENT** —
  the only multi-rank prod run (`run_wp_n162_L50_E700_mpi_propagate`, -np 2) uses
  wp_real_space_stats (safe), not plane_screen → no corrupted results to date.
  Logged in `docs/validation/inqkit-errors.md`. ⚠ first engine multi-rank test
  = highest harness risk; prove trivial 2-rank ctest before invariance asserts.

- **Wavepacket physics** (T04/T28/T29/T31/T26) — DONE. T04 RESOLVED by reading
  (wp_momentum_stats uses `|FFT(ψ)|²` correctly). T28+T04 discriminating twin
  test (BL-wp-1: identical-envelope WPs k₀=0 vs k₀≠0 → ⟨p⟩=0 vs k₀, + cod-slope
  group-velocity cross-check). T29 orthogonalisation (reuse BL-dens-1 w/ WP
  overlapping He). T31 data-driven from T29 residual (single vs iterated GS).
  T26 = oneoff experiment, out of CI scope.

Phase-0 baseline registry: BL-coord-1a/1b, BL-dens-1 (also yields the
parallel-invariance np=1 golden + T29 ortho report), BL-wp-1 — ~4 distinct
tiny configs, all in `docs/validation/inqkit-tests.md`.

- **Θ-vector + Θ-naming + singles** (T07/T09/T10/T11/T12/T13/T15/T18/T19/T20) —
  DONE. T07/T09 → small **pure** `inqkit::detail::Vec3` (keeps COD+observables
  pure). T10/T13/T19 renames inline under characterization. T11 (`0.0L`=long
  double literal, not a unit), T12, T20 = doc only. T15 (δn vs fixed t₀ ref —
  unit test), T18 (write_every NOT inherited, two composing gates — unit test).
- **Deferred to restructure phase:** T08, T23, T24, T26, T30, T14, T16.

## GRILLING PHASE COMPLETE (2026-06-10)
All 31 TODOs triaged into tests / resolved-by-reading / deferred. Full design in
`docs/validation/inqkit-tests.md`; confirmed bug in
`docs/validation/inqkit-errors.md` (E01).

## IMPLEMENTATION STARTED — pure tier (2026-06-10)
- **Pure harness stood up:** `inq-stack/tests/cpp/CMakeLists.txt` (Catch2 v3
  reused from `inq/build/_deps/catch2-src`, no network; `pure` ctest label;
  `inqkit_add_pure_test` helper). `tests/cpp/` was empty (greenfield).
- **Smoke test** `test_grid_layout.cpp` (flatten_index, step_suffix) — green;
  proved the pipeline before touching source.
- **T01/T05 DONE via red→green:** `fft_shift_index` → `detail/grid_layout.hpp`;
  density.hpp (6) + orbital.hpp (3) call sites re-pointed; orbital.hpp dropped
  the density.hpp include. **Test A** `test_fft_shift.cpp` (4 cases: even table,
  odd parity, origin→0, bijection) green. `ctest -L pure` = 2/2, 0.03s.
- **⚠ UNVERIFIED:** density.hpp/orbital.hpp engine compile against INQ (pure
  relocation, low risk, but not compiled in a real run.cpp yet). Rebuild one
  production run before next launch.
- **Source files changed (first inqkit edits this effort):**
  `inq-stack/include/inqkit/detail/grid_layout.hpp` (+fft_shift_index),
  `inq-stack/include/inqkit/fields/density.hpp`,
  `inq-stack/include/inqkit/fields/orbital.hpp`. NOT committed.

## Pure tier — more done (2026-06-10, cont.)
- **T15 DONE:** `tests/cpp/test_density_delta.cpp` (6 cases) — proves fixed-t₀
  base (discriminator L2=36 not 16), lazy capture, set_reference, dV scaling,
  grid-mismatch throw. Pure (emit flags off → no I/O). Green.
- **T07/T09 partial:** `detail/vec3.hpp` (new pure type) + `test_vec3.cpp`
  green; **COD characterization** `test_center_of_density.cpp` (5 cases, incl.
  anti-transposition (1,2,3), validates T11 Bohr units + T12 total_weight/guard)
  green. **Struct-swap DEFERRED** (engine callers read `cod.x_bohr`).
- **Pure tier total: 5 test files, all green, `ctest -L pure` 0.05s.**
- New source (uncommitted): `detail/vec3.hpp` (only new addition; no engine
  caller touched). Earlier: grid_layout.hpp/density.hpp/orbital.hpp (fft_shift).

## Engine harness — IN PROGRESS (2026-06-10)
- **GPU confirmed working**: CUDA runtime sees 2× A30 (sm_80), cudaMalloc OK.
  `nvidia-smi` "Driver/library version mismatch" is COSMETIC here — ignore it;
  build/run GPU-enabled.
- Engine tier lives in `inq-stack/tests/cpp/engine/` (separate cmake project so
  the pure tier stays INQ-free). Mirrors `inq-run`: `add_subdirectory(${INQ_
  SOURCE} inq EXCLUDE_FROM_ALL)` provides the `inq` target AND `Catch2::Catch2`
  (INQ FetchContents it — version-matched). Custom main
  (`test_engine_main.cpp`) copies INQ's `unit_tests_main.cpp`:
  `inq::input::environment::global()` then `Catch::Session().run()`.
- Configure flags (from `shared/config.sh`): `-DENABLE_CUDA=ON
  -DCMAKE_CUDA_ARCHITECTURES=80 -DCMAKE_CUDA_COMPILER=/lsc/opt/cuda-12.6.2/bin/
  nvcc -DPython_EXECUTABLE=<venv>/python3`. Test env: INQ_SHARE_PATH +
  PSEUDOPOD_SHARE_PATH (set via ctest ENVIRONMENT prop).
- **First engine test** `test_density_total_engine.cpp`: He atom, initial_guess,
  `density::total`, assert dims==basis + integral≈2. Doubles as the **fft_shift
  engine-compile verification** (density.hpp now uses grid_layout::fft_shift).
- **✅ DONE (2026-06-10):** configure 33.6s; build exit 0; `ctest -L engine` =
  **1/1 pass, 5.12s** on GPU. Integral ≈ 2.0 (He). **fft_shift engine-compile
  caveat RESOLVED** — density.hpp verified against INQ on GPU. Engine harness is
  now the template for all engine tests.

## Parallel workflow (2026-06-10, user directive)
Locked+implemented tests build/run in the BACKGROUND (one shared `engine/build/`
→ serialize builds, don't launch concurrent ones); we lock MORE tests in the
foreground after checks.

**Engine tests status:**
- `test_density_total_engine` — ✅ green (He, integral 2.0).
- `test_coord_mapping_engine` (Test B, off-centre He even+odd) — written;
  BUILDING (bg task b2jf2ewoi at time of writing).
- `test_density_semantics_engine` (T02/T29: He+WP, ∫total=3, bath=2, ortho
  report) — written + in CMake; QUEUED for next batch build.

**Locks added:** #2 stats testability refactor — extract `compute()→Moments`
from WPMomentumStats/WPRealSpaceStats (characterize CSV unchanged first), then
unit-test compute() for T28/T04 ⟨p⟩≈k₀ and parallel-invariance via struct
compare. See inqkit-tests.md.

**WavePacket API (verified):** `.center(x,y,z).sigma(s).k0(kx,ky,kz)
.orthogonalise_against_occupied(e[,tol]).inject_into_last_extra_state(e,occ)` →
`InjectionReport{state_index,norm_before,norm_after,max_overlap,passed_tolerance}`.
Injection step 5 sets `electrons.occupations()[0][ist_wp]=occ`. Needs
`options::electrons{}.extra_states(N)`. Grid spacing via `.spacing(0.5_bohr)`;
even/odd via L=10.0/10.5.

## ⚠ MAJOR FINDINGS from engine tests (2026-06-10) — see docs/validation/inqkit-errors.md
- **E02 (high):** `electrons.density()` returns a CACHED `spin_density_`, NOT
  refreshed by manual WP injection. Test confirmed: pre-propagation
  `∫density::total = 2` (WP excluded) → `total_excluding_orbital` double-subtracts
  (bath=1). After `real_time::propagate` refreshes it, includes WP (=3).
  **Overturns my earlier resolved-by-reading T02 guess.** **Suspected production
  bug:** `_wf` runs compute `full0 = density::total` at t=0 BEFORE propagate
  (run.cpp:179 vs 264) → t=0 bath/total frames likely over-subtracted/mislabeled;
  t>0 fine. USER ADJUDICATES the production fix (their bath definition).
- **E03 (med):** single-pass orthogonalisation gives `passed_tolerance=false` for
  a strong-overlap WP (0.97). Confirms T31: needs iterated (2nd) GS pass.
- Test `test_density_semantics_engine.cpp` REWRITTEN to encode the verified
  reality (pre/post-propagation density; E03 single-pass limitation) — should be
  green. Building now (bg bwiadplr9) with test_orbital_engine.

## Engine tests status
- ✅ test_density_total_engine, test_coord_mapping_engine (Test B both parities).
- building: test_density_semantics_engine (revised), test_orbital_engine.

## LOCKED (designs, not yet implemented)
- #2 stats testability refactor (compute()→Moments); #3 E01 parallel-invariance
  via **cross-rank agreement** (self-contained, mpirun -np 2, separate target);
  #4 T18; deferred Vec3 struct-swap + doc-comment batch.

## PLANNING PHASE COMPLETE (2026-06-10) — full plan locked
Engine suite now 4/4 green (incl. E02 verification: post-propagation density=3).
Pure 5/5. **9 test files green.**

**Locked global decisions (this round):**
- **Document all findings, FIX NONE** (E01/E02/E03 characterized, not fixed).
- **Tests only — NO source changes** this round (fft_shift move already done).
  Defer: 3 bug fixes, Vec3 struct-swap, stats compute() extraction, doc-comment
  batch (comments = source), restructure stubs (T08/T23/T24/T26/T30).
- **Full public-API coverage** — every non-trivial component (see the FULL
  COVERAGE TEST MATRIX in docs/validation/inqkit-tests.md).

**Locked test forms:**
- T28/T04 = **two-route ⟨p⟩ cross-validation** (real-space ∇ vs reciprocal ∫k|ψ̃|²,
  must agree ≈ k₀) — user-locked. + CSV class-smoke for the stats class.
- E01 = cross-rank agreement, `[!shouldfail]` (documents bug till fixed),
  separate target under mpirun -np 2.
- config/tsubonoya = SKIP.
- Known bugs characterized via `[!shouldfail]` or assert-current.

## FIX ROUND — DONE & VERIFIED (2026-06-10): E01, E03, Vec3
- **E01 FIXED + verified:** 2 all_reduce (basis+set comm) in plane_screen::extract;
  `test_plane_screen_parallel_engine` (mpirun -np 2) now PASSES (ranks agree,
  no [!shouldfail]). NB: WP injection throws multi-rank (separate limitation) —
  test uses He alone.
- **E03 FIXED + verified:** iterated 2-pass GS + residual measurement in
  wavepacket.hpp (`max_overlap`=pre-ortho, residual gates `passed_tolerance`);
  T29 case flipped to CHECK(passed_tolerance), passes.
- **Vec3 swap DONE + verified:** CenterOfDensityResult → {Vec3 center_bohr;
  total_weight}. Migrated: test_center_of_density (pure, passes), production
  run_template.hpp:433 + run_wp_e1000.../run.cpp:309 (cod.center_bohr.x).
  ⚠ production .cpp compile UNVERIFIED (mechanical change; verify next inq-run).
- **stats compute()→Moments DONE (BOTH wp_momentum_stats + wp_real_space_stats):**
  extracted `compute(electrons) const` returning a Moments struct; `accumulate()`
  = compute() + same CSV row (unchanged). Direct unit tests
  `test_wp_momentum_compute_engine` (⟨p⟩=k₀; N>0 = raw reciprocal Parseval,
  convention-dependent, NOT ≈1) + `test_wp_real_space_compute_engine` (⟨r⟩=centre;
  N≈1 real-space norm; Var≈σ²/2; node convention so no E04 offset).
- **Suite: 24 green** (pure 9 + engine 15).
- **STILL TODO (last 2, precisely scoped):**
  1. **momentum_distribution** test — MOCK Viewables {electrons()/iter()/time()};
     ctor `MomentumDistribution(csv, wp_idx, l_bohr, cfg{n_bins=64})`; inject k₀
     WP; snapshot; parse binned CSV (rows: step, k_bohr-bin, n_total, n_wp);
     assert n_wp peaks in the |k₀| bin. (Lower priority — two-route validated units.)
  2. **state_energy_writer** test — REAL short `real_time::propagate` callback
     (Viewables needs `data.ham()`); CSV cols step,time_au,kpoint,state,weight,
     occupation,E_expect_ha[,E_variance_ha2]; assert per-state E_expect_ha finite
     + row count == n_states.
- DEFERRED (document only): E02 (production t=0 frame; **user-proposed proper fix
  = orbital-based `total_from_orbitals`, see errors doc**), E04 (COD half-cell).
- ⚠ production run_template.hpp / run_wp_e1000.../run.cpp Vec3 migration compile
  UNVERIFIED (mechanical; verify on next inq-run).

## FIX ROUND PLAN (locked 2026-06-10) — "fix none" lifted
APPLY now (user-selected): **E01** (add 2 all_reduce to plane_screen::extract →
multi-rank test flips, remove [!shouldfail]); **E03** (iterated 2nd GS pass in
wavepacket.hpp → T29 flips CHECK_FALSE→CHECK(passed_tolerance)); **Vec3
struct-swap** (CenterOfDensityResult→{Vec3 center_bohr; total_weight},
current/dipole→Vec3; migrate callers cod.x_bohr→cod.center_bohr.x in
run_template.hpp + run.cpp; COD char test updated); **stats compute()→Moments**
(extract from WPMomentumStats/WPRealSpaceStats, CSV unchanged; add direct
compute() unit test). DEFER (document only): **E02** (production t=0 frame),
**E04** (COD half-cell — changes production WP trajectories). Also finish the 2
remaining characterization tests: momentum_distribution (mock), state_energy
(real RT).
**USER RULE:** every new component/feature must ship with a test.

## ⚠ E04 — half-cell coordinate bug (2026-06-10, NEW, high value)
`test_wp_real_space_replica_engine` caught it: WP injected at (1.0,−1.0,0.5)
recovered at (1.25,−0.75,0.75) = **+dx/2 on every axis**. INQ grids are
NODE-centred (`real_space.hpp:127-129` `rvector=symmetric_coord·dx`), so correct
field coord = `origin+ix·dx`. **`center_of_density.hpp:64` uses `(ix+0.5)·dx`
→ WP centroid biased by +dx/2 (~0.25 Bohr).** `wp_real_space_stats` uses INQ
rvector (node) → correct; the two centroid codes disagree. Test B's 1-cell
tolerance + the self-consistent COD pure test missed it. See E04 in errors doc.
Fix deferred ('fix none'); test corrected to node convention.

## IMPLEMENTATION PROGRESS (2026-06-10) — 18+ test files green
**Engine green (9):** density_total, coord_mapping (Test B), density_semantics
(T02/E02), orbital, observables_writer, plane_screen extract, eigenvalue_dump,
wp_momentum **two-route ⟨p⟩** (T28/T04: reciprocal=1.0000=k₀ exact, FD agrees
within discretization), wp_real_space (⟨r⟩/var, node convention).
Building (bxp9yy2z1): occupations_writer (mock Viewables), orbital_overlap (t=0
identity). Remaining: momentum_distribution (mock), state_energy (real RT), leed,
E01 multi-rank ([!shouldfail]).

## Findings tally: E01 (plane_screen Allreduce), E02 (cached density t=0),
## E03 (single-pass GS), E04 (COD half-cell). All documented, none fixed.

## OLD progress notes:
**PURE tier — COMPLETE, 9/9 green:** grid_layout, fft_shift (Test A), vec3,
center_of_density, density_delta (T15), real_field_writer (round-trip),
complex_field_writer (round-trip), vti_writer (T06 reorder), jellium_shells
(162 closure). `coarse_grain` SKIPPED (static-private, binary VTI).
**ENGINE tier — green: density_total, coord_mapping (Test B), density_semantics
(T02/E02), orbital.** Building now (bg6i1woi2): observables_writer (StepContext
POD CSV), plane_screen extract (single-rank). Prepared/queued:
eigenvalue_dump (written).

### Remaining ENGINE tests (API notes for delay-free continuation)
- **eigenvalue_dump** — WRITTEN; `dump_eigenvalues(electrons,csv)` free fn; CSV
  cols kpoint_index,kx..kz,weight,state_index,eigenvalue_ha,_ev,occupation.
- **wp_momentum_stats (T28/T04 two-route, USER-LOCKED, high value)** — inject WP
  k₀≠0; Route1 real-space `Im∫ψ*∂_dψ dV / ∫|ψ|²` via `operations::gradient`;
  Route2 reciprocal `∫k|ψ̃|²/∫|ψ̃|²` via `operations::transform::to_fourier`
  (copy k-grid loop from wp_momentum_stats.hpp:140-200); assert both ≈k₀ & agree.
- **wp_real_space_stats** — inject WP; assert ⟨r⟩≈WP centre, ⟨r²⟩ width≈σ
  (replica or short-RT CSV).
- **momentum_distribution** — injected plane-wave; |ψ̃(k)|² peak at k₀.
- **occupations_writer / state_energy_writer** — `snapshot(Viewables)` needs RT
  data → drive via short `real_time::propagate` callback; parse CSV.
- **orbital_overlap** — `OrbitalOverlapMatrix(electrons,n_ref=wp_idx,dir)`;
  snapshot writes overlap_XXXXXX.csv; t=0 identity: diag≈1, off≈0, WP col≈0.
- **leed_pattern_accumulator** — `accumulate(electrons,dt)` ×2 steps; pattern =
  Σ slice·dt (thin wrapper over plane_screen).
- **plane_screen E01 multi-rank** — separate target, `mpirun -np 2`, cross-rank
  agreement via world-comm signature (sum+sumsq, all_reduce min==max),
  `[!shouldfail]` (documents E01 till fixed). HARDEST (first multi-rank).

## Build cadence note
Engine builds are SLOW (each target links INQ, ~3-5 min) + SHARED build dir →
run as background tasks, reconfigure after adding targets, serialize.

## IMPLEMENTATION QUEUE (write tests one after another, no source changes)
PURE: io/real_field_3d_writer (round-trip), complex_field_3d_writer, vti
(T06 reorder), observables_writer (CSV cols), density_delta coarse_grain,
jellium/shells (magic numbers), leed_pattern_accumulator.
ENGINE: eigenvalue_dump, occupations_writer, state_energy_writer,
momentum_distribution, orbital_overlap (t=0 identity), wp_momentum_stats
(two-route ⟨p⟩), wp_real_space_stats (⟨r⟩,⟨r²⟩), plane_screen extract
(single-rank) + E01 multi-rank ([!shouldfail]).

## Build note
Engine build dir is SHARED (`engine/build/`) → serialize builds; **reconfigure
after adding any new target** before building (a build without reconfigure fails
with 'No rule to make target'). Add pure tests to the pure CMakeLists.

## Deferred to a later implementation round (source changes)
E01/E02/E03 fixes; Vec3 struct-swap (+cod.x_bohr caller migration); stats
compute()→Moments; iterated GS (T31); doc-comment batch; restructure stubs.
2. **Phase 0 baseline runs** (user-gated GPU): BL-coord-1a/1b, BL-dens-1,
   BL-wp-1 — freeze goldens.
3. Engine tests (Test B, T02, T28/T04, T29, parallel-invariance) + **E01 fix**
   (plane_screen Allreduce, red→green).
4. Deferred source changes once engine compile verifiable: Vec3 struct-swap;
   doc-comment batch (T11/T12/T20/T13 + T15 answer + T03 keep-rationale).

**GOVERNING RULE (user, 2026-06-10):** every baseline/reference run a test needs
is executed and frozen in **Phase 0, before any source change**. Registry of
Phase-0 runs lives in `docs/validation/inqkit-tests.md`.

Subtask 1 (mapping) substantially done earlier — `understand-anything` graph +
`docs/inqkit_map.md`. No tests written yet; no source touched.

## Subtask 1 progress (2026-06-08)
- Installed + built the `understand-anything` plugin (corepack pnpm; fixed the
  pnpm-11 ignored-builds gate by building core via direct `tsc` and using the
  tree-sitter prebuilds — see Known issues).
- Ran the 7-phase pipeline over `inq-stack/include/inqkit/`:
  37 files → knowledge graph **88 nodes / 112 edges / 10 layers / 14 tour
  steps**, inline validation **0 issues**. Artifacts in
  `inq-stack/include/inqkit/.understand-anything/` (gitignored via `.*`).
  Explore with `/understand-dashboard`.
- Authored `docs/inqkit_map.md` (per-header role/API/tier/formula/test
  approach) + logged findings to `docs/notes/inqkit-rejuvenation-ideas.md`.
- **Key refinement:** I/O writers + `center_of_density` + `density_delta` are
  `pure` (POD inputs), not engine — verified by signature. Net surface:
  13 pure + 11 engine targets + 13 deferred stubs.

## What changed
- Resolved the full design tree: scope, test tiers, C++ harness, CI topology,
  verification-agent model, mapping approach, ideas-file location, error-log
  location, flow gate + definition of done.
- Discovered: 11 inqkit headers are **0-byte/TODO placeholders** never
  implemented on any branch (git history checked). `screens/*` is real code.
  CLAUDE.md is drifted (advertises empty `jellium/analytics` + `core/*`).
- `understand-anything` plugin installed this session (skills now live).

## Files touched (all new docs; no code)
- `/local/data/public/skcb2/tddft/CONTEXT.md` — glossary.
- `/local/data/public/skcb2/tddft/docs/adr/0001-inqkit-test-harness.md`
- `/local/data/public/skcb2/tddft/docs/adr/0002-ci-topology-local-first.md`
- `/local/data/public/skcb2/tddft/docs/plans/inqkit-rejuvenation.md` — locked plan.
- `/local/data/public/skcb2/tddft/docs/notes/inqkit-rejuvenation-ideas.md` — running ideas.

## Commands run
- Inventory: `find`/`wc`/`git log --follow`/`git rev-list --all` +
  `git cat-file -s` to prove the placeholder headers were never non-empty.
- No build, no test, no source edits.

## Tests and validation
- None run yet. Harness not built. `inq-stack/tests/cpp/` is empty;
  `tests/python/` has 2 pre-existing tests.

## Trusted sources used
- INQ's own Catch2 test idiom (`inq/src/observables/density.hpp`) as the
  framework precedent.

## Attribution notes
- C++ harness mirrors INQ's Catch2 usage but externalises tests (ADR 0001).

## Known issues / blockers
- Building an `engine`-tier C++ test (links INQ) is unproven — first engine
  test (fields tier) must validate the CMake INQ discovery before the rest.
- Python `engine`/VTK marker scheme is recommended but not yet confirmed by
  the user (confirm at Python phase).

## Assumptions still in play
- `config/tsubonoya_2014_coronene.hpp` is pure data (treated as `pure` tier) —
  verify when mapping.
- `screens/*` completed behaviour is testable; unfinished TODO path is not.

## Exact next steps
1. **Subtask 1**: run understand plugin over `inq-stack/include/inqkit/`;
   assemble `docs/inqkit_map.md` (role | public API | deps | pure/engine |
   formula-bearing? | testability) for the ~24 real headers; append ideas to
   `docs/notes/inqkit-rejuvenation-ideas.md`. Get user LOCK on the map.
2. **Subtask 2**: per header (bottom-up, pure tier first) interview the
   behaviour/expected/failure/tier; user locks each test plan.
3. Stand up `inq-stack/tests/cpp/CMakeLists.txt` with the `pure` label and the
   first pure test (`detail/grid_layout`) to prove the harness before engine.
