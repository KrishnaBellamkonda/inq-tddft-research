# inqkit test design & baseline-run registry

Rolling design doc for the inqkit unit/integration test effort (the TODO-driven
rejuvenation). Companion to:
- `docs/code-revitalisation/inqkit-review-and-next-steps.md` — user's review.
- `docs/code-revitalisation/inqkit-todo-catalogue.md` — T01–T31 catalogue.
- `docs/plans/inqkit-rejuvenation.md` — locked plan.
- ADR 0001 (test harness), ADR 0002 (CI topology).

---

## FULL COVERAGE TEST MATRIX (locked 2026-06-10)

Scope: **full public-API coverage**. Round = **tests only, no source changes**
(fft_shift move already done). Known bugs characterized via `[!shouldfail]` or
assert-current-behaviour. Source fixes (E01/E02/E03) + refactors (Vec3 swap,
stats `compute()`) deferred to a later round.

**Skipped:** 11 placeholder stubs (config/simulation_config, core/*, detail/
filesystem|text_io|validation, ground_state/*, io/manifest|text_summary,
jellium/analytics) — restructure phase. Trivial PODs (injection_report,
real_field_3d, complex_field_3d, step_context) — covered indirectly.

### Pure tier
| Component | Status | Test idea / key assertions |
|---|---|---|
| detail/grid_layout | ✅ done | flatten_index, fft_shift (Test A), step_suffix |
| detail/vec3 | ✅ done | dot/norm/arithmetic |
| observables/center_of_density | ✅ done | centroid, Bohr units (T11), weight+guard (T12), anti-transpose |
| observables/density_delta (compute) | ✅ done | fixed-t₀ base (T15), dV, mismatch throw |
| io/real_field_3d_writer | ✅ done | `.raw`+`.meta` round-trip (values, dims, origin, spacing) |
| io/complex_field_3d_writer | ✅ done | `_real.raw`+`_imag.raw` round-trip |
| io/vti_image_data_writer | ✅ done | **T06 pure half:** ascii VTI, per-axis-distinct values → x-fastest reorder verified (no transpose) + origin/spacing |
| io/observables_writer | ✅ done (ENGINE) | **reclassified engine** (StepContext pulls INQ): header+row reflect selection; hand-filled POD, no GPU |
| observables/density_delta (coarse_grain) | ⊘ SKIP | static-private helper writing BINARY VTI (base64) — impractical to parse; viz-only, low value. Revisit if made public. |
| jellium/shells | ✅ done | shell table + 162-electron closure + partial-shell truncation |
| screens/leed_pattern_accumulator | ☐ TODO (engine) | takes `electrons`; accumulate 2 steps → pattern = Σ slice·dt |
| config/tsubonoya_2014_coronene | ☐ maybe | assert geometry constants present (low value; confirm if wanted) |

### Engine tier
| Component | Status | Test idea / key assertions |
|---|---|---|
| fields/density::total | ✅ done | dims==basis, integral==N |
| fields/orbital::wavefunction | ✅ done | normalised complex field |
| fields/density (T02/E02) | ✅ done | cache stale→refresh, bath subtraction |
| wavepacket (T29/E03) | ✅ done | ortho report, single-pass limitation |
| coord_mapping (Test B) | ✅ done | off-centre He, both parities |
| observables/eigenvalue_dump | ✅ done | He GS → CSV: state count, ascending eigenvalues, occ (2,0,…), Ha→eV |
| observables/occupations_writer | ✅ done | MOCK Viewables → CSV occupations: state 0 = 2, WP slot = 1, others 0 |
| observables/state_energy_writer | ☐ TODO | **REAL short propagate** (Viewables needs `ham()`) → CSV per-state energy rows finite + count |
| observables/momentum_distribution | ☐ TODO | **MOCK Viewables** + injected k₀ WP → distribution CSV peaks at k₀ (lower priority — two-route already validates units) |
| observables/orbital_overlap | ✅ done | direct snapshot(electrons) at t=0 → identity block diag≈1/off≈0 + WP column valid |
| observables/wp_momentum_stats compute() | ✅ done | refactor extracted compute()→Moments; direct unit test ⟨p⟩=k₀, N≈1, ekin sanity |
| screens/leed_pattern_accumulator | ✅ done | pattern = Σ slice·dt over 2 steps |
| screens/plane_screen (extract single-rank) | ✅ done | He z=0 slice: shape, non-negative, non-trivial |
| screens/plane_screen (E01 multi-rank) | ✅ FIXED+done | mpirun -np 2 cross-rank agreement; E01 all_reduce applied, test green |

**Viewables driving (locked 2026-06-10):** RT-callback observables take a
templated `Viewables` with `.electrons()/.iter()/.time()` (+`.ham()` for
state_energy). HYBRID: tiny mock struct duck-types it for occupations +
momentum_distribution (no propagation, exercises the real templated code);
state_energy uses a real short `real_time::propagate`; orbital_overlap takes
electrons directly. Mock example:
`struct MockView { systems::electrons& e; int it; double t; auto& electrons() const {return e;} int iter() const {return it;} double time() const {return t;} };`
| observables/wp_momentum_stats (T28/T04) | ☐ TODO | **CROSS-VALIDATE TWO ⟨p⟩ ROUTES (user-locked):** on the same orthonormalised injected WP (k₀≠0), compute ⟨p⟩ via (1) real-space gradient `Re⟨ψ\|−i∇\|ψ⟩` and (2) reciprocal-space `∫k\|ψ̃\|²/∫\|ψ̃\|²` (wp_momentum_stats' method, via `to_fourier`). **Assert the two routes AGREE and both ≈ k₀** (units Bohr⁻¹ + complex-ψ + correctness). Both in-test, no source change. PLUS a class-level CSV smoke (drive the real class via short RT, assert columns/finite). |
| observables/wp_real_space_stats | ☐ TODO | injected WP → real-space ⟨r⟩≈WP centre, ⟨r²⟩ width ≈ σ (replica or short-RT CSV) |
| screens/plane_screen (extract, single-rank) | ☐ TODO | He+WP → extract z=0 slice → assert slice peak/sum consistent with density on that plane |
| screens/plane_screen (E01 multi-rank) | ☐ TODO | **cross-rank agreement, `[!shouldfail]`:** mpirun -np 2, assert all ranks' slice signatures equal — currently RED (documents E01); flip when fixed |

## Governing workflow rule — baseline before any change

**Phase 0 (baseline capture) runs before any source restructuring or change.**

Any reference/baseline simulation a characterization (golden-master) test needs
**must be executed and its outputs frozen first**, while the code is still in its
current state. Only once every required baseline is captured do we begin the
refactor/fix phase. The characterization tests then prove "results unchanged" by
comparing post-change output against the frozen Phase-0 baseline. A change to
source is permitted **only** if a test fails *and* the failure has no test-side
explanation (i.e. the test itself is confirmed correct).

Practical consequences:
- All "simple, fast runs needed for tests to work" are collected up front in the
  **Baseline-run registry** below and executed in one Phase-0 pass.
- Runs are small and GPU-cheap (single light atom, tiny box) — see registry.
- Each baseline's frozen artefacts get a fixed path + a recorded checksum/shape
  so drift is detectable.

---

## Test tiers (recap, from ADR 0001)

- **pure** — no INQ link; C++17 only; runs every CI.
- **engine** — links INQ; GPU; runs locally / occasionally.
- **char** — characterization/golden-master: freezes current output to detect
  refactor drift (may be pure or engine).

---

## Cluster Θ-coord — coordinate convention & fft_shift  (T01, T05, T06, #12)

### ✅ DONE (2026-06-10): T01 + T05 move + Test A green
- `fft_shift_index` moved to `detail/grid_layout.hpp` (pure). Removed from
  `density.hpp` (redirect note left); call sites in `density.hpp` (6) and
  `orbital.hpp` (3) re-pointed to `grid_layout::fft_shift_index`. `orbital.hpp`
  no longer includes `density.hpp` (T05 resolved; remaining `density::` mentions
  are comments only).
- **Test A** (`tests/cpp/test_fft_shift.cpp`) written RED-first (referenced the
  destination symbol → compile error), GREEN after the move. 4 TEST_CASEs:
  documented even table, odd parity, origin→0 all parities, bijection. Pure
  harness (`tests/cpp/CMakeLists.txt`, Catch2 reused from INQ build) + smoke
  (`test_grid_layout.cpp`) also green. `ctest -L pure` = 2/2 pass, 0.03s.
- **✅ ENGINE COMPILE VERIFIED (2026-06-10):** engine harness stood up at
  `tests/cpp/engine/` (GPU); `test_density_total_engine` (He, initial_guess,
  `density::total`) builds + links + runs, integral ≈ 2.0. So `density.hpp` with
  the relocated `grid_layout::fft_shift_index` is confirmed against INQ on GPU.
  `orbital.hpp` shares the same header path (compiles transitively; a dedicated
  orbital engine test will confirm directly).

### Architecture change (T01 + T05)
Move `fft_shift_index(int output_idx, int size)` **out of**
`inqkit::fields::density` (in `fields/density.hpp`) **into**
`inqkit::detail::grid_layout` (in `detail/grid_layout.hpp`), alongside
`flatten_index`. That header is already the pure index-convention file
(includes only `<string>`), so the helper stays pure.

Knock-on edits:
- `fields/density.hpp` — call `grid_layout::fft_shift_index(...)`; drop the local def.
- `fields/orbital.hpp` — call `grid_layout::fft_shift_index(...)`; **remove**
  `#include <inqkit/fields/density.hpp> // fft_shift_index` (resolves T05).

This change is gated by **Test A** (must stay bit-identical).

### Cross-reference finding (corroborates the code's physical claim)
INQ `basis::grid::to_symmetric_range` (`inq/src/basis/grid.hpp:78-84`) subtracts
`sizes` once `ix >= (size+1)/2`, so **array index 0 ↔ symmetric coord 0 ↔ cell
centre**, and `fft_shift_index = (output_idx + (size+1)/2) % size` is exactly the
inverse (`from_symmetric_range`). The documented table in `density.hpp:48`
matches. → Flagship #12 is a *regression guard + end-to-end validation*, not a
bug hunt. Two unknowns remain that only a real run settles: a possible
**half-cell offset** (peak at `c·dx` vs `(c+½)·dx`) and **even/odd parity**.

### Test A — `fft_shift_index` permutation  (pure · unit · gates the refactor)
Catch2 test over `detail/grid_layout.hpp` only. Cases:
1. **Documented even table** (size 6): the 6 mappings in `density.hpp:48`.
2. **Odd-size parity** (sizes 5 and 7): `(size+1)/2` rounding keeps
   physical-origin → FFT index 0.
3. **Permutation invariant**: for a sweep of sizes (both parities) the outputs
   are a bijection of `[0,size)` (no collisions). No inverse fn added.

### ✅ Test B DONE (2026-06-10): `tests/cpp/engine/test_coord_mapping_engine.cpp`
Off-centre He at (1.5,−2.0,1.0) Bohr, EVEN (L=10.0) + ODD (L=10.5) grids,
dx=0.5; `density::total` argmax recovered within one cell on every axis (three
distinct signed coords → catches transposition/sign-flip). **Both parities
green, 2.67s.** Coordinate mapping validated end-to-end (C++ side). VTI→Python
read-back (T06 Python half) deferred to the inqview phase.

### Test B — off-centre single-atom coordinate round-trip  (engine · integration · char)
Oracle: **a single light atom whose electron density peaks at the nucleus**,
placed at an **asymmetric off-centre position with three distinct signed
coordinates** (so an x↔z swap or a sign flip is caught — a symmetric position
would not be). Pipeline exercised end-to-end:

  GS density (`electrons.density()` → `fields::density::total`)
    → 3D RealField write (+ meta Origin = −L/2)
    → assert `argmax(field)` cell == ion cell (within half-cell)
    → VTI write → Python read-back → assert same coordinate.

This single test exercises Θ-coord across C++ → VTI → Python (T06 + #12) and
validates the half-cell/parity unknowns.

**T06 coverage:** the VTI writer already converts its z-fastest C buffer into
VTK's x-fastest PointData order on the fly (`vti_image_data_writer.hpp:38-51`,
iterating `iz` outermost → `ix` innermost). Test B's Python read-back asserts
that reorder introduces no axis transposition — so T06 needs no separate test.

**Assertions:**
- C++ side: `argmax(field)` cell index == the cell the He nucleus was placed on,
  **exact cell match** (nucleus sits on a cell centre by construction).
- Python side: VTI loaded by inqview resolves the He peak to the same physical
  `(x,y,z)` coordinate (within half a cell).

**Baseline-validation checkpoint:** if the Phase-0 He run shows density that does
*not* cleanly peak at the nucleus (e.g. pseudopotential smoothing splits/flattens
the peak), swap the oracle to a synthetic injected Gaussian *before* freezing the
baseline. This is exactly why the run happens in Phase 0.

### Concrete run spec — BL-coord-1a / BL-coord-1b
- Atom: **He**, default LDA, default pseudopotential.
- Cubic box, spacing **dx = 0.5 Bohr**.
  - **1a (even):** 20 pts/axis, L = 10.0 Bohr.
  - **1b (odd):** 21 pts/axis, L = 10.5 Bohr.
- Nucleus at grid offset **(+3, −4, +2)** cells from box centre (three distinct
  signed coords) → physical ≈ **(+1.5, −2.0, +1.0) Bohr**, on a cell centre.
- Frozen artefacts per run: the written `.raw` + `.meta.txt` (RealField3D) and
  the `.vti`, with recorded shape + checksum + the argmax cell index.

---

## Cluster Θ-density-semantics — total vs WP-excluded density  (T02, T03, T04)

### T02 — does `electrons.density()` include the WP?  (engine · integration)
**Established (code + project history):** `density::total()` → `electrons.density()`
sums `|ψ|²·occ` over **all** states in the kpin set; the WP is injected as an
extra occupied orbital (`wavepacket.hpp::inject_into_last_extra_state`), so
`total()` is the **total** density (WP included). `total_excluding_orbital`
recovers the bath. Review asks to *confirm empirically*, not assert.

**Test (direct WP-inclusion):** baseline **BL-dens-1** — He bath (N=2) + one WP
orbital injected (occ=1) via the production path. Assert:
- `∫ total()·dV ≈ 3.0` (bath 2 + WP 1) → proves the WP is included.
- `∫ total_excluding_orbital(wp_index)·dV ≈ 2.0` → proves the bath subtraction.

This one test answers T02 **and** validates T03's function on the same path the
jellium `_wf` runs use. Physical oddity (WP in a He box) is irrelevant — it is a
state-counting test.

**⚠ EMPIRICAL RESULT (2026-06-10) — overturned the resolved-by-reading guess.**
`test_density_semantics_engine` ran: pre-propagation, `∫total = 2` (WP NOT
included), `∫total_excluding_orbital = 1` (double-subtraction). `electrons.
density()` is a **cached** field (`spin_density_`), stale after manual injection.
→ **E02** in `inqkit-errors.md` (incl. suspected `_wf` t=0-frame production bug).
T02's real answer: WP-inclusion depends on refresh state (stale-bath pre-prop;
full post-propagate). Test expectation must be revised after the user
adjudicates E02; verification = propagate ≥1 step and recheck.

### T03 — is `total_excluding_orbital` used?  (RESOLVED by grep — keep + document)
**Used, load-bearing.** Computes the canonical bath density `n_total − n_wp` in
jellium `_wf` production: `shared/cpp/run_template.hpp:284` (bath0, t=0), `:421`
(bath_f, final), and every `run_wp_*_wf/run.cpp:181,237`. The two overloads are
recompute-from-electrons vs reuse-precomputed-fields (a perf split for per-step
callbacks). **Decision: keep; replace the TODO with a doc comment** pointing at
the canonical-bath definition. BL-dens-1 doubles as its characterization.

### T04 — is the full complex ψ used in observables?  (partly resolved by grep)
`orbital_overlap.hpp` **does** consume the full complex field (builds
`ComplexField3D` ref + evolved wfns via `orbital::wavefunction`). But
`wp_momentum_stats` does **not** route through `ComplexField3D` — risk that a
momentum distribution is built from `FFT(|ψ|²)` instead of `|FFT(ψ)|²`.
**Test (unit):** inject a known plane wave `ψ = e^{ik·r}` → assert the momentum
observable peaks at **k**, not 0. Folded with T28 (k-units) in Θ-parallel/
wavepacket grilling.

---

## Cluster Θ-parallel — GPU+MPI reduction correctness  (T21, T22, T25, T27)

### Reference-correct pattern (already in the repo)
`wp_real_space_stats.hpp:162-214`: local `gpu::run`/`gpu::reduce` (7 partial
sums: n, r·n ×3, r²·n ×3) → then **two** `all_reduce_in_place_n(buf, 7, +)`:
one over `basis.comm()` (FFT-grid decomposition) and one over `phi.set_comm()`
(state decomposition). This is the template every reducing observable should
follow. T21/T22 = write a test that **confirms** this is correct (expected GREEN
from the start = characterization).

### T25 — `plane_screen` missing Allreduce  (CONFIRMED BUG → fix via red→green)
`plane_screen.hpp:104-126`: loops `ist < phi.set_part().local_size()` (local
states only) with **no** `all_reduce` afterwards, and **no guard** (unlike
`density.hpp` which throws on multi-rank). → silent wrong slice under state
parallelism. **Fix:** add the same two `all_reduce_in_place_n` calls as
`wp_real_space_stats`. Logged in `docs/validation/inqkit-errors.md`.

### Parallel-invariance test harness (covers T21, T22, T25, T27)
Invariant: **result(1 rank) == result(N ranks)** elementwise (within tol),
across BOTH decomposition axes. Implementation:
- Phase 0 freezes the **np=1** plane_screen slice + wp_real_space_stats 7-vector
  (+ wavepacket post-injection norm/overlap reduces, T27) as golden artifacts.
- The **np=2** engine test recomputes and asserts elementwise equality vs frozen.
- `plane_screen`: RED pre-fix (drops states), GREEN post-fix.
- `wp_real_space_stats` (T21/T22) + `wavepacket` reduce (T27): GREEN from start
  (already do the all_reduces) — characterization guards against regression.

**⚠ Harness risk (highest):** this is the first **engine multi-rank** test.
Standing up the CMake INQ link + `mpirun -np 2` ctest path (`INQ_EXEC_ENV` per
CLAUDE.md) is unproven. Prove a trivial 2-rank engine test compiles+runs before
wiring the invariance assertions.

### System: reuse BL-dens-1 (He + WP = 3 states, splits 2/1 across 2 ranks).

---

## Cluster Wavepacket physics  (T04, T28, T29, T31, T26)

### T04 — momentum uses the full complex ψ?  (RESOLVED by reading — correct)
`wp_momentum_stats.hpp:140,162`: `operations::transform::to_fourier(kpin[0])`
then `|ψ̃(k)|²` = the correct `|FFT(ψ)|²` (not `FFT(|ψ|²)`), with the two
`all_reduce`s (211-214). No bug. Proven by the T28 twin test below.

### LOCKED (2026-06-10): stats testability refactor (governs T28/T04 + parallel-invariance)
`WPMomentumStats`/`WPRealSpaceStats` currently compute moments then only WRITE
CSV (no return), consuming an RT `Viewables`. **Refactor:** extract a public
`compute(...) → Moments` (px,py,pz, σ², ekin, N); `accumulate()` = `compute()` +
CSV write. Behaviour-preserving — **characterize the CSV output first** (freeze a
golden from a baseline run), refactor, re-run, assert CSV byte-identical, THEN
unit-test `compute()` directly. Makes both observables unit-testable and the
parallel-invariance test (#3) a struct comparison (no CSV parsing).
Sequence: Phase-0 CSV golden → extract compute → verify CSV unchanged → new tests.

### T28 + T04 — ⟨p⟩ two-route cross-validation  (engine · integration) — USER-LOCKED 2026-06-10
**Supersedes the earlier twin-WP design.** On the same orthonormalised injected
WP with known `k₀ ≠ 0`, compute the mean momentum two independent ways:
- **Route 1 (real space):** `⟨p⟩ = Re ⟨ψ| −i∇ |ψ⟩` via the real-space gradient
  on the WP orbital.
- **Route 2 (reciprocal):** `⟨p⟩ = ∫ k |ψ̃(k)|² dk / ∫ |ψ̃|² dk` via `to_fourier`
  — the exact method `wp_momentum_stats` uses.
**Assert:** the two routes AGREE (within tol) AND both ≈ `k₀`. Agreement of two
independent methods validates units (Bohr⁻¹), complex-ψ usage (T04), and the
correctness of the reciprocal-space momentum machinery in one shot. Both routes
computed in-test → no source change. Optional cross-check: short free
propagation, cod slope ≈ k₀ (group velocity). Baseline **BL-wp-1**.
(The before/after-ortho momentum-shift comparison the user also described is a
separate *physics* experiment, parked with T26.)

### T29 — orthogonalisation rigor  (engine · integration)
`orthogonalise_against_occupied` does modified Gram-Schmidt vs occupied states.
Reuse **BL-dens-1** (He bath + WP) — **but place the WP overlapping the He
density** so pre-orthogonalisation overlap is non-trivial (else the test passes
vacuously). Assert:
- post-ortho `max|⟨ψ_wp|ψ_occ⟩| < 1e-6` (from `InjectionReport.max_overlap`),
- `norm_after ∈ [0.97, 1.03]` (dev-feedback-loop WP norm band),
- record `norm_before`/pre-ortho overlap to prove the subtraction did work.

### T31 — GS protocol limitations  (data-driven from T29)
Code self-answers (wavepacket.hpp:262-264): KS orbitals mutually orthonormal →
single-pass GS doesn't reintroduce cross-overlaps in exact arithmetic; only
finite-precision residual remains. **Decision deferred to T29's measured
residual:** if residual ~1e-12, single pass suffices; if ~1e-6+, add an iterated
(second) GS pass. No pre-commitment.

### T26 — momentum-space Gram-Schmidt A/B  (oneoff · later)
Experiment, not a CI test: build a k-space GS variant, compare observables vs
the real-space GS, decide. Out of the test scope; tracked for the architecture
proposal. Links overall #10.

---

## Cluster Θ-vector + Θ-naming + doc singles  (mechanical, no new runs)

### T07/T09 — vector unit  (REFAC under characterization)
Add **`inqkit::detail::Vec3 {double x,y,z}`** (pure POD, header-only, `.dot()`,
`.norm()`, operator overloads mirroring `inq::vector3`'s API). Consolidate:
`CenterOfDensityResult{x,y,z}`, `current_{x,y,z}`, `dipole_{x,y,z}`. **Keeps
center_of_density + observables_writer PURE.** Characterization net: CSV/field
output byte-identical before/after the type swap.

**✅ Partial DONE (2026-06-10):**
- `detail/vec3.hpp` created + `tests/cpp/test_vec3.cpp` (4 cases: dot/norm,
  arithmetic, k̂-projection) green.
- **COD characterization** `tests/cpp/test_center_of_density.cpp` (5 cases)
  green — locks the CURRENT `CenterOfDensityResult{x_bohr,…}` behaviour incl. an
  anti-transposition case (centroid (1,2,3) all-distinct), cell-centre +0.5
  convention, Bohr units (**T11 validated**), `total_weight=∫f dV` + w>0 guard
  (**T12 validated**).
- **⚠ DEFERRED (engine-caller change):** the actual struct-swap
  `CenterOfDensityResult → {Vec3 center_bohr; double total_weight}` +
  `current/dipole → Vec3` touches production callers
  (`run_template.hpp:433`, `run_wp_e1000_L40x40x150/run.cpp:309` read
  `cod.x_bohr`). Do under the COD char test once engine compile is verifiable;
  callers move to `cod.center_bohr.x`.

### T10/T13/T19 — naming  (REFAC under characterization)
`f→field`, `w→weight`, `mx/my/mz→moment_{x,y,z}` (yes, `m`=moment). `function_`
trailing underscore = **C++ member-variable convention** (Python analogue
`self._x`); document the cross-language mapping. Apply renames inline **when each
file gets its characterization test**, so behaviour-unchanged is proven by the
same frozen output. Codify the convention in the plan/coding-standard note.

### Resolved-by-reading singles (doc only, no behaviour change)
- **T11** — `0.0L` is a C++ long-double literal (precision accumulator), not a
  unit suffix. Coordinates already Bohr. → clarifying comment.
- **T12** — `total_weight` = `∫f dV`, a normalisation diagnostic; the "condition"
  is the divide-by-zero guard when the field integrates to ~0. → comment.
- **T20** — RAII destructor flushes + closes the CSV stream when the writer goes
  out of scope. → comment explaining the idiom.

### Deferred to the restructure phase (REFAC/feature, tracked not tested now)
- **T08** richer meta schema (ions/electrons/GS-RT). **T23** generalise
  plane_screen to all axes. **T24** time-averaged density. **T26** k-space GS
  (oneoff). **T30** multi-kpoint injection (currently gamma-only throws).

### density_delta + write_every singles (RESOLVED by reading)
- **T15** — δn uses a **single fixed reference t₀ for all timesteps**
  (`delta = current − ref_`, ref set once / lazy-captured), NOT a rolling
  `t→t+1` diff. ✅ **DONE (2026-06-10):** `tests/cpp/test_density_delta.cpp`,
  6 cases incl. the discriminator (`L2(C)=36` fixed-t₀, ≠ `16` rolling) + lazy
  capture, explicit set_reference, dV scaling, grid-mismatch throw. Pure (emit
  flags off → no file I/O). Green. (TODO comment at density_delta.hpp:135 to be
  replaced with the answer in the later doc-comment batch with T11/T12/T20.)
- **T18** — `write_every` is **NOT inherited** from the RT session. Two
  independent gates: `RealTimeSession.write_every_` (task dispatch) and
  `wp_real_space_stats cfg_.write_every` (accumulation); run_template.hpp:324
  sets them differently (10×WRITE_EVERY vs WRITE_EVERY) → they **compose**
  (footgun). → doc clarification + **unit test** on the `iter % write_every` gate.
- **T14/T16** — choice of reference (t=0 vs t=dt) + step-by-step δn = oneoff
  experiments, out of CI scope.

---

## Baseline-run registry (Phase 0 — execute before any change)

| ID | Purpose | System | Status |
|---|---|---|---|
| BL-coord-1a | Test B oracle (even grid): off-centre He GS density | He, 20³ pts, L=10.0, dx=0.5, ion at (+1.5,−2.0,+1.0) Bohr | defined, not run |
| BL-coord-1b | Test B oracle (odd grid): off-centre He GS density | He, 21³ pts, L=10.5, dx=0.5, ion at (+1.5,−2.0,+1.0) Bohr | defined, not run |
| BL-dens-1 | T02/T03 (counting) + T29 (orthogonalisation) | He bath (N=2), 20³ pts L=10.0 dx=0.5, 1 WP injected occ=1 **placed overlapping the He density**, orthogonalisation ON; record ∫total, ∫total_excl_orbital, WP orbital density, InjectionReport (norm_before/after, max_overlap) | defined, not run |
| BL-parallel-1 | T21/T22/T25/T27 oracle: np=1 reduces | SAME config as BL-dens-1, launched np=1; freeze plane_screen slice + wp_real_space_stats 7-vector + wavepacket norm/overlap reduces as golden artifacts | defined, not run |
| BL-wp-1 | T28/T04: k-units + complex-ψ | vacuum + WP, TWIN runs k₀=0 and k₀≠0 (same b,σ), short free propagation; record ⟨p⟩(t=0), cod(t) slope | defined, not run |

(More rows added as later clusters are grilled.)
