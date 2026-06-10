# inqkit error register

Structured log of confirmed/suspected defects found during the inqkit
rejuvenation. Per workflow: a defect is **documented here first**; source is
changed only via a test that fails on the defect and passes after the fix
(red→green), with the change recorded in `docs/validation/inqkit-tests.md`.

Status: `suspected` · `confirmed` · `fix-pending` · `fixed`

---

## E02 — electrons.density() is a cached field; stale after WP injection
- **TODO:** T02 / T03 (`fields/density.hpp`)
- **Status:** confirmed (caching) · **suspected production impact (t=0 frame)**
- **Severity:** high (suspected wrong t=0 bath/total frames in jellium `_wf` runs)
- **Confirmed by test** (`test_density_semantics_engine`, 2026-06-10): after a
  GS + `WavePacket::inject_into_last_extra_state(occ=1)`, **with no propagation**,
  `density::total` integrates to **2 (He bath only), not 3** — the WP is NOT
  included. `total_excluding_orbital` then **double-subtracts**: `bath = 2 − 1 =
  1` instead of 2.
- **Root cause:** INQ `electrons.density()` (`inq/src/systems/electrons.hpp:444`)
  returns `observables::density::total(spin_density_)` — a **cached member
  field** computed during SCF/propagation. Manual injection writes the orbital +
  sets `occupations()[0][ist_wp]` but does **not** refresh `spin_density_`, so
  `density()` returns the stale post-SCF (bath) density.
- **Corrected T02 answer:** "does `density()` include the WP?" depends on the
  refresh state. **Right after injection (pre-propagation): NO** (stale bath).
  During/after `real_time::propagate` (which rebuilds `spin_density_` from all
  occupied states incl. the WP at occ=1): **YES** (full). My earlier
  resolved-by-reading guess (always-included) was WRONG — the test overturned it.
- **Suspected production impact:** the `_wf` runs' **t=0 block** computes
  `full0 = density::total(electrons)` (line ~179) BEFORE `real_time::propagate`
  (line ~264). If `spin_density_` is stale there (as the test shows), then
  `full0` is the bath, `bath0 = full0 − wp0` is over-subtracted, and the t=0
  `total`/`system` frames are wrong. **t>0 frames refreshed by propagate are
  fine.** NEEDS confirmation: propagate ≥1 step in the test and recheck (does
  density::total become 3?); and check whether the RT callback also emits a
  (correct) iter-0 frame that supersedes the standalone t=0 block.
- **Do NOT fix yet** — this is a production-correctness question for the user to
  adjudicate (the `_wf` bath definition is theirs). Document first.
- **PROPER FIX (user-proposed 2026-06-10):** compute density directly from the
  orbitals — `ρ(r) = Σ_i occ_i |ψ_i(r)|²` over `electrons.kpin()[0]` +
  `electrons.occupations()[0]` — instead of reading the cached
  `electrons.density()`. This always reflects the current orbitals (incl. an
  injected WP), so it is WP-aware WITHOUT needing a propagation step to refresh
  the cache. `plane_screen::extract` already uses this approach for exactly this
  reason. Trade-off: O(n_states) sum per call vs one cached field read.
  Recommended as a new `fields::density::total_from_orbitals(electrons)` (with
  its own integral==N test), leaving the cached `total()` available; production
  WP-density paths switch to the orbital-based one.

## E03 — single-pass orthogonalisation misses tolerance for strong overlap
- **TODO:** T29 / T31 (`wavepacket/wavepacket.hpp`)
- **Status:** FIXED 2026-06-10 (iterated 2-pass GS + residual measurement)
- **Fix:** wrapped the modified-Gram-Schmidt projection loop in 2 passes; the
  pre-ortho overlap (first pass) is reported as `max_overlap`, the FINAL-pass
  residual gates `passed_tolerance`. Root cause was partly measurement: the old
  `max_ov` captured the pre-subtraction overlap, so `passed_tolerance` could
  never be true for an overlapping WP. Test `test_density_semantics_engine`
  T29 case flipped CHECK_FALSE→CHECK(passed_tolerance).
- **Severity:** medium (WP not rigorously orthogonal to occupied states)
- **Confirmed by test:** a WP centred on the He (sigma 1) has pre-ortho
  `max_overlap = 0.966` (nearly the He 1s). After
  `orthogonalise_against_occupied`, `norm_after ≈ 1.0` (✓) but
  **`passed_tolerance == false`** — the post-ortho overlap exceeds the 1e-6
  tolerance. The single-pass modified Gram-Schmidt does not reach tolerance when
  the WP is strongly inside the occupied subspace.
- **This is exactly T31's question** (wavepacket.hpp:262): single pass leaves a
  residual; an **iterated (second) GS pass** is the standard fix. The test gives
  the data the T31 decision needed: at strong overlap, one pass is insufficient.
- **Caveat:** an overlap of 0.97 is a pathological stress case (WP ≈ He orbital).
  Realistic WPs (a fast Gaussian far in k-space from occupied states) overlap
  weakly; one pass may suffice there. Recommend: test BOTH a weak-overlap case
  (should pass single-pass) and this strong case (needs iterated) to bound it.

## E04 — center_of_density: half-cell (+dx/2) coordinate offset
- **TODO:** Θ-coord cluster / overall #12 (`observables/center_of_density.hpp`)
- **Status:** confirmed
- **Severity:** medium-high (WP centroid/trajectory biased by half a grid cell)
- **Confirmed by test** (`test_wp_real_space_replica_engine`, 2026-06-10): a WP
  injected at centre (1.0,−1.0,0.5) Bohr (via INQ `rvector_cartesian`) is
  recovered at (1.25,−0.75,0.75) when the field index is mapped with
  `origin+(ix+0.5)·dx` — **uniformly +0.25 = +dx/2** (dx=0.5) on every axis.
- **Root cause:** INQ grid points are NODE-centred —
  `real_space.hpp:127-129` `rvector(ix)=to_symmetric_range(ix)·rspacing`, i.e.
  physical = symmetric_coord·dx (point 0 at exactly 0). The correct inqkit field
  coordinate is `origin + ix·dx`. `center_of_density.hpp:64` uses
  `origin + (ix+0.5)·dx` (cell-centre) → a **+dx/2 systematic offset** in the
  computed centroid.
- **Inconsistency:** `wp_real_space_stats` computes ⟨r⟩ via
  `point_op.rvector_cartesian` (INQ node convention) → CORRECT.
  `center_of_density` uses the +0.5 form → OFFSET. Two centroid codes, two
  conventions.
- **Why earlier tests missed it:** the `center_of_density` PURE test is
  self-consistent (builds its toy field with the same +0.5 convention, never
  compares to INQ) so it passes. `Test B` used a one-cell (dx) argmax tolerance,
  which absorbs a dx/2 error. Only a test comparing against the INQ-placed WP
  (this one) exposes it.
- **Fix (deferred, per 'fix none'):** drop the `+0.5` in
  `center_of_density.hpp` (use `origin + ix·dx`). Until then production WP
  trajectories from COD carry a +dx/2 bias (≈0.25 Bohr at dx=0.5).

## E01 — plane_screen: missing MPI Allreduce under state parallelism
- **TODO:** T25 (`screens/plane_screen.hpp`)
- **Status:** FIXED 2026-06-10 (two all_reduce calls added to extract())
- **Fix:** after the slice loop, flatten and `all_reduce_in_place_n(+)` over
  `phi.basis().comm()` (domain) and `phi.set_comm()` (state), then unflatten —
  every rank returns the complete slice. Verified by
  `test_plane_screen_parallel_engine` (mpirun -np 2, cross-rank agreement;
  was `[!shouldfail]`, now a normal green).
- **Severity:** high (silent wrong results, no guard)
- **Symptom:** `accumulate`/slice loop runs `for ist < phi.set_part().local_size()`
  (local states only) and performs **no `all_reduce`** afterwards. Unlike
  `fields/density.hpp` (which throws on multi-rank), `plane_screen` has **no
  guard**, so a state-parallel multi-rank run silently returns a slice missing
  every non-local state's contribution.
- **Evidence:** the correct pattern exists 100 lines away in
  `observables/wp_real_space_stats.hpp:209-214` — two `all_reduce_in_place_n`
  over `basis.comm()` and `phi.set_comm()` after the local `gpu::reduce`.
- **Fix (planned):** add the same two `all_reduce_in_place_n` calls to the
  plane_screen slice before returning.
- **Test:** parallel-invariance harness — np=1 frozen slice vs np=2 recompute,
  elementwise. RED pre-fix, GREEN post-fix. Baseline BL-parallel-1.
- **Manifests in production? — LATENT, not yet triggered (verified 2026-06-10).**
  Only one production run is genuinely multi-rank:
  `ResearchProject/systems/jellium/run_wp_n162_L50_E700_mpi_propagate` (`mpirun
  -np 2`; INQ may split states(101) or basis(168)). That run uses
  `wp_real_space_stats` + `wp_momentum_stats` (run.cpp:99-100), **not**
  `plane_screen`, and `wp_real_space_stats` does both `all_reduce`s → its results
  are safe. `plane_screen` is wired into `shared/cpp/run_template.hpp` (jellium +
  coronene) but has only ever been exercised single-rank.
  → **No corrupted results to date. Urgency = "fix before the next multi-rank
  run through run_template.hpp", not retroactive.**
