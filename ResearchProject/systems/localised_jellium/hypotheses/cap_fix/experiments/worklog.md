# cap_fix worklog — CAP energy-artifact removal

Session start: 2026-07-13. Goal: CAP setup with monotone reported E_total (rise ≤ 0.1 eV)
that still absorbs (absorbed_e ≥ 0.5). Baseline artifact: two-sided η=−0.2, 700 steps →
rise +23.49 eV. See `../autoresearch.md` for the full contract.

Data summary ("become one with the workload"): the artifact is the post-minimum rise of
E_total, driven by the norm-divided kinetic term as the CAP eats bound/slow density
(diagnosis campaign, conf 0.90). η=−1 is clean at 700 steps but rises by t=36 → winners
must pass a 950-step η=−1 confirmation segment. Runs are deterministic; noise floor
0.1 eV = 3× the worst clean-run rise (capon_reach +0.033 eV).

## Setup log

### Harness build + smoke (pre-run-1)
- Timestamp: 2026-07-13 (session)
- What: new `scripts/cap_fix/run.cpp` (ablation clone + EM_CAP_MODE two|wrap +
  charge.csv ∫n dV) + `inqkit/perturbations/absorbing_wrap.hpp` (cos² bump peaking at
  the periodic boundary; width 30 ≡ equal-integral twin of the two-sided default).
  Locked harness: autoresearch.sh / run_metrics.py / checks.sh.
- Sanity gate (from existing diagnosis data, same physics/binary lineage): known-bad
  scores 23.49 eV, known-clean scores ≤ 0.033 eV → the metric separates them. ✓
- Smoke: 10-step wrap-mode run on GPU 1 (build + new-code path + charge.csv sanity).

## Key Insights

- (from diagnosis) rise requires W·(bound/slow density) overlap — predicts footprint
  and strength arms beat topology arm; wrap-around fixes the W=0-at-boundary gap but
  not the overlap. To be TESTED, not assumed.

## Next Ideas

- If all setup arms fail → absorbed-energy accumulator ledger term (separate campaign).
- If wrap wins unexpectedly → boundary-crossing density (W=0 plane) was the real feeder;
  check phase-0 VTIs for boundary wrap-around flux.

### Setup fix: build-race + exit-code masking (pre-run-1, harness not yet locked)
- Timestamp: 2026-07-13 14:40
- What: launching runs 1+2 concurrently via `inq-run` raced in the shared build/
  dir → run 1's device link corrupted (nvlink undefined refs) while run 2's link
  won and ran. ALSO: piping the launcher into `tail` masked the non-zero exit
  (the exact trap the skill README warns about).
- Fix (deliberate, documented, pre-baseline): autoresearch.sh now REQUIRES a
  pre-built `./run` (freshness-checked vs run.cpp, exit 3 if stale) and executes
  it directly with the inq env — no per-run cmake at all. Launch logs go to
  `results/<name>.log`, never through a pipe.
- Harness LOCKED as of this fix.

### Wrap-profile test PASS (pre-run-1)
- `hypotheses/cap_fix/tests/test_wrap_profile.py`: wrap peaks (=|η|) at the
  boundary and is smooth across it; two-sided W=0 exactly there (the gap);
  equal integrals 15.0000η; equal |z|>25 footprint. Catalogue rows added.
- Smoke (10-step wrap, GPU 1): N_total=53.000 in-propagator (52 slab + 1 WP →
  absorbed_e metric measures WP absorption); stale step-0 row (52.0) kept in
  charge.csv as a WP-refresh diagnostic; run_metrics dedups keep-last.

### Runs 1–2 IN FLIGHT
- run01_baseline_two_eta0p2 (GPU 0, relaunched 14:40): two-sided η=−0.2 —
  expect rise ≈ +23.5 eV (reproduction gate for the new binary).
- run02_wrap_eta0p2_w30 (GPU 1, 14:32): wrap topology, equal-integral twin —
  the user's hypothesis-1 draft. ETA ~55 min each.
