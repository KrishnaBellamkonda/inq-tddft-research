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

### Run 1: baseline two-sided η=−0.2 — artifact_rise_eV=23.4947 (KEEP, baseline)
- Timestamp: 2026-07-13 15:26
- What changed: nothing (baseline on the new binary + charge.csv).
- Result: rise 23.4947 eV — BIT-IDENTICAL to the diagnosis capon_weak_partial;
  excursion 0.11, drain −23.38, absorbed 0.871e (87% of the WP), t_min 21.6.
- Insight: metric pipeline fully validated end-to-end; absorbed_e finally
  quantifies absorption (the diagnosis-era density_l2 could not).
- Next: strong-η and footprint arms.

### Run 2: wrap-around equal-integral twin — artifact_rise_eV=35.3384 (DISCARD)
- Timestamp: 2026-07-13 15:30
- What changed: EM_CAP_MODE=wrap, width 30 (same ∫W dz, same |z|>25 footprint).
- Result: rise 35.34 eV (+50% vs baseline), excursion 1.90 eV (17× worse),
  absorbed 0.876e (same), drain −33.4.
- Insight: **topology arm REFUTED.** Peaking W at the boundary — where the
  slowest, most diffuse density pools — eats MORE below-average-KE density:
  worse artifact at identical absorbing strength. Direct, controlled support
  for the covariance/filter mechanism. The two-sided W=0-at-boundary hole is
  NOT the artifact's feeder.
- Next: runs 3 (η=−2.0) + 4 (centre 35, width 10) in flight.

## Meta-review (after runs 1–2)
Topology class: loses. Mechanism-aligned classes (strength, footprint) now carry
the campaign; if both fail → setup-only fix falsified → absorbed-energy ledger
term is the answer (separate campaign).

### Run 3: strong-CAP two-sided η=−2.0 — artifact_rise_eV=0.0000 (KEEP — best)
- Timestamp: 2026-07-13 16:25
- What changed: EM_CAP_ETA −0.2 → −2.0 (one lever; geometry unchanged).
- Result: rise 0.0000 eV (monotone at every sample), excursion 0.0000, drain
  −153.2 eV, absorbed_e 0.9998 — the ENTIRE wavepacket removed by t=28.
- Insight: strong absorption finishes the job before the slow-nibble era —
  consistent with the mechanism (the artifact needs a long weak-absorption
  tail-eating phase). CAUTION: t_min = end-of-window; η=−1 looked identical
  at t=28 and rose +31 eV by t=36 (phase-0). Monotone-at-700 is necessary,
  not sufficient.
- Next: run05_confirm_eta2p0_950 (GPU 0, IN FLIGHT) — same config, 950 steps
  (t=38): the decisive late-rise test. Planned positive control: η=−1 @950
  (must REPRODUCE the +31 eV rise, proving the window catches late risers).

### Run 4: pushed-out footprint η=−0.2 c35 w10 — artifact_rise_eV=24.3767 (DISCARD)
- Timestamp: 2026-07-13 16:32
- What changed: CAP centre 32.5→35, width 15→10 (footprint 30–40 Bohr), η fixed −0.2.
- Result: rise 24.38 eV (~baseline), excursion 2.63 eV (24× baseline's 0.11 —
  WORSE), absorbed 0.755 (down from 0.871), drain −21.7.
- Insight: geometry does NOT beat the artifact at weak η — the slow spill
  reaches the CAP wherever it sits; weak absorption then nibbles it for the
  whole run. Less absorbed WP → shallower drain → the (unchanged) rise crosses
  zero HIGHER. The lever is absorption completeness/speed (η), not footprint.
- Next: run06_poscontrol_eta1p0_950 (GPU 1, IN FLIGHT) — positive control:
  η=−1 @950 must reproduce phase-0's late +31 eV rise, validating that the
  950-step window catches late risers (else "clean @950" means nothing).

## Meta-review (after runs 1–4, all four draft arms in)
Classes: topology LOSES (worse). geometry-at-weak-η LOSES (same rise, worse
excursion). strength WINS at 700 steps (η=−2: rise 0.0000, full absorption).
Open: does strength survive t=38? (run 5 in flight; run 6 = window validity.)

### Run 5: η=−2.0 @950 (confirmation) — artifact_rise_eV=3.5001 (KEEP, segment-1 baseline)
- Timestamp: 2026-07-13 18:45
- What changed: winner config unchanged; EM_N_STEPS 700→950 (t=38) [segment 1].
- Result: excursion 0.000 — never climbs above the t=0 reference (the
  user-visible pathology is ABSENT); residual late rise 3.50 eV from min at
  t=33.8 (9× smaller than η=−1's +31 eV); drain −167 eV; absorbed 1.0008e
  (whole WP + slab spill — the spill IS the residual-rise feeder).
- Insight: strong η doesn't eliminate the mechanism, it out-runs it: absorption
  completes before slow-density nibbling can accumulate — but the excited slab
  keeps evaporating spill into the CAP, giving a small late rise. Setup tuning
  may floor out above the 0.1 eV target; the excursion criterion is already met.
- Next: η=−4 ladder (run 7) + strong-η×pushed-footprint cross (run 8).

### Run 6: η=−1.0 @950 (positive control) — rise=169.32, excursion=+31.27 eV (DISCARD by design)
- Timestamp: 2026-07-13 18:45
- What changed: η=−1, 950 steps [control, op:debug].
- Result: reproduces phase-0's +31.3 eV above-zero excursion at t=36 —
  **confirmation window VALIDATED** (it catches late risers). Run TIMED OUT
  (exit 124) at step 904/950: per-step wall time degraded 4.7→21 s late in the
  run (cause unknown — GPU 1 shows free now; possibly propagator iteration
  growth at strong late-time absorption). Partial CSV (t=36) fully sufficient.
- Insight: partial-run metrics work (the locked extractor reads whatever rows
  exist); watch for recurring slowdowns → budget margin at 950 steps is thin.

## Meta-review (after 6 runs)
Validated instrument + validated mechanism. Ladder now optimizes a RESIDUAL
(3.5 eV rise, already below-zero everywhere). If η=−4 and the cross both floor
out above 0.1 eV, the honest verdict is: setup tuning fixes the EXCURSION
(user's visible artifact) but only bookkeeping (absorbed-energy term) yields a
strictly monotone reported ledger.

### Run 8: η=−2.0 × pushed footprint (c35 w10) @950 — artifact_rise_eV=24.6554 (DISCARD)
- Timestamp: 2026-07-13 20:05
- What changed: vs run 5, CAP centre 32.5→35, width 15→10 (one compound
  geometry lever, footprint 30–40 Bohr).
- Result: excursion 0.000 (still never above zero) BUT residual late rise
  24.66 eV — 7× worse than run 5 (3.50); absorption complete (1.000e).
- Insight: at strong η the footprint change HURTS: a narrower, steeper, more
  distant CAP leaves 5 Bohr of free slosh space and likely reflects more slow
  spill (steeper W = more reflective for low k), feeding a longer nibble era.
  Geometry loses again — strength + standard footprint remains the front-runner.
- Next: run09 wrap×η=−2 @950 (GPU 1) — last unexplored cross; run 7 (η=−4)
  still in flight on GPU 0 → 3-point η ladder verdict next.

### Run 7: η=−4.0 @950 — artifact_rise_eV=20.2110 (DISCARD)
- Timestamp: 2026-07-13 20:20
- What changed: vs run 5, η −2.0 → −4.0 (one lever).
- Result: excursion 0.000, but residual rise 20.21 eV — 6× WORSE than η=−2
  (3.50). Absorption complete (1.0009e), drain −164.7, t_min 32.2.
- Insight: the η ladder is NON-MONOTONE: η=−1 → 169.3, η=−2 → 3.50,
  η=−4 → 20.2. Classic CAP reflection trade-off (Riss & Meyer): too weak
  transmits (long nibble era), too strong REFLECTS slow spill back into the
  box where it lingers and is slowly re-eaten. An optimal absorption window
  exists near η ≈ −2 for this geometry/velocity.
- Next: run10 η=−3.0 @950 (GPU 0, launched) to bracket the right flank;
  η=−1.5 queued for the left flank after run 9 frees GPU 1.

## Meta-review (after 8 runs)
The campaign has found: (i) ALL strong-η configs (η ≤ −2) kill the above-zero
excursion — the user-visible artifact is fixed by strength alone; (ii) the
residual below-zero rebound has an η-optimum near −2 set by the
transmit-vs-reflect trade-off; geometry/topology changes only hurt. Remaining
question: how deep is the optimum (can any setup reach ≤0.1 eV strict
monotonicity, or does the floor sit at a few eV → bookkeeping needed for
strictness).

### Run 9: wrap × η=−2.0 (w30) @950 — artifact_rise_eV=0.3240 (KEEP — new best)
- Timestamp: 2026-07-13 21:35
- What changed: vs run 5, EM_CAP_MODE two→wrap (equal integral/footprint).
- Result: rise 0.324 eV (10× better than run 5's 3.50), excursion 0.000,
  absorbed 1.0005e, drain −168.3, t_min 36.4 (nearly end-of-window).
- Insight: **interaction effect between the user's two arms.** At weak η the
  wrap was worse (more slow-density eating at the boundary pool); at strong η
  the fast WP is fully absorbed either way, and the residual is dominated by
  slow spill near the BOUNDARY — exactly where the two-sided profile has its
  W=0 hole (spill leaks across the periodic plane unabsorbed and lingers).
  The wrap peak covers that plane; the rebound drops 10×. Neither arm alone
  sufficed; the cross wins. 0.324 eV is 3× above the strict 0.1 eV target.
- Next: run11 wrap η=−2 WIDTH 40 (gentler ramp, footprint |z|>20 — less
  reflection of slow spill) on GPU 1; run10 (two-sided η=−3) completing the
  two-sided ladder on GPU 0. Then possibly wrap η ladder (−2.5, −3).

### Run 10: two-sided η=−3.0 @950 — artifact_rise_eV=11.8320 (DISCARD)
- Timestamp: 2026-07-13 21:50
- What changed: η −2→−3 (right-flank bracket of the two-sided ladder).
- Result: 11.83 eV — sits between −2 (3.50) and −4 (20.2): the two-sided η
  ladder is convex with optimum ≈ −2, whose 3.50 eV floor is 10× the wrap's.
- Insight: the two-sided family's floor is set by its boundary hole — no η
  can fix a topology gap. Arm CLOSED.
- Next: wrap η ladder — run12 wrap η=−1.5 w30 (GPU 0, launched; hypothesis:
  wrap needs LESS strength since boundary coverage removes the transmit-leak
  penalty → optimum may sit at weaker η with less reflection). run11 (wrap
  w40) still on GPU 1.

### Run 11: wrap η=−2.0 w40 @950 — artifact_rise_eV=0.000000 (KEEP — TARGET REACHED)
- Timestamp: 2026-07-13 23:00
- What changed: vs run 9, wrap width 30→40 Bohr (gentler ramp; footprint |z|>20).
- Result: rise 0.000000 eV — monotone at EVERY sample, t_min = 38.0 = the final
  step; excursion 0.000; drain −178.3 eV; absorbed 1.015e.
- Insight: the gentler, earlier-starting cos² ramp absorbs the slow spill
  without reflecting it — the residual rebound vanishes below sampling
  precision. This config satisfies the campaign success criterion (≤0.1 eV,
  window validated by the run-6 positive control, 950⊃700 supersedes the
  screening segment). CAVEAT: absorbed_e = 1.015 — the |z|>20 footprint eats
  ~0.015e of static slab tail over t=38 (slow real-charge leak, 0.03% of 53e);
  flag in the production recommendation; run 12 (η=−1.5) may trade this off.
- Next: run13 = winner @1200 steps (t=48) hardening (GPU 1, in flight);
  run12 wrap η=−1.5 w30 (GPU 0) completes the wrap ladder → final
  recommendation between w40/η−2 and any equally-clean gentler config.

### Run 12: wrap η=−1.5 w30 @950 — artifact_rise_eV=8.5946 (DISCARD)
- Timestamp: 2026-07-14 00:20
- What changed: vs run 9, η −2.0 → −1.5 (wrap ladder left flank).
- Result: 8.59 eV — 26× worse than wrap η=−2 w30. The wrap family ALSO needs
  ≥ −2 strength; boundary coverage does not relax the transmit penalty.
- Insight: strength and topology are independently necessary: η≥2 to finish
  the fast era, wrap+wide ramp to absorb the slow spill without a hole or
  reflection. Winner (wrap η=−2 w40, rise 0.000000) unchallenged.
- Next: run14 = winner at HALF WP speed (K0 2.8465, one lever) @950 (GPU 0) —
  transferability: the user wants ALL runs clean, so probe a slower projectile
  (longer weak-absorption era = harsher test). run13 (winner @1200) on GPU 1.

### Run 13: winner @1200 steps (t=48) — artifact_rise_eV=0.000000 (KEEP — hardened)
- Timestamp: 2026-07-14 01:30
- What changed: vs run 11, EM_N_STEPS 950→1200 (t=48; hardening only).
- Result: monotone at EVERY sample to the final step; excursion 0.000;
  drain −192.0 eV; absorbed 1.020e.
- Insight: "monotone so far" now extends 12 a.u. beyond the era where η=−1
  rose +31 eV — the winner is not merely delaying the rise. The slab-tail
  nibble persists at ~5e-4 e/a.u. (0.005e over the extra window): a slow,
  steady, MONOTONE drain — it feeds no rebound, but long production runs
  should budget for it (or pull the wrap width back toward 35 if it matters).
- Next: run14 (winner at half WP speed) decides transferability → then close:
  study notebook + recommendation.

### Run 14: winner @ half WP speed — artifact_rise_eV=0.000000 (KEEP — transfers)
- Timestamp: 2026-07-14 02:40
- What changed: vs run 11, EM_K0 5.693→2.8465 (half projectile speed).
- Result: monotone to the end; excursion 0.0156 eV (early bound-tail nibble,
  6× below the noise floor); absorbed 0.796e (slower transit).
- Insight: the fix is not tuned to one velocity — the slow-projectile case
  (longer weak-absorption era, the harsher test) stays clean.

## FINAL SUMMARY (loop closed, 14 runs, target met)
Winner: **unified wrap-around CAP (absorbing_wrap), η=−2.0, width 40 Bohr** —
rise 0.000000 at t=38 AND t=48 AND at half speed; excursion ≤0.016 eV
everywhere; whole WP absorbed. Both user arms necessary, neither sufficient
(interaction effect). Falsified: topology@weak-η, geometry (both strengths),
η beyond the convex optimum, gentler wrap strength. Caveats: ~5e-4 e/a.u.
static-tail drain; reported E_total under ANY CAP remains
bookkeeping-incomplete (absorbed-energy accumulator = future campaign).
Study notebook: cap_fix_study.ipynb (executed, 0 errors).
