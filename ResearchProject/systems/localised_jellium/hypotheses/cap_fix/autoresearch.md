# Autoresearch: remove the CAP energy-ledger artifact (drain-then-rise) from localised-jellium runs

> The metric is the sole arbiter. "This looks better" never drives a keep/discard —
> only a measured improvement that clears the noise floor does.

## Objective

Find a CAP setup (topology + parameters) in which the reported `E_total(t)` of the
localised-jellium witness run shows **clean monotone decay** — no post-minimum rise, no
excursion above the t=0 reference — while the CAP **still absorbs** the outgoing
wavepacket. The diagnosed mechanism (energy_oscillation_diagnosis campaign, conf 0.90):
the reported ledger books no absorbed energy, so absorbing *bound / below-average-KE*
density RAISES the reported total (dominant channel: the norm-divided kinetic term,
`inq/src/hamiltonian/energy.hpp:55`). Candidate fixes under test (user-specified):
1. **Unified wrap-around CAP** — one smooth cos² bump peaking AT the periodic boundary
   (`inqkit::perturbations::absorbing_wrap`), vs the standard two sin² bumps which fall
   to W=0 exactly at the boundary. Equal-integral twin: wrap width 30 Bohr ≡ two-sided
   η, footprint |z|>25, ∫W dz = 15η.
2. **Parameter tuning** — η (strength), CAP centre/width (footprint), until the rise
   disappears. Mechanism prediction: pushing the CAP inner edge OUT (less overlap with
   slab tail / slow spill) and/or strong η (fast full absorption before the slow-nibble
   era) kill the rise; weak η in overlap is the worst case (= baseline).

Workload knowledge ("become one with the data", from the diagnosis campaign):
- Witness: slab_n52 GS + σ1 WP (m_eff 2.10, K0 5.693), ETRS, dt 0.04, **700 steps
  (t=28)**, η=−0.2 two-sided → drain −23.4 eV (min at t=21.6) then rise +23.5 eV,
  crossing to +0.11 eV. This IS the artifact.
- η=−1 two-sided, same window: monotone −138 eV (rise appears only after t=28; phase-0
  ran to t=36 and rose +31 eV → confirmation runs need 950 steps).
- CAP-off: conserved to −0.015 eV. Clean-run rise floor ≤ 0.033 eV.
- One 700-step run ≈ 45 min on one GPU (A30-class, both free at session start).

## Metrics

- **Primary**: `artifact_rise_eV` = E_total(final) − E_total(min) (eV, **lower** is
  better; 0 = monotone = fixed)
- **Secondary**: `excursion_eV` (max climb above t=0 ref), `drain_eV` (absorption
  context), `absorbed_e` (electrons removed, from charge.csv), `t_min_au`
- **Noise floor**: 0.1 eV (3× the worst clean-run rise: capon_reach +0.033 eV;
  deterministic propagator — no seed variance)

## Budget

- maxRuns: 24 | maxSeconds: none | targetMetric: 0.1 (rise at/below noise floor with
  checks passing)
- Per-experiment wall-clock cap: 2 h (`timeout 7200` in autoresearch.sh)

## How to Run

`cd ResearchProject/systems/localised_jellium/scripts/cap_fix && EM_...=... ./autoresearch.sh <run_name> <gpu_id>`
— outputs `METRIC name=value` lines. Deterministic (no RNG); no SEED argument.
Then `./checks.sh <run_name>` (correctness gate: the absorber must still absorb).
Two GPUs → run two experiments concurrently (one per GPU), still one atomic
config change per experiment.

## Files in Scope

- Experiment configs = **env vars only**: `EM_CAP_MODE` (two|wrap), `EM_CAP_ETA`,
  `EM_CAP_CENTER_BOHR`, `EM_CAP_WIDTH_BOHR`, `EM_WRAP_WIDTH_BOHR`. (`EM_N_STEPS` is
  fixed at 700 for the screening segment; a confirmation segment at 950 steps/η=−1
  re-inits with a new config header.)
- `hypotheses/cap_fix/autoresearch.jsonl`, `experiments/worklog.md`,
  `autoresearch-dashboard.md` — state (this loop's bookkeeping).

## Off Limits

- `scripts/cap_fix/autoresearch.sh`, `run_metrics.py`, `checks.sh` — **the eval
  harness is locked.**
- `scripts/cap_fix/run.cpp` and `inq-stack/include/inqkit/perturbations/absorbing_wrap.hpp`
  — the physics binary is part of the harness for this campaign (experiments are
  configs, not code edits). A genuine binary fix = deliberate re-baseline, noted.
- `inq/` (immutable, always), the GS checkpoint, the diagnosis campaign's results.
- Deviations from the stock skill (documented): no `git add -A` and no branch switch —
  the repo carries unrelated uncommitted work; commits are scoped to campaign files
  only. No `git checkout -- .` reverts — env-config experiments leave no code to revert.

## Constraints

- checks.sh must pass for any `keep` (absorbed_e ≥ 0.5, fallback drain ≤ −10 eV).
- GPU etiquette: probe `cudaMemGetInfo` before launching; never touch a GPU running
  another user's job (both free at session start, 2026-07-13).
- A config with rise ≤ noise floor at 700 steps is only PROVISIONALLY fixed — it must
  survive the 950-step η=−1 confirmation segment before being declared a fix.

## What's Been Tried

(see worklog; updated as runs accumulate)
- Pre-session (diagnosis campaign, reused as reference points): two-sided η=−0.2 →
  rise 23.49 eV (the artifact); two-sided η=−1 @700 → rise 0.03 (but rises +31 eV by
  t=36 in phase-0); CAP-off → rise 0.00.
