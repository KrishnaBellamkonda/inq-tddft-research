---
id: lj-cap-fix-experimentation
area: localised_jellium
title: CAP energy-artifact removal — wrap-around topology + parameter tuning
status: running
hypothesis: "A CAP setup exists (unified wrap-around topology and/or tuned eta/footprint) in which the witness run's reported E_total decays monotonically (artifact_rise_eV <= 0.1 eV noise floor) while still absorbing the wavepacket (absorbed_e >= 0.5) — and the winning setup stays clean in the eta=-1, 950-step confirmation regime."
handover: docs/handovers/cap-fix-experimentation.md
tasks:
  - { name: "harness: cap_fix binary (EM_CAP_MODE + charge.csv) + wrap perturbation + smoke", done: true }
  - { name: "baseline re-run on new binary (two-sided eta=-0.2, 700 steps) reproduces rise ~23.5 eV", done: true }
  - { name: "draft experiments: wrap-around twin, strong-eta, pushed-out footprint", done: true }
  - { name: "improve loop: refine best draft until rise <= 0.1 eV with checks passing", done: true }
  - { name: "confirmation segment: winner(s) at eta=-1, 950 steps (phase-0 regime) stay clean", done: false }
  - { name: "study notebook + ledger update + recommended production CAP config", done: false }
blocked_reason: ""
---

# CAP energy-artifact removal (cap_fix)

**User decision (2026-07-13, locked in-conversation):** run a fix campaign for the
diagnosed CAP energy-ledger artifact. Candidate fixes named by the user: (1) a unified
CAP that wraps around the periodic boundary (instead of the current two bumps), and
(2) tuning η, CAP length, and the other CAP parameters until the effect disappears.
Method: the `autoresearch` skill's locked-harness experiment loop (installed from
github.com/drivelineresearch/autoresearch-claude-code, MIT, see
`.claude/skills/autoresearch/PROVENANCE.md`); small runs, iterate, learn; start
immediately on the free GPUs. The stage-gated grill was superseded by these explicit
user locks — recorded here instead of a Stage-1..5 transcript.

## Question

Which CAP setup removes the drain-then-rise artifact from the reported E_total while
keeping the absorber functional — topology (wrap-around vs two-sided), strength (η),
or footprint (centre/width)? And does the winner generalise to the η=−1 long-window
regime (phase-0's +31 eV riser)?

## Prior knowledge (diagnosis campaign, conf 0.90 — do not re-derive)

- The rise is a CAP-gated bookkeeping artifact: the ledger books no absorbed energy;
  dominant channel = norm-divided kinetic filtering (`inq/src/hamiltonian/energy.hpp:55`;
  Graefe et al. 2010 Eq. 8). Witness numbers: two-sided η=−0.2, 700 steps → drain −23.4,
  rise +23.5 eV; η=−1 @700 monotone −138 eV but +31 eV by t=36; CAP-off conserved.
- The two-sided `perturbations::absorbing` bumps fall to W=0 exactly at the periodic
  boundary (sin² profile peaking at ±32.5 Bohr; fractional coords,
  `inq/src/perturbations/absorbing.hpp:44`) — the wrap-around hypothesis targets this.
- Mechanism predictions (test, don't assume): topology alone does NOT remove the rise
  (the rise needs W·(bound/slow density) overlap, which both topologies have); pushed-out
  footprint and/or strong η DO reduce it; weak-η-in-overlap is worst. A fixed reported
  ledger may ultimately need the absorbed-energy accumulator (separate campaign if this
  one falsifies the setup-only fix).

## Harness (locked; ADR-0007 placement)

- Binary: `ResearchProject/systems/localised_jellium/scripts/cap_fix/run.cpp` — clone of
  the diagnosis ablation binary + `EM_CAP_MODE=two|wrap` + `charge.csv` (∫n dV per write
  step; closes the diagnosis Part-IV gap). Wrap topology:
  `inq-stack/include/inqkit/perturbations/absorbing_wrap.hpp` (cos² bump peaking at the
  boundary plane; width 30 Bohr ⇒ same footprint |z|>25 and same ∫W dz as the two-sided
  default — topology is the only difference). Engine: inq-study (mass fork). `inq/`
  untouched.
- Loop protocol: autoresearch skill. State in `hypotheses/cap_fix/` (autoresearch.md,
  autoresearch.jsonl, experiments/worklog.md, autoresearch-dashboard.md). Harness
  scripts in `scripts/cap_fix/` (autoresearch.sh, run_metrics.py, checks.sh) — locked.
- Primary metric `artifact_rise_eV` (lower better), noise floor 0.1 eV, correctness
  gate absorbed_e ≥ 0.5 (fallback drain ≤ −10 eV). Budget: maxRuns 24 screening
  (700 steps, ~45 min each, 2 GPUs), then a re-init confirmation segment (950 steps,
  η=−1) for winners.

## Experiment matrix (draft-first, then improve)

| # | op | Config (one atomic change vs baseline) | Tests |
|---|----|----------------------------------------|-------|
| 1 | baseline | two-sided, η=−0.2, centre 32.5, width 15 | reproduces rise 23.5 eV on the new binary + charge.csv sanity |
| 2 | draft | **wrap**, η=−0.2, width 30 (equal-integral twin) | user hypothesis 1: topology |
| 3 | draft | two-sided, **η=−2.0** | strong-CAP arm (phase-0 hint: no rise in-window) |
| 4 | draft | two-sided, η=−0.2, **centre 35, width 10** (footprint 30→40) | pushed-out footprint arm |
| 5+ | improve | refine the best arm (η ladder; width/centre ladder; wrap×strong-η cross) | one lever per run |
| conf | confirm | winners at η=−1 equivalent strength, **950 steps** | phase-0 regime stays clean |

Success: a config with rise ≤ 0.1 eV in BOTH segments, checks passing. Failure of all
setup-only arms = evidence the fix must be the absorbed-energy ledger term (next
campaign), itself a valuable falsification.

## Guard rails

- GPU probe via `cudaMemGetInfo` (NVML broken, harmless) before each launch; both GPUs
  free at session start; never touch another user's job.
- Abort a run on NaN/complex energy in observables.csv; `timeout 7200` per run.
- Cutoff/aliasing guard: setup identical to the diagnosis ablation binary (same GS, WP,
  spacing, dt) whose guard PASSED — carried over, not re-run.
- No pilot gate on v-drift (light-projectile rule; not S(v) extraction anyway).
- Checkpointing: not added — runs are 45–70 min; a kill loses one run, acceptable per
  checkpoint-dont-block (its ~200-step-ckpt precondition targets multi-hour runs).

## Output contract

- Per-run: `scripts/cap_fix/results/<run_name>/` (observables.csv, charge.csv,
  run_summary.txt) — gitignored; provenance in the worklog + dashboard + JSONL.
- Study notebook `hypotheses/cap_fix/cap_fix_study.ipynb` at campaign end: per-arm
  ΔE_total(t) overlays (shared axes), the winner's confirmation run, the keep/discard
  tree, and the recommended production CAP config (or the falsification verdict).
- Handover: `docs/handovers/cap-fix-experimentation.md` (rolling).
