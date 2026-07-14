---
id: lj-pbc-open-z-oscillation
area: localised_jellium
title: Does electrostatic z-periodicity drive the energy oscillation? (Arm B)
status: done
hypothesis: "The post-minimum rise / period-lengthening of the reported E_total in slab_n52 witness runs is driven (at least in part) by the ELECTROSTATIC z-periodicity channel: with slab-truncated Poisson (periodicity 2, open-z electrostatics, matched p2 GS) the oscillation's t_min / rise amplitude / excursion differ from their periodicity-3 twins beyond the 0.1 eV noise floor. User suspicion: results WILL differ."
handover: docs/handovers/pbc-open-z-oscillation.md
tasks:
  - { name: "p2 GS for slab_n52 (task 0) + n(z) sanity gate", done: true }
  - { name: "EM_PERIODICITY knob in cap_fix binary + p2 smoke", done: true }
  - { name: "p2 twins: two eta=-0.2@700, two eta=-1@950, wrap eta=-2 w40@950", done: true }
  - { name: "p2-vs-p3 comparison table + verdict + handover", done: true }
blocked_reason: ""
---

# PBC vs open-z: the electrostatic-periodicity channel (Arm B)

**User decision (2026-07-14, locked in this grill):** after re-reading the cap_fix
results as *period-lengthening rather than artifact-removal* (t_min drifts 21.6 →
27.8 → ~33 → 36.4 → >48 along the "improvement" ladder; every sufficiently long
two-sided run turns up; runs with no turn have <1 period of data), the user
hypothesises the periodic boundary plays a causal role. The grill split "PBC vs
open-z" into two channels: **Arm A (density recirculation)** — the FFT wavefunction
always wraps in z; only absorbers stop it; testable only indirectly via L_z-scaling
of the period — **understood, NOT run**. **Arm B (electrostatic periodicity)** —
Hartree images along z + charged-cell G=0 convention; directly switchable via INQ
`periodicity(2)` (slab-truncated Poisson, `inq/src/systems/cell.hpp:103`,
`solvers/poisson.hpp`) — **this campaign. One GPU only (GPU 1; GPU 0 occupied by
another user at launch).**

## Fact base (verified in-session)

- ALL prior slab_n52 oscillation runs (diagnosis + cap_fix, 14 runs) are
  periodicity 3 (`.periodic()`, cap_fix/run.cpp + effmass_sigma1/wp/run.cpp).
  Arm B has never been tested on this lineage.
- Prior p2 usage exists elsewhere with correct GS hygiene (campaign_autorun N=82:
  GS converged at matching periodicity; matched pair in shared_gs/…L50x50x111_pbc/
  _per2) — but no oscillation-witness energy data (results dirs hold build smokes).
- Absolute energies are convention-dependent across p2/p3 (charged-cell G=0):
  compare ONLY ΔE_total(t) shapes (t_min, drain, rise, excursion), never absolute.

## Design (locked)

- **GS (task 0):** fresh p2 GS for slab_n52 (`scripts/pbc_open_z/gs/run.cpp`,
  clone of the p3 producer effmass_sigma1/gs with the periodicity knob), saved to
  `shared_gs/slab_n52_L40x40x80_dx0p333_per2`. Gate before any RT: interior n(z)
  within 30% of n0, tail < 1e-4 at the CAP footprint (|z|>32.5), converged.
  Loading the p3 GS into p2 RT is FORBIDDEN (not an eigenstate of the p2
  Hamiltonian → spurious t=0 kick contaminating exactly the signal under test).
- **Binary:** cap_fix binary + `EM_PERIODICITY` env knob (default 3 →
  byte-compatible with all prior runs; deliberate documented harness re-baseline).
  Same locked metrics (run_metrics.py): t_min, drain, rise, excursion; noise floor
  0.1 eV.
- **Runs (p2 twins of the three bracketing p3 witnesses, serialized on one GPU):**

| p2 run | p3 twin | why this witness |
|---|---|---|
| p2_two_eta0p2_700 | run01 (rise 23.49, t_min 21.6) | fastest full-turn oscillation signal |
| p2_two_eta1p0_950 | run06 (excursion +31.27, t_min 27.8) | the big above-zero riser |
| p2_wrap_eta2p0_w40_950 | run11 (rise 0.000 to t=38) | does the winner's cleanliness survive p2? |

- **Discriminator:** for each pair, differences in t_min / rise / excursion beyond
  0.1 eV (and qualitative period change) ⇒ Arm B implicated; shape-identical pairs
  ⇒ Arm B refuted for this system, pointing back to recirculation (Arm A) or the
  ledger bookkeeping as the oscillation's clock.
- **Executor:** `scripts/pbc_open_z/orchestrate.py` (Python, idempotent resume,
  per-step traceback, correctness gates block, cost never blocks), all steps
  pinned to one GPU via `PBC_GPU`.

## Guard rails

- GPU probe before launch (cudaMemGetInfo; NVML mismatch harmless); GPU 0 occupied
  by another user at launch → everything on GPU 1; never touch GPU 0.
- Same cutoff/aliasing situation as cap_fix (identical grid/WP) — carried over.
- Abort on NaN/complex energy; timeout 2 h per RT run (7200 s, in the harness).
- 45–75 min runs: no checkpointing needed (checkpoint-dont-block threshold).

## Output contract

- `hypotheses/pbc_open_z/comparison.md` (auto-written by step 4) + study addendum
  notebook after verdict; handover `docs/handovers/pbc-open-z-oscillation.md`.
- Runs in `scripts/cap_fix/results/p2_*` (gitignored; provenance via summaries).
