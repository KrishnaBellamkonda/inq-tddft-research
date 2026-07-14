# Handover: PBC vs open-z energy-oscillation test (Arm B)

Campaign: `docs/campaigns/localised_jellium/pbc-open-z-oscillation.md`
(`id: lj-pbc-open-z-oscillation`, status: running). Executor:
`ResearchProject/systems/localised_jellium/scripts/pbc_open_z/orchestrate.py`
(log: `scripts/pbc_open_z/orchestrate.log`).

## 2026-07-14 — designed via grill, launched on GPU 1

**User hypothesis (locked):** electrostatic z-periodicity (Arm B) plays a causal
role in the energy oscillation; expects p2 (open-z Poisson) results to differ from
p3. Arm A (recirculation) understood, NOT run. All experimentation on ONE GPU.
Preceded by the user's period-lengthening re-reading of cap_fix (verified: t_min
drifts 21.6→27.8→~33→36.4→>48 along the ladder; every long two-sided run turns up;
"no rise" runs have <1 period of data) — recorded in CONTEXT.md glossary
("Period-lengthening reading", "PBC-vs-open-z channels").

**Key facts verified before design:** INQ periodicity 0/2/3 supported
(cell.hpp:103; slab Poisson in solvers/poisson.hpp); ALL slab_n52 oscillation runs
to date are p3; campaign_autorun (N=82) used p2 with correct matched-GS hygiene but
holds no oscillation-witness energy data; wavefunction always wraps on the FFT grid
(p2 switches electrostatics only).

**Built (this session):**
- `scripts/pbc_open_z/gs/run.cpp` — p2-switchable slab_n52 GS (clone of
  effmass_sigma1/gs), checkpoint `shared_gs/slab_n52_L40x40x80_dx0p333_per2`,
  dumps GS density VTI for the gate.
- `scripts/cap_fix/run.cpp` — added `EM_PERIODICITY` (default 3, byte-compatible;
  documented harness re-baseline). run_summary now records periodicity.
- `scripts/pbc_open_z/orchestrate.py` — serial one-GPU chain, idempotent resume:
  GS → n(z) gate (interior ±30% of n0; tail <1e-4 at |z|>32.5; p3-GS-into-p2-RT
  forbidden) → rebuild+5-step p2 smoke → three p2 twins via the locked cap_fix
  harness (p2_two_eta0p2_700 / p2_two_eta1p0_950 / p2_wrap_eta2p0_w40_950, twins
  of run01/run06/run11) → auto comparison table
  `hypotheses/pbc_open_z/comparison.md`.
- Campaign doc + INDEX (35 campaigns) + CONTEXT.md glossary entries.

**Discriminator:** per-pair Δ(t_min, drain, rise, excursion) beyond the 0.1 eV
noise floor ⇒ Arm B implicated; shape-identical ⇒ Arm B refuted → suspicion moves
to recirculation (Arm A, L_z-scaling test) or ledger bookkeeping as the clock.
Compare ΔE_total(t) shapes ONLY (absolute energies are convention-dependent).

**State at handover:** orchestrator IN FLIGHT on GPU 1 (launched ~05:14; GPU 0
occupied by another user — untouched). ETA ≈ 5 h (GS ~0.5–1 h + smoke + 55/75/75
min runs). NOT yet done: comparison verdict, study addendum, task flag flips,
commit of pbc_open_z files (pending first results).

## 2026-07-14 (08:03) — CHAIN COMPLETE, VERDICT: Arm B REFUTED as cause/clock

Full chain ran clean on GPU 1 (05:14→08:03): p2 GS converged + gate PASS (interior
n(z) on n0; tail 2.2e-11 at |z|>32.5 — 6 orders tighter than p3), knob+smoke OK,
three p2 twins + auto comparison (`hypotheses/pbc_open_z/comparison.md`).

**Result:** the oscillation is FULLY REPRODUCED at periodicity 2 with identical
clock — t_min shifts only 0.2 a.u. in both two-sided witnesses (21.8/21.6 and
28.0/27.8); same drain-then-rise morphology; the wrap winner stays clean at both.
Amplitudes are slightly LARGER at p2 (rise +23.74/+174.1 vs +23.49/+169.3 eV;
excursion +1.35/+37.6 vs +0.11/+31.3 eV) — beyond the 0.1 eV floor, so open-z
electrostatics MODULATES the artifact's size (worse, not better) but does not
cause it or set its period. User's suspicion ("results would be different"):
amplitudes differ, mechanism/clock does not.

**Suspicion moves to:** Arm A (recirculation; period-vs-L_z scaling test, designed
but not run) and/or the diagnosis bookkeeping mechanism's own clock (slow-density
arrival at the CAP — consistent with t_min tracking CAP config, not Poisson
convention). Campaign `status: done` 4/4; INDEX regenerated. Uncommitted at this
milestone: pbc_open_z scripts, campaign, comparison, CONTEXT.md entries — commit
next.
