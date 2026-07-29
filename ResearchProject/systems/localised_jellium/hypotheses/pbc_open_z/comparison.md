# Arm B: p2 (open-z electrostatics) vs p3 (fully periodic)

| pair | conv | run | t_min | drain (eV) | rise (eV) | excursion (eV) |
|---|---|---|---|---|---|---|
| p2_two_eta0p2_700 | p2 | p2_two_eta0p2_700 | 21.8 | -22.39 | +23.739 | +1.346 |
| p2_two_eta0p2_700 | p3 | run01_baseline_two_eta0p2 | 21.6 | -23.38 | +23.495 | +0.110 |
| p2_two_eta1p0_950 | p2 | p2_two_eta1p0_950 | 28.0 | -137.57 | +174.131 | +37.623 |
| p2_two_eta1p0_950 | p3 | run06_poscontrol_eta1p0_950 | 27.8 | -138.05 | +169.324 | +31.272 |
| p2_wrap_eta2p0_w40_950 | p2 | p2_wrap_eta2p0_w40_950 | 38.0 | -176.31 | +0.000 | +0.002 |
| p2_wrap_eta2p0_w40_950 | p3 | run11_wrap_eta2p0_w40_950 | 38.0 | -178.31 | +0.000 | +0.000 |

## Verdict (2026-07-14)

**Arm B REFUTED as the cause/clock of the oscillation.** With open-z electrostatics
(p2, matched p2 GS):

1. **The oscillation is fully reproduced** — same drain-then-rise morphology in both
   two-sided witnesses, and the wrap winner stays clean (no turn to t=38) at both
   periodicities.
2. **The clock is periodicity-independent**: t_min shifts by only 0.2 a.u.
   (21.8 vs 21.6; 28.0 vs 27.8) — far below any period-relevant scale. Whatever
   sets the oscillation period, it is NOT the z-image electrostatics.
3. **Amplitudes are slightly LARGER at p2** (beyond the 0.1 eV floor): rise
   +23.74 vs +23.49 eV (weak η), +174.1 vs +169.3 eV (η=−1); excursion above zero
   +1.35 vs +0.11 eV and +37.6 vs +31.3 eV. Open-z electrostatics *modulates* the
   artifact's size (direction: worse, plausibly via the tighter p2 GS density and
   the physically-escaping field of the net-charged cell) but does not create it.

**Where suspicion moves:** Arm A (density recirculation through the always-periodic
FFT wavefunction — testable via period-vs-L_z scaling) or the intrinsic slab
spill/slosh timescale interacting with the CAP (the diagnosis's bookkeeping
mechanism, whose clock is the slow-density arrival time, consistent with t_min
drifting with CAP config but not with Poisson convention).

Provenance: p2 GS `shared_gs/slab_n52_L40x40x80_dx0p333_per2` (gate: interior n(z)
on n0, tail 2.2e-11 at |z|>32.5 — note: ~6 orders tighter than p3, yet the
oscillation persists, further weakening any GS-tail-electrostatics story). Runs:
`scripts/cap_fix/results/p2_*` (periodicity recorded in each run_summary.txt).
