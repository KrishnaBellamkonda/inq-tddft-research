# Autoresearch Dashboard: cap-fix-artifact-removal

**Runs:** 11 / 24 (+2 in flight) | **Kept:** 4 | **Discarded:** 7 | **Crashed:** 0 (1 timeout, partial data used) | **Checks-failed:** 0
**Baseline:** artifact_rise_eV: 23.49 eV (#1, two-sided η=−0.2 @700)
**Best:** artifact_rise_eV: **0.000000 eV** (#11, wrap η=−2.0 w40 @950 — **TARGET REACHED**)
**Noise floor:** ±0.1 eV | **Target:** ≤ 0.1 eV with absorbed_e ≥ 0.5 → **MET (pending 1200-step hardening)**

## Segment 0 — screening (700 steps, t=28)

| # | artifact_rise_eV | status | description |
|---|------------------|--------|-------------|
| 1 | 23.49 | keep | baseline two-sided η=−0.2 (bit-identical reproduction) |
| 2 | 35.34 | discard | wrap equal-integral twin @weak η — topology alone refuted |
| 3 | 0.00 | keep | two-sided η=−2.0 — strength arm wins screening |
| 4 | 24.38 | discard | pushed-out footprint @weak η — geometry alone refuted |

## Segment 1 — confirmation window (950 steps, t=38; validated by #6)

| # | artifact_rise_eV | excursion | status | description |
|---|------------------|-----------|--------|-------------|
| 5 | 3.50 | 0.000 | keep | two-sided η=−2.0 — excursion gone, small rebound remains |
| 6 | 169.32 | **+31.27** | discard (control) | η=−1 positive control — reproduces phase-0; **window validated** |
| 7 | 20.21 | 0.000 | discard | two-sided η=−4 — over-strong reflects |
| 8 | 24.66 | 0.000 | discard | η=−2 × pushed footprint — geometry hurts at strong η too |
| 9 | 0.324 | 0.000 | keep | **wrap × η=−2 (w30) — the arms compose; 10× better** |
| 10 | 11.83 | 0.000 | discard | two-sided η=−3 — ladder convex, optimum ≈ −2; two-sided arm CLOSED |
| 11 | **0.000000** | 0.000 | **keep — TARGET** | **wrap η=−2 w40: strictly monotone to t=38** (caveat: absorbed 1.015e, small slab-tail nibble) |
| 12 | *running (GPU 0)* | | | wrap η=−1.5 w30 — wrap ladder left flank (possible gentler equal winner) |
| 13 | *running (GPU 1)* | | | winner @1200 steps (t=48) — hardening |

## The story in one paragraph

The artifact tracks HOW LONG slow density is being absorbed, not where or by what shape
(runs 2, 4, 8). Strength alone kills the above-zero excursion but leaves a rebound with a
convex η-optimum (transmit-vs-reflect, runs 5/7/10). The user's wrap-around topology —
refuted alone at weak η — is the missing half at strong η: it plugs the two-sided
profile's W=0 hole at the periodic boundary where slow spill leaks and lingers (run 9),
and with a gentler w=40 ramp the reported E_total becomes strictly monotone to t=38
(run 11). **Production recommendation (pending runs 12–13): unified wrap-around CAP,
η=−2.0, width 40 Bohr.**
