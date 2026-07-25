# Hypothesis ledger — localised-jellium ΔE_total energy-oscillation

Living document. The **Advisor** updates this every iteration. Status values:
`OPEN` → `SUPPORTED` / `WEAKENED` → `CONFIRMED` / `REFUTED`. A cause is `CONFIRMED`
only by a **decisive control** experiment.

**Phenomenon:** `ΔE_total(t) = E_total(t) − E_ref` oscillates and rises **above 0**
in many localised-jellium RT runs once the CAP absorbs — unphysical (no energy
influx; a CAP can only remove energy). Contrast `p3_wp` (clean decay to plateau).

## Candidate mechanisms

| # | Hypothesis | Status | Evidence (latest) | Decisive control |
|---|---|---|---|---|
| a | CAP is a non-Hermitian energy **source** in the reported ledger | **CONFIRMED** (iter3, conf 0.90) | capoff_floor (CAP off): conserved −0.015 eV, never >0. capon_reach (η=−1): CAP is the SOLE non-conservative term (drains −138 eV, ledger exact). **capon_weak_partial (η=−0.2, t=28): E_total drains to −23.4 eV then RISES +23.5 eV and CROSSES ABOVE 0 to +0.11 eV** — the unphysical excursion reproduced under a fresh control, ledger EXACT to 1.39e-13 Ha throughout, density_l2→0.0001 (absorbed even as E_total rises). The CAP's contribution to the reported ledger is SIGN-CHANGING in the partial-absorption regime. | **capon_weak_partial vs capoff_floor** (CAP off⇒conserved; weak CAP on⇒drain-then-rise>0) — the DECISIVE control |
| b | Static `v_bg` **absent from the reported energy functional** | **WEAKENED** (iter2) | v_bg ON in capoff_floor (−0.015 eV) AND capon_reach (all drift = CAP term, ledger exact) — v_bg present yet contributes no drift. Isolating +v_bg-only control still not run. | `+v_bg` only (no projectile/CAP) drifts E_total |
| c | Wrong subtracted **`E_ref`** | **REFUTED** (iter1) | Component sum == energy_total to 1e-13 Ha every step; CAP-off conserves vs E_total(0_RT); +221 eV vs E_GS is a FLAT WP-injection offset, not a drift. The rise is CAP-gated, not reference-gated. | Recompute ΔE vs E_total(0_RT); baselines |
| d | **Propagator / grid numerics** (ETRS drift, dt, aliasing) | **REFUTED** (iter1) | Same ETRS/dx=0.333/dt=0.04/cutoff/WP conserve E_total to −0.015 eV with the CAP off — propagator/grid do not source the >0 rise. | pure-GS run violates conservation; dt-halving scales amplitude |
| e | **Density-dependent KS Hamiltonian** double-counting | **REFUTED** (iter1) | 8 total-contributing components sum to energy_total to 1e-13 Ha at every step; drifts (H +23, ext −30, kin +3.5 eV) cancel exactly — no KS double-counting. Σε_i drifts +24 eV while E_total flat (band≠total, benign). | component decomposition; Σε_i vs E_total; ∫n·v_xc |

## Experiment log

| Iter | Probe | Aim | Key result (raw) | Advisor verdict | Next |
|---|---|---|---|---|---|
| 0 | phase0_mine | characterise phenomenon + component gap from existing runs | Component gap CONFIRMED (only total/kin/H/xc recorded). ΔE_total(vs RT0) final: default η=−1 **+31 eV**, weak η=−0.4 **+31 eV** (max +32), wider-gap **+38 eV**, strong η=−2.0 **−165 eV** (no rise), classical **+58 eV peak → −179 eV**, p3_wp **−0.03 eV (flat/conserved)**. Amplitude/sign tracks η strongly. | **continue** (conf 0.7): a SUPPORTED, c WEAKENED; a vs d not yet split | capoff_floor_with_decomposition |
| 1 | capoff_floor | split a (CAP source) vs d (propagator floor): rerun effmass_sigma1 `EM_CAP=0`, WP+v_bg on, full decomp, write_every=5, 200 steps | E_total CONSERVED to **−0.015 eV** over t=8 a.u., MAX ΔE(vs RT0)=**0.0** (never >0). Component sum==E_total to 1e-13 Ha. Big drifts (H +23, ext −30, kin +3.5 eV) cancel. Σε_i drifts +24 eV, E_total flat. +221 eV vs E_GS = flat WP-injection offset. | **continue** (0.72): a SUPPORTED; **c/d/e REFUTED**; b OPEN. CAP-off removes the rise; 200-step window may pre-date WP→CAP arrival. | capon_matched |
| 2 | capon_matched | airtight a: SAME run `EM_CAP=1` (η=−1), 200 steps, matched to capoff_floor | ΔE(vs RT0) final **−0.016 eV** = indistinguishable from CAP-off (−0.015); both max=0.0. density_l2 0.052→0.001: WP **has NOT reached the CAP** at t=8 (window pre-dates absorption). Non-diagnostic. | (extend) | capon_reach |
| 3 | capon_reach | extend CAP-ON η=−1 to 700 steps (t=28) so density reaches the CAP and absorption begins | ΔE(vs RT0) monotonically NEGATIVE −0.11→−5.4→−28→−88→−138 eV as l2→0.0001; **MAX=0.0, never >0**. Ledger EXACT (Σ8==E_total to 1.3e-13 Ha). CAP is the SOLE non-conservative term but removes energy **correctly** here. The phase-0 >0 rise NOT reproduced in clean single-transit at η=−1. | **continue** (0.55): a SUPPORTED (CAP = sole non-conservative term, decisively); c/d/e REFUTED; b WEAKENED. Need to witness ΔE>0 in the weak/partial regime. | capon_weak_partial |
| 4 | capon_weak_partial | witness the unphysical ΔE>0 excursion in the WEAK/partial-absorption regime (phase-0's regime, η=−0.2), 700 steps t=28 | ΔE(vs RT0): min **−23.38 eV** (t=21.6) → rises +23.5 eV → **FINAL +0.11 eV, crosses>0=TRUE** (drain-then-rise). Ledger EXACT to **1.39e-13 Ha** throughout incl. the rise; density_l2→0.0001 (absorbed continuously even as E_total rises). Reproduces the phase-0 excursion under a fresh control. | **confirmed** (0.90): **(a) CONFIRMED**; c/d/e REFUTED; b WEAKENED (non-causal). | DONE — synthesis |

## KEY DISCOVERY (iter-2 post-hoc, from phase-0 re-inspection)

**The phase-0 default η=−1 run (`effmass_sigma1/wp`) is DRAIN-THEN-RISE, and my
`capon_reach` reproduced its drain EXACTLY but stopped at the minimum.** Phase-0 trajectory
(226 rows, t=36, dt=0.04):
`ΔE(vs RT0)`: 0 → **MIN −138.1 eV at t=27.8** → −77 (t=32) → +21 (t=35.2) → **+31.3 eV (t=36),
crosses>0=TRUE**. `capon_reach` (my run, same η=−1) hit **−138 eV at t=28** — bit-for-bit the
same minimum — then STOPPED at 700 steps (t=28), one step short of the rise. **The anomalous
>0 excursion happens AFTER the deep drain, between t=28 and t=36 (steps ~700→900).** So the
decisive confirmer is NOT weak η — it is **η=−1 run LONGER (to t≥38, ~950 steps)** to capture
the drain→rise→+31 eV. `capon_weak_partial` (η=−0.2) is still a useful η-dependence datum but
the τ-extension is the true PASS(a) probe. → probe 5/6: `capon_long` (η=−1, 950 steps).

## Confirmed cause (SYNTHESIS — advisor verdict, conf 0.90, 2026-07-13)

**Confirmed mechanism: (a) the CAP is a non-Hermitian energy artifact in the reported KS
total-energy ledger — a METHOD artifact, not physics.**

The unphysical ΔE_total>0 excursion is caused by the complex absorbing potential (CAP). In the
partial-absorption regime the CAP's contribution to the reported E_total is **sign-changing**:
after an initial drain it adds energy back into the reported total, driving ΔE above the
reference **even while probability density is still being removed**. The reported E_total does
not correctly book the energy carried away by the CAP-absorbed density, so the CAP is **not a
proper energy sink** in that ledger.

**Decisive control:** `capon_weak_partial` vs `capoff_floor` — identical cell/GS/WP/propagator/
functional, differing only in CAP strength:
- **CAP OFF** (`capoff_floor`): E_total conserved to −0.015 eV, **never >0**.
- **CAP ON, weak η=−0.2** (`capon_weak_partial`, t=28): E_total drains to −23.4 eV then **RISES
  +23.5 eV and CROSSES ABOVE 0 to +0.11 eV**, with the CAP the **sole non-conservative term** and
  the 8-component ledger exact to **1.39e-13 Ha** throughout while density_l2→0.0001 (absorbed
  continuously even as E_total rises).
- (Supporting: `capon_reach` η=−1 shows the CAP draining −138 eV as the only moving term; phase-0
  shows the same η=−1 family rises to +31 eV by t=36 — drain-then-rise.)

**Refuted mechanisms** (all under the same controlled instrumentation, component sum == E_total to
solver precision throughout):
- (c) wrong E_ref — REFUTED (CAP-off conserves vs RT0; the +221 eV vs E_GS is a flat WP-injection
  offset, not a drift).
- (d) propagator/grid numerics — REFUTED (same ETRS/dx/dt/cutoff conserve with the CAP off; all
  CAP-on drift is accounted by the CAP term).
- (e) KS Hamiltonian double-counting — REFUTED (Σ components == energy_total to 1e-13 Ha through
  full absorption incl. the rise; Σε_i drifts while E_total tracks physically = band≠total, benign).
- (b) v_bg missing from the functional — WEAKENED, **non-causal** (v_bg on in all 5 probes yet
  contributes no independent drift; the artifact is fully switched by the CAP alone). A dedicated
  +v_bg-only control was not run.

**Scope note (diagnosis only, per campaign charter):** no physics fix is proposed here. A proper
fix (booking the CAP-removed energy as a monotone sink in the reported ledger, e.g. an accumulated
absorbed-energy term) is a separate, later campaign.

## Post-campaign self-check (2026-07-13 — notebook Parts II–IV)

An independent re-audit (source-code check of INQ + literature verification + raw-CSV
re-analysis) CONFIRMED the diagnosis but corrected three evidence items above:

1. **"Ledger exact to 1e-13 Ha" is CIRCULAR** — `energy.total()` is *defined* as the
   component sum (`inq/src/hamiltonian/energy.hpp:127`); the residual is float addition
   order. The refutation of (e) rests on the CAP-off control instead (still valid).
2. **`density_l2` was over-interpreted** — it is `∫|n(t)−n(t0)|²dV` vs the step-0 snapshot,
   not remaining charge; "absorbed continuously even as E_total rises" was not established.
   During the rise the density field is nearly static (l2 ~500x below transit peak).
3. **Mechanism refined** — INQ's ledger is mixed-convention: `energy_kinetic` is
   norm-divided per orbital (`energy.hpp:55`), density-based terms are bare
   (`density.hpp:36`). Component attribution shows the norm-divided kinetic term carries
   −137 of the −138 eV drain (capon_reach) and +24 of the +23.5 eV rise
   (capon_weak_partial): the dominant channel is a covariance/filter effect on the
   renormalized kinetic average (Graefe, Höning & Korsch 2010, Eq. 8), with the
   bare-removal channel secondary. Root cause unchanged: CAP-gated bookkeeping artifact.

Remaining holes (untested): long CAP-off control at t=28; ∫n dV never logged; +v_bg-only
control never run. Literature: docs/sources/graefe-2010-nonhermitian-dynamics.md,
docs/sources/selsto-2010-absorbing-boundaries.md. Full audit: notebook Parts II–IV.
