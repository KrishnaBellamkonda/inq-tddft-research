# 03 — Results (recorded outcomes per task)

Neutral transcription of the recorded results. Numbers are copied from the
`artifacts/*.json` result files (also in `data/`); figures are in `figures/`.
Where the record itself notes a flag (e.g. a panel review, a provenance
measurement), it is reported as a **recorded fact**, not as a judgment. All
interpretation is left to the reader. Numbers rounded to 2–3 s.f.

---

## Era 1 — Signature gates (T1–T7)

**T1 — kernel pre-gate.** POD, DMD, form-factor formula-validation: **CONFIRM**
(×3, independent). Code-tests: **12/12 pass**. `F_ONCV(q) ≈ 1` within 5 % for
`q ≤ 1.9 /Bohr` (within 2 % for `q ≤ 1.2`) → the T2 prediction reduces to
`exp(−q²σ_pot²/2)` inside that window. (`foncv.json`, `T1_foncv.png`.)

**T2 — form-factor cut. Verdict: CONFIRM** (2/3 held-out agree). Held-out
(`T2_result.json`, agree = ≥50 % of q-window within ±20 %):

| σ_WP | σ_pot | frac within ±20 % | median rel. | σ_eff (fit) | agrees |
|---|---|---|---|---|---|
| 0.5 | 0.35 | 0.89 | 0.14 | 0.42 | yes |
| 3.0 | 2.1 | 0.20 | 0.52 | 3.2 | no |
| 8.0 | — | (third held-out) | — | — | — |

**T3 — wake gate. Verdict: INCONCLUSIVE** (1/3 agree). Held-out (`T3_result.json`,
agree = `|ω_DMD − ω_p|/ω_p ≤ 0.20` and Nyquist ok):

| E (eV) | v | ω_DMD (eV) | ω_p (eV) | λ_DMD | λ_theory | agrees | note |
|---|---|---|---|---|---|---|---|
| 25 | 1.4 | 12.0 | 3.5 | 19 | 67 | no | total_wp_included (caveat) |
| 100 | 2.7 | 3.5 | 3.5 | 130 | (≫L) | yes | — |
| 600 | 6.6 | — | 3.5 | — | — | — | — |

**T4 — localised slab. Verdict: INCONCLUSIVE.**

**T5 — dynamics. Status: DONE** (2-mode latent ODE fitted; cell E=100,
`total_wp_included` caveat; no forward-prediction test at this stage).
(`T5_result.json`, `T5_dynamics.png`.)

**T6 — exchange/diffraction. Status: EXPLORATORY** (carries SIE ≈ 7 eV caveat).

**T7 — synthesis roll-up** (`T7_result.json`): T2 CONFIRM (2/3), T3 INCONCLUSIVE
(1/3), T4 INCONCLUSIVE, T5 DONE, T6 EXPLORATORY.

---

## Era 2 — Governing-PDE discovery (T8–T14)

**T9 — PDE-FIND kernel validation:** recovered known synthetic PDEs
(heat/advection/wave); **6/6 tests pass**.

**T11 — `PDE_classical`. Verdict: REFUTE** (0/3 held-out validated).
Frozen config: order 2, threshold 0.06, poly 2, deriv_order 3, POD rank 6.
Held-out (`PDE_T11_result.json`):

| E (eV) | v | admitted equation | forward rel-L2 |
|---|---|---|---|
| 25 | 1.36 | `u_tt = 1.4·u_xx + 0.35·u·u_x + 0.12·u·u_xx + 0.21·u·u_xxx + …` | 0.91 |
| 100 | 2.71 | `u_tt = 0` (all terms thresholded) | 46.9 |
| 600 | 6.64 | (see JSON) | — |

**T12 — `PDE_WP`. Verdict: CONFIRM** (3/3 held-out validated).
Frozen config: order 2, threshold 0.04, poly 2, deriv_order 2, POD rank 8.
Held-out (`PDE_T12_result.json`):

| E (eV) | v | admitted equation | forward rel-L2 | validated |
|---|---|---|---|---|
| 25 | 1.36 | `u_tt = 1.4·u_xx` | 0.26 | yes |
| 100 | 2.71 | `u_tt = 7·u_xx` | 0.49 | yes |
| 600 | 6.64 | `u_tt = 43·u_xx` | 0.05 | yes |

Recorded numeric relation (coefficient vs velocity): the discovered coefficient
`c²` in `u_tt = c²·u_xx` is `{1.4, 7, 43}` at `v = {1.36, 2.71, 6.64}`, i.e.
`√c² = {1.18, 2.65, 6.56}` vs `v = {1.36, 2.71, 6.64}`.

**T13 / T14 — compare + synthesis** (`PDE_T13_result.json`, `PDE_T14_result.json`,
`T13_compare.png`).

**Panel review + re-analysis (recorded facts, 2026-07-04):**
- A 4-expert scientific panel reviewed the T12 result and flagged it (transcript
  referenced in the handover).
- Independent provenance re-measurement: `run_wp_n162_L50_E100` `density_system`
  = **163 e** (WP-included); matched `_wf` run `total − wp` = **162 e**
  (bath-only). The T11/T12 discovery used the axial reduction of these fields.
- Re-analysis on the clean `total − wp` bath (E100 v=2.71; `run_base…E1p5`
  v=0.33): **no low-order PDE validated** in either testable regime. On-axis line
  cut peak-to-peak ≈ 1.2 × 10⁻³ (≈30× the transverse-mean amplitude).

(These are recorded measurements; the reader draws the conclusion.)

---

## Era 3 — POD/DMD bath-structure sweep (the "true blob-free bath")

**σ-sweep, fixed v = 2.71** (`bathstruct_summary.json`, `bathstruct_sigma_sweep.png`):

| σ_WP | POD rank(90 %) | leading-mode energy | dominant DMD (eV) |
|---|---|---|---|
| 0 (classical point) | 4 | 0.62 | 210 |
| 0.5 | 1 | 0.94 | 9.7 |
| 3 | 1 | 0.93 | 7.3 |
| 8 | 1 | 0.99 | 11.4 |

**Classical velocity sweep** (`bathstruct_velocity_sweep.png`):

| E (eV) | v | POD rank(90 %) | leading energy | dominant DMD (eV) |
|---|---|---|---|---|
| 20 | 1.21 | 3 | 0.62 | 225 |
| 25 | 1.36 | 3 | 0.63 | 250 |
| 50 | 1.92 | 3 | 0.67 | 322 |
| 100 | 2.71 | 4 | 0.62 | 210 |
| 600 | 6.64 | 4 | 0.67 | 613 |

**WP velocity points (mixed σ/grid — recorded, caveated):**

| run | σ_WP | v | POD rank(90 %) | leading energy | dominant DMD (eV) |
|---|---|---|---|---|---|
| `base…E1p5` | 5 | 0.33 | 1 | 0.92 | 1.6 |
| `plasmon…varyv` | 3 | 0.50 | 12 | 0.25 | 476 |
| `wp…E100_σ3` | 3 | 2.71 | 1 | 0.93 | 7.3 |

Recorded caveats (from the summary): the runs are non-stationary → DMD frequencies
are approximate (POD rank / coherence is the robust descriptor); DMD growth rates
> 0; single density (`r_s = 5.69`); the `varyv` point uses a coarse `dt = 4.0`.

---

## Era 4 — Linear-response residual test

**σ = 0.5 (recorded as the only SNR-adequate pair)** — radial channel
(`linres_residual_sigma0.5.json`, `linres_residual_sigma0.5.png`):
- Fitted exponent `a = 0.192`, `σ_fit = 0.62`, fit R² = 0.63; `matches_sigma_wp =
  true` (σ_fit near σ_WP = 0.5 rather than σ_pot = 0.35).
- High-q excess over noise ≈ 0.15σ.
- Temporal flatness of `|R(q,t)|`: full window **0.43** (`t_flatness`); early half
  **0.19**, late half **0.19**. (Recorded flag: `snr_adequate = true`,
  `decel = true` in the summary.)
- Overlap window `T ≈ 4.8 a.u.`, 59 common frames.

**σ = 3, σ = 8 (recorded as SNR-dead, excluded):** form factor e-folds within
~1–4 shells (`q_efold = √2/σ = 0.47 / 0.18`); descending-arm fit hits the
blob-subtraction noise floor → `σ_fit = 1.46 / 1.26` (flagged not adequate).

**Fork A (√2 trap):** recorded as **INCONCLUSIVE on existing data** — only one
SNR-adequate σ (need ≥2); the available point leans σ_WP.
(`linres_residual_summary.json`, `linres_forkA_collapse.png`.)

---

## Held-for-launch (recorded, not executed)

The panel's proposed decisive new run (spec'd, **not launched** — expensive-sim
gate, user owns launches): a matched classical+WP pair at `v ≈ 1–2`,
`L ≥ 100 Bohr`, `≥ 3` plasma periods (≈150 a.u.), writing both `_wf` bath and
`density_wp`. This is recorded as a next step, with no result.
