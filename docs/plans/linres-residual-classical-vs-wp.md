# Plan: linear-response residual test — classical vs WP induced density

Status: **executed** (2026-07-06) — Tier-1 (existing-data) done; Tier-2 new run
held for user GO. Supersedes the governing-PDE line (dead: blob artifact) and
extends the POD/DMD bath-structure contrast. Chosen by the scientific panel
(opus×9, 2026-07-06) over Floquet/Koopman, HAVOK, optimal transport,
wavelet+transfer-entropy.

## RESULT (2026-07-06)

- σ=0.5 (only SNR-adequate pair): |R(q)| ~ Gaussian of width σ_fit≈0.62 → **near
  σ_WP=0.5, not σ_pot=0.35**; high-q excess ≈0.15σ. But |R(q,t)| NOT flat over 4.8
  a.u. (t_flatness=0.43) → static-linear-filter null REJECTED. Early-vs-late split
  (both halves flat ≈0.19, full 0.43) → monotone drift = **DECELERATION** (light WP
  slows; classical ghost holds v).
- σ=3, σ=8: SNR-DEAD (form factor e-folds in ~1–4 shells; fit hits the total−wp
  blob floor) → excluded (gate a≥0.15σ²). **Fork A INCONCLUSIVE** on existing data
  (1 usable σ), leans σ_WP.
- Answer: instantaneously the WP induced density = classical through a Gaussian
  low-pass F(q)≈exp(−q²σ_WP²/2); dynamically they diverge via WP deceleration.
  Artifacts: `artifacts/linres_residual_summary.json` + 4 PNGs. 7/7 kernel tests.

## Question

How, specifically, does the induced electron-bath density of a **classical point
charge** differ from that of a **quantum Gaussian wavepacket** in matched-velocity
bulk jellium — deeper than POD/DMD's "coherent vs incoherent", and immune to the
d'Alembert kinematic trap?

## Reduced model (the null the test challenges)

Linear response: `n_ind(q,ω) = χ(q,ω)·V_ext(q,ω)`. χ is a property of the medium
(the HEG) → **identical for both projectiles**. Only the drive differs:

- point charge: `V_ext(q) = 4π/q²`
- Gaussian WP:  `V_ext(q) = 4π/q² · exp(−q²σ²/2)`

So the WP is a **low-pass-filtered point charge**, `F(q) = exp(−q²σ²/2)`. In
linear response, frame-by-frame:  `n_WP(q,t) = F(q)·n_cl(q,t)`. The ratio cancels
χ **in the time domain** — no ω-binning, which is essential because the matched
runs are only T ≈ 4.8–17 a.u. (0.1–0.5 plasma periods) → Δω ≈ 10–36 eV ≫ ω_p =
3.5 eV, so **no frequency-resolved technique is extractable** on existing data.

Two discriminants make it deeper than POD/DMD and d'Alembert-safe:
1. Dividing by V_ext removes any rigid `f(z−vt)` → translation kinematics cancel
   (in magnitude the rigid phase cancels exactly).
2. `|R(q,t)|` must be **flat in t**. Any t-drift ⇒ the equal-trajectory/linear
   premise breaks (differential deceleration, WP spreading, nonlinearity),
   localized in q.

## Outcomes (pre-registered)

- **t-flat collapse onto F(q)** ⇒ the whole classical↔WP difference is one linear
  filter (clean null; POD/DMD demoted to a q-filter corollary). Still answers the
  campaign decision.
- **t-flat high-q excess at σ=0.5** ⇒ genuine nonlinear / quantum-projectile
  fingerprint (near-field nonlinear screening of the point charge at r_s ~ 5–6,
  or WP momentum spread).
- **t-drift** ⇒ deceleration mismatch is the story (re-test in the early v≈v₀
  window per the light-projectile rule).

## Fork A resolved empirically (the √2 trap)

Do NOT hardcode F(q). Fit the Gaussian exponent `a` in `|R(q)| ~ exp(−a q²)` from
data, then test `a(σ)` across σ=0.5/3/8: slope 0.5 in σ² ⇒ physical width is
`σ_WP`; slope 0.25 ⇒ `σ_pot = σ_WP/√2`. The data selects the width. (Panel Fork A.)

## Data (verified on disk, all L=50, 125³, dq=0.13 a.u.)

| pair (σ) | classical | WP (_wf, blob-free) | overlap T | binding |
|---|---|---|---|---|
| 0.5 | run_classical_n162_L50_E100_v2 (16.6 a.u., 431f) | run_wp_..._sigma0p5_wf (4.8 a.u., 241f) | 4.8 a.u. | **primary** |
| 3   | same classical | run_wp_..._sigma3_wf (12.9 a.u., 324f) | 12.9 a.u. | cross-check |
| 8   | same classical | run_wp_..._sigma8_wf (6.6 a.u., 333f) | 6.6 a.u. | SNR-dead (collapse only) |

bath: classical = density_total − GS; WP = density_total − density_wp − GS.
EXCLUDE the dt=4.0 a.u. varyv run (aliased). All at v=2.71 a.u.

## Method (two lenses, existing data, CPU, no new run)

1. **3-D radial (headline).** Per frame, GS-subtract, `rfftn`, bin |q| shells →
   `n(|q|,t)` = shell-mean magnitude; `noise(|q|)` = shell-std/√N_shell (real SNR,
   the reason 3-D beats 1-D). Isotropic |q| is the physically correct form-factor
   axis. Stream frame-by-frame (never hold the whole (T,125³) stack).
2. **1-D axial (cross-check, settles panel Q#4).** transverse-mean → u(z,t) →
   rfft_z → `n(q_z,t)` on the q_⊥=0 line. Cheap; F(q_z)=exp(−q_z²σ²/2) still holds.
3. Restrict to `[0, T_overlap]`; resample both onto a common time grid (linear,
   per-q) using REAL frame times (`t = linspace(0, total_time_au, n_frames)`).
4. `|R(q,t)| = |n_WP|/|n_cl|` where `|n_cl| > 3·noise`. Report:
   δ(q,t)/noise(q) (normalized residual), σ_fit vs σ_WP/σ_pot, t-flatness
   (std_t/mean_t per q), high-q excess after best-fit F.

## Deliverables

- `kernels/formfactor_residual.py` (pure numpy) + `tests/test_formfactor_residual.py`
  (synthetic known-case: exact-F recovery, injected high-q excess, injected t-drift).
- `linres_residual_test.py` runner: idempotent per-pair JSON, figures, 4-part
  email, per-pair try/except, autonomous. CPU-only.
- Catalogue rows in `docs/validation/test-catalogue.md`; handover milestone.

## Held for user go (expensive-sim gate — user owns launches)

The panel's single *decisive* new run: matched classical+WP pair at **v≈1–2**,
**L≥100 Bohr** (dq≈0.063, puts q*=ω_p/v on-grid), **≥3 plasma periods ≈150 a.u.**
(Δω≈1.1 eV), writing both `_wf` bath and `density_wp`. Alone answers
plasmon-vs-Doppler (is the WP's ~7–11 eV mode real plasmon coupling at 3.6 eV, or
the kinematic co-moving Doppler line ω=q_min·v≈9.3 eV). Spec'd, NOT launched.
```
```
