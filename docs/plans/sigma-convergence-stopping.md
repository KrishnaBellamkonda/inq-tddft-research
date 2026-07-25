# Plan: σ-convergence of rt-TDDFT stopping power → point-charge Lindhard

**Branch:** `convergence-gaussian-electron-previous-scheme`
**Created:** 2026-06-14 (grill-with-docs session)
**Status:** design locked; awaiting go-ahead to build + launch.

---

## 1. Goal & hypothesis

Test whether the rt-TDDFT electronic stopping power S(v) of a classical
erf-smoothed Gaussian electron in r_s=5.69 jellium **converges toward the single
analytical point-charge Lindhard reference as the projectile width σ shrinks**.

The σ=0.5 runs sit *below* the point-charge curve near the stopping peak
(sim/LR ≈ 0.59 at v≈0.6) but *on* it at high v (≈1.01 at v≈3). A finite width
suppresses high-q coupling (form factor e^{−q²σ²}); shrinking σ should lift the
near-peak points toward the reference. **Hypothesis:** S(v;σ) → S_LR^point(v) as
σ→0, monotonically in σ.

## 2. The single reference (DONE)

- **One** analytical curve everywhere: point-charge (σ→0) Lindhard,
  `inqview.analysis.lindhard_elf.stopping_power_point(v, kF)`, infinite system,
  **no finite-box correction**. Never a per-σ form-factor family; never σ=0.2.
- Verified converged 0.00 % (qmax margin 2→8, n_q 4k→16k); f-sum <1e-3.
  Natural kinematic cutoff `qmax=2(kF+v)+margin` (no 1/σ blow-up).
- Reference values (Ha/Bohr): v=0.20→0.0144, 0.44→0.0392, 0.62→0.0603 (peak),
  1.00→0.0416, 1.94→0.0170, 2.98→0.0088.
- Tests: `test_lindhard_elf.py::test_stopping_point_converged`,
  `::test_stopping_point_above_finite_sigma` (16/16 pass).
- All plots updated to this reference (notebook + analyse_sv.py +
  make_sv_comparison.py); Barkas language removed everywhere.

## 3. Sanity checks settled (DONE)

1. **Form-factor convention — PASS.** The sim UPF on disk is exactly
   `C·erf(r/(σ√2))/r` (rel err 5e-11); σ is the charge-cloud RMS width, and the
   Lindhard form factor e^{−q²σ²} pairs with V(q)∝e^{−q²σ²/2}. No stray √2.
2. **q-grid — fixed for the point reference** via `stopping_power_point`
   (the legacy `stopping_power_sigma` under-resolves at small σ; not used for
   the reference).

## 4. Locked design decisions

| Decision | Value | Rationale |
|---|---|---|
| New σ ladder | **{0.15, 0.25, 0.35}** | three new; reuse σ=0.5 as soft anchor. UPFs generated + verified. |
| Reference | single **point-charge** Lindhard | §2; target is a point particle, not σ=0.2. |
| Velocities (per σ) | **v₀ ∈ {3.0, 2.0, 1.3, 0.8, 0.6, 0.2}** | match the σ=0.5 set for direct comparability + new low-v 0.2 anchor. |
| Projectile mass | **m_e for ALL runs, no exceptions** | direct comparability (user, firm). v=0.2 free-decelerates; window-mean-v point in the low-v friction band — handled by Method A. |
| Bath / grid | L=50 cubic, N=162, dx=0.40, LDA, GS `checkpoints/gs_L50_cubic_N162_dx0p40` | identical setup across all runs; reuse validated GS. |
| Observable set | full **ADR-0006 jellium_classical** minimum set + COD + manifest | "thorough analysability" + derived observables. |
| Density VTI cadence | **6 frames / run** (total/system/delta) | user: don't overcrowd memory. ~36 MB/run. |
| Cheap scalars cadence | **every step** (energy/current/dipole/density_l2/COD) | dense Method-A sampling; cheap. |
| Orbital overlap | **t=0 and t=end only** | user. |
| Momentum distribution | N/A (WP-only observable) | classical projectile. |
| Run layout | **3 run dirs, one per σ**, each `results/<vtag>/raw/...` | user (ADR 0007). |
| Analysis home | `systems/jellium/hypotheses/06_sigma_convergence/` | user (ADR 0007); executed .ipynb + README + figures. |
| Schedule | 18 runs, 2 GPUs, **~2.5 nights** | user accepted full 6-velocity sweep. |

### Consequence accepted
6 density frames ⇒ **no spectral loss-function / density_fourier** (needs a dense
time series). S(v) is unaffected (energy + track only). A cheap per-step axial
n_q(t) CSV writer could restore the loss function later if requested — NOT in scope now.

## 5. Pseudopotentials (DONE)

`shared/pseudopotentials/electron_gaussian_sigma0p{15,25,35}.upf` generated via
`inqview.io.gaussian_psp.generate_gaussian_psp` from `electron-ONCV-1.2.upf`
(Rydberg template, C=2). Verified to `C·erf(r/(σ√2))/r` at 5e-11; V(0)=√(2/π)/σ·C
= 10.64 / 6.38 / 4.56 Ha. Existing σ=0.4, 0.5 reused.

## 6. Run machinery (TO BUILD)

Per-σ dirs: `run_classical_n162_L50_sv_sigma0p15`, `_sigma0p25`, `_sigma0p35`.

- **run.cpp**: derive from the FULL classical template
  `run_classical_n162_L50_E100/run.cpp` (writes the ADR-0006 minimum set +
  manifest + COD + overlap). **Modify** to decouple cadences:
  - `VTI_EVERY` = N_STEPS/5 (→ 6 density frames) for density_total/system/delta.
  - scalar `observables.csv` every step.
  - overlap_full snapshots at step 0 and N_STEPS only (drop mid).
  - emit `observables_manifest.json` (RunType::jellium_classical).
  - σ via `SV_PSEUDO`, v₀ via `PROJ_V0`, N_STEPS via `SV_N_STEPS`, out subdir via
    `SV_OUT_SUBDIR` — one build per σ dir serves its 6-velocity ladder.
- **config header**: `shared/configs/sv_ladder_L50_sigma_sweep.hpp` (or extend
  `sv_ladder_L50_sigma0p5.hpp`) with the three σ structs (psp path only differs).
- **dispatcher**: per-σ, loops the 6 velocities on a chosen GPU (reuse
  `dispatch_ladder.sh` pattern: `v0:nsteps:subdir`).

### N_STEPS per velocity (dt=0.020)
| v₀ | 3.0 | 2.0 | 1.3 | 0.8 | 0.6 | 0.2 |
|---|---|---|---|---|---|---|
| N_STEPS | 300 | 450 | 700 | 700 | 700 | 1000 |
| VTI_EVERY (6 frames) | 60 | 90 | 140 | 140 | 140 | 200 |

(v=3..0.6 match the σ=0.5 ladder exactly; v=0.2 sized for ~3–4 Bohr decelerating
path. All stay inside the box from launch z=−20 — no wrap.)

## 7. S(v) extraction (per `docs/handovers/stopping-power-measurement.md`)

- **Method A (primary):** ΔE_sys(t)=E_total(t)−E_total(t₀); discard first 20 % of
  sim time; linregress ΔE_sys vs path s; slope = S (positive), stderr = uncertainty.
- **Method B (cross-check):** projectile KE loss / speed-integrated path over the
  same window.
- **Confidence gate:** Method A vs B within **10 %** is a good sign (current σ=0.5
  set achieves 1–7 %).
- The per-run `analyse.py` `stopping` + `bath_energy` phases implement these;
  the study notebook re-derives them uniformly across all 18 runs.

## 8. Deliverables

- Per-(σ,v) `analysis/REPORT.md` + full classical pipeline figures (minus
  density_fourier).
- `hypotheses/06_sigma_convergence/`:
  - executed `.ipynb` assembling all 18 S(v) points (+ existing σ=0.5/0.4) onto
    the single point-charge Lindhard reference; one panel per σ or one overlay
    with σ-coloured series + the reference.
  - method A/B cross-check parity plot.
  - `README.md` + figures (canonical theme).

## 9. Validation status

- DONE: `stopping_power_point` + tests (16/16); UPF verification; reference
  convergence/f-sum; figures regenerated.
- PENDING: smoke test of the new full-observable run.cpp at one (σ,v) before the
  full launch (verify manifest, 6 VTI frames, energy drift <1 mHa, Method A/B
  agreement). Per validation-gates: launch only after the smoke test passes.
- Record test-catalogue row for `stopping_power_point`.

## 10. Risks / open items

- **Grid floor:** dx=0.40 ⇒ g_Nyquist=7.85; the grid cannot resolve σ<~0.15, so
  σ=0.15 is borderline (form factor 0.86 at Nyquist) — it behaves like the
  grid's near-point charge. This is the *intended* small-σ extreme; the sim S at
  σ=0.15 should land closest to the point reference. (σ=0.05 was rejected:
  unrepresentable on this grid.)
- v=0.2 m_e free-deceleration ⇒ representative point at window-mean v (~0.12–0.15)
  with a visible velocity band; this is physical, not a defect.
- ~2.5-night compute on 2 GPUs; the user's idle `weather-climate` Jupyter kernel
  (PID 1816189) holds a GPU context — consider killing to free memory.

## 11. Exact next steps

1. Build config header + the three per-σ run dirs from the full classical
   template with decoupled VTI cadence + manifest.
2. `inq-run` build once per σ dir (GPU).
3. **Smoke test** one (σ=0.25, v=1.0) short run; validate manifest + 6 frames +
   energy drift + Method A/B.
4. On pass: launch the 6-velocity ladders across 2 GPUs (per-σ dispatchers),
   ~2.5 nights.
5. Per-run `analyse.py`; then the `hypotheses/06_sigma_convergence/` study notebook.

---

## 12. Extension 2026-06-15 — large-σ probe (σ=3) + plot revisions

The σ∈{0.15,0.25,0.35,0.5} sweep finished (2026-06-15 07:46Z): S(v) is essentially
σ-independent at low/mid v and only fans out at high v — the near-peak shortfall vs
point-charge Lindhard does **not** close as σ→0. Next objective (user): probe the
*other* direction — a σ large enough to deviate strongly — and refine the peak.

### 12a. Plot revisions (grill-with-docs 2026-06-15, applied to `sigma_sweep_report.py`)
- **x-position = nominal launch v₀** (not the window-mean). Decided despite the
  large low-v₀ deceleration (v₀=0.8→window-mean 0.61; v₀=0.2→0.09) — the point is
  a mild upper bound on the velocity at which S was sampled; noted in the caption.
- **Error bars: vertical only** = linregress stderr on the slope S. Remove the
  horizontal (vlo→vhi) bars entirely.
- **★ peak marker → vertical dashed line** annotated "Lindhard peak" (at the
  point-charge peak velocity). The existing dashed line annotated **"k_F"**.
- **Companion energy figure** `sv_convergence_energy.png`: x-axis in **eV**
  (log scale), E₀ = 0.5·v₀²·27.211; Lindhard reference + both dashed lines mapped
  to energy. Velocity figure `sv_convergence.png` retained as primary.

### 12b. σ=3.0 Bohr run set
- **Why σ=3:** form factor e^{−q²σ²} suppresses all q≳0.3 Bohr⁻¹ ⇒ expected large
  downward deviation from point-charge Lindhard. V(0)·σ≈1.596 Ha ⇒ V(0)≈0.532 Ha
  (weak, diffuse, repulsive). Grid floor irrelevant (very smooth; dx=0.40 fine).
- **Geometry:** STANDARD boundary rule (σ=3 feasible in L=50):
  launch_z = −25+4·3 = **−13**, stop_z = +22, traversal 35 Bohr. NOT the existing
  binary (launch −20 would push the σ=3 trailing 4σ tail 9 Bohr past the −25 face)
  ⇒ new Cfg `SV_Ladder_L50_sigma3p0` (launch −13) + new run dir + own build.
- **UPF:** `shared/pseudopotentials/electron_gaussian_sigma3p0.upf`, generated like
  the others (`inqview.io.gaussian_psp`, V(r)=C·erf(r/(σ√2))/r, C=2 Ry).
- **Ladder (7 runs):** v₀ ∈ {3.0, 2.0, 1.3, **1.0 (new, peak refinement)**, 0.8,
  0.6, 0.2}. N_STEPS reused from the existing map {3.0:300, 2.0:450, 1.3:700,
  0.8:700, 0.6:700, 0.2:1000}; v₀=1.0→700 (14 a.u. window, same as neighbours).
  Max centroid reach (v₀=3, const v) = −13+18 = +5 ≪ stop +22 — all within bounds.
- **Cadences:** identical to the sweep (6 density VTI frames; dense scalars; COD +
  L2; overlap t=0/end). One shared build serves the 7 velocities (PROJ_V0 /
  SV_N_STEPS / SV_OUT_ROOT env).
- **Plot:** σ=3 added to `sigma_sweep_report.py` SIGMAS as a 5th series (nested
  layout) on both the velocity and energy figures.
- **Run:** launch now across GPU0/GPU1 (both free); email both figures + the
  Method A/B + sim/LR table to chiddukanna@gmail.com on completion.

### 12c. Steps
1. Plot revisions in `sigma_sweep_report.py`; regenerate both figures on existing
   23 points; verify.
2. Generate σ=3 UPF; verify analytic form.
3. New Cfg `SV_Ladder_L50_sigma3p0` (launch −13, σ=3 UPF).
4. New run dir `run_classical_n162_L50_sv_sigma3p0/` + run.cpp + `inq-run` build.
5. Smoke test (σ=3, v=1.0, ~30 steps): manifest + 6 frames + Method A/B parser.
6. Orchestrate the 7 velocities across 2 GPUs; email both figures + table at end.
