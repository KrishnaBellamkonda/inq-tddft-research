# 02 — Tasks performed (what was done)

Neutral, factual record of every task in the `ml-patterns` campaign, in
chronological order across four eras. Each entry gives the **input data**, the
**method/algorithm** (see `01_algorithm_papers.md`), the **code path**, and the
**output artifacts**. Verdicts and numbers are in `03_results.md`. No
interpretation is offered here — that is left to the reader.

Common context (all eras):
- **System:** pure bulk jellium, homogeneous electron gas, `r_s ≈ 5.69`, box
  `L = 50` Bohr, grid `125³` (`dx ≈ 0.4`), plasma frequency `ω_p ≈ 3.5 eV`
  (plasma period ≈ 49 a.u.). rt-TDDFT / Ehrenfest runs, produced by INQ+inqkit
  elsewhere; this campaign only post-processes their field output.
- **Two projectiles being contrasted:** a **classical** point/Gaussian charge
  and a **quantum wavepacket (WP)** of width `σ_WP`.
- **Induced density (the object of study):** the bath response
  `n_bath = n_total − n_wp`, GS-subtracted to give the induced field
  `δn(r,t) = n_bath(r,t) − n_bath^GS`. Fields loaded via `inqview.load_vti`
  (physical order; never `fftshift`).
- **All post-processing is CPU, no new INQ runs**, campaign-local under
  `docs/campaigns/ml-patterns/`.

---

## Era 0 — Data inventory (T0)

- **T0 — run database.** Built a filterable inventory of every available TDDFT
  run (581 runs × 137 columns) and independently validated it.
  - Code: `build_run_database.py`, `validate_run_database.py`
  - Output: run DB CSV/JSON + `run_database_validation.md`

## Era 1 — Original signature gates (T1–T7, 2026-07-01)

Autonomous orchestrator (`orchestrate.py`) ran T1→T7; each phase idempotent, wrote
`artifacts/T*_result.json` + a notebook, emailed a 4-part summary.

- **T1 — pre-gate the kernels.** Built and validated the numerical kernels before
  any science use: POD (`pod.py`), DMD (`dmd.py`), form factor (`formfactor.py`),
  VTI normaliser/subtraction ladder (`normaliser.py`), cell resolver (`celldb.py`),
  pipeline (`pipeline.py`).
  - Validation: 3 independent `formula-validation` agents CONFIRM (POD, DMD,
    form-factor); 12/12 known-case code-tests pass (`tests/test_kernels.py`).
  - Sub-task: computed `F_ONCV(q)` from the actual `electron-ONCV-1.2.upf` local
    potential (Coulomb-tail-subtracted radial FT) to establish the q-window where
    the point projectile ≈ unity form factor.
  - Output: `artifacts/foncv.json`, `artifacts/T1_foncv.png`.

- **T2 — form-factor cut (bulk).** Metric: the q-space ratio of induced bath
  densities `R(q) = |δn_WP(q)| / |δn_classical(q)|` compared to the known
  `F_WP(q)/F_ONCV(q)`. Fixed `E = 100 eV`, σ-sweep. Pinned split (ADR 0011):
  calibration `σ_WP ∈ {1, 5}`, held-out `σ_WP ∈ {0.5, 3, 8}`. Held-out rests on
  clean `_wf` bath-only runs.
  - Algorithm: form factor + radial q-spectrum. Code: `pipeline.py` (T2 path).
  - Output: `artifacts/T2_result.json`, `T2_heldout.png`, `T2_pod_delta.png`.

- **T3 — wake gate (bulk).** Windowed DMD on the `σ_WP = 5` velocity sweep;
  dominant DMD frequency vs `ω_p`, wavelength vs `λ = 2πv/ω_p`. Split: calib =
  even-velocity-index E, held-out = odd.
  - Algorithm: exact windowed DMD. Code: `pipeline.py` (T3 path).
  - Output: `artifacts/T3_result.json`, `T3_heldout.png`.

- **T4 — localised slab (Rung 1b).** Transferred the frozen pipeline to the
  `σ_WP = 0.5` `sigma_matched_gauss` slab geometry.
  - Output: notebook `rung1b_slab.ipynb` (verdict in `03_results.md`).

- **T5 — dynamics (Rung 2).** DMD/Koopman + SINDy on the bulk induced density;
  latent-mode ODE (2-mode) fitted in POD coordinates; mode spectrum vs Bohm-Gross.
  - Output: `artifacts/T5_result.json`, `T5_dynamics.png`, `rung2_dynamics.ipynb`.

- **T6 — exploratory exchange/diffraction.** On the vacuum-WP-subtracted field
  (labelled exploratory only, carries the externally-unverified SIE ≈ 7 eV caveat).
  - Output: `artifacts/T6_exploratory.png`, `exploratory_exchange_diffraction.ipynb`.

- **T7 — synthesis.** Cross-task synthesis notebook + verdict roll-up.
  - Output: `artifacts/T7_result.json`, `T7_synthesis.png`, `synthesis.ipynb`.

## Era 2 — Governing-PDE discovery redo (T8–T14, 2026-07-03/04)

Bulk-only redo (designed via grill-with-docs; ADR 0012). Two-track: **Track A**
re-runs the T2/T3 gates clean on bulk; **Track B** discovers a governing field PDE
via PDE-FIND. Ran through `orchestrate_pde.py` (idempotent, 12 h wall-cap,
per-phase Gmail).

- **T8 — scope + cell pin.** Fixed the bulk-jellium cells from the validated run
  DB (form-factor σ-sweep + classical/WP velocity sweeps). Pinned Track-B split:
  calibration `E ∈ {20, 50, 300}`, held-out `E ∈ {25, 100, 600}`.

- **T9 — pre-gate the PDE-FIND kernel.** Built `kernels/pdefind.py` (STRidge) +
  `kernels/discovery.py` (cell → axial field `n(z,t)` → PDE). Validation: recover
  KNOWN PDEs (heat/advection/wave) from synthetic data.
  - Validation: 6/6 tests pass (`tests/test_pdefind.py`); catalogue rows added.

- **T10 — Track-A gates (bulk).** Re-ran form-factor + wake held-out verdicts clean
  on bulk-only.

- **T11 — discover `PDE_classical`.** PDE-FIND on the classical velocity-sweep
  axial field `n(z,t)`; broad agnostic library, target order m = 2; three
  validation walls (calib/held-out cell split; forward-integration prediction;
  bootstrap stability).
  - Output: `artifacts/PDE_T11_result.json`.

- **T12 — discover `PDE_WP`.** Same procedure on the WP velocity-sweep
  (`σ_WP = 5` primary). Discovered on the run's `density_system`.
  - Output: `artifacts/PDE_T12_result.json`, `T12_heldout.png`.

- **T13 — compare + latent ODE + interpret.** Cross-compared the classical vs WP
  discovered terms; latent-ODE cross-check; post-hoc physics naming.
  - Output: `artifacts/PDE_T13_result.json`, `T13_compare.png`.

- **T14 — synthesis judge.** Aggregated the two tracks.
  - Output: `artifacts/PDE_T14_result.json`.

- **Scientific-panel review + re-analysis (2026-07-04).** A 4-expert panel
  reviewed the T12 WP result; a subsequent re-analysis re-measured the input data
  provenance. Recorded facts:
  - velocity-sweep run `run_wp_n162_L50_E100` `density_system` = **163 e**
    (WP-included); the matched `_wf` run `total − wp` = **162 e** (bath-only).
  - The T11/T12 PDE discovery used the axial reduction of these fields.
  These measurements are reported in `03_results.md` without a verdict.

## Era 3 — POD/DMD bath-structure sweep (2026-07-04)

Applied POD + DMD to the **blob-free bath** (`n_total − n_wp − GS` for WP runs;
`n_total − GS` for classical) across a σ-sweep and a velocity sweep, tabulating
POD rank(90 %), leading-mode energy fraction, and dominant DMD frequency.

- Code: `bath_structure_sweep.py`
- Runs used (σ-sweep, fixed `v = 2.71`): classical point (`E100_v2`) + WP
  `σ ∈ {0.5, 3, 8}` (`_wf` bath-only). Velocity sweep (classical): `E ∈ {20, 25,
  50, 100, 600}`. WP velocity points (mixed σ/grid): `E1p5` (v=0.33, σ5),
  `varyv` (v=0.5, σ3), `E100 σ3` (v=2.71).
- Output: `artifacts/bathstruct_{sigma,velocity}_sweep.png`,
  `bathstruct_summary.json`, per-run `bathstruct_run_*.json`,
  `bath_pod_dmd_compare.png`.

## Era 4 — Linear-response residual test (2026-07-06)

Panel-chosen time-domain technique. Null: since `χ` is a medium property (same for
both projectiles), the WP induced density is a Gaussian low-pass filter of the
classical one, so the frame-by-frame q-space ratio `|R(q,t)| = |n_WP(q,t)| /
|n_cl(q,t)|` should equal `F(q) = exp(−q²σ²/2)` and be flat in time. Two
discriminants (both d'Alembert-safe): `|R(q)|` shape vs `F(q)`, and `|R(q,t)|`
temporal flatness. Fork A resolves the √2 trap empirically by fitting the exponent
`a` in `|R(q)| ~ exp(−a q²)` across σ.

- Code: `kernels/formfactor_residual.py` + runner `linres_residual_test.py`
- Validation: 7/7 known-case tests pass (`tests/test_formfactor_residual.py`);
  catalogue rows added.
- Matched classical+WP pairs analysed at `E = 100`, `σ ∈ {0.5, 3, 8}`, `maxf = 100`.
- Plan: `docs/plans/linres-residual-classical-vs-wp.md`.
- Output: `artifacts/linres_residual_summary.json`, per-σ
  `linres_residual_sigma*.json`, figures `linres_residual_sigma*.png`,
  `linres_forkA_collapse.png`.
