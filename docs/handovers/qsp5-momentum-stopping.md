# Handover: momentum-KE-loss stopping, orbital-free (qsp5 sweep + σ=1 plateau pair)

---

## Update: 2026-07-28 (later) — methods #3 and #4: Ehrenfest drag + momentum-space notebooks (v1p3)

Status: done. Two new method notebooks on the clean run, both executed 0 errors:
- `hypotheses/qsp_phase5/p5_wp_v1p3_ehrenfest_drag.ipynb` (builder
  `build_v1p3_ehrenfest_drag_report.py`, summary `ehrenfest_drag_summary.json`).
  Method #3, LOCAL: S = −d(⟨p_z⟩²/2)/dz_c from per-step orbital moments + COD.
  KEY PHYSICS: electron packet is ACCELERATED into jellium (image attraction,
  pz 1.30→1.45 by t≈13) → the light-projectile "v≥0.85·v0 early window" rule
  FAILS for attractive projectiles (returns ~0.03); headline window = POST-PEAK
  (v≥0.85·v_peak): S(v0) = 1.2 eV/Bohr at v̄=1.37. Full-deceleration interior
  slope = 3.0 ± 0.55 at v̄=0.47 (whole sweep; not a v0 number). Ehrenfest
  identity v_c=⟨p_z⟩ breaks at t≈5 (CAP non-unitarity) — selection pushes pz up
  → drag values are lower bounds. Rule-update candidate for
  light-projectile-stopping.md (attractive-projectile post-peak window).
- `hypotheses/qsp_phase5/p5_wp_v1p3_momentum_space.ipynb` (builder
  `build_v1p3_momentum_space_report.py`, summary `momentum_space_summary.json`).
  Method #4, ASYMPTOTIC-spectral: signed P(k_z) marginal from complex
  wavefunction_wp VTIs (FFT_z + transverse Parseval); launch gate exact
  (⟨k_z⟩=1.2975 = recorded pz_mean to the digit; Var=2.003); transmitted
  ensemble = corridor-windowed envelope over t* scan (coverage 0.75);
  rank-matched S = 0.38 ± 0.27 at u≈1.67 — AGREES with TOF (0.30 ± 0.42):
  spectral-ψ and continuity-n machineries cross-validate.
Method disagreement to resolve: local drag (~1.2 at v̄1.37) vs asymptotic
methods (~0.3-0.4 at u~1.7) — candidate explanations: elastic/reflected momentum
return, CAP-selection bias, u-dependence (S falls with u per TOF band), exit
deceleration outside drag window. Discuss with user.

## Update: 2026-07-28 — combined low-velocity S(v) added to the qsp5 notebook

Status: done. New section "Combined low-velocity S(v)" in
`qsp_phase5_momentum_stopping.ipynb` (42 cells, 0 errors): flux-weighted average
of the grade-A/B S(u) bands only (v1p3 + v3p0; aliased v4p0–v6p0 excluded),
weight w_r(u) = −dN_in/du. Output `momentum_stopping_lowv.json`:
S(1.5)=0.25, S(2.0)=0.22, S(2.5)=1.4, S(3.0)=1.5, S(3.5)=0.82 eV/Bohr
(± 0.41 syst). Caveat: the v1p3/v3p0 bands barely overlap — the u≈2.0→2.5 step
is a band seam (each run's band edge abuts its excluded fast-tail zone), not a
physical threshold; same for the S(3.5) droop (v3p0 upper edge). Summary json
shape kept back-compatible (per-run top level; low-v in separate file).

## Update: 2026-07-27 (evening) — third notebook: snapshot kinematics of the surviving projectile (v1p3)

Status: done.
`/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/p5_wp_v1p3_snapshot_kinematics.ipynb`
(21 cells, 0 errors; builder `build_v1p3_snapshot_report.py`; summary
`snapshot_kinematics_summary.json`). User-driven snapshot picture: projectile ≡
excess density in the vacuum corridors at t* ∈ {30,40,50,60}; Madelung split
T_full = T_W[n] (exact 3D) + T_v (flow; longitudinal from continuity-J, transverse
free-dispersion estimate). Results: S_snap = 0.16/0.53/0.53/−0.95(junk, coverage
0.04) eV/Bohr across the scan, bracketing TOF S_drift = 0.30 ± 0.42, vs deposit
2.37. Baseline gate: N=1.000, <p_z>=1.297/1.3, T_drift=0.840/0.845 Ha. Key
numerical finding: T_W(0)=5.44 Ha vs analytic 3.00 is PURE discretisation
(σ_z=0.7·dx; analytic-Gaussian-on-grid control reproduces 5.442); from frame 1
cached T_W tracks the free-dispersion law to ≲3% (0.6402 vs 0.6402 at t=0.48).
T_v,z(0) unusable (one-sided ∂t on under-resolved profile) — replaced by the
flat-phase identity T_v(0)=T_drift(0). Next: user questions on the notebooks.

---

## Milestone: 2026-07-27 (later) — plateau-pair notebook built: term-by-term plateau dissection (σ=1 cap/nocap)

### Current status
DONE. Second notebook
`/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/wp_cap_energy_plateau/wp_cap_energy_plateau_momentum_stopping.ipynb`
(30 cells, executed, 0 errors) applies the validated TOF/rank-matching method to
the clean-launch σ=1, E=100 eV cap/nocap pair and — the centrepiece — dissects
the CAP energy plateau term by term using the pairwise `interactions.csv` (eV)
ledger, directly testing the naive estimator S from (E(0)−E_plateau)−T_loc.

### Key numbers (plateau_dissection_summary.json, 2 s.f.)
- Arrival ledger E(0)−E_GS = 124 eV = 100 (drift) + 20 (T_loc) + 12 (e_pp WP
  self-Hartree) + launch-interaction/self terms.
- CAP removed R = 86 eV; kept by slab D = 38 eV; N_wp(end) = 0.004 (fully drained).
- D ≈ T_loc(20.4) + e_pp(11.9) + ~5 eV remainder → the plateau "absorption" is
  mostly zero-point + self-energy, NOT stopping.
- S_naive = [E_drift0 − (R − T_loc0)]/L = 1.4 eV/Bohr; S_deposit = D/L = 1.5;
  both close to the σ=0.5 qsp5 S_drift(u≈3) = 1.3 ± 0.56.
- TOF S_drift here: −0.91 ± 1.9 (cap), −1.5 ± 1.8 (nocap) — cap/nocap agree, but
  **this geometry is resolution-limited**: the entrance flight is only 4 Bohr
  (launch −20.5, plane −16.5), and the geometry-specific null test prices the
  systematic at ±1.6 eV/Bohr. The TOF cannot discriminate here; the ledger can.
- nocap conservation control: |ΔE_total| = 0.0017 eV over 100 a.u. (numerics
  excellent); nocap TOF valid only to t = 19 a.u. (no CAP → edge wrap).

### Files added
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/wp_cap_energy_plateau/plateau_kinematics.py` — cache extractor (cap+nocap, 251 frames each)
- `.../hypotheses/wp_cap_energy_plateau/build_plateau_momentum_report.py` — two-pass builder
- `.../hypotheses/wp_cap_energy_plateau/wp_cap_energy_plateau_momentum_stopping.ipynb` — executed notebook (0.7 MB, GIFs path-referenced from results/<run>/report/)
- `.../hypotheses/wp_cap_energy_plateau/plateau_dissection_summary.json`, `cache/`, `cache_extract.log`

### Data facts verified for this pair (will bite again)
- `energies.csv` in Ha; `interactions.csv` in **eV** (e_ss+e_pp+e_ps = hartree×27.211 exactly).
- n_total(0) = n_gs + n_wp(0) to machine precision (∫ = 102 + 1 = 103).
- CAP: η = −0.7 Ha, region 60 < |z| < 70 (cap_mid_frac ±0.4643 · 140 = ±65).
- E_GS(slab) = −830.0242258 Ha (shared_gs run_summary).

### Next steps (this milestone)
1. User review: does the plateau ledger (D ≈ T_loc + self-Hartree + ~5 eV) settle
   the "plateau too high" question for the report?
2. If a sharper σ=1 TOF number is wanted: re-run with launch further back
   (z0 ≈ −45 in this 140-Bohr box) to lengthen the entrance flight — the only
   fix for the ±1.6 systematic.
3. Charged-cell caveat stands: pairwise differences are trend-level once the CAP
   changes the net charge; the total-energy plateau itself is convention-safe.

---

## Milestone: 2026-07-27 — S_drift(u) extracted orbital-free from qsp_phase5; notebook built & executed (0 errors)

### Current status
DONE. `qsp_phase5_momentum_stopping.ipynb` (40 cells, executed, 0 errors) delivers
a second, KS-orbital-free stopping estimate for the qsp_phase5 WP velocity sweep
(k0 = 1.3/3/4/5/6), using only density-level data. Method validated in-session by
a three-part battery (see below). Headline: the deposit-based S over-counts by
~2–8×; for the clean run (v1p3) S_drift = 0.30 ± 0.42 eV/Bohr vs S_deposit =
2.37 eV/Bohr — the gap is the WP localisation energy (41 eV/electron, verified
three independent ways) + capture binding, which projectile-KE stopping must
exclude. Conceptual chain that led here (drift vs localisation KE, TOF planes,
rank matching) lives in the 2026-07-27 conversation and is restated in the
notebook's method sections.

### Method (as finally locked — differs from first plan draft)
1. ρ(z,t) = ∫∫(n_total − n_gs)dxdy from `density_total` VTIs (NOT `density_delta`
   — that is n(t)−n(0); verified numerically). Cached per run.
2. J(z,t) via 1D continuity + CAP sink 2W(z)ρ, W = 0.7·sin²(π(|z|−35)/10) on
   35<|z|<45 (shape from `inq-study/src/perturbations/absorbing.hpp`; sink
   REGION verified empirically from band-resolved norm decay; closure test
   −dN/dt = A₊+A₋ passes at few-% level). **Side-adaptive integration**:
   entrance plane from left box edge, exit plane from right — short paths avoid
   accumulated gradient noise (right-integrated entrance was badly biased).
3. TOF detector planes z_in=−15.5 (left-int), z_out=+22 (right-int); u(t)=J/ρ.
   +15.5/+18/+20 exit planes are contaminated by bath-polarisation oscillation
   rectification (N_trans>N_in observed); ≥+22 converges.
4. **Exceedance-matched (rank-transport) S(u)**: incident vs transmitted
   N(>u) curves matched at equal exceedance q (top-down; robust to capture
   truncating the slow end); S(u_in(q)) = ½(u_in²−u_out²)/25 · 27.211 eV/Bohr.
   Naive per-electron ⟨K⟩_in−⟨K⟩_out is INVALID at low v (capture selects fast).
5. **Null test**: free analytic Gaussian pushed through the identical pipeline
   (grid dx=0.5, cadence 0.48 au, finite-diff, extended domain ±1100 Bohr so
   nothing escapes — the real runs' CAP guarantees J→0 at edges) must give S≡0.
   Trusted rank window q ∈ [0.30, 0.90]·q_top; systematic ±0.41 eV/Bohr.
   Top ~25% of ranks (fastest tail) excluded — entrance under-resolution.

### Key results (momentum_stopping_summary.json, 2 s.f.)
| run | u_ref | S_drift [eV/Bohr] | S_deposit | ratio | grade |
|---|---|---|---|---|---|
| v1p3 | 1.6 | 0.30 ± 0.42 | 2.4 | 7.9× | A |
| v3p0 | 3.0 | 1.3 ± 0.56 | 2.6 | 1.9× | B |
| v4p0 | 4.0 | 2.3 ± 1.7 | 4.5 | 2.0× | C |
| v5p0 | 4.8 | 4.1 ± 1.0 | 9.8 | 2.4× | F |
| v6p0 | 5.4 | 9.9 ± 2.3 | 19 | 1.9× | F |

QC-1 (forensic aliasing): dx=0.5 → k_Ny=6.28; σ_p=1.41 ⇒ measured ⟨p_z⟩(0)/k0 =
1.00/0.95/0.82/0.53/0.10 — v4p0+ launched corrupted (predates the mandatory
cutoff-aliasing guard). All S points stand at MEASURED u_ref, never nominal k0.
t=0 ledger closes three ways: orbital internal energy = lobe T_W = 3/(8σ_r²) =
1.5 Ha = 41 eV.

### Files touched
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/qsp5_momentum_kinematics.py` — cache extractor (one 3D pass/run; ~15 min parallel; caches in `cache/*_kinematics.npz` with z-profiles ρ, ρ_wp, lobe T_W/N/Z1)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/build_momentum_stopping_report.py` — two-pass nbformat builder (pass 1 executes → summary json; pass 2 injects computed takeaway numbers)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/qsp_phase5_momentum_stopping.ipynb` — executed notebook (1.3 MB; GIFs path-referenced from per-run notebook fig dirs)
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/momentum_stopping_summary.json` — headline numbers
- `/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/hypotheses/qsp_phase5/cache_extract.log` — extraction log
- `/local/data/public/skcb2/tddft/docs/plans/qsp5-momentum-stopping-notebook.md` — plan (all steps checked)

### Commands run
```bash
# cache extraction (all 5 runs, parallel, nohup-detached)
cd .../hypotheses/qsp_phase5 && PYTHONPATH=.../inq-stack/python nohup \
  .../venv/bin/python3 qsp5_momentum_kinematics.py > cache_extract.log 2>&1 &
# builder (two passes, executes notebook)
PYTHONPATH=.../inq-stack/python .../venv/bin/python3 build_momentum_stopping_report.py
```

### Tests and validation
- Run: (1) estimator exactness on analytic ρ,J free Gaussian — exact (58.5 vs
  58.2 eV truth, cadence-insensitive); (2) full-pipeline null test — |S| ≤ 0.41
  eV/Bohr in trusted window (also a cell in the notebook); (3) CAP-sink closure
  −dN/dt vs A₊+A₋ — few-% residual, all runs (notebook cell); (4) notebook
  executes 0 errors, all figure refs resolve.
- Remaining gaps: rank-matching ordering assumption unbounded (needs current
  frames); transverse channel invisible to z-planes.

### Trusted sources used
- Muga, Palao, Navarro, Egusquiza, Phys. Rep. 395, 357 (2004) — CAP removes
  density at local rate 2Wρ (sink term in continuity).
- von Weizsäcker kinetic functional T_W = ∫|∇n|²/8n (standard; Parr & Yang).

### Known issues / data caveats (verified, will bite again)
- `density_delta` VTIs are n(t)−n(0), NOT n−n_gs.
- `density_total` INCLUDES the WP here (∫=83) — config-dependent per memory.
- CAP mid/width are passed as FRACTIONS (40/90, 10/90) to absorbing.hpp which
  compares against Bohr coordinates — empirically the sink IS at |z|∈[35,45]
  as run_summary claims; do not "fix" without re-verifying against data.
- qsp_phase5 v4p0/v5p0/v6p0 launches are momentum-aliased (grades C/F).
- Bash-tool background tasks die with the session: extraction was nohup-detached.

### Assumptions still in play
- Rank matching assumes the slab preserves velocity ordering (no overtaking) and
  capture removes only the slow end. Plausible, unverified.
- Bath polarisation contribution at the detector planes is negligible beyond the
  chosen buffers (3 Bohr entrance / 9.5 Bohr exit).
- S normalised by L=25 Bohr (slab thickness); buffers are vacuum → no stopping.

### Exact next steps
1. User reviews the notebook (`qsp_phase5_momentum_stopping.ipynb`) — especially
   the S(u) collapse figure and whether S_drift vs S_deposit tells the intended
   localisation-energy story for the report.
2. If the method is accepted: apply it to the σ=1 cap/nocap pair
   (`wp_cap_energy_plateau`, 11005 VTIs, `interactions.csv` pairwise ledger) for
   the full plateau dissection (drift / localisation / capture attribution).
3. Optional hardening: re-run one velocity with saved current frames to close the
   rank-ordering and transverse-channel gaps.
4. Consider a clean high-v re-sweep at dx≤0.35 (or σ larger) to replace grades C/F.
