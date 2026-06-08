# RESUME NOTE — presentation evidence batch (2026-05-31)

**Tools (Bash + Read) are intermittently returning EMPTY or CORRUPTED output this
session.** Writes/edits succeed; stdout capture and file reads are unreliable.
**START A FRESH SESSION.** This note is self-contained so you can resume without
re-deriving anything. Read this first, then `mphil_midterm_presentation.md` milestone D.

---

## VERIFIED FACTS (do not re-derive)

### Density fields (run_wp_n162_L50_E20_sigma1_v2, by integration at t=0)
- `density_rt_system` = 163 e  (INCLUDES WP)
- `density_rt_total`  = 163 e  (IDENTICAL to system; total−system=0)
- `density_wp`        = 1 e
- `total − wp`        = 162 e  = N_electrons = **TRUE BATH**
- Code reason: run.cpp writes `sys_f = density::total(electrons)` which already
  includes the WP orbital (occ 1.0); `total_wr` writes the same `sys_f` (the
  `add_real_fields(sys_f, wp_f)` "total" line is NOT used for the saved total here —
  both system_wr and total_wr get sys_f). So both saved fields include the WP.
- **CLASSICAL runs:** no WP orbital in e-density → `density_rt_total`/`density_rt_delta`
  IS the pure bath wake; no subtraction needed.

### Which energy-sweep runs to use (σ=1, v2 — all have density_wp):
run_wp_n162_L50_E{20,25,50,100,200,300}_sigma1_v2  (from stopping_power_data.py
get_L50_wp_sigma1_runs). YES these are the 1-Bohr-sigma runs. ✓

### σ-sweep runs (E=100): σ=0.5/1/3/8 — DO NOT have density_wp saved.
Their run.cpp lacks the WF_WRITE_EVERY block. Their density_system = 163e (WP-incl).
User DECISION: **re-run the σ-sweep with density_wp + high-cadence wavefunction.**

### Loss-function plot uses SymLogNorm (spectral_weight.py:375 use_symlog=True)
→ compresses dynamic range → makes all (q,ω) channels look equally excited.
That is the cause of the user's "all channels equally excited?" concern.
Linear-scale script written: run_wp_n162_L50_E20_sigma1_v2/loss_function_linear.py
(reuses spectral_weight internals; fixed grid=125 from VTI). It RAN (3 PNGs written:
loss_function_linear.png, loss_function_1d_cuts.png, spectral_weight_response_linear.png)
but the 1D-cut diagnostic printed L≈0 everywhere — **UNVERIFIED whether physical
(weak response) or a clip/index bug. RE-CHECK in fresh session.**

---

## 2026-06-01 — #9 LINEAR-LOSS L≈0 RESOLVED (root cause = threshold bug, NOT symlog)

**Diagnosis (systematic-debugging, instrumented & confirmed):** the "L≈0 everywhere"
was a genuine BUG, not weak physics — L was *exactly* 0 (max=min=−0.000e+00).
Root cause: the deconvolution threshold `thr = 1e-3 * max|V_ext|` was poisoned by the
q=0 (DC) row. At qz=0, `q²` is clamped to 1e-10, so `V_ext = -(4π/q²)·δn_wp` is
inflated ~10³× at q=0 (measured: max|V_ext|=5.95e8 at q=0 vs 4.56e5 for all physical
modes). thr then exceeded every physical mode → 0/154 survivors → χ≡0 → L≡0.
Response WAS present (max|δn_resp_w|=2.1e2, max|W_resp|=4.5e4) — refuting "no response".

**Same bug existed in the PIPELINE** `inqview/postprocess/spectral_weight.py:303`, so
`loss_function_qz_omega.png` (the M8 figure) was ALSO L≡0 → a uniform map. **That is
the likely true cause of the user's "all channels equally excited" — NOT SymLogNorm.**
(Previous session's symlog attribution is superseded. User owns the verdict.)

**FIX (both scripts):** zero the unphysical q=0 V_ext row before the threshold max.
Physically grounded: L=-(4π/q²)Im[χ] is undefined at q=0 (zero momentum transfer; the
stopping integral ∫dq/q starts at qmin>0).
- `inq-stack/python/inqview/postprocess/spectral_weight.py` (~line 301-305)
- `ResearchProject/.../run_wp_n162_L50_E20_sigma1_v2/loss_function_linear.py` (~line 75)

**Known-case test PASS (mechanical):** survivors 0→154/154, max|L| 0→2.98, dyn range 9.8×.
3 PNGs regenerated in `run_wp_n162_L50_E20_sigma1_v2/results/analysis/observables/`.

**⚠️ NEW PHYSICS CAVEAT (the fixed L is noise-dominated — user must judge before M9):**
- 48.1% of L values are NEGATIVE (a passive-medium loss fn must be ≥0).
- 65% of Σ|L| weight sits at ω>8 eV, peaks pinned at the 20 eV window edge — NOT at the
  plasmon ω_p=3.47 eV. Classic deconvolution-noise amplification at high ω (V_ext decays
  → χ=δn_resp/V_ext blows up noise). **The raw extracted L is NOT a clean spectral
  density.** The cleaner observable is `W_resp = |δn_resp|²` (no deconvolution) —
  `spectral_weight_response_linear.png`.
- **M9 IMPACT:** ω-weighted ∫ωL dω stopping integral would be DOMINATED by this high-ω
  noise → loss-function stopping number unreliable without band-limiting / L≥0 clamp /
  stronger regularisation. **M9 path is now a user decision (see grill checkpoint).**

## 2026-06-01 — FUNDAMENTAL LIMIT: loss-function spectroscopy is time-window-starved

While "investing in cleaner extraction" (user choice) for M9, found the binding limit is
NOT the deconvolution algorithm but the **propagation time**:
- E20 v2 ran T = 11.40 a.u. = 0.276 fs. DFT freq resolution Δω = 2π/T = **15.0 eV**
  (independent of frame count / zero-padding).
- The 0–20 eV L(q,ω) window holds only **1.33 INDEPENDENT frequency bins**; ω_p=3.47 eV
  is at 0.23 bins; the run spans **0.23 of one plasmon period** (plasmon period 49.3 a.u.).
- ⇒ The loss-function map is ~1 real frequency point × 16× zero-pad interpolation. The
  apparent (q,ω) structure (incl. the η-recovered "plasmon at 3.74 eV") is interpolation
  of a single bin + deconvolution noise. No regularisation creates missing resolution.
- Even the 19.8 a.u. base-spec runs give Δω ≈ 8.6 eV > ω_p — still sub-plasmon-period.

**Prototype evidence** (`run_wp_n162_L50_E20_sigma1_v2/loss_function_clean_prototype.py`):
η-damping + Tikhonov/Wiener inversion (χ=δn_resp·conj(Vext)/(|Vext|²+λ²)) DID move the
small-q peak onto ω_p (η_factor=6, reg=0.01–0.03 → peak 3.74 eV) and cut %neg 50→32%,
but residual high-ω weight (~64%) and ~32% negative L persist — consistent with the
1-bin resolution limit + imperfect free-WP subtraction (WP self-repels, occ 1.0; the
run also saved actual density_wp at coarse 58-frame cadence — usable for a cleaner
subtraction but does not fix the resolution wall).

**Method note (grounded):** RT-TDDFT absorption/loss spectroscopy needs long propagation
(tens of fs) for sharp spectra; freq resolution = 2π/T (Yabana-Bertsch convention).
The WP runs were built for short single-pass transit, so they are resolution-starved for
spectroscopy. Stopping power does NOT require the freq-resolved loss fn — it is available
time-domain via bath energy absorbed / drift-energy / force (resolution-independent).

**OPEN — user decision (grill checkpoint, 2026-06-01):** M8/M9 loss-function route is
challenged. Options surfaced: (A) reframe M9 to time-domain stopping (bath_energy/drift,
no loss fn); (B) a dedicated long kick-response GS run for a real loss fn (new sim);
(C) present the resolution limit itself as the M8/M9 finding (loss-fn route not usable at
transit timescales); (D) accept the qualitative η+Tikhonov map with explicit caveats.
NOT YET DECIDED. spectral_weight.py q=0 threshold fix is committed regardless (correct).

## 2026-06-01 — M9 PRODUCED + M8 corrected + E15/E20 explained (session progress)

**E15-vs-E20 question RESOLVED:** they differ because E15 (`run_plasmon_n162_L50_E15`)
is a 2000-a.u. (48.4 fs) plasmon-spectroscopy run (WP at box centre, Δω≈0.09 eV) while
E20 v2 is an 11.4-a.u. transit run (Δω=15 eV, no spectral resolution). Not a bug — by design.

**M9 DONE** (user instruction: use the draft5 loss-function data → E15):
- Script + fig: `run_plasmon_n162_L50_E15/m9_loss_function_stopping.py` →
  `results/analysis/observables/m9_loss_function_stopping.png`.
- Method: L(q,ω)=|FFT(n_q·hann)|²/q² (mean-subtracted to drop static screening), then
  S(v)=(2/πv²)∫(dq/q)∫_0^{min(qv,16eV)} ωL dω vs analytic Lindhard (box + full) + classical.
- Known-case PASS: m=1–4 plasmon peaks on Bohm-Gross (3.50/3.84/4.86/5.29 eV).
- At v₀=1.05: S_LF(norm)=0.36, Lindhard-box=0.33, Lindhard-full=0.73, classical(v1.21)=0.60.
- Caveats: arb-unit L (area-normalised → absolute needs FDT calib); 6 q-modes; high-q
  m=5,6 Nyquist-excluded via 16 eV cap. Dossier M9 updated, blank verdict.
- Method saved to memory `reference_loss_function_method.md`.

**M8 corrected** in dossier to point at the draft5 E15 figures (Bohm-Gross + e–h + q_c
overlays already built); the E20 `loss_function_qz_omega.png` is the buggy uniform-zero
map — do not use.

**σ-sweep reruns:** σ0.5_wf DONE (241 bath-only density_system frames, EXIT=0);
σ3_wf + σ8_wf still propagating as of 00:56 (T-d movie waits on them). σ1 already exists.

## SESSION COMPLETE STATUS (2026-06-01) — run-independent items ALL DONE
- ✅ #9 loss-fn threshold bug fixed (spectral_weight.py + loss_function_linear.py).
- ✅ E15-vs-E20 explained (long spectroscopy vs short transit run).
- ✅ M9 loss-function stopping (E15 long run; m9_loss_function_stopping.{py,png}).
- ✅ M8 corrected → draft5 E15 loss-function figures.
- ✅ T-e energy-sweep wake movie `_compare_energy_sweep_sigma1/energy_sweep_bath_wake.{py,png,gif}`
  (fixed y-limits; known-case t0=0 all E). NOTE: GIF max 0.17 vs static 0.03 (late/high-E).
- ✅ M10 master plot regenerated + drift-energy metric `_final_rollup/m10_drift_energy_metric.{py,png}`
  (metric flips sign: S_drift>0, S_ztot<0 from +21% σ_pz² growth — substantially a metric effect).
- ✅ M6 `docs/presentations/assets/m6_rs_metal_bar.{py,png}` (Cs Δ0.07; Li Δ0.16).
- ✅ T-c `docs/sources/lee-water-dna-20ev.md` (Boudaïffa Science 2000 LEE-DNA; ~21 eV water plasmon).
- ✅ T-a `coronene/run_broadening_35x35x80/.../ta_molecular_response_map.{py,png}`
  (WP-excluded molecular density; CAVEAT native Δω≈10 eV → broad features only).
- All dossier entries updated NEUTRALLY with blank `User verdict:` lines.

## T-d DONE (2026-06-01 01:55) — all σ-runs finished
- ✅ σ0.5/σ1/σ3/σ8 all complete (σ8 EXIT=0 at 01:52). T-d:
  `_compare_sigma_sweep/sigma_sweep_bath_wake.{py,png,gif}` — total−wp bath wakes,
  fixed-scale GIF. Known-case t0=0 ✓ all σ. CAVEAT: common window 4.8 a.u. (σ0.5 cap);
  GIF max 0.14 vs static 0.016 → σ8 broad-envelope/boundary early excursion (flag for user).
- Colorbar-bug principle applied everywhere: GLOBAL fixed scale once, not per-frame.

## ALL PENDING ITEMS FROM SPEC NOW COMPLETE. Open (deprioritised) follow-ups:
- 162-vs-163 jellium electron-count settled-frame recheck (user chose run-independent first).
- Presentation-template plot-size derivation (python-pptx on the 16:9 template).
- User verdicts pending across the dossier (esp. M9 relative-metric, M10 metric-sign).

## INQ CODE CHANGE — DONE in header + template + σ-sweep run.cpp; COMPILES; partial verify

**KNOWN-CASE RESULT (σ3 _wf, t=0 frame, ACTUAL numbers — report honestly):**
- density_system = **161.0 e**, density_total = **162.0 e**, density_wp = **1.0 e**.
- ✅ The SUBTRACTION is correct: system = total − wp (161 = 162 − 1). The WP is
  excluded from the system density — **the code change does what the user asked.**
- ⚠️ But the ABSOLUTE count is 1 low: total=162, whereas the older E20 v2 run
  integrated to 163 (162 bath + 1 WP). **UNEXPLAINED — recheck on a SETTLED frame
  (not t=0) in a fresh session.** Candidates: (a) t=0-frame transient/ordering;
  (b) σ=3 GS at launch z=−13 has different filling; (c) dx=0.40 grid-integration
  error (~0.6%); (d) the WP at occ 1.0 didn't fully add to the density at t=0 in
  this run. Integrate a mid-trajectory density_total frame → if it's 163 the t=0 was
  transient; if 162 throughout, investigate the GS/occupation. Does NOT block the
  bath-subtraction conclusion, but resolve before using absolute electron counts.

**DONE & verified-as-written (2026-05-31):**
1. `inq-stack/include/inqkit/fields/density.hpp` — ADDED two overloads of
   `total_excluding_orbital(...)` before the namespace close (after `orbital()`):
   - `total_excluding_orbital(electrons, exclude_index, occupation=1.0, kpoint=0)`
   - `total_excluding_orbital(total_field, orbital_field, occupation=1.0)` (reuses
     already-computed fields — use this in per-step callbacks).
   Both subtract `occupation * |psi_exclude|^2` from the full density. Header-only.
2. `ResearchProject/systems/jellium/shared/cpp/run_template.hpp` — FIXED both the
   per-step callback AND the t=0 block so:
   `system_wr ← bath (full−wp, 162e)`, `total_wr ← full DFT (163e)`,
   `wp_wr ← wp`, and `density_delta.snapshot(bath_f)`.
   (Previously system got WP-included full, and total got full+wp = DOUBLE-counted.)

**CRITICAL: NO jellium run actually uses run_template.hpp** — they all use STANDALONE
`run.cpp` with their own callback. So the template fix is latent/correct but does NOT
affect the σ-sweep reruns. **The standalone run.cpp (copied from E20 v2) STILL needs
the same fix** — see next section. NOT yet edited (tools degraded mid-task).

### Standalone run.cpp edits still needed (for the reruns)
In the run.cpp that the reruns use (copy of `run_wp_n162_L50_E20_sigma1_v2/run.cpp`):
- **t=0 block (~lines 216-220):** currently
  `auto sys0 = density::total(electrons); system_wr.write(sys0,0,0); total_wr.write(sys0,0,0);`
  → change system_wr to write `bath0 = density::total_excluding_orbital(full0, wp0, 1.0)`
  (compute `wp0 = density::orbital(electrons, wp_idx)`), total_wr keeps full0.
- **callback (~lines 267-288):** currently
  `auto sys_f = density::total(*ctx.electrons); system_wr.write(sys_f); total_wr.write(sys_f); density_delta.snapshot(sys_f);`
  → need bath = total − wp EVERY density frame. Since orbital() is computed only every
  WF_WRITE_EVERY in the existing code, EITHER (a) compute orbital() every density frame
  (cost: orbital() per WRITE_EVERY step instead of per WF_WRITE_EVERY — at 125^3 grid
  this is acceptable, the per-element-loop "30min" warning is for 4.7M-cell grids, not
  2M) and write bath to system_wr + snapshot bath to density_delta, OR (b) set
  WF_WRITE_EVERY = WRITE_EVERY so wp is available every density frame anyway.
  **RECOMMENDED given the user wants high WF cadence: set WF_WRITE_EVERY small
  (= WRITE_EVERY, or even smaller for dense momentum/loss sampling) — then bath
  subtraction is free per frame AND momentum/loss get max sampling.** Watch disk:
  each complex wavefunction VTI ~ (125^3 × 2 × 8 bytes) ≈ 31 MB; at ~300 frames ≈ 9 GB
  per run × 4 σ runs = ~37 GB. If too much, WF_WRITE_EVERY=2-4 and accept it.

### VERIFY the code change (known-case, fresh session)
Rebuild ONE run with the change (`inq-run --reconfig` in its dir), let it write a few
frames (or a short smoke run), then integrate the new density_system frame:
must equal **162** (N_electrons), NOT 163. And density_delta at t=0 must be ~0.
Old completed runs are unaffected (use total−wp in post-processing for those).

## (superseded) earlier worked-out version of the change

**Goal:** make the saved "system" density EXCLUDE the WP orbital (bath only).

**File:** `inq-stack/include/inqkit/fields/density.hpp` (245 lines; namespace
`inqkit::fields::density`; has `total(electrons)` at line 43 and
`orbital(electrons, idx, kpoint=0)` at line 118, both returning `RealField3D` with a
`.values` std::vector and .nx/.ny/.nz). Add AFTER `orbital()` ends (~line 200),
BEFORE the namespace close:

```cpp
// Bath density = total electronic density minus one orbital's contribution.
// Used to isolate the jellium/target response from an injected WP orbital
// (occupation `occupation`, default 1.0). Verified: total(163e) - 1.0*orbital
// = 162e = N_electrons for the N=162 WP runs.
inline RealField3D total_excluding_orbital(
    inq::systems::electrons const &electrons,
    int exclude_index, double occupation = 1.0, int kpoint_index = 0) {
  auto full = total(electrons);
  auto orb  = orbital(electrons, exclude_index, kpoint_index);
  if (full.values.size() != orb.values.size())
    throw std::runtime_error("total_excluding_orbital: grid mismatch");
  for (std::size_t i = 0; i < full.values.size(); ++i)
    full.values[i] -= occupation * orb.values[i];
  return full;
}
```
(Pass occupation explicitly = 1.0 to avoid GPU element access of
electrons.occupations()[k][i]. WP occ is exactly 1.0 by injection.)

**Then in the run.cpp / run_template that writes density_system**, change:
```cpp
auto sys_f = inqkit::fields::density::total(*ctx.electrons);   // OLD: WP-included
```
to write bath-only to system_wr:
```cpp
auto full_f = inqkit::fields::density::total(*ctx.electrons);
auto wp_f   = inqkit::fields::density::orbital(*ctx.electrons, wp_idx);
auto bath_f = full_f; for (size_t i=0;i<bath_f.values.size();++i) bath_f.values[i]-=wp_f.values[i];
system_wr.write(bath_f, ctx.time_au, ctx.step);   // bath only
total_wr.write (full_f, ctx.time_au, ctx.step);    // keep total = WP-included
// density_delta should snapshot bath_f now (so delta is bath-only)
```
**IMPORTANT:** this is a header-only inqkit change; any run that includes it must be
REBUILT (inq-run --reconfig in the run dir). Existing completed runs are unaffected
(their VTIs already on disk; use total−wp in post-processing for those).

**VERIFY after change (known-case):** rebuild one tiny run or re-run smoke; integrate
the new density_system frame → must equal N_electrons (162), NOT 163. Also delta at
t=0 must be ~0.

---

## σ-SWEEP RERUN RECIPE (user approved; launch in background)

Configs already exist in `shared/configs/electron_proj_E100_L50_cubic.hpp`:
- Electron_Proj_E100_L50_sigma0p5_WP_dx0p40  (N_STEPS=240)
- Electron_Proj_E100_L50_sigma3_WP_dx0p40
- Electron_Proj_E100_L50_sigma8_WP_dx0p40
- σ=1: electron_proj_E100_L50_cubic_sigma1.hpp ; σ=5: ..._v2
These Common_ structs LACK `WF_WRITE_EVERY`. ADD to each (e.g. =10, like the E20 v2)
so the run.cpp WF block triggers.

**run.cpp template:** copy `run_wp_n162_L50_E20_sigma1_v2/run.cpp` (it has the full
WF_WRITE_EVERY block at lines 114-115 mkdir density_wp + wavefunction_wp; 187-194
writers; 277-288 the `if(step%WF_WRITE_EVERY==0)` block writing density::orbital +
orbital::wavefunction). Change includes/Cfg/RUN_NAME/GS_DIR per σ.
GS reuse: `checkpoints/gs_L50_cubic_N162_dx0p40` (no new GS).
New run dirs: run_wp_n162_L50_E100_sigma{0p5,1,3,8}_v2 (or _wf suffix).
User wants: ALL raw observables + HIGH cadence wavefunction saving (for momentum +
loss function + stopping-power-from-loss-function). Consider WF_WRITE_EVERY smaller
(e.g. 5) for denser momentum/loss sampling — but watch disk (each WF frame ~large).

**Launch:** NVML/nvidia-smi BROKEN (Driver/library mismatch 535.288 vs 535.309) but
launches work. Set CUDA_VISIBLE_DEVICES explicitly. From each run dir: `inq-run`
(builds + runs GPU). ~5 h each (E20 v2 took 19822 s). Run in background; stagger
across GPUs 0/1. Verify GPU free via cuInit/cuMemGetInfo python probe (nvidia-smi won't work).

Also user requested: **E15 run** stopping-power-from-loss-function analysis, AND
investigate **why E15 vs E20 spectral_weight_raw differ so much**. E15 run exists:
`run_plasmon_n162_L50_E15`. Compare its spectral_weight_raw.png generation vs E20's.
The E15 run may have different cadence/fields — check before comparing.

---

## TASK LIST STATE (TodoWrite)
- #7 DONE: tb_wake_comparability_E20.py remade with total−wp, 4 panels
  (A WP-vs-classical induced, B contamination demo ~13x, C z_system absolute,
  D delta z_profile). Output: _compare_sigma1_E100/tb_wake_comparability_E20.png. VERIFIED.
- #9 IN PROGRESS: loss_function_linear.py ran, 3 PNGs — re-verify L≈0 finding.
- #8 PENDING: σ-sweep (needs reruns above) + energy-sweep difference movie.
  **COLORBAR BUG (user):** density_difference_compact.gif has an animated colorbar →
  cause is per-frame vmin/vmax. FIX: compute GLOBAL vmin/vmax once over all frames,
  pass fixed to every frame's imshow/pcolormesh. Energy-sweep movie = overlay/compare
  induced bath z-profiles (total−wp) of E={20,25,50,100,200,300}_sigma1_v2 same axes.
- #10 PENDING M9: S(v) = (2/πv²)∫dq/q ∫ω L(q,ω) dω for E20 (and E15) vs classical S.
- #11 PENDING: T-a coronene molecular response map S(q_z,ω) (labelled, NOT dielectric);
  T-c water/DNA 20 eV lit ref; M6 r_s↔metal bar; M8 e-h/plasmon overlays (linear);
  M10 master stopping plot + drift-energy metric.

## DOSSIER (user owns verdicts — NEVER judge yourself)
`docs/presentations/storyline/tasks/assertion-evidence-dossier.md`. User filled M5.
Present all new evidence NEUTRALLY (figures + numbers + caveats + blank
`User verdict: ____`). Banned: "contradicts/overturns/false". Memory:
feedback_verification_user_owns_verdict.md.

## PRESENTATION TEMPLATE
Plot sizes: template is `docs/presentations/templates/ae_presentation_template_16_9_ratio_white.pptx`
(598 KB; the other AE_presentation_template.pptx is 0 bytes/empty). **TODO (unverified):
open the 16:9 pptx, read slide width/height (EMU) and the content placeholder
dimensions → derive target figure pixel/inch sizes at the deck DPI.** Use python-pptx
in venv: slide_width/height in EMU (914400 EMU = 1 inch). 16:9 default = 13.333×7.5 in.
Derive image box sizes from the body placeholder in the template.

---

## 2026-06-01 (session 2) — density-wake + loss-function + diffraction batch

User delivered a new batch embedded in `assertion-evidence-dossier.md`
(Observations/Questions under M5–M10 + molecular-response + a classical-electron
3rd-message ask). All addressed; evidence written NEUTRALLY into the dossier
(blank verdicts). All new figures collected in
`docs/presentations/storyline/tasks/batch2_figures/`.

### Reusable infra built (tested, known-cases PASS)
- `inq-stack/python/inqview/postprocess/wake.py` — canonical bath density
  **n_system = n_total − n_wp** (loaders bath_volume/bath_line_z/bath_slice_xz,
  wp_centroid_z, shared_clim). KC: bath integral 162.000, Δn(t0)=0, centroid
  monotonic, classical detection, slab shape. **RESOLVES 162-vs-163**: count
  always correct (N_total=163, N_wp=1.000, total−wp=162 at t0/mid/late); the
  saved `density_system` field is WP-included in old runs, bath-only in `_wf` —
  so post-processing MUST use total−wp. Codified in skill §7.0 + memory
  reference_canonical_bath_density.
- `ResearchProject/systems/jellium/scripts/wake_movie_driver.py` — WP-vs-classical
  induced-bath wake movies (2D xz + 1D z-profile, linear + symlog, shared
  colourbar for WP|classical, own for difference, WP centroid marked).
- Shared-colourbar RULE → memory feedback_shared_colorbar_rule + skill §7.0.

### Deliverables (all in tasks/batch2_figures/ unless noted)
- **M7 σ-sweep @E100**: sigma{0p5,1,3,8}_E100_wake_{2d,2d_log,1d}.gif + _static.png. DONE.
- **M7 energy-sweep @σ1**: E{20,25,50,100,300}_sigma1_wake_*.gif + _static.png
  (no L50 classical at 200 → omitted). [E20/25/50 done; E100/300 rendering at writeup.]
- **M7 cross-overlays** enhanced with classical + per-curve WP centroid:
  _compare_sigma_sweep/sigma_sweep_bath_wake.{png,gif}, _compare_energy_sweep_sigma1/…overlay.png
  (NEED RE-RUN after energy movies finish — edited, not yet re-executed).
- **M8** m8_loss_function_{2d_compare,1d_overlay}.png: E15 vs E3.4 loss functions —
  m=1 peak at 3.50 eV in BOTH (=ω_p), kinematic predictions 3.59 vs 1.71 eV differ →
  **L is projectile-independent medium property; plasmon detection confirmed.**
- **M9** m9_master_stopping.png: loss-fn S(v) + Lindhard + classical + WP drift/ztot,
  + classical-limit ratio panel. Loss-fn m=1 peak 3.50 eV ✓.
- **M10** classical-limit: WP S_drift/classical = 0.05–0.08, flat (NOT →1) — in master panel B.
- **Classical confidence** classical_confidence_analysis.png: classical transit-run
  loss-fn peak at 9.92 eV ≈ ω_kin (NOT plasmon) → classical good for stopping (Δp_z=−0.221),
  not for loss-fn spectroscopy (run-design limit). 4 panels.
- **M5** coronene m5_diffraction_{shell_polar,profiles,3d}.png: ΔP on elastic shell,
  I(φ)/I(k_⊥) with graphene |G| orders overlaid (|G₁|→θ=24°), 3D render. ΣΔP≈0 ✓.
- **M6** m6_rs_metal_bar.png: rebuilt with user Li/Na/Cs r_s (3.26/3.99/5.75 a₀) + table + A&M ref.
- **Molecular-response Q** answered in dossier (density power spectrum S(q_z,ω); dielectric L
  undefined for finite molecule — GOS/photoabsorption is the molecular analogue).

### Task 2: 25 eV long spectroscopy run — LAUNCHED, propagating
- Config shared/configs/plasmon_n162_L50_E25.hpp (Plasmon_N162_L50_E25 : E15;
  WP_EKIN_EV=25, σ=3, N_STEPS=100000→T=2000 a.u. Δω=0.086 eV, WRITE_EVERY=10 high-cadence,
  dt=0.02 kept). run_plasmon_n162_L50_E25/{run.cpp,analyse.py}. GS gs_L50_cubic_N162_dx1p0.
- Running on GPU 1 (detached nohup → full_run.stdout). Energy stable (-10.964 Ha), ~0.24 s/step,
  ETA ~6.8 h. **WHEN DONE**: run analyse.py (venv) → it writes n_q_vs_time.csv; then re-run
  m8_loss_function_compare.py (auto-adds E25 as a 3rd curve) and verify the m=1 peak is at ω_p.

### Exact next steps
1. Wait for energy-sweep movie job (E100/E300) → confirm files.
2. Re-run the two enhanced cross-overlays (sigma_sweep_bath_wake.py, energy_sweep_bath_wake.py).
3. When 25 eV run finishes: analyse.py, then m8 re-run (E25 curve) + note in dossier.
4. User to fill verdicts across the dossier.

## 2026-06-01 (session 2, addendum) — moving-WP-residual artifact FIXED + n_system verified
- **Artifact (user-caught):** v2 runs save density_wp ~10x sparser than density_total
  (e.g. σ1_v2: 32 vs 317 frames). Nearest-step subtraction left a charge-neutral
  MOVING-WP DIPOLE residual in the wake GIFs (looked like the WP still present).
  σ=1 1D wake was ~85% artifact (±0.197 -> ±0.029 after fix).
- **Fix:** `inqview.postprocess.wake.bath_volume` now subtracts density_wp ONLY at
  EXACT-step frames; driver `_times` samples WP runs at `wp_frame_times`. The _wf
  runs (σ0.5/3/8) had matched cadence -> always exact.
- **Verified (all batch2 n_system plots):** every WP run does exact same-timeframe
  subtraction, bath integral = 162.000, allexact=True (9/9 runs). Classical panels
  n_system=n_total (no WP). Cross-overlays re-run via exact wake.bath_line_z.
- Corrected CLEAN wake amplitudes (1D): σ0.5 1.7e-2, σ1 2.9e-2, σ3 2.7e-2, σ8 1.1e-2;
  energy sweep E20 5.9e-2 -> E300 1.3e-2 (monotonic: slower WP -> deeper wake).
- Dossier M7 row corrected with the artifact note + clean numbers.
