# Smoke-test & launch-readiness — QKE + loss-function campaigns

**Date:** 2026-06-22. **Purpose:** make both GPU campaigns *launch-ready* (build,
GS-load, propagate, no-NaN validated; one-command launch) so they start with no
fuss when their gates lift. Tier-A smoke only (seconds–minutes); the expensive
production runs were NOT run. `check-stopping-power` has no GPU runs → N/A.

**GPU:** A30 ×2; `cudaMemGetInfo` probe used (NVML/nvidia-smi broken on this host).
Both free at end of session.

---

## Campaign 1 — quantum-kick-extension (omega_peak(v) drift, 5 new multi-k runs)

**Status: LAUNCH-READY (blocked on the fourier-analysis skill, per campaign).**

Infrastructure added/fixed:
- 5 Cfg structs `Li_54_v0p0375/100/175/250/375` in
  `QuantumKickExtension/inq-codebase/Li/shared/configs/base_li_54.hpp`.
- 5 run wrappers `run_propagate_v{0p0375,0p100,0p175,0p250,0p375}/run.cpp`
  (thin wrappers over `shared/cpp/run_template.hpp` + the Cfg struct).
- `QKE_N_STEPS` env override added to `run_template.hpp` (defaults to
  `Cfg::N_STEPS`; lets any run be smoke-tested in N steps).
- **BUG FIXED (caught by the smoke):** `run_template.hpp:249` used the old
  `CenterOfDensityResult` API `cod.x_bohr` — the struct now nests the centroid in
  `cod.center_bohr.x`. Without this fix the whole sweep would NOT compile against
  the current inqkit.
- Launch script `scripts/launch_omega_peak_sweep.sh` (new Cfg-struct pattern; the
  legacy `run_velocity_sweep.sh` uses the old monolithic pattern and is NOT used).

**Tier-A smoke (v0p175, GPU 1, `QKE_N_STEPS=20`): PASS**
- Build clean (after the fix); GS `li_54_2x2x2_T200_xyz` loaded;
  `num_states=101, num_electrons=162`; `energy_total(step0)=−389.3143 Ha` (= GS);
  20 steps, **0 NaN/inf**, all rows finite. Smoke output cleared so production runs fresh.

**To launch (after the gate lifts):**
```bash
cd QuantumKickExtension/inq-codebase/Li
bash scripts/launch_omega_peak_sweep.sh 1            # 5 runs, GPU 1, ~30 h each
# smoke any single run:  QKE_N_STEPS=20 bash scripts/launch_omega_peak_sweep.sh 1
```

---

## Campaign 2 — cap-jellium-loss-function (9 runs: classical/wp/kick × E15/E20/E30)

**Status: LAUNCH-READY.** (Simulations may run before the Fourier gate; only the
spectral analysis, tasks 5–7, is gated.)

Infrastructure built:
- `scripts/cap_loss_function/run.cpp` — env-driven (`CAP_MODE`, `CAP_V0`,
  `CAP_N_STEPS`, `CAP_WRITE_EVERY`, `CAP_OUT_SUBDIR`). classical/wp are
  byte-identical to the proven `cap_baselines` b2/b3; the **kick** mode is new
  (`perturbations::kick{cell,{0,0,K}}`, velocity gauge, CAP OFF). Built ONCE
  against **inq-study** (binary at `scripts/cap_loss_function/run`).
- `scripts/cap_loss_function/dispatch.py` — 9-run matrix, `--smoke`, GPU
  auto-pick via `cudaMemGetInfo`, per-run email hook.
- **BUG FIXED (caught by the smoke):** the velocity-gauge kick makes INQ's
  `data.iter()` unreliable *inside the observer* — the kick run wrote only step 0
  (1 VTI frame) while classical/wp wrote every `WRITE_EVERY`. Fixed by driving the
  write cadence from a self-incremented counter (`wcount`) in `step_fn` rather
  than `data.iter()`. Without this the kick (q=0 plasmon reference) runs would
  have produced no n_q time series.
- **BUG FIXED:** `dispatch.py` binary path (`run`, not `build/run` — `inq-run`
  emits the binary in the run dir).

**Tier-A smoke (all 3 modes @ E15, `CAP_N_STEPS=20`): PASS**
| mode | obs rows | VTI frames | NaN |
|---|---|---|---|
| classical | 12 (steps 0,2,…,20) | 11 | none |
| wp | 12 | 11 | none |
| kick | 12 | 11 | none |

Kick KE jumped by ½·162·1.05²≈89 Ha (velocity-gauge boost of the electron sea, as
intended); energy conserved post-kick (CAP off). Smoke outputs cleared.

**NOT validated by the smoke (deferred to the campaign's gated pilot):** "CAP
absorbs the projectile" (needs ≫20 steps) and "kick rings at ω_p" (needs the
gated Fourier analysis). The campaign's pilot-first numeric gate still applies
before the long runs.

**To launch:**
```bash
cd ResearchProject/systems/jellium/scripts/cap_loss_function
python3 dispatch.py --smoke              # 20-step readiness re-check, all 9
python3 dispatch.py --gpu 1              # PRODUCTION: T~2000 a.u. (~100k steps), HEAVY (~12 h each)
# rebuild if needed:
#   INQ_SOURCE=…/inq-study INQ_SHARE_PATH=…/inq/install/share \
#   PSEUDOPOD_SHARE_PATH=…/inq/install/share/pseudopod inq-run run.cpp
```

---

## Value of doing this now
The smoke tests caught **three real breakages** that would each have blocked a
"come back later and launch" attempt: the QKE `center_of_density` API drift
(compile failure), the kick observer-cadence bug (kick runs would yield no n_q
data), and the dispatcher binary-path bug. All fixed and re-validated.
