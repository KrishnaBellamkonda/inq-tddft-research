# Handover: per-state mass fork (inq-study) — muon / band-structure

Rolling handover. Task: add a tunable per-state mass to inq-study so INQ can
simulate a muon projectile, all-muon jellium, and mass-tuned band structure.
Plan: `docs/plans/muon-mass-fork-implementation.md`. Design + call chain:
`docs/campaigns/muon_projectile/inq_study_engine_notes.md`. User understanding +
feedback: `docs/notes/muon_plan_understanding.md`.

## 2026-07-11 (evening) — CAP study VERDICT: current CAP validated, weak-η is the only bad one

All 3 variants + baseline complete; 3 run-notebooks built. Cross-CAP metrics
(N from ∫density_total; reflection = backward WP density past launch z<−16.5 in
the late half; absorption = WP orbital norm end):

| run | CAP | N_min | WP-absorbed | refl%(late) | E_end (Ha) |
|---|---|---|---|---|---|
| baseline | η=−1.0 [25,40] gap12.5 | 51.999 | 100.0% | 0.13 | −27.67 |
| R1 | η=−1.0 [32,40] **gap19.5** | 52.003 | 99.7% | 1.15 | −27.42 |
| R2 | **η=−0.4** [25,40] | 52.012 | 98.75% | **4.14** | −27.67 |
| R3 | **η=−2.0** [25,40] | 51.999 | 100.0% | **0.001** | −34.88 |

**VERDICT (decision metric = total N(t), per user):**
1. **No CAP eats bath charge.** N_min ∈ [51.999, 52.012] for ALL four → the CAP
   absorbs exactly the 1-e⁻ WP and nothing else, in every configuration. The
   suspected bath/wake-eating problem is **NOT present**.
2. **Current CAP (η=−1.0) is well-chosen**: complete WP absorption, only 0.13%
   reflection, no bath loss. KEEP IT.
3. **Weak η=−0.4 is the only defective one**: 4.1% edge reflection (leaky branch —
   packet reaches box edge and reflects) + 1.25% un-absorbed WP residue. DO NOT
   weaken below η≈−1.0.
4. **Strong η=−2.0**: even cleaner reflection (0.001%), complete absorption, N=52
   → η can safely be raised, no need. BUT its **E_end=−34.9 Ha vs −27.7** elsewhere
   confirms the CAP removes strength-dependent energy → **the WP energy-method S is
   CAP-dependent (a bound), NOT a clean stopping observable**. Use the classical twin
   / near-field WP momentum for S, never E_total after absorption.
5. **Wider gap (R1)** reduces wake-under-CAP (0.33→0.22 e⁻) but its narrower 8-Bohr
   width is less adiabatic → reflection up to 1.15%, 0.3% residue. Marginal; only
   adopt if a future wake analysis shows clipping at |z|=25 matters.

Notebooks: `hypotheses/muon_mass_fork/effmass_sigma1_cap_{gap19p5,eta0p4,eta2p0}_wp_run.ipynb`
(45 cells each, 0 errors, linear|log GIFs + ΔE|N(t) panel). Runs:
`scripts/muon_mass_fork/effmass_sigma1/wp/results/cap_{gap19p5,eta0p4,eta2p0}`.
NOTE reflection metric bug avoided: WP launches at z=−16.5 (behind slab), so an
early-frame z<−12.5 window mis-reads the launch as reflection; use z<−16.5 in the
LATE half of the run.

## 2026-07-11 (later) — CAP-parameter study LAUNCHED + run-notebook rules enforced

**Run-notebook skill changes (enforced, 2026-07-11):**
1. **Log alongside linear everywhere.** `density_gifs.py::_save_gif` now renders
   every battery GIF **LINEAR | LOG** side by side — density: linear + `LogNorm`;
   Δn (delta0/dstep): linear + `SymLogNorm` (linthresh=vmax/100). SKILL.md + builder
   captions updated.
2. **N(t) beside ΔE_total.** `delta_total_energy_fig` is now a 2-panel
   [ΔE_total(t) | N(t)=∫n dV]; N from `electron_number.csv` (classical) else ∫
   `density_total` VTIs (WP). Old standalone `energy_and_number_fig` removed.
   Both rules validated on real data + logged in `docs/validation/test-catalogue.md`.

**CAP-parameter study (user-requested; Fable-5 advisor design, Opus-approved).**
Suspected CAP problem in the σ=1 WP run. Measured on the current run (η=−1.0,
region [25,40], gap 12.5): WP norm→0, **reflected fraction ≤0.13%**, **total N(t)
53→52, never <52** (bath preserved), induced wake under |z|≥25 peaks 0.33 e⁻.
Reflection is already tiny → the real suspects are **wake clipping** + η-robustness.
**Decision metric = total N(t)** (user: judge absorption on whole-cell electron
count, NOT WP orbital norm alone) — the contamination signature is N(t) dipping
below 52.

Three variants of the SAME run (parametrised `wp/run.cpp` — CAP now env-overridable
`EM_CAP_ETA`/`EM_CAP_CENTER_BOHR`/`EM_CAP_WIDTH_BOHR`, defaults reproduce the
original; binary rebuilt + validated):
- **R1 cap_gap19p5** η=−1.0, centre 36, width 8 → region [32,40], **gap 19.5** (protect wake). cap_inner=32.
- **R2 cap_eta0p4** η=−0.4, centre 32.5, width 15 → [25,40], gap 12.5 (weak-η branch, OD≈4.4).
- **R3 cap_eta2p0** η=−2.0, centre 32.5, width 15 → [25,40], gap 12.5 (strong-η branch, OD≈22).

Autonomous detached orchestrator:
`scripts/muon_mass_fork/effmass_sigma1/cap_study/orchestrate.sh` (setsid, survives
logout). R1(gpu0)+R2(gpu1) concurrent, R3 on first freed GPU; then builds a
run-notebook per variant (`effmass_sigma1_cap_{gap19p5,eta0p4,eta2p0}_wp_run.ipynb`).
Launched 2026-07-11 01:14; ~75 min/run → ~2.3 h total. **GPUs shared with another
user's Li/Frenkel jobs** — lean config has headroom, runs started at normal speed.
Log: `cap_study/orchestrate.log`. Decision rule: if R1 matches current within noise
AND R2/R3 bracket it flat (and all keep N(t)=52), current CAP validated; else adopt
R1 geometry as production.

## 2026-07-11 — σ=1 LEAN pair COMPLETE + both run-notebooks built

Both `effmass_sigma1` runs finished cleanly (both GPUs now free) and each has a
deep single-run notebook in `hypotheses/muon_mass_fork/` (skill-local builder,
`--gif-seconds 17`; args: cap-inner 25, rs 5.684, proj-sigma 0.7071, launch-z
−16.5, v0 2.7111, l-slab 25; WP also `--e-gs-ha -36.9404590471`).

- **WP (chirped focus)** `effmass_sigma1/wp/results/sigma1`, wall **4155 s = 1.15 h**,
  900 steps, run_completed=true. Notebook `effmass_sigma1_wp_run.ipynb`
  (**45 cells, type=wp, 0 exec errors**; 9 density GIFs + 12 PNGs).
  Physics: focusing **σ_z 0.864→0.729 Bohr waist at t=1.28 au** (target 0.707 at
  the slab face) → arrives compact, comparable to the rigid classical; **complete
  absorption norm 1.000→0.000** (clean CAP, no reflection); z_cod −16.5→+14.7.
  WP energy-method `S=[E_total(t_f)−E_GS]/L_z` prints ~**10 eV/Bohr but is a BOUND**
  — CAP fully absorbs the WP (E_total also books CAP energy removal) and the WP
  makes the cell net −1 charged (convention-dependent Hartree, see
  [[reference_charged_cell_hartree_convention]]). NOT directly comparable to the
  classical channel; reconciliation = the WP−classical analysis, deferred.
- **Classical twin (m=2.10, σ_pot=0.707)** `effmass_sigma1/classical/results/classical`,
  wall **4527 s = 1.26 h**, park_step=426 (z=26.27), run_completed=true. Notebook
  `effmass_sigma1_classical_run.ipynb` (**42 cells, type=classical, 0 exec errors**;
  3 induced-wake GIFs). **S = 0.208 eV/Bohr** (equal-potential slab faces ±12.5,
  ΔKE≈0.20 Ha / 25 Bohr; v_z 2.57→2.53 across the slab, ≥0.85·v₀ → clean
  initial-drag window per light-projectile rule). The raw-file final v_z→0 / KE→0
  is the **park-at-CAP, NOT physical stopping**.

**Library gap noted (pre-existing, harmless to the deliverable):** the pipeline
`stopping` phase `[fail]`s on any classical run whose `electron_track.csv` lacks an
`fz` column (only `fullsuite/p5` ever logged force; p3/p4/wide_wp/this do not).
The builder computes the PRIMARY classical stopping panel itself from `ke_ion_ha`
(`classical_transport.png`, S=0.208), so the notebook is complete; the failed phase
only makes the un-embedded secondary `dE_kinetic_vs_z.png`. Fix (if wanted): make
`inqview/pipeline/stopping.py` select track columns dynamically and guard the
force plot on `"fz" in columns` — does NOT change the run-notebook content.

**Next:** WP−classical comparison (isolate purely-quantum Δ at matched
width/mass/velocity/charge) — the actual deliverable the user is after; requires
reconciling the WP energy channel vs the classical KE channel (CAP energy
accounting + charged-cell Hartree convention).

## 2026-07-09 — σ=1 LEAN re-run AUTONOMOUS (detached orchestrator, user may log off)

First σ=1 launch (50×50×101 @ dx=0.333, 61 states) hit the 24 GB VRAM wall →
104 s/step → 34 h ETA → KILLED at step ~36 (physics clean; in-medium focusing
CONFIRMED: σ_z 0.864→0.729 @ t=1.28 toward 0.707 waist). Lean re-plan (user
geometry): `docs/plans/effmass-sigma1-lean-rerun.md`.

**Lean config:** 40×40×80 @ dx=0.333 (120×120×240 = 3.46 M pts), density-matched
**N=52** (r_s=5.679), 26 occ + 10 extra = **36 states**, N=900 steps (3.9×
traversal), dt=0.04, **3 checkpoints** (300/600/900-final, both runs,
extension-ready), CAP η=−1.0 region [±25,±40] (15-Bohr width kept), launch
−16.5, `focus_z(4.0, 2.10)`, k0=5.693, m=2.10. Empirical cost model (3
calibration points, log-log fit): **5–8 s/step ≈ 1.5–2 h**; GS ~35–45 min.

**AUTONOMOUS CHAIN — running without any session:**
`scripts/muon_mass_fork/effmass_sigma1/orchestrate.sh`, launched setsid/nohup
(pid 3985254, PPID=1, own SID — survives SSH logout). It: (1) waits for /
relaunches the lean GS → `shared_gs/slab_n52_L40x40x80_dx0p333`; (2) launches
WP (GPU 0, `wp/`, out `results/sigma1`) + classical twin (GPU 1, `classical/`,
out `results/classical`) CONCURRENTLY; (3) logs every transition to
`effmass_sigma1/ORCHESTRATOR_STATUS.txt` (check this file first on resume).

- GS was ~10 min into SCF at orchestrator launch (23:47 start; GPU 0 at 62 % =
  NOT memory-walled ✓ key prediction holding).
- Classical binary PRE-BUILT + link-verified (fail-fast exit 2 on missing GS).
  Classical run.cpp is NEW: chunked Ehrenfest park+remove (park |z|≥25), ckpt
  saves electrons + ion state, **fixed KE ledger 0.5·M_ME·v² (template used
  0.5·1.0·v² — would under-log KE 2.1× at m=2.10)**; resume = fresh-chunk
  continuation restoring ion state, density_delta re-anchored to GS (UNTESTED).
- Classical UPF `electron_gaussian_wpsigma1p0.upf` (σ_pot=0.707) generated +
  DATA-VERIFIED by V(r) (erf model ratio 1.000 at all r).
- WP resume smoke-tested earlier (bit-faithful, effmass_12h). focus_z vacuum
  test PASS at dx=0.333 (test-catalogue row added).
- NEXT on completion: verify pilot rate vs 5–8 s/step prediction, run-notebooks
  for both runs, WP−classical comparison.

## 2026-07-08 (later) — σ-comparison twin runs: σ≈2 (rigid) + σ≈1 (concentrated, chirped)

**Goal (user):** two mass-fork WP runs on the localised jellium slab, each comparable
to a rigid classical Gaussian-charge projectile, to isolate purely-quantum effects.
Launched on two GPUs. effmass_12h (m=2.51) is DISTRUSTED → demoted to a later phase;
"figure out the mass correctly." Uses a **Fable-5 advisor / Opus decision-maker** loop
(user-requested).

### Run 1 — σ≈2 "important" (GPU 1) — LAUNCHING
- `scripts/muon_mass_fork/effmass_sigma2/{gs,wp}/run.cpp`. GS **done**
  (`shared_gs/slab_n82_L50x50x90_dx0p40`, r_s=5.667, E=−72.50 Ha, 61 states).
- **σ_WP=2, v=2.711, dx=0.40, k0=6.793, m=2.506 (inv_mass=0.3991), E=251 eV**,
  box **50×50×90 N=82** (no transverse shrink → no wake-wrap, unlike effmass_12h),
  launch z0=−19.21 (6.71 Bohr standoff = 3% spread), CAP **η=−1.0** region [±35,±45],
  dt=0.05, **N=819** (3.5× traversal τ=11.7), **CKPT_EVERY=234** (4 checkpoints).
  thru-slab spread ~54% (partial rigidity; σ≥3 needed for true rigidity).
- WP run building against **inq-study** on GPU 1. ETA ~6 h (3.5M pts, 61 states, ~30 s/step est).

### Run 2 — σ≈1 "concentrated" (GPU 0) — PENDING CHIRP FEATURE
Fable-5 deliberation (advisor) + Opus decision. Key findings:
- **σ=1 cannot be rigid through the slab** (theorem: needs dx≤0.04). It is the
  *dispersion-dominated end-member* of a σ-family {0.5,1,2,3.5}. WP−classical ΔS at
  σ=1 is dispersion-dominated, NOT purely quantum. Defensible deliverables instead:
  entry-window S, **in-medium dispersion measurement** (a quantum observable), σ-family trend.
- **<3% centroid spread at impact is incompatible with clean 4σ injection** (centroid
  enters at σ=1.22 = 22%). **User chose Option B: CHIRPED (focusing) launch** — a
  quadratic phase exp(iα(z−z0)²) so the WP waist (σ=1.0, zero spread-rate) sits AT the
  slab face → literally <3% at impact with a clean 4σ standoff.
- Config (Fable): σ=1, v=2.711, **m=2.10, k0=5.69**, dx=0.40, box **50×50×101**
  (reuse existing GS `shared_gs/slab_n82_L50x50x101_h0p40`), CAP **η=−1.0, 15 Bohr/side**
  inner ±35.5, launch z0=−16.5 (4σ), dt=0.05, **N=1400** (3.6× traversal), 4 checkpoints,
  ~8 h. Classical twin **mass 2.10** (MUST match WP eff-mass), σ_pot=0.707,
  UPF `electron_gaussian_wpsigma1p0.upf` (generate + verify by V(r)).

### NEXT (Run 2 blockers)
1. **Add chirp/focus to `inqkit/wavepacket/wavepacket.hpp`** — `.focus(focal_z)` or
   `.chirp(alpha)` applying exp(iα(z−z0)²); α set so waist lands at the slab face.
   Ships a **vacuum focus test** (known-case: WP converges to σ=1.0 at focal point,
   free-particle law) per code-test. inqkit change (allowed) — record in test-catalogue.
2. Generate `electron_gaussian_wpsigma1p0.upf` (σ_pot=0.707) for the classical twin.
3. Build σ=1 wp/run.cpp (adapt effmass_sigma2, add chirp), launch on GPU 0.
- Fable also recommends a cheap **vacuum-WP pilot** (same σ,m,v, no slab) to validate
  the fork + bound SIE before science claims.

## 2026-07-08 — effmass run RE-PLANNED for ≤12 h (N3/A) + LAUNCHED with checkpointing

### Context
The first effmass run (`effmass_pair`, dx=0.333, m=3.09, 82 e⁻, 61 states) measured
**~200 s/step** (contended) → ~4–7 day ETA. Old cost model (calibrated on qsp_phase4)
was ~15× optimistic. User: stop it, re-plan for **≤12 h**, relax spread <1%→**<2%**,
smaller cell + coarser grid, **add RT checkpointing**, defer the classical twin to a
draft phase. Full spec: `docs/plans/muon-effmass-12h-run.md`.

### Locked spec (N3/A) — user-selected via grill-with-docs
- Cell **36×36×80**, dx=**0.40**, **N=42** (r_s=5.689), 31 states. New GS + config.
- σ_WP=2, k0=**6.7933**, **m=2.506** (inv_mass=0.39907), E=**251 eV**, v=2.711 (=100 eV e⁻).
- launch_z=−16.389, impact spread 1.0%, **dt=0.05**, **846 steps** (T=42 au=3.5×traversal),
  CAP η=−0.7 retuned for Lz=80.
- Velocity v=2.711 LOCKED (S(v) anchor). Transverse 50→36 → **wake-wrap** validation flag.

### New files
- `shared/configs/slab_n42_L36x36x80.hpp`
- `scripts/muon_mass_fork/effmass_12h/{gs,quantum}/run.cpp`
- GS checkpoint: `shared_gs/slab_n42_L36x36x80_dx0p40` (E=−30.128 Ha, 31 states, r_s=5.689).

### CHECKPOINT / RESUME (user request) — VALIDATED
Quantum run.cpp uses INQ's **native `start_step`**: every `EM_CKPT_EVERY` (100) steps
`electrons.save()` + `rt_state.txt`; `EM_RESUME=1` reloads, **re-applies inverse_mass**
(save/load does NOT persist it — verified in inq-study source), resumes via `start_step`.
Static jellium ions + static bg/CAP ⇒ bit-faithful. Smoke test PASSED: continuous
step-6 E vs (0→3, save, resume 3→6) step-6 E **identical, abs diff 0.000e+00 Ha**.
dt=0.05 stable (H·dt=1.54; no NaN).

### Build note (important)
Quantum MUST build against **inq-study** (`INQ_SOURCE=…/inq-study`) — `inverse_mass`
lives there with a modified propagator. Pristine `inq/` has no such member (first
compile failed on exactly that). inq-study propagate HAS `start_step`; save/load is
byte-identical to inq/ (GS built on inq/ loads fine).

### STATUS: RUNNING (2026-07-08 09:33)
PID 2314600, GPU 0, `results/quantum/`. **~6.9 s/step → ETA ~1.6 h** (31 states + freer
GPU beat the ~11.5 h conservative estimate). Energy stable (−20.84 Ha = GS −30.13 +
WP kinetic +9.2 Ha). Next on completion: build run + comparison notebooks, email.
Then scaffold the **classical draft** (free Ehrenfest, m=2.506, electron_gaussian_wpsigma2p0.upf).

## 2026-07-06 — Implementation complete, compilation PENDING

### DONE (code written, NOT yet compiled/tested)
Strategy A implemented in `inq-study` (immutable `inq/` untouched). Mechanism:
per-state inverse mass on `electrons`, opted into the Hamiltonian via a
`set_inverse_mass()` setter, applied by distinct `_states` kernels. An **empty-
factor guard** routes all-mass-1 (electron) cases through the ORIGINAL scalar
path → intended bit-for-bit Tier-0 invariant.

| File | Change | State |
|---|---|---|
| `src/systems/electrons.hpp` | `inverse_mass_` member + accessors + reextent/fill(1.0) at both alloc sites | written |
| `src/operations/laplacian.hpp` | new `laplacian_states`, `laplacian_add_states`, `laplacian_expectation_value_states` (kernel `fac[ist]`) | written |
| `src/hamiltonian/ks_hamiltonian.hpp` | `kinetic_factor_` member + `set_inverse_mass()` setter + empty-guard branch at the 3 kinetic call sites | written |
| `src/real_time/propagate.hpp` | `ham.set_inverse_mass(electrons.inverse_mass()[0])` after ctor | written |
| `src/ground_state/calculator.hpp` | `ham_.set_inverse_mass(...)` in ctor body | written |
| `src/ground_state/initial_guess.hpp` | NOTE only — no change (overlap-only Hamiltonian) | done |

Constructor signatures UNCHANGED (setter approach) → all other call sites
(`initial_guess`, `paw`) source-compatible, untouched.

### Key design decisions (survive compaction)
- **Empty `kinetic_factor_` ⇒ scalar path.** `set_inverse_mass` computes the
  deviation `Σ|im-1|`; if 0, leaves `kinetic_factor_` empty so the unforked scalar
  `laplacian(-0.5)` code runs verbatim (bit-for-bit). Only non-trivial mass builds
  `kinetic_factor_[ist] = -0.5*im[ist]`.
- **Distinct `_states` names, not overloads** — protects existing callers.
- **Alignment safety:** per-state mass only matters in RT on the full set (muon in
  electron jellium), where `kinetic_factor_` aligns with `phi.set_part()`. GS mass
  is uniform (all-1 or global muon) so any block indexing is harmless.
- **Gamma-only:** setter uses `inverse_mass()[0]` (kpin 0). Multi-kpin = future.
- **Spinor:** implemented for `spinor_dim==1` (jellium). The `_states` asserts
  `factors.size()==local set size`; a spinor run would trip it → revisit if needed.
- **Vector-potential caveat:** the `+A²/2m` term is NOT mass-scaled → valid at
  gamma / zero vector potential (the WP regime). Documented in laplacian.hpp.

### NOT DONE / PENDING (do next)
1. **COMPILE.** CPU configure running in background → `inq-study/build-cpu`
   (`/tmp/inqstudy_cpu_configure.log`). Next: build the `operations/laplacian` and
   `hamiltonian/ks_hamiltonian` unit-test targets; fix any errors. **No correctness
   claim until this passes.**
2. **Tier 0** (bit-for-bit inert-when-off), **Tier 1** (plane-wave kinetic
   eigenvalue `k²/2m`, electrons-unchanged, ledger consistency) — cheap, run first.
3. **Tier 2/3/4** GPU sims (σ(t) spreading, group velocity, muon-vs-electron
   jellium, MPI partition, GPU-vs-CPU) — `simulation-validation`, user-approved.
4. Catalogue rows in `docs/validation/test-catalogue.md`; formula-validation
   subagent for the spreading law + plane-wave oracles.

### Compile log (2026-07-06)
- CPU build `inq-study/build-cpu` (ENABLE_CUDA=OFF) configured OK (deps on disk,
  no network). Standalone Tier-1 test `tests/muon_mass_fork.cpp` added (plane-wave
  oracle + electrons-unchanged), built via `ctest -R muon_mass_fork`.
- **Compile error #1 (FIXED):** `inverse_mass_.fill(1.0)` on a 2D `gpu::array` is
  ill-formed in this `boost::multi` (treats the scalar as a per-row range →
  `begin(double)`). Fix: value-constructor
  `inverse_mass_ = gpu::array<double,2>(occupations_.extensions(), 1.0)` at both
  alloc sites in electrons.hpp. Rebuild in flight.

- **Rebuild after fix: BUILD OK.** `tests/muon_mass_fork.cpp` compiles and links.

### Validation status (2026-07-06)
- **Tier-1 PASSED** (CPU, `ctest -R muon_mass_fork`, 0.65 s):
  - `laplacian_states`: per-state factor `factor[ist]·(-|k|²)·ψ` applied exactly
    (electron −0.5, muon −0.5/206.77) — diff < 1e-8.
  - `laplacian_add_states` (fourier in-place, ks_ham:235 path): 2× accumulation OK.
  - **electrons-unchanged**: electron states bit-identical with/without the muon
    slot present — no leakage.
- **STILL UNVERIFIED:** Tier 0 (bit-for-bit inert-when-off — needs a full GS/RT
  compare), Tier 1.2/1.5 (expectation-value + ledger consistency — needs an
  orbital_set / ks_hamiltonian engine test), all Tier 2/3/4 GPU sims.
- **NOT built with CUDA yet** — CPU only. A GPU build is required before any
  physics run (validation-gates: GPU default). The `_states` kernels use the same
  `gpu::run`/`GPU_LAMBDA`/`begin(array)` idioms as the scalar versions, so GPU
  compilation is expected to work, but is UNCONFIRMED.

### Next steps (priority order)
1. GPU build of inq-study (mirror inq/build: cuda-12.5 nvcc, arch 80) → confirm
   `_states` kernels compile under nvcc; rerun Tier-1 on GPU.
2. Tier-1.2/1.5 engine test (ks_hamiltonian `kinetic_expectation_value` on an
   orbital_set with a muon slot → ⟨T⟩=k²/2m per state; ledger).
3. Tier-0 bit-for-bit: an electron GS + short RT with the fork present but all
   mass=1 vs a pristine `inq/` build → identical energies.
4. Tier-2/3 GPU physics (σ(t) spreading, muon-vs-electron jellium) — user-approved.

## 2026-07-06 — Campaign authored + orchestrator verified
- Campaign: `docs/campaigns/muon_mass_fork/muon_mass_fork.md` (status: draft; 7
  phase-gated tasks). INDEX regenerated (27 campaigns).
- Orchestrator: `ResearchProject/systems/localised_jellium/scripts/muon_mass_fork/orchestrate.py`.
  **Engine VERIFIED** (no GPU build / no real emails): parses; nvcc path
  `/lsc/opt/cuda-12.5/bin/nvcc` present; `inqview.email.send_run_email` signature
  matches; GPU probe runs (real `cudaMemGetInfo` OK despite NVML breakage). Logic
  tests (stubbed) all pass: all-pass, strict-stop-on-fail, checkpoint-pause,
  resume-past-checkpoint-via `muon_xc_pick.json`, resume-skip-done, single-phase.
- **Bug fixed during check:** Phase-4 checkpoint would loop forever → now gated on
  `muon_xc_pick.json` (user writes it to resume into Phase 5).
- **GPU occupancy (flag):** both GPUs occupied by other users — GPU 0 ~4 GB free,
  GPU 1 ~0 GB. A real run should target GPU 0 (GPU=0) and may be tight.
- **Executable scope:** Phase 1 body (GPU build + `muon_mass_fork` ctest) is
  complete and runnable — but the GPU compile of the `_states` kernels is
  UNVERIFIED (CPU-only so far). Phases 2/3/3b/5 correctly BLOCK (missing run
  scripts: vacuum-WP `run.cpp`, regression harness, muon-jellium `run.cpp`) — the
  next build-out. So the orchestrator runs cleanly and blocks-by-design, NOT a
  full end-to-end autonomous run yet.

## 2026-07-07 — GPU build attempted; kernel bug found + fixed
- Orchestrator launched autonomously (auto-picked GPU 0, 23 GB free). Phase 1 GPU
  build ran; the nvcc compile of `muon_mass_fork.cpp` FAILED.
- **Diagnosis (real, GPU-only bug):** `laplacian_expectation_value_states`'s reduce
  lambda captures a device array-iterator (`fac = begin(factors)`) and dereferences
  it (`fac[ist]`) in the RETURNED expression. The crude 4-D reduce path
  (`external_libs/gpurun/include/gpu/reduce.hpp:376`) wraps the kernel and deduces
  its return type via `decltype(kernel(...))` in HOST code → trips nvcc's
  "extended __device__ lambda return type queried in host code" static assertion.
  The scalar `laplacian_expectation_value` escapes this only because all its
  captures are plain doubles. (NOT a missing CUB fix — inq/inq-study reduce.hpp are
  byte-identical and both carry `cuda::proclaim_return_type`.)
- **Fix applied:** explicit `-> double` on the reduce lambda in
  `laplacian_expectation_value_states` (+ `double lapl`). Rebuild + GPU Tier-1 in
  flight. The apply variants (`laplacian_states`, `laplacian_add_states`) are plain
  element-wise kernels (no reduce, no return-type query) → unaffected.
- **Orchestrator diagnostics fixed:** phase1 now merges stderr into stdout
  (compiler errors go to stderr; the first launch logged only stdout and hid the
  real error).
- CPU Tier-1 still green; this was a GPU-compile-only issue — exactly the kind the
  Phase-1 GPU gate exists to catch before any physics run.
- **`-> double` was insufficient.** Deeper cause: `ks_hamiltonian::kinetic_expectation_value`
  compiles BOTH branches (runtime `if`, not `if constexpr`), so EVERY energy
  computation instantiates `laplacian_expectation_value_states` → the broken reduce
  broke the ENTIRE engine's GPU build (confirmed: existing test `bo`, dipole,
  singularity_correction all fail the same assertion — not just my test).
- **Better fix (applied):** kinetic energy is LINEAR in the prefactor, so
  `laplacian_expectation_value_states` now calls the GPU-proven SCALAR
  `laplacian_expectation_value(ff,-0.5,..)` then rescales each per-state value by
  `-2*factors[i]` (= inverse_mass[i]) with a PLAIN element-wise kernel (no reduce,
  no device-iterator return-type query). Rebuild of `bo` + `muon_mass_fork` pending.

## 2026-07-07 — ROOT CAUSE: wrong CUDA version (12.5 vs 12.6.2)
- After the scale-based `_states` fix, `laplacian.hpp` is GONE from the GPU error
  trace — my fork code is clean. Remaining failures were in UNMODIFIED INQ
  (`observables/dipole.hpp:27`, `hamiltonian/singularity_correction.hpp:74`).
- inq vs inq-study for those files + `gpu/run.hpp` + `gpu/reduce.hpp`: **byte-
  identical** (no drift). inq/install binary dated 2026-03-27.
- **THE BUG WAS MINE (orchestrator config):** NVCC was `/lsc/opt/cuda-12.5/bin/nvcc`
  — copied from the STALE, never-completed `inq/build` cache. The version that
  actually compiles every working GPU run is **`/lsc/opt/cuda-12.6.2/bin/nvcc`**
  (`shared/config.sh` `INQ_CUDA_COMPILER`, used by `inq-run`). 12.5 trips the
  nvcc extended-__device__-lambda return-type assertion on INQ's reduce lambdas;
  12.6.2 does not (every dipole/energy run compiles under it).
- **Fix:** orchestrator `NVCC` → 12.6.2. Wiped `build-gpu` (12.5-configured),
  reconfiguring + building `bo` + `muon_mass_fork` under 12.6.2 (fresh, ~20 min).
- **Two real fixes stand regardless of CUDA version:** (1) `laplacian_expectation_value_states`
  rewritten to reuse the scalar reduce + elementwise rescale (no device-iterator
  reduce); (2) orchestrator phase1 now logs stderr; (3) orchestrator GPU auto-select.

## 2026-07-07 — PHASE 1 GREEN ON GPU
- Fresh `build-gpu` under cuda-12.6.2: `build exit 0`, **0 static assertions**.
  `bo` (unmodified INQ) BUILT → engine GPU-compiles; my `_states` rewrite didn't
  break the engine. `muon_mass_fork` BUILT + **GPU Tier-1 PASSED** (1.65 s) — the
  per-state mass fork is validated on real hardware (nvcc/sm_80), not just CPU.
- Orchestrator `state.json`: `done: [phase1]`. GPU build dir: `inq-study/build-gpu`
  (12.6.2). Reusable — Phase-1 rebuilds are now incremental/fast.
- **Phase-1 status: COMPLETE (kernel-level).** Still-open Phase-1 TODOs (nice-to-have
  before physics, per plan): expectation/ledger engine test, wrong-slot, GPU-vs-CPU,
  MPI-partition — the orchestrator's phase1 has a TODO marker for these.
- **NEXT (Phase 2):** write the vacuum-WP `run.cpp` (HERE/vacuum_wp/run) — the σ(t)
  spreading sim — now that it can link against the GPU-built engine. Then phase2
  dispatch is unblocked.

## 2026-07-07 — DUAL-GPU AUTONOMOUS RUN STARTED (Phase 2 build)
- User: "two free GPUs now … peruse both … start running all phases one after
  the other." Confirmed both idle via the cudaMemGetInfo probe: GPU 0 ≈ 23.1 GB,
  GPU 1 ≈ 23.8 GB free (NVML/nvidia-smi still broken — probe is the source of truth).
- **Orchestrator made dual-GPU** (`orchestrate.py`): `resolve_gpu()` now builds a
  POOL `GPUS` of every device ≥4 GB free (freest first) instead of one; new
  `run_sims_parallel(jobs, label)` dispatches independent runs across the pool,
  ≤1 job per GPU, idempotent-skip, backfill as GPUs free. `run_sim` gained a `gpu=`
  arg. Single-GPU path preserved when `GPU=<n>` is pinned.
- **Phase-2 vacuum-WP run written**: `scripts/muon_mass_fork/vacuum_wp/run.cpp`.
  Empty cubic box, one injected Gaussian WP in the single extra_state, per-state
  `inverse_mass` set on the WP orbital, propagated under `non_interacting`
  (H = −∇²/2m; no self-Hartree/XC → exact free-Gaussian law). Records
  `wp_real_space_stats` (σ_ρ(t), centroid), `energy_kinetic` (⟨T⟩), norm; boundary
  guard prints when 4σ_ρ ≥ L/2. Env: WP_OUT/L/SPACING/DT/TSTEPS/SIGMA/K0/INV_MASS/
  WRITE_EVERY. Builds against inq-study via inq-run (INQ_SOURCE=inq-study,
  cuda-12.6.2); binary lands at `vacuum_wp/run` (orchestrator's expected path).
- **Oracle checker written**: `vacuum_wp/check_oracle.py` — validates σ_z2(0)=σ_WP²/2,
  parabola fit σ_ρ²=σ_ρ0²+b·t² → m_fit (tol 5%, relaxed 25% for barely-spreading
  muon over a short window), v_group=k0/m (2%), norm drift (<1e-3), ⟨T⟩ drift (<5e-3).
- **Smoke test in flight**: electron σ_WP=0.5, 150 steps (t=3, stays in box). First
  build compiles the full inq-study subset (~15–20 min, one-time); recompiles fast.
  NOT yet validated — awaiting build+run completion before claiming Phase 2 works.
- **STILL TODO before phase2 launch:** add optional density-VTI output to run.cpp
  (WP_EMIT_VTI) for the xz-density-vs-σ deliverable; finalise phase2 job list
  (spread elec/muon, vgroup elec/muon, mass-dial m=10, xz σ∈{0.5,1,2,4}); wire
  `run_sims_parallel` + `check_oracle.py` into `phase2()`.

## 2026-07-07 — PHASE 2 VALIDATED + LAUNCHED ON BOTH GPUs
- **Vacuum run.cpp fixed** (two bugs found by running): (1) explicit
  `input::environment{}` double-inits MPI → removed (INQ lazily inits via
  `environment::global()`); use `electrons.root()`. (2) `extra_electrons(0)` +
  no ions → INQ throws "system does not have any electrons" (electrons.hpp:239);
  fix = **1 spectator electron** + `ground_state::initial_guess`. Under
  `non_interacting` the spectator is a decoupled free particle (no Hartree/XC
  coupling) → does NOT affect the WP; every WP observable reads only wp_idx.
- **KE oracle sourced from `wp_momentum_stats`** (⟨k²⟩ conservation, mass-
  independent) not total energy (which would include the spectator). Added
  `wp_momentum_stats` + optional density VTI (WP_EMIT_VTI) for the xz deliverable.
- **BOTH FORK PATHS VALIDATED to ~1e-5 (GPU, cuda-12.6.2):**
  - **Electron** (mass 1 → scalar guard path): mass_fit=0.9999983 (rel 1.7e-6),
    v_group=0.5 exact, σ_z²(0)=0.125, ⟨k²⟩ & norm conserved. (also a mini
    bit-for-bit of the mass-1 path.)
  - **Muon** (mass 206.77 → `_states` fork path): **mass_fit=206.773 (rel 1.3e-5)**,
    v_group=0.002418=k0/206.77, ⟨k²⟩ drift=0.0, norm conserved. The fork scales
    the kinetic operator by 1/m per state, correctly, on real hardware.
  - Oracle checker: `vacuum_wp/check_oracle.py` (σ_z2(0), parabola m-fit, v_group,
    norm, ⟨k²⟩, ⟨k_z⟩). Both runs PASS all oracles.
- **Phase 2 running autonomously across GPU 0+1** (10:00): orchestrator
  `run_sims_parallel` dispatches 9 vacuum runs (spread elec/m10/muon, vgroup
  elec/muon, xz σ∈{0.5,1,2,4}) across both GPUs, then gates on check_oracle.
  Muon runs use L=24 (60³, ~0.07 s/step) to stay cheap at 6000 steps. spread_elec
  done in 94 s; muon runs (~7 min each) in flight. Log:
  `scripts/muon_mass_fork/orchestrate_phase2.log`.
- **Phase 3 regression harness written + BUILDING** (parallel to phase2):
  `regression/run.cpp` = He-atom LDA GS + kicked short RT, NEVER touches
  inverse_mass (mass-1 → scalar path). Built against inq-study (`fork/`) AND
  pristine inq (`pristine/`); `compare_regression.py` diffs GS+RT energies +
  density (etol 1e-9). Both first-builds ~15–20 min, in flight.

## 2026-07-07 — PHASE 3 BIT-FOR-BIT PASS + PHASE 4 RESEARCH DONE
- **Phase 3 HARD TRUST GATE PASSED.** He-atom LDA GS+kicked RT built against the
  fork (inq-study, mass untouched=1) vs pristine inq: GS energies (total/kinetic/
  hartree/external/xc) **exactly identical to 14 digits (|d|=0.0)**; RT trace (41
  steps) total/external agree to 9.99e-15, kinetic/hartree/xc exactly 0.0; GS
  density (5832 pts) |d|=0.0. `compare_regression.py` → BIT-FOR-BIT PASS. The fork
  is provably inert when off → muon physics is attributable to the mass, not an
  engine edit. Harness: `regression/{run.cpp, fork/, pristine/, compare_regression.py}`.
- **Phase 3b wired** (orchestrator): muon WP + electron control under FULL LDA
  (interacting) via the new `WP_THEORY=lda` toggle in vacuum_wp/run.cpp — tests the
  forked kinetic path inside the Hartree+XC propagator (Phase 2 was
  non_interacting). Pass = no NaN/inf in the WP momentum trace. NEEDS a vacuum_wp
  rebuild (source edited; rebuild only AFTER phase2 releases the binary).
- **Phase 4 research COMPLETE (checkpoint pending user pick).** 3 grounded source
  notes: `docs/sources/{heg-mass-scaling-xc, kreibich-gross-multicomponent-dft,
  car-parrinello-fictitious-mass}.md`. Candidate summary +
  recommendation: `docs/campaigns/muon_mass_fork/phase4_muon_xc_candidates.md`.
  - **Physics:** one-component all-muon jellium is EXACTLY mass-scalable — exchange
    is mass-independent, correlation carries the mass as ε_xc(M·r_s). At physical
    r_s=5.69 the muon effective r_s≈1177 (near-Wigner), where electron-LDA is a long
    extrapolation → correlation-dominated, LDA-suspect regime (the campaign's point).
  - **Candidates:** A = mass-rescaled LDA (exact for the muon HEG; wrapper feeding
    r_s^μ=m·r_s and ×m energy) — RECOMMENDED; B = naive electron-LDA at physical r_s
    (control/baseline, INQ default); C = multicomponent/NEO (only if muon is a
    distinct species — deferred). CP fictitious mass explicitly EXCLUDED (numerical,
    not XC).
  - **Recommended Phase 5:** run all-muon r_s=5.69 twice — B (naive LDA) vs A
    (mass-rescaled). Pick via `muon_xc_pick.json`.

## 2026-07-07 — PHASE 2 PASSED (all oracles) + NOTEBOOKS BUILT
- **Phase 2 fully green.** All 9 vacuum runs completed across both GPUs (10:00–10:28);
  all 5 physics runs pass EVERY analytic oracle (re-run after fixing the oracle
  path bug `results/<lab>/` — runs skip idempotently):
  - mass_fit: spread_elec **0.999998 (2.1e-6)**, spread_m10 **9.99995 (4.6e-6)**,
    spread_muon **206.769 (4.9e-6)** — mass-dial continuous over 2 decades.
  - v_group: vgroup_elec 0.5 (8.6e-6), vgroup_muon 0.0024181 (4.5e-6) = k0/206.77.
  - ⟨k²⟩ drift ≤ 2.5e-11 (muon exactly 0), norm 1.2e-6, σ_z²(0)=0.125 all runs.
  - Success email sent; state.json: phase1+phase2 done.
- **Oracle path bug fixed** in orchestrate.py phase2: check_oracle rd is
  `RUNS/phase2/<lab>/results/<lab>` (run.cpp writes to results/<WP_OUT>/).
- **Notebooks built + executed** (venv kernel, canonical theme) under
  `hypotheses/muon_mass_fork/`: `index.ipynb` (read-order guide),
  `phase2_physics.ipynb` (σ(t) overlays + mass-dial + v_group + xz-density-vs-σ,
  4 PNGs), `phase3_regression.ipynb` (bit-for-bit table), `phase4_xc_research.ipynb`
  (candidate summary). Builder: `hypotheses/muon_mass_fork/build_notebooks.py`.
- **vacuum_wp rebuilt** with the WP_THEORY toggle (for phase3b) — in flight.

## 2026-07-07 — AUTONOMOUS HORIZON COMPLETE (Phases 1→4-checkpoint), PAUSED
- Full `orchestrate.py` ran phase1/2 (skip, done) → **phase3 BIT-FOR-BIT PASS** →
  **phase3b 2/2 PASS** (muon + electron WP under full LDA, no NaN — the forked
  kinetic path is stable in the interacting Hartree+XC propagator; the sanity Phase
  2's non_interacting runs could not test) → **phase4 CHECKPOINT** (emailed the user
  to pick the muon-XC functional; paused cleanly). Both GPUs used throughout.
- **state.json:** done=[phase1,phase2,phase3]; phase3b passed (its runs are the gate,
  no separate done-flag needed as phase4 is the pause). Campaign frontmatter: phases
  1,2,3,3b + index = done; status=paused; blocked_reason = awaiting muon-XC pick.
- **TO RESUME PHASE 5:** copy `scripts/muon_mass_fork/muon_xc_pick.example.json`
  → `muon_xc_pick.json` (edit if desired) and re-run `orchestrate.py`. Recommended
  pick: A = mass-rescaled LDA vs B = naive electron-LDA baseline (see
  `docs/campaigns/muon_mass_fork/phase4_muon_xc_candidates.md`).
- **PHASE 5 NOT YET BUILT (correctly gated on the pick).** Needs: (1) all-muon GS
  run.cpp for the r_s=5.69 bath (N=162, L=50 cubic, dx=0.40) with global
  inverse_mass=1/206.77 via the calculator.hpp mass path; (2) incident muon WP run
  under B (stock LDA) AND — if A is picked — a NEW inqkit mass-rescaled-LDA wrapper
  (feed r_s^μ=m·r_s, ×m energy); (3) S/wake/ledger extraction vs the electron
  r_s=5.69 reference, SIE floor from a vacuum-WP control; (4) phase5 notebook.

## REMAINING autonomous flow (this session)
1. Phase 2 finishing (m10 long pole ~120³); orchestrator gates on check_oracle +
   emails. Already MANUALLY validated: electron mass_fit rel 1.7e-6, muon 1.3e-5.
2. Rebuild vacuum_wp/run (WP_THEORY toggle) once phase2 releases the binary.
3. Run full `orchestrate.py` → phase3 (re-run skips, compares → PASS), phase3b
   (muon-LDA sanity), phase4 (emails checkpoint, writes blocked, PAUSES for pick).
4. Build notebooks: phase1–4 + index (phase5 after user pick).
5. Phase 5 (gated on muon_xc_pick.json): all-muon GS + incident muon WP, A vs B.

## Open design questions (grill, unresolved)
- Q3: muon-in-electron-jellium distinguishability — accept+bound SIE (a) vs own
  spin channel (b). Recorded, not decided.
- Q4: all-muon jellium regime — fixed r_s (rescaled/validation) vs fixed physical
  density (strongly-correlated, genuine new physics). User wants it as a real
  physics phase; rescaling caveat (S∝m²) still to be reconciled. See engine notes
  §3 "CRITICAL rescaling fact".

## 2026-07-07 — NEW Phase 4 designed: momentum-matched effective-mass projectile

New phase requested BEFORE the XC-research phase (which stays a subagent + user
pick). Goal: an LDA RT run of a projectile WP on the r_s=5.69 localised jellium
slab, DIRECTLY comparable to the existing 100 eV electron run, WITHOUT a 10×
finer grid and WITHOUT appreciable spreading.

**Key physics settled (all in the spec notebook + broadening plot):**
- The 3 projectile knobs are σ_WP, E, m. Spreading time τ_s = m·σ_WP². Momentum
  k₀ = √(2mE) = m·v; grid needs k₀+3σ_p ≤ π/dx, σ_p = 1/(√2 σ_WP).
- **Match INITIAL MOMENTUM (not energy, not velocity)** to a reference electron:
  same k₀ ⇒ same de Broglie wavelength ⇒ same grid. Matched-velocity needs a 330×
  finer grid; matched-energy (the 300 eV muon) needs 11× — both rejected.
- **Spreading-per-DISTANCE is mass-independent at fixed k₀** (mass cancels). Mass's
  real role: decouple k₀ (grid+spread) from v,E (stopping physics+comparability).
- **σ_WP=0.5 (old runs) CANNOT give <1% spread at impact** at any mass (needs
  D_launch<0.18 Bohr). Strict <1% needs σ_WP≥2.
- **The literal muon (m=206.77) is RULED OUT** by (T<80 au + grid ≤1.5×): it needs
  k₀=561 (60× grid) or ~5500 au. Constraints select an effective mass ≈ 3.

**SELECTED config (user 2026-07-07, "more concentrated packet"):** σ_WP=1.0,
m=2.7 mₑ, v=2.71 a.u. (=100 eV electron), **E=269 eV**, k₀=7.3, **dx=0.333 (1.5×,
grid cap)**, reuse slab_n82_L50x50x90 (r_s=5.665, N=82 even), launch z0=−14.6,
dt=0.02, T=77 au (3.5× traversal), LDA/ETRS/γ-only. Spread: **4.1% at impact /
285% (→σ_ρ≈2.7 Bohr) at slab exit** — relaxes the <1% target for a tighter probe;
σ=1 and σ=2 curves CONVERGE inside the slab (~1.5 Bohr), so little lost. Read
stopping from n(k,t) coherent peak (robust to real-space spread).
- **Est. GPU wall: ~15 h (1 GPU, dt=0.02, ~3850 steps)**; calibrated to measured
  qsp_phase4 (3.75 s/step @ dx=0.50, 61 states) × 3.66 (6.08 M-pt grid). ~8 h on
  2 GPUs or dt=0.04. CONFIRM with a ~100-step pilot.

**Artefacts (hypotheses/muon_mass_fork/):** `muon_effmass_spec.ipynb` (scan +
spec + spreading + GPU time), `broadening_estimate.py`/`muon_wp_broadening.png`
(earlier 300 eV muon variant), `build_spec_notebook.py`. Concept note user-owned:
`docs/misc/thoughts/mass-as-a-knob-in-simulations.md`.

**Grid budget relaxed to 1.7× (user 2026-07-07)** → dx_min=0.294, k_max=10.68.
This pulls the clean <1%-at-impact packet down to σ_WP=1.5 (was 2.0 at 1.5×).
Notebook re-scanned; four contenders (all v=2.71, r_s=5.665, dx=0.294, <1% impact):
- A tight: σ_WP=1.25, m=3.3, 331 eV, impact 0.79% BUT 2.3% t=0 slab overlap (not clean).
- **B concentrated ★ (RECOMMENDED): σ_WP=1.5, m=3.4, 342 eV, launch 2.75σ (z0=-15.4),
  impact 0.97%, 0.30% overlap, slab-growth 65% (σ_ρ 1.06→1.75 Bohr).**
- C clean: σ_WP=1.75, m=3.5, 349 eV, impact 0.82%, growth 40%.
- D balanced: σ_WP=2.0, m=3.6, 355 eV, impact 0.61%, growth 25%.
- **GPU: ~22 h (1 GPU, dt=0.02, ~3900 steps) / ~11 h (dt=0.04 or 2 GPUs)** — 1.7× grid
  is ~9 M pts, ≈2× the 1.5× cost. Read stopping from n(k,t) (through-slab <1% impossible).

**LOCKED: contender B** (σ_WP=1.5, m=3.4, 342 eV, v=2.71, dx=0.294, launch z0=-15.4,
r_s=5.665 N=82 slab, LDA/ETRS/γ). Notebook re-centered on B: widths ledger +5 panels
(σ_ρ(t), packet-vs-slab real space, momentum-vs-Nyquist [3σ_p edge exactly at k_max,
0.13% aliased], geometry schematic, B-vs-old-electron). Through-slab growth 65%
(σ_ρ 1.06→1.77 Bohr) → read stopping from n(k,t) coherent peak.

**dt SMOKE-TESTED (2026-07-07) — answer: dt=0.02 for B.** Ran the built regression probe
(He/LDA/ETRS, same propagator) at the production cutoff, swept dt, checked energy trace:
ETRS cliff is H·dt=E_cut·dt ≈ 2.2, consistent across grids. On B's 1.7× grid (E_cut=57 Ha):
dt=0.08 NaN@step3, **dt=0.04 NaN@step9** (my 2.8 estimate was too generous — real cliff ~2.2),
dt=0.03 stable-but-edge (drift 8e-3), dt=0.025 stable (4e-3), **dt=0.02 STABLE safe (1.5e-3
Ha/300 steps)**. Control: dt=0.08 @coarse-20Ha stable (confirms grid sets it). Also measured
the 1.5× grid (D, E_cut=44 Ha): dt=0.05 NaN, dt=0.04 stable.
- **User ladder (0.08→0.04→0.02): both 0.08 & 0.04 fail → dt=0.02.** Wall (1 GPU): B@0.02
  ~22h (~12h on 2 GPUs); B@0.025 ~18h (edge-ish); **D@0.04 ~7h**.
- **TRADEOFF SURFACED:** the concentrated packet B is ~3× the wall-time of D — the 1.7× grid
  both adds points AND forces dt=0.02, vs D's 1.5× grid allowing dt=0.04. Awaiting user call:
  stick with B@0.02 (~22h) or drop to D@0.04 (~7h, σ_ρ,0 1.41 vs 1.06, still <1% impact).
Smoke-test scratch cleaned; regression probe reused (no new binary). Probe is a proxy (He,
same E_cut/propagator) — the production GS run's first steps re-confirm before committing.

## 2026-07-07 — effmass_pair BUILT + LAUNCHED (quantum WP vs classical)

User locked **contender D** (σ_WP=2) at **dt=0.04**, and asked to run BOTH the quantum
effective-mass WP and its CLASSICAL analogue autonomously on the two GPUs, with run + a
model-comparison notebook auto-built. dt=0.04 forces the **1.5× grid** (dt=0.04 diverges on
the 1.7× grid), so the actual config is:
- **σ_WP=2.0, dx=0.333 (E_cut≈44 Ha=88.8 Ry), k0=8.364, m_eff=3.0852 m_e, v=2.7111 a.u.
  (=100 eV e-), E=309 eV, inverse_mass=0.324127, launch z0=-16.743 (3σ), dt=0.04, N_steps=2000
  (T≈80 au), CAP η=-0.7 at ±35..±45, LDA/ETRS/γ, r_s=5.665 N=82 slab.** Bath = mass-1 electrons.
- Momentum-matched: classical particle m=3.0852, v=2.7111 → same k0=8.36 as the WP.

**Files (all NEW, under scripts/muon_mass_fork/effmass_pair/):**
- `gs/run.cpp` — bare-jellium GS at dx=0.333 → saves shared_gs/slab_n82_L50x50x90_dx0p333.
- `quantum/run.cpp` — WP (σ=2, inverse_mass fork) injected into last extra state; WP real-space
  + momentum (n(k,t)) stats, density VTIs, energy ledger. Loads the shared GS.
- `classical/run.cpp` — copied from qsp_phase4/classical, retargeted: UPF wpsigma2p0, mass
  3.0852/1822.8885 amu, v=2.7111, spacing 0.333, launch -16.743, dt 0.04, GS dir dx0p333.
- `orchestrate.py` — polls for GS → build-gates both binaries → launches quantum(GPU0) +
  classical(GPU1 if ≥18GB free else SEQUENTIAL on GPU0) → build_notebooks.py → emails milestones.
- `build_notebooks.py` — comparison ("phase") notebook effmass_pair_comparison.ipynb +
  per-run notebooks + effmass_pair_stopping.png, into hypotheses/muon_mass_fork/. Defensive.
- **UPF `shared/pseudopotentials/electron_gaussian_wpsigma2p0.upf`** — generated via
  generate_gaussian_psp(2.0), DATA-VERIFIED: σ_charge=1.414=σ_WP/√2, V(0)=0.564 Ha, tail→+1/r.

**RUNNING NOW (2026-07-07 ~16:38):** GS pid 1180798 (GPU0, dx=0.333, SCF converging ~step60).
Orchestrator pid 1217967 (polling GS checkpoint, up to 4h). On GS-ready it builds+launches.

**Compile fixes already applied:** quantum needed `io/observables_writer.hpp` (not observables.hpp)
+ `observables/density_delta.hpp` + `real_time/step_context.hpp`; and its CAP `pert` ternary had
mismatched types → rewritten as single sum with η=0 when CAP off. Classical include paths deepened
+1 level (moved dir). **Both RT binaries' compiles are GATE-VERIFIED by the orchestrator**
(build_only → email+STOP if a binary fails to build) — NOT yet independently confirmed by me.

**KNOWN RISKS / watch:** (1) GPU1 had only ~5.3 GB free (another user) — orchestrator will run the
pair SEQUENTIALLY on GPU0 (~14h) unless GPU1 frees to ≥18 GB. (2) dt=0.04 on the 1.5× grid drifts
~7e-3 Ha/200 steps in the He probe (bounded) — watch the quantum energy trace. (3) classical compile
+ GS-load-with-refined-grid unverified until the orchestrator runs. (4) WP carries a ~few-eV SIE the
classical lacks — bound with the vacuum-WP control before over-reading small quantum-vs-classical gaps.

**Earlier open item (300 eV electron twin) deferred** — the quantum/classical pair is the comparison now.

The notebook `hypotheses/muon_mass_fork/muon_effmass_spec.ipynb` (spec, scan, dt smoke tests,
GPU time) still reflects the σ-scan + B/D analysis; its contender table shows D at the 1.7× grid
(k0=9.62) — the RUN correctly uses the 1.5×-grid D (k0=8.36, dt=0.04). Minor doc inconsistency noted.

---

## 2026-07-09 — σ=1 mass-only run + vacuum CAP calibration (NEW sweep `sigma1_massonly`)

**Why:** user rejected all increased-mass σ=2 runs (effmass_sigma2, effmass_12h) — CAP
**reflection** + packet **too wide**. New ask: a CONCENTRATED σ=1 projectile, mass-only,
100 eV, ≤1–2 h; second GPU runs a vacuum calibration of CAP + spreading.

**Sanity-check verdict (grounded, presented to user):** a σ=1 packet CANNOT be held rigid
across the 25-Bohr slab by mass tuning — fractional spread `√(1+(d/k₀σ²)²)` improves only as
1/√m, and m is capped by aliasing (≤5.5 at dx=0.4) AND cost (heavier=slower=more steps). At the
cost-feasible mass it is ~3.4× wider at the slab centre than the classical Gaussian. **User chose
(knowingly) the mass-only dispersive end-member** over chirp-focus/σ≈1.5, and **vacuum calibration**
for the 2nd GPU. cutoff_guard (fork-corrected E_eff=E·m=345 eV): PASS, tail 0.00%.

**Locked config (slab run):** σ_WP=1, E=100 eV, **m=3.45 (INV_MASS=0.289855)**, **k₀=5.0356**
(=√(2·E·m), fork-corrected), v=1.460, dx=0.40, box **50×50×64** (4σ standoffs shrink L_z 90→64),
slab 25 Bohr (N=82, r_s≈5.67), launch z₀=−16.5, **CAP η=−0.6** (gentler than σ=2's −1.0) sin² band
peak ±26.25 width 11.5 → band [±20.5,±32], dt=0.05, **N=1192** (3× traversal), write_every=4,
ckpt_every=300 (3 interior ckpts + resume machinery, cloned from effmass_sigma2/wp).

**Files created:**
- `shared/configs/slab_n82_L50x50x64.hpp` (clone of L90, LZ=64, dx=0.40)
- `scripts/muon_mass_fork/sigma1_massonly/gs/run.cpp` → `shared_gs/slab_n82_L50x50x64_dx0p40`
- `scripts/muon_mass_fork/sigma1_massonly/wp/run.cpp` (slab WP, checkpoint/resume, graded CAP)
- `scripts/muon_mass_fork/sigma1_massonly/vacuum/run.cpp` (CAP+spreading calibration, no slab)
- `docs/plans/localised-jellium-wp-sigma1-100eV-massonly.md`

**GOTCHA recorded:** `cutoff_guard.py` computes p0=√(2E) assuming **m=1** — for the mass fork feed
**E_eff = E·m** so p0=√(2·E·m)=k₀. Fork-uncorrected guard hugely under-reports aliasing.

**STATUS (2026-07-09 ~12:35):** both compiled clean against inq-study; all gates GREEN.
- **Vacuum calibration** (GPU0, `sigma1_massonly/vacuum`, non_interacting): **COMPLETE — PASS.**
  Fork VERIFIED (WP v=1.459=k₀/m ✓). (1) CAP absorption CLEAN: norm 1.000→**0.00016 (99.98%)**,
  monotonic, no plateau; reflected norm (z<0) peaks ≤**0.25%** then decreases → not a coherent
  escaping lobe. η=−0.6 fixes the σ=2 reflection complaint. (2) Free spreading matches the oracle
  to **ratio 1.000** (σ_z 0.707→3.154 at t=0→15, norm=1) — validates m=3.45 in the DYNAMICS and
  confirms the dispersion (4.5× by t=15). Wall 23 min.
- **GS** (GPU1): COMPLETE. E=−39.398 Ha (vs L90 −72.498 — the ~33 Ha gap is the harmless
  box-size electrostatic offset; cancels in ΔE). Density ISOLATED: n/n0=0.0008 at CAP inner edge
  (±20.5), 0.0000 at box edge → CAP won't drain the GS. Saved `shared_gs/slab_n82_L50x50x64_dx0p40`.
- **Slab run** (GPU1, `sigma1_massonly/wp` → `results/quantum`): **COMPLETE** 14:01.
  wall **5351 s = 1.49 h** (within budget ✓), N=1192, 3 ckpts (step 300/600/900). WP norm
  1.0→**0.0007 (99.93% absorbed)** — clean, matching the vacuum calibration. σ=1, m=3.45,
  E=100.0 eV, v=1.460, CAP η=−0.6. Run-notebook building → `hypotheses/muon_mass_fork/
  sigma1_massonly_wp_run.ipynb`.

**Both GPUs free (14:01).** Matched classical twin still offered; **awaiting user go** (held).

**Open / next:** matched classical rigid-Gaussian twin (σ_pot=0.71, m=3.45, UPF
`electron_gaussian_wpsigma1p0.upf` — NOT yet generated) on GPU0 for the WP-vs-classical pair;
run-notebook for the slab run once complete (`hypotheses/muon_mass_fork/`, run-notebook skill).
Light-projectile S-extraction (initial-drag window, packet decelerates) + phantom-absorbed-orbital
caveat apply. NOTE: slab-face→CAP gap is only 8 Bohr — watch the wake analysis for CAP draining the
dynamic screening response (GS density there is ~0, but plasmon flux may reach it).
