# Handover — wavepacket twin of the high-density classical S(v) benchmark

Plan: `/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/docs/plans/wavepacket-highdensity-sv-twin.md`
Parent campaign: `docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md`
(id `classical-highdensity-sv`), handover `docs/handovers/classical-highdensity-sv-benchmark.md`.
Branch `quantum-stopping-power`. Machine: **CSD3**, `ampere` partition,
account `mphil-nikiforakis-skcb2-sl2-gpu`. Started 2026-07-30.

## Goal (one line)

Re-run the classical high-density S(v) benchmark with an electron **wavepacket**
projectile instead of a classical Gaussian charge, every other physical parameter
held fixed, and extract S with the KS-orbital definitions
(`docs/plans/bulk-jellium-ks-stopping.md` §4).

---

## STATUS 2026-07-31 — σ = 2 and σ = 3 campaigns launched (σ sweep)

**User instruction (2026-07-31).** Run the same four velocities at σ_WP = 2 and
3 Bohr, chained (one campaign's four points in parallel, then the next), and
produce an S(v) plot by the localised-jellium **deposit** definition
`S = (E_total(t_final) − E_GS)/L_slab_z` carrying traces for σ = 0.5, 2 and 3.

**Why.** At σ_WP = 0.5 the packet disperses at 1/(√2σ) = 1.414 Bohr/a.u. Launched
at z = −24 with the slab face at −12.5, it is **4.7–8.1 Bohr wide on arrival** and
**14.8–25.8 Bohr on departure** — at v = 2.0 it leaves wider than the 25 Bohr slab
it just crossed. It is not a localised projectile while it does the physics. At
σ = 2 / 3 the rate falls to 0.354 / 0.236 Bohr/a.u.; the packet arrives at
1.8–2.5 / 2.3–2.5 Bohr. See `sigma_sweep_widths.png`.

**Ceiling on σ, and why it is 2–3 and not 4.** The minimum-final-width choice is
σ\* = √t_out ≈ 3.2–4.3, which would give exactly √2 growth over the whole path.
It is NOT reachable with the launch point and CAPs held fixed: the −z CAP inner
edge is z = −30, so clearance is 6 Bohr = 4.2 density-std at σ = 2, 2.8 at σ = 3,
2.1 at σ = 4 (1.8 % of the packet inside the absorber at t = 0). σ = 3 is the
practical ceiling and already starts with ≈0.23 % inside the CAP — a one-off loss,
reproduced by `vac_s3p0_*`, flagged in the notebooks.

### Chain (submitted 2026-07-31 03:30, `shared/bin/submit-wp-hd-sigma-sweep.sh`)

| # | stage | job | state at write |
|---|---|---|---|
| 1 | smoke σ=2 (t=0 gates) | 32439180 | **DONE — ALL GATES PASSED** |
| 2 | sweep σ=2, array 0–3 | 32439181 | **RUNNING** (4 GPUs, parallel) |
| 3 | vac σ=2 (CAP baselines) | 32439182 | pending `afterany` 2 |
| 4 | smoke σ=3 | 32439183 | pending `afterany` 3 |
| 5 | sweep σ=3, array 0–3 | 32439184 | pending `afterok` 4 |
| 6 | vac σ=3 | 32439185 | pending `afterany` 5 |
| 7 | notebooks + σ sweep | 32439186 | pending `afterany` 6 |

Expected ~7 h end to end (2.75 s/step measured at σ=2; v=2.0 is the long pole at
3623 steps ≈ 2.8 h per campaign). No GS rebuild: σ does not enter the bath, so
both campaigns load the existing `shared_gs/slab_n100_L35x35x85_dx0p4_per2`.

### σ = 2 smoke gate results (job 32439180) — verified, not assumed

| gate | σ=0.5 (previous) | σ=2.0 |
|---|---|---|
| σ_pz² deviation | ~1.6 % (1.25 grid pts/σ) | **2.8e-5 %** (5 grid pts/σ) |
| z-momentum weight past Nyquist | measurable | **7e-60 %** |
| T₁−T₂ | 81.63 eV | **5.10 eV** (= 3/(4σ²), as predicted) |
| max overlap w/ occupied manifold | — | **9.1e-5** (want <1e-3) |

That last row was a genuine risk — a 4× broader packet has far more room to
overlap the occupied bath states — and it is fine.

### What changed in code (no engine changes; σ was already a knob)

`LJ_SIGMA` (`wp/run.cpp:167`) and `WPC_SIGMA` (`cap_check/run.cpp:112`) already
existed and every t=0 gate is derived from `SIGMA_WP`, so **no C++ edit and no
re-validation of the binary was needed**. Changes are dispatcher + analysis:

- `shared/bin/run-wp-hd-wp.slurm`, `run-wp-hd-vac.slurm` — accept `LJ_SIGMA`,
  prefix run names `s2p0_` / `s3p0_`. σ = 0.5 keeps BARE names, so the completed
  runs, their notebooks and `wp_S_summary.csv` resolve unchanged.
- `shared/bin/submit-wp-hd-sigma-sweep.sh` — NEW, the chain above.
- `run-wp-hd-notebooks.slurm` — walltime 4 h → 12 h (12 run notebooks now).
- `hypotheses/.../wp_hd_stopping.py` — `set_campaign(σ)`, `sigma_tag`,
  `vac_name_for`, `has_campaign`, `transverse_overlap_time`, `aliasing_bias_pct`.
- `hypotheses/.../build_run_notebooks.py` — σ-aware run/synthesis builders plus
  **`build_sigma_sweep()`**, the requested cross-σ deposit figure.

### Two things established by measurement, worth not re-deriving

1. **The fit window is σ-dependent and widens hugely.** Transverse periodic images
   overlap when 6·σ_d(t) = L_xy, i.e. t = σ·√(2(L_xy/6)² − σ²): **4.12 a.u. at
   σ=0.5, 16.0 at σ=2, 23.1 at σ=3**. The σ=0.5 campaign's 4 a.u. window was a
   consequence of the width choice, not a fixed property of the box.
2. **Aliasing vanishes.** σ_p = 1/(√2σ) = 0.354 / 0.236 against k_Nyq = 7.85, so
   the σ=2/3 campaigns have zero moment bias to machine precision at every
   velocity **up to v = 4.5**. v = 4.0/4.5 would be recoverable in these campaigns
   (not requested; the four-point grid is held for comparability).
   `aliasing_bias_pct()` reproduces the recorded σ=0.5 dx=0.4 table exactly; against
   the one MEASURED aliased point it is good to ~7 %/~13 % — enough to sort
   "negligible" from "fatal", NOT a calibrated correction.

### Bug caught live (worth keeping)

`deposit_stopping()` read a **still-propagating** run and returned a perfectly
plausible-looking S: 86 of 3623 steps gave `S_deposit = 2.35 eV/Bohr` with
`norm_final = 1.000`. It now returns `steps_done` / `steps_target` / `complete`,
and the σ-sweep notebook **excludes** incomplete points from the figure while
still listing them. Any future consumer must check `complete`.

### Verified so far / NOT yet verified

- VERIFIED: σ=2 t=0 gates; dispatcher name resolution for all three campaigns;
  `aliasing_bias_pct` against the recorded table; the whole notebook path executed
  end to end (`sigma_sweep.ipynb`, 0 errors, correctly excluding in-flight points).
- NOT yet verified: σ=3 gates (chain stage 4); any σ=2/3 production physics; the
  σ=3 CAP-clearance loss magnitude (predicted 0.23 %, must be read off `vac_s3p0_*`).

### Open, not done

- No classical twin exists at σ = 2 or 3 (benchmark is σ_pot = 0.354 only), so the
  classical curve is orientation-only on those traces. New classical runs would
  need `electron_gaussian_sigma1p41.upf` / `sigma2p12.upf`.
- The bath-internal deposit (ΔE_SS + ΔE_SB from `interactions.csv`), which would be
  genuinely comparable to the classical estimator, is still not built.

---

## STATUS 2026-07-30 — chain launched, fully autonomous, nothing to babysit

SLURM chain submitted by `shared/bin/submit-wp-highdensity-sv.sh`:

| # | stage | job | state at write |
|---|---|---|---|
| 1 | GS dx=0.50 (fidelity check) | 32418534 | **DONE — PASSED** |
| 2 | GS dx=0.40 (production) | 32418535 | running |
| 3 | smoke (build + t=0 gates) | 32418536 | pending `afterok` 2 |
| 4 | sweep, array 0–3 | 32418537 | pending `afterok` 3 |
| 5 | vacuum CAP controls | 32418538 | pending `afterok` 3 |
| 6 | notebooks + synthesis | 32419008 | pending `afterany` 4,5 |

Stage 6 makes the campaign autonomous end to end: notebooks are built AND executed
without further input. Stage 6 uses `afterany`, so one failed velocity still yields
notebooks for the rest.

---

## User decisions (2026-07-30) — all locked

1. **σ_WP = 0.5**, the exact width match (classical σ_pot = 0.35355 = 0.5/√2).
2. **CAPs added** — 12.5 Bohr per z face, |η| = 1 Ha. Departs from the campaign's
   CAP-free rule, necessarily (see "engine facts" below). Applied as **η = −1.0 Ha**:
   the user said "eta is 1" meaning strength; +1.0 would be an exponentially
   *growing* gain potential.
3. **Classical half NOT re-run** — compare against published `S_summary.csv` only.
4. **1.5× the classical step count**, same dt = 0.04.
5. **dx = 0.40** (not the classical 0.50) **and the velocity grid cut to four
   points** (v = 2.0, 2.5, 3.0, 3.5). v = 4.0 and 4.5 dropped — see aliasing below.
6. Vacuum CAP controls added by me as a direct consequence of finding 4 below;
   user was informed, not asked (cheap, and follows from the finding).

---

## Engine facts established (verified, file:line)

- **`periodicity(2)` does NOT make z open for orbitals.** It is consulted only by
  `inq/src/solvers/poisson.hpp:189,206` (Rozzi slab kernel),
  `ionic/periodic_replicas.hpp:39`, `ionic/interaction.hpp:282-312`,
  `perturbations/kick.hpp`, `operations/spatial_partitions.hpp`. The wavefunction
  basis and kinetic operator are a plain 3-D FFT over all three axes
  (`basis/fourier_space.hpp:60-151`, `hamiltonian/ks_hamiltonian.hpp:200-204`).
  **A KS orbital travelling +z wraps and re-enters at −z.** Independently
  confirmed by a subagent and by `docs/handovers/pbc-open-z-oscillation.md:20`.
  This is *why* the classical CAP-free design does not transfer.
- **`perturbations::absorbing` takes FRACTIONAL cell coordinates**, not Bohr. It
  compares `point_op.rvector()[2]`, which uses the point_operator's
  **contravariant** spacing (`basis/real_space.hpp:105,129`) and lies in
  [−0.5, 0.5). Note `real_space.hpp` has TWO `rspacing_` members — the basis-level
  one (line 28) IS Bohr; the point_operator one (line 105) is not. Passing Bohr
  would put the CAP at z ≈ 0.4 Bohr, through the slab centre.
  For 12.5 Bohr per face on Lz = 85: `CAP_WIDTH_FRAC = 0.147058823529`,
  `CAP_MID_FRAC = 0.426470588235` (= 36.25 Bohr).
- **A CAP requires `inq-study`, not stock `inq`** — stock keeps the scalar
  potential real so `vk[...] += complex(0.0,...)` does not compile; inq-study
  complexifies it (`self_consistency.hpp:176`). inq-study IS built on CSD3.
- **`compute_coulomb_wp` + `orbital_density_field` already exist** in
  `inqkit/jellium/interaction_energies.hpp` and implement the WP ledger closure
  (`E_hartree = E_SS+E_PS+E_PP`, `E_external = E_SB+E_PB`) with built-in check
  columns. No new kernel was needed.

---

## Findings that changed the design

**1. GS fidelity — PASSED (this licenses the whole comparison).**
The classical GS was lost with the `/local/data/public` machine. Recomputed here:
`E_GS = 207.183239622 Ha` vs recorded `207.183221561` → **Δ = 1.8e-5 Ha**,
`num_states = 74`, `∫n dV = 100.000`, `r_s = 4.18147`. Checkpoint at
`ResearchProject/systems/localised_jellium/shared_gs/slab_n100_L35x35x85_dx0p5_per2`.

**2. Classical raw data did NOT survive the migration.** Only
`hypotheses/classical_highdensity_sv/sv_sweep/S_summary.csv` (6 rows) and per-run
`REPORT.md`/`result.json`. No per-step overlays or WP−classical difference GIFs
are possible.

**3. CAP validated (figure: `hypotheses/wp_highdensity_sv/cap_check/cap_validation.png`).**
Geometry confirmed (`+z band [30,42.5]`, `−z [−42.5,−30]`). The CAP-off control
**wraps** (circular centroid −24 → +41.4 → −28.2, norm conserved 0.998). With the
CAP on, norm 1 → 0.226 over 48 a.u. and `min ⟨p_z⟩ = +0.61` — **no reflection**.

**4. ⚠️ THE CAP ITSELF DECELERATES THE PACKET.** In vacuum with *no bath and no
forces*, the CAP alone drags ⟨p_z⟩ from 2.00 → 0.61 over 48 a.u. (control holds
1.985 flat). Cause: σ_WP = 0.5 spreads at 1.414 Bohr/a.u., so the *leading* edge
reaches the +z band first and is preferentially removed. **Any S fitted where the
CAP is active measures the CAP.** Hence the vacuum controls (stage 5) and the
short fit window.

**5. FIT WINDOW = t ∈ [0.5, 4.12] a.u. (steps ~12–100).** Two independent limits
agree: transverse periodic images overlap at 6σ_d = L_xy → t = 4.12 a.u.; CAP
attrition stays < 0.3 % of norm until ~4 a.u. Everything later is recorded (GIFs,
absorption physics) but is not slope data.

**6. Momentum aliasing — the reason the velocity grid was cut.** σ_p = 1/(√2σ_WP)
= 1.414 is fixed by the width match, so the k-distribution folds at k_Nyq = π/dx.
MEASURED at dx = 0.5, v = 4.5: ⟨p_z⟩ = 3.44 vs 4.5 (−24 %), σ_pz² = 9.05 vs 2.0
(+353 %). A fold model reproduces this and gives, at dx = 0.40:

| v | 2.0 | 2.5 | 3.0 | 3.5 | 4.0 | 4.5 |
|---|---|---|---|---|---|---|
| σ_pz² error | +0.05 % | +0.26 % | +1.24 % | +5.06 % | +17.9 % | +55.1 % |

v = 4.0/4.5 excluded rather than caveated. **Recoverable at dx = 0.30** (≤0.11 %
everywhere) — the two missing points are the only gap vs the classical 6-point
curve. This is a WP-specific problem; a classical Gaussian *charge* has no
momentum content of its own.

**7. Gate tolerances must be RELATIVE.** σ_WP = 0.5 on dx = 0.5 is one grid point
per σ, giving a real ~1.6 % discretisation error that aborted the first cap_check
(job 32416846) against a 0.02 *absolute* bound. Now percent-level relative bounds
(⟨p_z⟩ 2 %, σ_pz² 10 %, T₁ 3 %, T₁−T₂ 5 %) — loose enough for the grid, tight
enough to catch a factor-2 blunder (which is exactly what caught v = 4.5).

---

## Files written (all new)

**Run machinery** — `ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/`
- `gs/run.cpp` — GS, env `GS_SPACING`/`GS_TAG`/`GS_DIR`; E_GS gate applies only at
  dx = 0.5 (a finer grid legitimately shifts it).
- `wp/run.cpp` — production. WP injection, two CAPs, full energy decomposition,
  `WPMomentumStats` + `WPRealSpaceStats` every step, pairwise ledger every step,
  VTIs (total / WP orbital / induced Δn), complex WP wavefunction, t=0 gates that
  abort, **5 retained numbered checkpoints** (`ckpt_step<N>`) + rolling
  `checkpoint` + `rt_state.txt`, `LJ_RESUME` with segment-suffixed CSVs.
  **COMPILES CLEAN** (verified, job 32418545).
- `cap_check/run.cpp` — free-WP CAP validation replica; also produces the vacuum
  controls (`WPC_H` selects the grid).

**Dispatchers** — `shared/bin/`
- `submit-wp-highdensity-sv.sh` — the whole 6-stage chain.
- `run-wp-hd-gs.slurm <0.5|0.4>`, `run-wp-hd-wp.slurm <smoke|0-3>`,
  `run-wp-hd-vac.slurm`, `run-wp-hd-notebooks.slurm`, `run-wp-cap-check.slurm`.

**Analysis** — `ResearchProject/systems/localised_jellium/hypotheses/wp_highdensity_sv/`
- `wp_hd_stopping.py` — adapter over the existing
  `jellium/hypotheses/bulk_ks_stopping/ks_stopping.py` (which already implements
  T₁/T₂ × s₃/s₄), plus `cap_corrected()` = slab minus vacuum twin.
- `build_run_notebooks.py` — per-velocity notebooks + synthesis.
- `cap_check/build_cap_check_figure.py` + `cap_validation.png`.

---

## What to check when the chain finishes

1. `results/<v>/run_summary.txt` → `run_completed = true`, and the t=0 gate block
   in the job log (all PASS, plus the `[info] ALIASING:` line).
2. **Ledger closure**: `hartree_residual` / `external_residual` in the notebooks'
   section 8 must be ~0. This REPLACES energy conservation, which the CAP breaks.
3. **Ehrenfest residual** s₃ − s₄ inside the fit window — should be ~0; divergence
   means CAP non-unitarity or orbital norm leakage.
4. **Norm at t = 4.12 a.u.** — if much below ~0.99, the fit window is already
   CAP-contaminated and must be shortened.
5. The vacuum-control panel (section 6): the CAP-corrected Δ⟨p_z⟩ is the physical
   signal; the raw ⟨p_z⟩ drop is not.

## Resume / extend

Never recompute from step 0:
```
LJ_RESUME=1 sbatch shared/bin/run-wp-hd-wp.slurm <idx>    # with a larger N_STEPS
```
Segment CSVs (`observables.from<N>.csv` …) are concatenated by `ks_stopping._concat_segments`.

## Not done / open

- Notebooks not yet executed (stage 6 pending — runs need hours).
- **v = 4.0 and 4.5 have no quantum counterpart.** Re-run at dx = 0.30 to recover
  them; needs a third GS.
- No WP−classical difference GIFs (finding 2).
- Campaign frontmatter for `classical-highdensity-sv` NOT touched — this is a
  separate sweep; consider its own campaign file if it grows.

---

## Milestone 2026-07-30 22:30 — CAP kinetic-energy/norm question resolved

**User asked whether the "divide KE by norm" CAP energy fix is present in this
codebase. VERIFIED ANSWER: NO — it is not in the engine here.**

- `diff -q inq/src/hamiltonian/energy.hpp inq-study/src/hamiltonian/energy.hpp`
  → **IDENTICAL**. The norm division is live in both at `energy.hpp:50-55`:
  `occ[ip]*real(arr[ip])/real(nor[ip])`, used ONLY for `kinetic_` (`:83`).
  `eigenvalues_` (`:95`) uses the bare 2-arg overload; Hartree/external/xc are
  density-based and already extensive.
- The complete inq-study divergence from stock inq is: the muon per-state mass
  fork (`ks_hamiltonian`, `laplacian`, `electrons`, `propagate`, `calculator`,
  `initial_guess`), the CAP complexification (`self_consistency`), and the new
  `absorbing_monomial.hpp`. Nothing touches `energy.hpp`.
- The historical remedy was POST-PROCESSING only:
  `scripts/wp_cap_energy_plateau/wp_kinetic_normalization_fix.py` (present here,
  but its `BASE` still points at the dead `/local/data/public/...` path).

**COULD NOT PULL the user's new engine change.** `git fetch` →
`git@github.com: Permission denied (publickey)`. Key `~/.ssh/id_ed25519`
(`SHA256:YfCaJSnG6PeScr2g9xOrPpVmpu+fcloOjkx55TFR8Ck`, label `csd3-skcb2`) is
rejected by GitHub; no ssh-agent, no credential helper, and a global
`url.git@github.com:.insteadof=https://github.com/` rewrites HTTPS to SSH so that
route does not help either. Confirmed the change is NOT on disk: inq-study at
`8c59be9`, `energy.hpp` mtime Jul 29 15:35, nothing under `inq-stack/`,
`inq-study/src/` or `shared/` modified in the last 3 h.

**POST-PROCESSING CONFIRMED — validated on live run data**, so the runs in flight
need no re-run. Implemented as `wp_hd_stopping.wp_kinetic_norm_correction()`:

    E_total_corrected = E_total_reported - occ * T1 * (1 - norm)      (occ = 1)

since `e_kin_ha` (WPMomentumStats) IS the norm-divided <T>/norm. All three inputs
are written EVERY step by the current runs, so this is exact at full cadence —
better than the original script, which reconstructed <T> and norm from ~100 sparse
wavefunction VTIs. NOTE: use the REAL-SPACE `norm_check` (wp_real_space_stats,
~1 at t=0), NOT the momentum-space Parseval constant.

Measured on v2p0 (1757 steps, t = 0 -> 70.2 a.u.):

| t (a.u.) | norm_WP | T1 (eV) | correction (eV) | E_raw (eV) | E_corrected (eV) |
|---|---|---|---|---|---|
| 0.0 | 1.0000 | 136.06 | 0.00 | 5777.58 | 5777.58 |
| 23.4 | 0.5734 | 88.78 | 37.88 | 5726.57 | 5688.69 |
| 46.8 | 0.2232 | 66.04 | 51.30 | 5705.47 | 5654.17 |
| 70.2 | 0.1173 | 60.55 | 53.45 | 5700.35 | 5646.90 |

Correction is exactly 0 at t=0 (norm=1) — sanity gate passes. NON-CIRCULAR
cross-check: the WP's BARE kinetic content (T1*norm) falls 136.06 -> 7.10 eV, i.e.
**128.95 eV leaves with the absorbed packet**, against a **corrected** E_total
drift of **130.68 eV**. They agree to 1.7 eV. The RAW ledger shows only 77.23 eV
of drift — it hides ~53 eV because the norm division keeps reporting the
per-particle mean instead of letting absorbed kinetic energy leave.

Scope: only `energy_kinetic`/`energy_total` are affected. T1/T2 (inqkit
WPMomentumStats, independent), the pairwise Coulomb ledger, and the
`hartree_residual`/`external_residual` closure gates are all unaffected.

When the upstream real-time column lands, run one job with both and cross-check
the live column against this post-processed value.

---

## Milestone 2026-07-30 23:45 — first results; notebooks built for 3 of 4 runs

**Notebooks executed with 0 errors:** `run_v2p5.ipynb`, `run_v3p0.ipynb`,
`run_v3p5.ipynb` in `hypotheses/wp_highdensity_sv/` (~43-48 MB each; 9 density
GIFs per run embedded inline at the top per `.claude/rules/notebook-density-gif.md`).
v2p0 was still at 87 % at the time of writing.

**Two builder bugs fixed (both would have silently corrupted the analysis):**
1. `run_summary()` split on the FIRST `=`, but run_summary.txt packs several
   `key = value` pairs per line, so `save_every` came back as
   `"7  wf_every = 21  stats_every = 1  ckpt_every = 414"`. Now a regex. **This is
   the exact bug the classical campaign already hit and documented** — re-fixed,
   not re-learned.
2. The duck-typed `_Shim` for the ks_stopping layout failed (`load_wp_run` calls
   `Path(run_dir)` internally). Replaced with an explicit `_load_at(obs_dir)`.
   It also now uses the REAL-SPACE norm (`norm_check_r`), not the momentum-space
   Parseval constant that `ks_stopping` picks.

### RESULTS — S over the IN-SLAB TRANSIT (eV/Bohr)

| v | S_13 (full ⟨p²⟩/2m) | S_23 (drift ⟨p⟩²/2m) | classical benchmark | broadening = S_13−S_23 | norm at slab exit |
|---|---|---|---|---|---|
| 2.5 | 1.432 | **0.849** | 0.970 | 0.584 | 0.833 |
| 3.0 | 1.049 | **0.676** | 0.709 | 0.373 | 0.893 |
| 3.5 | 0.715 | **0.490** | 0.509 | 0.226 | 0.930 |

**The drift-only quantum stopping tracks the classical benchmark to 4–12 %**
(and improves with velocity), while the full-KE definition sits 30–70 % higher.
That gap IS the momentum-width (angular-scattering + localisation) channel — the
quantity this study was built to isolate. Both fall with v, Bethe-like.

### The window question, settled by the data

Two windows are now reported in every notebook, neither silently preferred:
- **in-slab transit** — centroid inside ±12.5 Bohr; t ∈ [3.83, 12.17] at v = 3.0.
  This is where a stopping power is defined, and it gives the table above.
- **localised window** t ∈ [0.5, 4.12] — gives NEGATIVE S (−0.19 to −0.33).
  With the campaign-matched launch at z = −24 there is 11.5 Bohr of vacuum
  standoff, so that window sees the packet ACCELERATING down the slab's
  attractive gradient before it ever reaches the medium. It is a diagnostic, not
  a stopping power.

**The earlier CAP worry does not bite in the transit window.** Norm at slab exit
is 0.83–0.93, and the vacuum control sits at essentially the same value there
(0.903 vs 0.893 at v = 3.0) — so CAP attrition during the transit is small AND
common to both, i.e. it cancels in the comparison. The dramatic CAP-only ⟨p_z⟩
collapse measured in the replica happens at t ≳ 20 a.u., long after the transit.
The launch position therefore did NOT need changing (user was right).

### Norm-corrected energy ledger — applied

`wp_kinetic_norm_correction()` is wired into every notebook (section 8) and the
engine fix is still ABSENT here (`inq-study` energy.hpp byte-identical to `inq`,
still at `8c59be9`; observables.csv carries no norm column). Per-run:

| v | corrected E_total drift | raw drift | hidden by the artefact |
|---|---|---|---|
| 3.0 | 203.6 eV | 146.7 eV | 57.0 eV |
| 3.5 | 248.7 eV | 191.1 eV | 57.6 eV |

Correction is exactly 0.00000000 eV at t = 0 in every run (norm = 1) — gate passes.

---

## FINAL 2026-07-31 01:00 — campaign COMPLETE, all 4 runs + 5 notebooks

Whole chain finished autonomously. Queue empty. Every SLURM stage COMPLETED:
GS×2, smoke, sweep array 0-3, vacuum controls, notebooks (32419008, 19 min).

**Wall times:** v2p0 3:09:39, v2p5 2:43:30, v3p0 2:12:42, v3p5 1:57:33 (~10 GPU-h).
All four ran to their full 1.5× step count with 5 retained checkpoints each.

**Artefacts** in `hypotheses/wp_highdensity_sv/`:
`run_v2p0/v2p5/v3p0/v3p5.ipynb` (43–51 MB, 25 cells, **0 errors**, 9 density GIFs
embedded inline at the top of each), `synthesis.ipynb`, `wp_S_summary.csv`,
`wp_vs_classical_Sv.png`, `cap_check/cap_validation.png`.

### HEADLINE RESULT — S over the in-slab transit (eV/Bohr)

| v | t_in | t_out | S_13 (full ⟨p²⟩/2m) | S_23 (drift ⟨p⟩²/2m) | classical | S_23 vs cl | broadening | norm at exit |
|---|---|---|---|---|---|---|---|---|
| 2.0 | 5.75 | 18.25 | 1.844 | **1.004** | 1.087 | −7.6 % | 0.840 | 0.745 |
| 2.5 | 4.60 | 14.60 | 1.432 | **0.849** | 0.970 | −12.5 % | 0.584 | 0.833 |
| 3.0 | 3.83 | 12.17 | 1.049 | **0.676** | 0.709 | −4.7 % | 0.373 | 0.893 |
| 3.5 | 3.29 | 10.43 | 0.715 | **0.490** | 0.509 | −3.8 % | 0.226 | 0.930 |

**The drift-only quantum stopping reproduces the classical benchmark to 4–13 %
across all four velocities**, sitting slightly below it throughout. The full-KE
definition lies 30–85 % higher; the gap (S_13 − S_23) is the momentum-width
channel — angular scattering + localisation — and it FALLS monotonically with
velocity (0.84 → 0.23), i.e. broadening matters most for the slowest packet,
which spends longest in the medium. This is the quantity the study was built to
isolate, and it is a clean, monotone signal.

### Bugs caught and fixed in the analysis (all would have shipped wrong numbers)

1. `run_summary()` split on the first `=` — same trap the classical campaign
   documented. Now regex.
2. `_Shim` layout hack failed (`load_wp_run` calls `Path()`); replaced by
   `_load_at(obs_dir)`, which also uses the REAL-SPACE norm not the Parseval one.
3. **The synthesis originally fitted the LOCALISED window** and wrote an
   all-negative `wp_S_summary.csv`. Now fits the in-slab transit and keeps the
   localised values as `*_localised` diagnostic columns.
4. Two string-escaping faults in generated cells (`\n` collapsed; `\rangle` →
   carriage return). LaTeX removed from the affected labels.
5. `run-wp-hd-notebooks.slurm` called `build_cap_check_figure.py` at the wrong
   path (it lives in `cap_check/`). Fixed.

### Still open

- **v = 4.0 and 4.5 have no quantum counterpart** (momentum aliasing at dx=0.40:
  σ_pz² +17.9 %/+55.1 %). Recover by re-running at dx = 0.30 (≤0.11 % everywhere);
  needs a third GS. This is the only gap vs the classical 6-point curve.
- **Engine norm fix never pulled** — GitHub rejects this machine's key
  (`SHA256:YfCaJSnG6PeScr2g9xOrPpVmpu+fcloOjkx55TFR8Ck`). `inq-study` energy.hpp
  still byte-identical to `inq`. The post-processing correction is applied in
  every notebook instead; when the engine column lands, run one job with both and
  cross-check.
- No WP−classical difference GIFs (classical raw data lost with the old machine).

---

## Addendum 2026-07-31 — deposit-based S added to the synthesis plot (user request)

Added `wp_hd_stopping.deposit_stopping()` and a fourth curve to
`wp_vs_classical_Sv.png`: the CLASSICAL Definition-2 estimator evaluated on the
wavepacket runs,

    S_deposit = (E_total(t_final) - E_GS) / L_slab,   L_slab = 25 Bohr

referenced to the **dx = 0.40 production GS, 207.18323030158 Ha** — NOT the
dx = 0.50 fidelity value the classical campaign used (they agree to 9e-6 Ha, but
the deposit must reference the GS the run was actually started from).

| v | S_deposit (corrected) | S_deposit (raw) | S_23 | classical | E_deposit corr (eV) | norm at t_f |
|---|---|---|---|---|---|---|
| 2.0 | 0.239 | 2.471 | 1.004 | 1.087 | 5.98 | 0.064 |
| 2.5 | 0.200 | 2.447 | 0.849 | 0.970 | 5.00 | 0.046 |
| 3.0 | 0.167 | 2.446 | 0.676 | 0.709 | 4.18 | 0.032 |
| 3.5 | 0.131 | 2.434 | 0.490 | 0.509 | 3.27 | 0.018 |

**The plotted curve uses the NORM-CORRECTED ledger.** The raw column is kept as a
diagnostic and is the sharpest demonstration yet of why the correction matters:
on the raw ledger this estimator is velocity-INDEPENDENT at ~2.44 eV/Bohr, which
is unphysical (S must fall with v in this regime). It is flat because the residual
is dominated by the norm-divided kinetic term of the ~2-6 % of packet that
survives, not by deposited energy. After correction it falls monotonically.

**Why it sits ~4-5x BELOW the classical curve — expected, not a discrepancy.**
Classically the projectile is an EXTERNAL perturbation: never in the electronic
ledger, and the CAP-free z-open box lets it leave without removing ledger energy,
so `plateau - E_GS` is purely the slab's gain. Here the WP IS part of the system
and the CAP REMOVES it, so `E_total(t_final) - E_GS` counts only what remains in
the box after everything the CAP absorbed is already gone. It is a **lower bound**
on the deposit, not the classical quantity. The two must not be read as the same
measurement — which is exactly why S_23 (drift momentum loss), not S_deposit, is
the estimator that tracks the classical benchmark to 4-13 %.

---

## Addendum 2026-08-01 — sweep observables versioned into the repo (user request)

The user wanted the non-VTI results of the WP S(v) sweep on GitHub (to read from
another device). The production outputs live under
`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/wp/results/`,
which `.gitignore` excludes wholesale (`**/results/`), so nothing of the sweep
was ever pushed.

**Done (verified):**
- New `hypotheses/wp_highdensity_sv/export_sweep_data.py` copies, for each of
  the 12 production points (sigma_WP in {0.5, 2, 3} x v in {2.0, 2.5, 3.0, 3.5}),
  every `.csv`/`.txt` in `raw/observables` (resume segments included) plus
  `run_summary.txt`/`rt_state.txt` into
  `hypotheses/wp_highdensity_sv/sweep_data/<run_name>/`, and writes
  `observables_corrected.csv` — the segment-concatenated INQ ledger merged with
  the CAP norm-correction columns (`norm_wp`, `correction_ev`,
  `e_total_raw_ev`, `e_total_corrected_ev`, `wp_kinetic_bare_ev`,
  `energy_total_corrected` in Ha):
  `E_total_corrected = E_total_raw - T1*(1 - norm_WP)` (see
  `wp_hd_stopping.wp_kinetic_norm_correction`).
- Self-check ran and PASSED for all 12 runs: the final-step corrected energy
  reproduces the published `S_deposit_corrected = (E(t_f) - E_GS)/25 Bohr` in
  `sigma_sweep_S_deposit.csv` (atol 5e-4). Full cadence confirmed (e.g. v2p0:
  3624 steps in both source and export).
- Committed alongside: the sweep summary CSVs (`sigma_sweep_S_deposit.csv`,
  `wp_S_summary*.csv`), the sweep figures (`sigma_sweep_*.png`,
  `wp_vs_classical_Sv*.png`, force-added past the global `*.png` ignore) and the
  four small synthesis notebooks (`sigma_sweep.ipynb`, `synthesis*.ipynb`,
  force-added past `*.ipynb`) — these render on GitHub. The 23–51 MB per-run
  notebooks stay untracked.

**Not done / unchanged:** vacuum CAP-control runs (`cap_check/results/vac_*`)
are not exported — the deposit estimator does not use them; re-run
`export_sweep_data.py` with an extended run list if they are wanted too.

---

## 2026-08-03 — disk cleanup: raw VTI frames of the 16 GIF'd runs purged

**What was deleted (user-approved policy: "remove VTIs where the density GIFs
have already been made using run notebooks").** The `raw/vti/` trees (~229 GB,
~19,600 files) of all 16 production runs under
`/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/systems/localised_jellium/scripts/wp_highdensity_sv/wp/results/`:
`{,nl_,s2p0_,s3p0_}v{2p0,2p5,3p0,3p5}`.

**Pre-deletion gate (verified per run, all 16 passed):** the matching
`hypotheses/wp_highdensity_sv/run_<run>.ipynb` contains 9 embedded `image/gif`
outputs (the full density-GIF battery) and its file mtime postdates the newest
VTI file of that run — the animations are base64-embedded and survive without
the sidecar files.

**What was KEPT (verified after deletion):**
- every run's rolling `checkpoint/` + `rt_state.txt` (`last_step` 2070–3623) —
  all runs remain extendable per the final-timestep-checkpoint rule;
- all `raw/observables/*.csv` (incl. segment-suffixed resume files), the
  exported `hypotheses/wp_highdensity_sv/sweep_data/`, summary CSVs, figures,
  and all notebooks;
- the four smoke runs' VTIs (110 files — no run notebooks, fail the GIF gate);
- everything in `sigma56_sv` (active campaign, other conversation), the
  cylindrical_jellium system + `vacuum/wp_selfinteraction` (active SIC
  campaign), `vacuum/wp_traversal_energy` (no notebooks/GIFs built), and the
  jellium sigma1 bulk sets (GIF battery never built).

**Irreversible consequence:** run notebooks and twin GIFs for these 16 runs can
no longer be REBUILT from raw fields; the embedded notebook outputs are the
record. Extending a run (resume) regenerates VTIs only for the new segment.
Filesystem after purge: 746G/1.0T used (279G free).

**Addendum (2026-08-03, same cleanup):** also deleted the five interior
`ckpt_step*` snapshots of `wp/results/s3p0_v2p0` (~10 GB, steps 724–3620); its
final state remains in the rolling `checkpoint/` (`last_step=3623`, verified) so
the run is still extendable.
