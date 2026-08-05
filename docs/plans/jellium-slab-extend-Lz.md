# Plan — Slab→bulk generalisation: L_slab sweep {15, 25, 35} at fixed per-σ geometry

Status: **MACHINERY BUILT + CHAIN SUBMITTED (2026-08-05), pilot-first.**
Branch `quantum-stopping-power`. Machine CSD3, `ampere`. Proposed sweep name: **`lz_bulk_sweep`**.
Parents: `docs/handovers/sigma56-sv-twin.md` (σ=5/6 twins, 105-box — geometry template
and analysis machinery), `docs/handovers/wavepacket-highdensity-sv-twin.md`
(σ=0.5/2/3 WP sweep, 85-box — the σ=0.5 anchor).

## Goal

Generalise the slab S(v) results to bulk jellium by sweeping slab thickness
L_slab ∈ {15, 25, 35} Bohr and extrapolating the corrected deposit stopping
power in 1/L:

    S(L) = S_bulk + c/L        (surface / per-traversal terms scale as 1/L)

Inference, to be tested: linearity in 1/L across three points. The extrapolation
also discriminates the open σ classical/WP gap (1.9–3.4× at σ=6): a
per-traversal excess (self-Hartree, surface) must close as 1/L→0; a genuine bulk
quantum excess persists.

## Locked decisions (user, 2026-08-05 interview)

| Decision | Value |
|---|---|
| L_slab | 15, 25 (existing anchors), 35 Bohr |
| σ_WP | 0.5 and 5, WP + classical twins (bracket the S range) |
| Velocities | **all four** (2.0, 2.5, 3.0, 3.5) at every (σ, L) |
| σ=0.5 @ L=25 | **NOT re-run** — existing 85-box data is the anchor |
| Geometry rule | **per-σ**: each σ's new runs replicate ITS anchor's launch standoff and face→CAP gap, so arrival width and CAP distances match across L within each σ |
| VTI cadence | coarsened ~3× vs sigma56 (disk: ~280 GB free) |
| Execution | fully autonomous SLURM chain (afterok-gated smokes, warn-don't-block, self-healing finalizer) once smoke-validated |
| **Pilot-first** (user, 2026-08-05) | **v = 3.0 runs FIRST at all four boxes, both halves (~20 GPU·h); an automated pilot-gate job checks completeness, ledger closure and S sanity, and only on PASS does SLURM release the remaining ~70 GPU·h** (production arrays sit afterok on the gate) |

## Geometry — two families (all: L_xy=35, dx=0.40, dt=0.04, r_s=4.18147, CAP width 12.5, η=−1.0 Ha)

**σ = 0.5 family** (85-box layout: standoff 11.5, face→CAP 17.5 ⇒ L_z = L_slab + 60):

| L_slab | box L_z | faces | launch z | CAP inner | N_e | ~states | config |
|---|---|---|---|---|---|---|---|
| 15 | 75 | ±7.5 | −19.0 | ±25 | 60 | ~44 | NEW `slab_n60_L35x35x75.hpp` |
| 25 | 85 | ±12.5 | −24.0 | ±30 | 100 | 74 | existing 85-box (anchor data only, no new runs) |
| 35 | 95 | ±17.5 | −29.0 | ±35 | 140 | ~104 | NEW `slab_n140_L35x35x95.hpp` |

**σ = 5 family** (105-box layout: standoff 15, face→CAP 27.5 ⇒ L_z = L_slab + 80):

| L_slab | box L_z | faces | launch z | CAP inner | N_e | ~states | config |
|---|---|---|---|---|---|---|---|
| 15 | 95 | ±7.5 | −22.5 | ±35 | 60 | ~44 | NEW `slab_n60_L35x35x95.hpp` |
| 25 | 105 | ±12.5 | −27.5 | ±40 | 100 | 74 | existing `slab_n100_L35x35x105.hpp` (sigma56 twins REUSED) |
| 35 | 115 | ±17.5 | −32.5 | ±45 | 140 | ~104 | NEW `slab_n140_L35x35x115.hpp` |

CAP fractions per box: width = 12.5/L_z, mid = (CAP_inner + 6.25)/L_z (compute
exactly in each config, same convention as `slab_n100_L35x35x105.hpp`).
Classical σ_pot = σ_WP/√2 (sigma-wp-convention rule); `dyn_direct` lineage (erf/r).

## Anchors at L = 25 (existing data, used as-is)

- **σ=5:** sigma56_sv twins, `s56_S_summary.csv` (complete, corrected).
- **σ=0.5 WP:** wp_highdensity_sv deposit values (norm-corrected,
  `sigma_sweep_S_deposit.csv`: 0.239/0.200/0.167/0.131 eV/Bohr at v=2.0–3.5).
- **σ=0.5 classical:** merged classical CAP sweep, **E_PS-corrected only**
  (e.g. 0.76 at v=2.0 — never the raw S_B, which is ~40 % monopole artefact).
  **VERIFY before the figure:** confirm its launch z (assumed −24, 85-box) from
  the merged CSV provenance; if unverifiable, state the assumption in the caption.

## Run matrix (new runs only)

- 4 new production sets: (σ=0.5, L=15), (σ=0.5, L=35), (σ=5, L=15), (σ=5, L=35)
  — each 4 v × {WP, classical} = 8 runs → **32 production runs**.
- 4 new ground states (one per new box: n60_L75, n140_L95, n60_L95, n140_L115),
  gated on ∫n dV = N_e and r_s = 4.183.
- Vacuum CAP baselines: 4 v per (σ, new box) = 16 cheap WP-only runs.
- Every run: final + interior checkpoints, LJ_RESUME, decomposed interactions.csv,
  segment-suffixed CSVs (project rules).

## Comparability instruments

1. Within each σ family, launch standoff and face→CAP gap are IDENTICAL across L
   ⇒ arrival width depends on (σ, v) only; CAP distances match.
2. Residual in-transit ⟨σ_d⟩ drift with L is analytic and reported per point:
   σ=5: ≤ ~10 %; σ=0.5: ×~1.4 from L=15→35 (⟨σ_d⟩ ≈ 0.71·(23+L)/v at standoff 11.5).
3. Every run carries its MEASURED effective width ⟨σ_r⟩ (1 %-norm-loss window,
   sigma56 machinery) in the summary CSV and figure labels; new same-⟨σ_r⟩
   (σ, L) pairs extend the sigma56 collapse test.
4. One corrected estimator everywhere:
   `S = [E_total(t_f) − E_GS − E_PS(t_f)]/L_slab`, WP norm-corrected,
   each L referencing its OWN box's E_GS.
5. Cross-σ comparisons are bracket-level only (the two families differ in vacuum
   length; corrected estimators are insensitive to it — inference, noted).

## Registered caveats (agreed)

- **The quantitative 1/L extrapolation rests on σ=5.** The σ=0.5 WP disperses
  (in-slab width grows with L; transverse images overlap before arrival at
  v ≤ 3) — its trace is qualitative context. σ=0.5 CLASSICAL is the clean upper
  bracket.
- The WP deposit is NOT monotone in σ (recorded: σ=0.5 sits below σ=3/6 because
  ⟨σ_r⟩ = 10.4 Bohr is widest and the norm correction dominates) — the bracket
  claim is for the classical family.
- **L=15 bulk-likeness gate:** Friedel period ≈ 6.9 Bohr ⇒ ~2 oscillations
  across a 15-Bohr slab; σ=5 packet ±3σ_d ≈ 21 Bohr > slab. GATE: run both
  L=15 GS first, require n(z=0) within ~1–2 % of n0 (profile compared to the
  L=25 GS) before TDDFT. If it fails: still plot L=15, quote the extrapolation
  from L=25/35, add contingent L=50.
- σ=5 transverse windows clean at all v (slab exit ≤ 25 a.u. < image overlap
  ~33 a.u.). σ=0.5 aliasing at dx=0.40 unchanged (σ_pz² bias ≤ +5.1 % at v=3.5).

## Step counts and cost (WARN-not-gate; re-measure at smoke)

    N_STEPS = round(4.36·(|z0| + L_z/2)/(v·dt)),  dt = 0.04   (sigma56 formula)

| set | box | N_STEPS v=2.0/2.5/3.0/3.5 | ~s/step (states×grid from measured 3.15 @ 105-box) | ~GPU·h (×2 halves) |
|---|---|---|---|---|
| σ=0.5, L=15 | 75 | 3079/2463/2053/1760 | ~1.3 | ~7 |
| σ=0.5, L=35 | 95 | 4170/3336/2780/2383 | ~4.0 | ~28 |
| σ=5, L=15 | 95 | 3815/3052/2543/2180 | ~1.7 | ~11 |
| σ=5, L=35 | 115 | 4905/3924/3270/2803 | ~4.8 | ~40 |

Projected total **~90 GPU·h** production + 4 GS + vacuum baselines ≈ ~100 GPU·h.
Disk: VTI cadence ~3× coarser than sigma56 ⇒ ~8 GB/run worst case; GIF batteries
unaffected. Checkpoint/resume per `final-timestep-checkpoint` rule; exclude
known-bad nodes (`--exclude=gpu-q-2,gpu-q-25`).

## Execution — the pilot-gated chain (AS BUILT, `submit-lz-bulk-sweep.sh`)

1. **GS ×4** (`run-lzb-gs.slurm <preset>`; first builds, rest exec afterok).
   Hard gates in the binary: ∫n dV = N_e, r_s = 4.183. Bulk-likeness of the
   interior (n(z=0) vs n0 from the GS density VTI) is a **WARN in the pilot
   gate**, not a block — a thin slab is an interpretation problem, not a
   correctness one (checkpoint-dont-block); the pilot itself is the GPU-budget
   protection.
2. **Smoke** (`run-lzb-smoke.slurm`): builds wp+classical binaries, runs 20-step
   t=0-gate smokes for all 4 presets × both halves; any failure afterok-blocks
   everything downstream.
3. **PILOT**: v = 3.0 at all four boxes × {wp, classical} (8 jobs, ~20 GPU·h).
4. **Pilot gate** (`pilot_gate.py`, afterany the pilot): HARD on completeness /
   ledger closure ≤ 1e-5 Ha / finite non-negative S; WARN on bulk-likeness,
   S(L)-vs-anchor ordering, and the measured-cost re-projection (emailed).
   Exit ≠ 0 ⇒ production arrays stay DependencyNeverSatisfied + PILOT_REPORT.md
   says why.
5. **Production** (afterok gate): wp + cl arrays 0-11%4 (v = 2.0/2.5/3.5,
   boxes cheap-first), vacuum baselines ×4 boxes.
6. **Finalize ×2** (afterany): status → repair (per-run timeout cap ~2× measured
   cost — the sigma56 flaw is fixed) → `build_lzb_figures.py` (S vs 1/L with
   1/L fits + S_bulk intercepts; every figure ALSO written to
   `docs/reports/report2/drafts/draft1/figures/jellium_slab/slab_*.png` at the
   house standard) → CAMPAIGN_REPORT.md → best-effort email.

Files: configs `shared/configs/lzb_boxes.hpp` (runtime presets, env LZB_CFG);
binaries `scripts/lz_bulk_sweep/{gs,wp,classical,vac}/run.cpp` (sigma56 clones,
one binary per type serves all boxes); dispatchers `shared/bin/lzb-params.sh` +
`run-lzb-{gs,smoke,wp,cl,vac,pilotgate,finalize}.slurm` +
`submit-lz-bulk-sweep.sh`; analysis
`hypotheses/lz_bulk_sweep/{lzb_stopping,pilot_gate,finalize,build_lzb_figures}.py`.
Run notebooks (density-GIF battery) are a follow-up after the pilot; the
finalizer already calls `build_run_notebooks.main()` if the builder appears.

## Open

- Contingent L=50 point — only if the L=15 bulk-likeness WARN fires badly or
  S(15) bends off the 1/L line.
- Run notebooks for the new runs (density-GIF rule) — follow-up once the pilot
  lands; wire `build_run_notebooks.py` into `hypotheses/lz_bulk_sweep/`.
- Gmail credentials absent on this machine — gate/finalizer report to disk;
  `python -m inqview.email setup` would enable email.

## Resolved

- **σ=0.5 classical anchor geometry VERIFIED (2026-08-05):** `S_of_v_cap.csv`
  (classical_highdensity_sv/dyn_direct) carries z_final/v_final consistent with
  the 85-box launch z = −24, and the monopole-corrected
  S = (E_absorbed − 100/z_final·Ha)/25 reproduces the recorded 0.760 eV/Bohr at
  v = 2.0 (cross-check: S_A_keloss = 0.763, an estimator sharing no machinery).
  `lzb_stopping.anchors()` implements exactly this.
