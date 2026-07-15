# Plan: σ=1 mass-pair decay runs (clean geometry, full + pairwise energy decomposition)

**Date:** 2026-07-15 · **Status:** approved (grill-with-docs interview complete)
**Goal:** two WP runs in the *clean* qsp_phase3 geometry (the one whose E_total
decays to a fixed value — reproduced bit-for-bit 2026-07-14), changing ONLY the
wavepacket: σ_WP=1, heavier mass, aliasing-safe, spreading ≤10% at the slab,
with checkpoints and the full energy ledger + pairwise Coulomb decomposition.
No classical twin, but observables follow the twin-run-generation contract so a
twin can be added later.

## Locked decisions (interview trail)

| decision | value | why |
|---|---|---|
| spreading target | ≤10 % at slab face (relaxed from 4 %) | 4 % needs k0 ≥ 19.7 at d=11.25 — impossible on dx=0.5 (k_Nyq=6.28). Mass cancels out of spreading-at-arrival: σ(d)=σ0·√(1+(d/(2·k0·σ0²))²) — see CONTEXT.md "Spreading-at-arrival law". |
| launch_z | **−16.5** (4σ = 4 Bohr from slab face −12.5) | 10 % at σ=1 needs d ≤ 4.4 Bohr at the guard ceiling; 4σ honours the boundary rule. CAP standoff unchanged. |
| k0 | **4.5** | spreading 9.4 %; aliased tail 0.58 % (guard WARN — fails strict 3σ, well under 2 % BLOCK). |
| masses | **m=2 and m=3** (pair) | at fixed k0, mass only apportions v and E: m=2 → v=2.25, E=138 eV (1.7 eV/e⁻); m=3 → v=1.5, E=92 eV (1.1 eV/e⁻, matches clean run's deposit). Same spreading — isolates deposit strength. m=1 rejected (E=276 eV, hotter than the oscillating family). |
| everything else | byte-identical to `p3_wp_m1_rerun` | box 50×50×90, dx 0.5, N=82 (61 states, wp_idx 60), GS `shared_gs/slab_n82_L50x50x90`, two-sided CAP η=−0.7 region ±35..±45, dt 0.04, 2500 steps (τ=100), WRITE_EVERY=WF_EVERY=8, periodicity 3, LDA, ETRS. |
| validation gate | cutoff_guard + 50-step smoke | smoke checks: finite energies, interactions closure ≤1e-9, early σ(t) vs formula, checkpoint+resume round-trip. |
| GPU | GPU 1 only, runs serial | GPU 0 occupied (standing constraint). |

Timing sanity: m=2 reaches far CAP inner face (+35) at t≈23, absorbed ~t≈30;
m=3 at t≈34, absorbed ~t≈45. Both leave ≥55 a.u. settle window before τ=100.
m≥4.5 rejected: absorption completes too late; m=10 never arrives (t≈114>τ).

## File placement (ADR 0007)

- `ResearchProject/systems/localised_jellium/scripts/sigma1_masspair/wp/run.cpp`
  — build-once, env knobs `LJ_OUT LJ_INV_MASS LJ_K0 LJ_SIGMA_WP LJ_LAUNCH_Z
  LJ_N_STEPS LJ_DT LJ_WRITE_EVERY LJ_WF_EVERY LJ_RESUME`.
- Runs: `…/sigma1_masspair/{wp_m2_k4p5, wp_m3_k4p5}/` (outputs; logs gitignored).
- Analysis: `…/hypotheses/sigma1_masspair/` (combined notebook, post-run).

## run.cpp composition (all pieces have working references)

1. **Base:** `scripts/qsp_phase3/wp/run.cpp` — geometry, CAP, GS load, WP
   injection, full observable suite (density/wavefunction VTIs, momentum
   distribution, wp_momentum_stats, **wp_real_space_stats** → verifies the
   9.4 % spreading prediction, overlaps, state energies, occupations,
   electron_number).
2. **Mass fork:** `electrons.inverse_mass()[0][wp_idx] = LJ_INV_MASS` after
   injection (reference `scripts/muon_mass_fork/sigma1_massonly/wp/run.cpp`;
   engine validated by the 2026-07-14 bit-identity rerun + code audit).
   k0 relation: v = k0/m, E = k0²/2m.
3. **Checkpoints + resume:** from `scripts/localised_jellium_dynamics/phase5_wp/run.cpp`
   — interior RT checkpoint every **200 steps** (overwrite one ckpt dir,
   ~1.8 GB) + final ckpt + `rt_state.txt` (`last_step/time_au/dt/wp_idx`) +
   `LJ_RESUME=1` branch (reload, re-apply inverse_mass, segment-suffixed CSVs,
   `propagate(..., START)`), per `.claude/rules/checkpoint-dont-block.md`.
4. **Pairwise decomposition:** per-step `interactions.csv` via
   `inqkit/jellium/interaction_energies.hpp::compute_coulomb_wp`
   (n_total, n_wp=orbital_density_field, φ₊ precomputed once; E_BB constant;
   columns `step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,e_hartree_check,
   e_external_check,norm_wp,norm_total`). 2 Poisson solves/step ≈ negligible
   vs the 4.2 s step (phase5_wp precedent emits every step).
5. **Ledger extension:** ObservableSelection adds `energy_external` +
   `energy_nonlocal` (twin-run-generation contract).

## Gates → launch sequence (serial, GPU 1)

1. `cutoff_guard.py --kind wp --sigma 1.0` with **effective energy 276 eV**
   (guard assumes m=1; pass E'=k0²/2 so p0=4.5). Expect WARN (tail 0.58 %).
2. Build once (`inq-run`, `TMPDIR=$PWD/build/tmp` — /tmp is full).
3. 50-step smoke (m=2 knobs): finite energies; closure ≤1e-9; σ_z(t) matches
   σ(t)=σ0√(1+(t/(2mσ0²))²); kill+`LJ_RESUME=1` round-trip reproduces the
   uninterrupted trace.
4. Production `wp_m2_k4p5` (2500 steps, ~4–5 h) → then `wp_m3_k4p5`.
5. Post-run: measured spreading at slab vs 9.4 % / arrival times vs 1.8 (m=2)
   & 2.7 a.u. (m=3, slab face); E_total decay shape vs the rerun; pairwise
   table narrative per twin-run-analysis rules (E_PP drop = dispersion, etc.).

## σ convention correction (smoke-gate finding, 2026-07-15)

The 50-step smoke at `sigma(1.0)` measured σ_z(0)=0.707 — the house
`WavePacket.sigma()` parameter is the WAVEFUNCTION width; the density std is
σ/√2 (consistent with the σ_WP labelling rule and the guard's σ_p=1/(√2·σ_WP)).
The interview's spreading arithmetic (9.4 % at the slab) used σ0=1 as the
DENSITY std, so at `sigma(1.0)` the ≤10 % target is physically unreachable
(min-uncertainty core spreads ~35 % by the face; measured second-moment ~73 %,
inflated by orthogonalisation tails — the p3 run shows the same birth
broadening at a mild 8 %). Actual aliased weight at k0=4.5 was fine (~0.01 %
wrapped; the 4.5 % Gaussian estimate used the tail-inflated moment σ).

**Resolution (autonomous, user away):** run with `sigma(√2)` → density std
1.0 Bohr (house label σ_WP=1.41). This is the packet ALL approved plan numbers
(9.4 % spreading, aliasing margins, launch −16.5, m=2/3, 138/92 eV) were
computed for; guard verdict improves WARN→strict PASS (tail 0.02 %). The
user's twice-negotiated spreading bound was prioritised over the literal
"σ=1" label; flagged for morning review.

## Watchdog + intervention (user requirement, 2026-07-15)

The 9 h window must stay utilised; an agent intervenes if a run stops midway.

- `scripts/sigma1_masspair/orchestrate.sh` runs the serial chain
  (guard → build → smoke → m2 → m3) in the background, logging per-run.
- **Stall watchdog:** alongside each production run, a monitor checks the run
  log every 5 min; no new step line for 15 min → kill the hung process, exit
  with code 42 (STALL). A crashed run exits with its own nonzero code.
- **Idempotent resume:** on start, each run checks its `rt_state.txt`; if
  `last_step < N_STEPS`, it relaunches with `LJ_RESUME=1` (checkpoint loses at
  most 200 steps ≈ 15 min). Completed runs (run_summary `run_completed=true`)
  are skipped.
- **Intervention:** any abnormal orchestrator exit notifies the interactive
  session, which spawns a **fable agent** briefed to: read the tail of the
  failing log + `rt_state.txt`, classify (OOM / CUDA error / stall / disk),
  fix what is safe (TMPDIR, stale lock), and relaunch `orchestrate.sh`
  (auto-resumes). Per `.claude/rules/checkpoint-dont-block.md`: never
  self-block on cost — resume and warn.

## Risks / watch-items

- **WARN-zone aliasing** (0.58 % tail): accepted by design; any k0/σ change
  re-runs the guard.
- **Deposit-strength risk (H2):** 138 eV (m=2) is 1.4× the clean run's 100 eV.
  If it oscillates while m=3 (92 eV) decays, that is itself an H1-vs-H2 datum —
  record either way in `docs/handovers/energy-oscillation-debugging.md`.
- CAP tuned at v=2.71; our v=2.25/1.5 are same order — absorption expected
  fine; verify via electron_number drain completing before τ−50.
- Guard/boundary conventions used: σ_p=1/(√2·σ_WP)=0.707, k_Nyq=π/0.5=6.283,
  4σ launch rule at σ=1.
