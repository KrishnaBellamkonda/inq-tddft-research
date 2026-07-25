# Plan: σ=1 concentrated (chirped) run — lean re-run after VRAM-wall abort

Date: 2026-07-09. Status: USER-STEERED GEOMETRY (40×40×80), awaiting GO.
Supersedes the first σ=1 launch (killed at step ~36 after ~1 h; ~104 s/step → 34 h ETA).
Parent effort: σ-comparison twin runs, `docs/handovers/muon-mass-fork.md`.

## Why the first launch was infeasible

The chirped focusing launch (`WavePacket::focus_z()`, vacuum-validated) requires
**dx = 0.333** — at dx = 0.40 the chirp's momentum tail (k0 + 3σ_p = 7.81) sits at
Nyquist (k_max = 7.85) and grid dispersion breaks the focus (measured: min σ_z
0.775 vs ideal 0.707). At dx = 0.333 the 50×50×101 box = 150×150×303 ≈ 6.9 M pts
× 61 states ≈ 6.7 GB per wavefunction set; INQ's several sets saturated the
24 GB card (0 MB free). Past the memory wall the per-step cost went super-linear:
**104 s/step** (vs 8.3 s/step for the σ=2 run at 3.5 M pts) → 34 h for 1200
steps. Physics was clean (E flat to 1e-5 Ha; in-medium focusing confirmed:
σ_z 0.864 → 0.729 at t = 1.28, heading to the 0.707 waist at t ≈ 1.48).

## User geometry directive (2026-07-09)

Keep dx = 0.333. Shorten z; x = y = 40 Bohr. Hold the slab DENSITY fixed
(adjust N with the transverse area). Produce an empirically-backed time
estimate from measured run costs.

## Locked physics (unchanged)

| Item | Value | Why |
|---|---|---|
| σ_WP | 1.0 | user requirement (concentrated) |
| v | 2.711 | S(v) anchor, locked |
| m | 2.10 (inv 0.476190), k0 = 5.693 | Fable-5 deliberation |
| Chirp | `focus_z(4.0, 2.10)` | waist (σ_ρ,z = 0.707) at slab face; vacuum-validated at dx = 0.333 |
| dx | 0.33333 | chirp Nyquist requirement; E_cut = 44.4 Ha |
| dt | 0.04 | H·dt = 1.78 < 2.2 empirical cliff (20 % margin) |
| CAP | η = −1.0, width 15 Bohr/side | reflectivity-curve tuned (σ-sweep) |
| Classical twin | m = 2.10, `electron_gaussian_wpsigma1p0.upf` (σ_pot = 0.707) | UPF generated + VERIFIED by V(r) (erf ratio 1.000 at all r) |

## New geometry + electronic structure

| Item | Value | Derivation |
|---|---|---|
| Box | **40 × 40 × 80** Bohr | user; z-budget below |
| Grid | 120 × 120 × 240 = **3.46 M pts** | dx = 0.333; FFT-friendly (2³·3·5 / 2⁴·3·5) |
| z-layout | slab ±12.5, launch −16.5 (4σ), **CAP inner ±25**, CAP [25, 40] | packet σ_ρ,z(exit) ≈ 3.2 → 3σ front at +22 < 25 ✓; back-tail −20 > −25 ✓ |
| N electrons | **52** | n0 = 1.312e-3 (σ=2 run) × 40·40·25 = 52.5 → 52 (even); **r_s = 5.679** (0.2 % from 5.667) |
| States | 26 occupied + **EXTRA_STATES = 10** = **36** | T = 100 K smearing; effmass_12h precedent; WP slot = last extra |
| N_STEPS | **900** (36 au) | slab traversal 9.2 au → 3.9× (user window 3–4×); flight+absorption done by ~19 au |
| CKPT_EVERY | **225** | 4 checkpoints, resumable beyond final time |
| Transverse clearance | ±20 = 3.8σ_ρ,x at slab exit, ~3.1σ at full absorption | τ_s = mσ² = 2.1 au; σ_ρ,x(t) = 0.707·√(1+(t/2.1)²) |

## Empirical cost model (calibrated on this campaign's measured runs)

| Run | Grid pts | States | pt·states | s/step (measured) | unit cost s/(pt·state) |
|---|---|---|---|---|---|
| effmass_12h | 1.62 M | 31 | 5.0e7 | 5.8 | 1.16e-7 |
| σ=2 (complete) | 3.52 M | 61 | 2.15e8 | 8.3 | 3.9e-8 |
| σ=1 aborted | 6.89 M | 61 | 4.2e8 | 104 | 2.5e-7 ← VRAM wall |

Unit cost FALLS with workload (GPU saturation) until VRAM saturates (6× jump).
Log-log interpolation through the two healthy points (slope −0.78) at the new
workload 3.46 M × 36 = 1.24e8 pt·states → **7.3 s/step**; naive linear scaling
from σ=2 gives 4.8 → adopt **5–8 s/step**.
Memory: 3.46 M × 36 × 16 B ≈ 2.0 GB/set → **~6–8 GB total** ≪ 24 GB ✓.

## Time estimate (empirically backed)

| Stage | Estimate | Basis |
|---|---|---|
| Lean GS (new) | ~35–45 min | 0.30× pt·states of the ~100-min aborted-config GS |
| WP run, 900 steps | **1.5–2.0 h** | 5–8 s/step model above |
| Classical twin (GPU 1, concurrent) | ~1–1.5 h | 35 states, same grid |
| **Total wall** | **~2.5–3 h** | GS + WP sequential on GPU 0; twin hidden |

## Execution sequence

1. Config header `shared/configs/slab_n52_L40x40x80.hpp` (n0 = 1.312e-3, N = 52,
   EXTRA_STATES = 10, slab half-width 12.5).
2. GS on GPU 0 → `shared_gs/slab_n52_L40x40x80_dx0p333`.
3. σ=1 chirped WP run on GPU 0 (inq-study build; EM_RESUME-capable; pilot gate below).
4. Classical twin on GPU 1 concurrently (free Ehrenfest, m = 2.10, verified UPF).
5. Run-notebooks + WP−classical comparison.

## Pilot gate (first ~30 steps, ~4 min)

- s/step ≤ 12 (else the model is wrong — stop and re-examine).
- E drift < 1e-4 Ha over 30 steps.
- σ_z(t) decreasing toward 0.707 by t ≈ 1.4 (chirp active in-medium).
- NO gate on velocity drift (light-projectile rule).

## Rejected alternatives (recorded)

- **dx = 0.40 + lower m (~1.8) to clear Nyquist:** similar cost, loses the
  validated chirp regime, needs fresh vacuum revalidation.
- **Let the 34 h run grind:** infeasible turnaround; killed at step 36.
- **Drop the chirp:** reverts to the rejected 22 %-spread-at-impact problem.
- **Hold N = 82 in the smaller box:** raises density (r_s 5.22), breaking
  σ-family comparability — user chose density-matched N = 52.
