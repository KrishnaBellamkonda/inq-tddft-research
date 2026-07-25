# Plan: effective-mass muon WP run, re-planned for ≤12 h (N3/A)

Status: **SPECIFICATION LOCKED, not yet launched** (2026-07-08).
Supersedes the aborted `effmass_pair` run (dx=0.333, ~200 s/step, ~4–7 day ETA).
Owner task: Phase 4 (muon effective-mass fork demonstration).

## Why we re-planned

The live `effmass_pair` quantum run measured **~200 s/step** at dx=0.333 (6.08 M
pts), vs the old cost model's **~14 s/step** prediction — a ~15× miss. Cause: the
old model was calibrated on a *different* run (`qsp_phase4`, 3.75 s/step @ dx=0.50,
61 states) and extrapolated with n·log n across a 3.4× grid jump + different state
count + GPU contention. The cost model is now **re-anchored on the measured
200 s/step** (`hypotheses/muon_mass_fork/replan_12h_scan.py`).

## Locked physics anchors (unchanged from earlier design)

- **Velocity v = 2.711 a.u.** = 100 eV-electron velocity ⇒ same S(v) comparison target.
- **Density r_s ≈ 5.69** localised slab, 25 Bohr thick (half_width 12.5, sharp Θ).
- **Effective-mass fork**: mass m = k0/v (want m meaningfully >1 for the demo).
- **T_total = 3.5 × traversal**; dt from the ETRS cliff (H·dt = E_cut·dt ≲ 2.2).
- Readout: **coherent momentum-peak n(k,t)** → initial-drag S(v0) (robust to real-space spread).

## Relaxed this round (user, 2026-07-08)

- Impact-spread threshold **<1% → <2%**.
- Grid + cell + mass are free cost levers to hit **≤12 h on one GPU**.

## THE RUN — N3/A (user-selected)

Quantum effective-mass muon WP on the shrunk r_s=5.69 slab.

| Parameter | Value | Note |
|---|---|---|
| Cell | **36 × 36 × 80 Bohr** | transverse shrunk 50→36 to fit 12 h; z 90→80 |
| Grid spacing dx | **0.40 Bohr** | k_max = π/dx = 7.854; 1.62 M pts |
| Electrons N | **42 (even)** | n0 = 1.296e-3, **r_s = 5.69** (held density) |
| Extra states | 10 | 21 occ + 10; fewer states than old run (61) ⇒ cheaper |
| σ_WP | **2.0 Bohr** | σ_ρ,0 = 1.414 |
| k0 | **6.793 /Bohr** | maxed to Nyquist: k_max − 3σ_p (σ_p=0.354) |
| Mass m | **2.506 m_e** | inverse_mass = **0.39907** (the fork) |
| Energy E | **251 eV** | ½ m v² |
| Launch z | **−16.39 Bohr** | 2.75 σ_ρ,0 before the slab face (−12.5) |
| Impact spread | **1.02 %** | launch→slab face; <2 % ✓ |
| Aliased tail | 0.14 % | 3σ_p edge at Nyquist; <2 % BLOCK ✓ |
| **dt** | **0.05 a.u. (SMOKE-TEST FIRST)** | E_cut=30.8, H·dt=1.54 (below 2.2 cliff) |
| N_steps | **846** | T = 42 a.u. = 3.5 × traversal (12.1 a.u.) |
| CAP | η=−0.7, mid |z|=32 (0.40·Lz), width 8 (0.10·Lz) | retuned for Lz=80 |
| **Est. wall** | **~11.5 h / 1 GPU** | conservative (grid-only); likely faster (½ the states) |
|  | ~6.4 h / 2 GPU | if the 2nd GPU frees |

### dt smoke-test procedure (mandatory — the last model was optimistic)
1. Reuse the new GS; run ~15 steps at **dt=0.05**. Pass = finite, energy stable.
2. If it diverges (NaN / runaway): fall back to **dt=0.04** (H·dt=1.23) AND trim
   T to 3.0× (N_steps≈725) to stay ≈12 h. dt=0.05 is expected to pass.

## Ground state (new — required by the transverse shrink)

Bare r_s=5.69 slab at the new cell/grid. Config header + GS binary:
- **New config**: `shared/configs/slab_n42_L36x36x80.hpp` (N=42, cell 36×36×80,
  SPACING 0.40, half_width 12.5, n0=1.296e-3, EDGE_WIDTH 0, EXTRA_STATES 10).
- **GS**: `scripts/muon_mass_fork/effmass_12h/gs/run.cpp` (adapt `effmass_pair/gs`,
  point at the new config; Broyden mix; `electrons.save(shared_gs/slab_n42_L36x36x80_dx0p40)`).

## Quantum run machinery

- **New sweep folder**: `scripts/muon_mass_fork/effmass_12h/quantum/` (adapt
  `effmass_pair/quantum/run.cpp`: include the new config; env defaults →
  EM_SPACING 0.40, EM_SIGMA_WP 2.0, EM_K0 6.7933, EM_INV_MASS 0.39907,
  EM_LAUNCH_Z −16.389, EM_DT 0.05, EM_N_STEPS 846, CAP mid/width for Lz=80,
  EM_GS_DIR shared_gs/slab_n42_L36x36x80_dx0p40).
- Outputs → `effmass_12h/quantum/results/` (logs gitignored).

## Classical twin — DEFERRED draft phase (built only after quantum succeeds)

User choice: **free Ehrenfest at the matched mass**, NOT constant-velocity.
- `scripts/muon_mass_fork/effmass_12h/classical/` (adapt `effmass_pair/classical`).
- Gaussian-electron UPF **`electron_gaussian_wpsigma2p0.upf`** (σ_pot=σ_WP/√2=1.414;
  ALREADY generated + DATA-verified).
- **mass = 2.506 m_e** (EM_MASS_ME=2.506 → /1822.8885 amu); v0 = 2.711; same GS.
- Projectile decelerates → **S(v0) from the early (vz ≥ 0.85 v0) window** per the
  light-projectile rule. park+remove after |z| ≥ CAP inner face.
- Marked **draft**: scaffold now, run after the quantum run + its analysis land.

## Validation notes

- **Transverse wake-wrap**: box shrunk 50→36 (periodic in x–y). At r_s≈5.7 screening
  is long-ranged; ±18 Bohr box vs ±8.3 Bohr exit packet gives ~9.7 Bohr clearance.
  Momentum-peak readout (projectile's own KE loss) is far less sensitive to wake
  wrap than a force-based readout. **Later check**: compare S(v0) at 36 vs a spot
  50-box point if the number looks off.
- **dt cliff** confirmed empirically last round (H·dt≈2.2 divergence); dt=0.05 sits
  at H·dt=1.54 with margin. Smoke-test still mandatory.

## Sequence (on user GO)

1. Write config header + GS/quantum run.cpp (adapt effmass_pair).
2. **Stop** the current effmass_pair quantum run (pid 1411675) + orchestrator
   (pid 1217967). ← irreversible; needs explicit GO.
3. Build + run new GS (bare slab, N=42, dx=0.40).
4. dt smoke-test (15 steps @ 0.05).
5. Launch quantum N3/A on GPU 0 (846 steps).
6. Auto-build run notebook + phase notebook on completion; email.
7. Scaffold (do not run) the classical draft phase.
