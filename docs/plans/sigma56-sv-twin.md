# Plan — σ_WP = 5 and 6 classical+wavepacket S(v) twins (localised jellium slab)

Sweep name: **`sigma56_sv`**
Branch `quantum-stopping-power`. Machine **CSD3**, `ampere`,
account `mphil-nikiforakis-skcb2-sl2-gpu`. Started 2026-08-02.

Parent work: `docs/handovers/wavepacket-highdensity-sv-twin.md` (σ = 0.5/2/3 WP
sweep), `docs/handovers/classical-highdensity-sv-benchmark.md` (the σ = 0.5
classical benchmark), `docs/handovers/effective-sigma-near-launch.md` (launch-vs-
arrival width).

---

## 1. Goal (one line)

Produce **matched classical + wavepacket twin pairs at σ_WP = 5 and 6 Bohr**
across the velocity grid, and combine them with the existing σ = 0.5/2/3 results
into one S(v) figure, to find the width at which the classical and quantum
projectiles stop being distinguishable.

## 2. Why σ = 5 and 6 — the width stops moving

A free Gaussian packet's density width is σ_d(t) = √(σ²/2 + t²/2σ²), so σ_WP sets
both the initial width **and** the spreading rate 1/(√2σ). Over the in-slab
transit at the new launch point:

| σ_WP | growth ×, v = 2.0 → 3.5 | ⟨σ_d⟩ over transit | σ_eq = √2⟨σ_d⟩ |
|---|---|---|---|
| 0.5 | ×3.2 | 17 → 9.7 | 24 → 14 |
| 2 | ×2.7 → 2.2 | 4.5 → 2.8 | 6.4 → 4.0 |
| 3 | ×1.9 → 1.4 | 3.6 → 2.7 | 5.1 → 3.8 |
| **5** | **×1.23 → 1.08** | 4.06 → 3.72 | **5.74 → 5.26** |
| **6** | **×1.12 → 1.04** | 4.56 → 4.35 | **6.45 → 6.15** |

At σ = 5 and 6 the packet is effectively constant-width and — crucially — the
**label agrees with the time-average** (σ_eq ≈ 6.2–6.5 for a σ = 6 packet),
whereas at σ = 2 the effective width is 4.0–6.4 depending on velocity, so a
single σ label there is meaningless. That is what makes a classical twin at a
*fixed* σ_pot a fair comparison at σ = 5/6 and only an approximation below.

**Free consistency test.** A σ = 6 packet sits at σ_eq ≈ 6.45 at v = 2.0 —
exactly where the existing σ = 2, v = 2.0 run already sits on the time-averaged
axis. If both give the same S, time-averaged σ is a valid collapse variable; if
not, the pattern in the reference figure is coincidental. This costs nothing
extra and is the single most informative comparison in the campaign.

## 3. Locked decisions (user, 2026-08-02)

| Decision | Value |
|---|---|
| σ_WP | **5 and 6 Bohr**, both halves |
| Classical binary | **`dyn_direct`** — direct erf/r potential (`moving_gaussian_projectile_potential`), NOT the Poisson perturbation. No charge in the cell ⇒ the projectile may sit outside the box and there is no clip/exit transient. |
| Classical σ | σ_pot = σ_WP/√2 derived **inside** the binary (`.claude/rules/sigma-wp-convention.md`). Runs are labelled σ_WP. σ_pot = 3.5355 (σ=5), 4.2426 (σ=6). |
| Velocity grid | **4 points: v = 2.0, 2.5, 3.0, 3.5** (matches the existing WP campaigns point-for-point) |
| CAP | **ON in BOTH halves**, η = −1 Ha, 12.5 Bohr per z face — so E_absorbed/L_slab is the same estimator on both. Plus one **CAP-free classical control at v = 3.0 per σ** to measure what the CAP costs. |
| Box | **L_z 85 → 105** (+20 Bohr of vacuum). L_x = L_y = 35, slab 25 Bohr, N = 100 unchanged ⇒ **r_s = 4.183 unchanged**. |
| Launch z | **−27.5, common to both σ and both halves** (WP and its classical twin must start identically) |
| Figure | new σ = 5/6 twins + existing σ = 0.5/2/3 WP traces as-is, L_z difference stated in the caption |

## 4. Geometry

Box z ∈ [−52.5, +52.5]; slab [−12.5, +12.5]; CAP bands [−52.5, −40] and
[+40, +52.5].

    CAP_WIDTH_FRAC = 12.5/105 = 0.119047619048
    CAP_MID_FRAC   = 0.5 - W/2 = 0.440476190476     (= 46.25 Bohr)

`perturbations::absorbing` takes **fractional** cell coordinates, not Bohr
(`docs/handovers/wavepacket-highdensity-sv-twin.md`, engine facts) — passing Bohr
would put the CAP through the slab centre.

Launch z = −27.5 gives 15 Bohr of standoff to the slab face and 12.5 Bohr of
clearance to the CAP. Both clearances are comfortable at both widths:

| σ_WP | σ_d(0) | to CAP | in σ_d | weight in CAP at t=0 | to slab | in σ_d | weight in slab at t=0 |
|---|---|---|---|---|---|---|---|
| 5 | 3.536 | 12.5 | 3.54 | **0.020 %** | 15.0 | 4.24 | **0.001 %** |
| 6 | 4.243 | 12.5 | 2.95 | **0.16 %** | 15.0 | 3.54 | **0.020 %** |

Both are at or below the 0.23 % t=0 CAP loss already accepted for the σ = 3
campaign. Transverse periodic images overlap at t_ov = 32.8 (σ=5) / 34.0 (σ=6)
a.u., and the in-slab transit ends by t = 20.0 a.u. at the slowest velocity — so
the whole transit is transversely clean, unlike σ = 0.5 where the two windows did
not intersect at all.

Momentum aliasing is irrelevant here: σ_p = 1/(√2σ) = 0.141 / 0.118 against
k_Nyq = π/0.40 = 7.85.

## 5. Run matrix (18 production runs + controls)

Per σ ∈ {5, 6}:

| half | runs | name |
|---|---|---|
| wavepacket | 4 velocities | `s5p0_v2p0` … `s6p0_v3p5` |
| classical (CAP on) | 4 velocities | `cl_s5p0_v2p0` … |
| classical (CAP off) | v = 3.0 only | `cl_nocap_s5p0_v3p0` |
| vacuum CAP control | 4 velocities | `vac_s5p0_v2p0` … |

Plus one ground state `shared_gs/slab_n100_L35x35x105_dx0p4_per2`.

### Step counts

Inherited convention, calibrated on the existing campaign (3623 steps at v = 2.0,
launch −24, L_z = 85):

    N_STEPS = round( 4.36 * (|launch_z| + L_z/2) / (v * dt) ),  dt = 0.04

(reproduces 3624 vs the recorded 3623). With launch −27.5 and L_z = 105:

| idx | v | N_steps | t (a.u.) | save/ | wf/ | ckpt/ |
|---|---|---|---|---|---|---|
| 0 | 2.0 | 4360 | 174.4 | 14 | 43 | 872 |
| 1 | 2.5 | 3488 | 139.5 | 12 | 35 | 698 |
| 2 | 3.0 | 2907 | 116.3 | 10 | 29 | 581 |
| 3 | 3.5 | 2491 | 99.6 | 8 | 25 | 498 |

The packet exits the slab by t = 20 a.u. and reaches the CAP by t ≈ 34 a.u. at
the slowest velocity, so the remaining ~140 a.u. is plateau time for the deposit
estimator. `ckpt/` gives 5 retained numbered checkpoints per run plus the rolling
`checkpoint` that `LJ_RESUME=1` loads.

### Cost projection — WARN, not a gate (`.claude/rules/checkpoint-dont-block.md`)

Grid grows 88×88×213 → 88×88×264 (×1.24), so ~3.4 s/step against the measured
2.75 s/step at L_z = 85. Estimated **≈ 60 GPU-h total** (WP ≈ 25 h, classical
≈ 21 h, vacuum controls ≈ 8 h, CAP-free controls ≈ 5 h, GS ≈ 1–2 h), i.e. ~15–20 h
wall clock with four GPUs per sweep stage. Every run checkpoints every ~N/5 steps
and supports `*_RESUME=1`, so a kill costs at most one checkpoint interval. This
is reported to the user, not used to block the launch.

## 6. Estimator

    S = [E_total(t_final) - E_GS] / L_slab_z ,   L_slab_z = 25 Bohr

referenced to the **L_z = 105 production GS** (not the L_z = 85 value — the
deposit must reference the GS the run actually started from).

For the WP half apply the norm correction `E_total - T1*(1 - norm_WP)`: INQ
reports the kinetic term as occ·⟨ψ|T|ψ⟩/⟨ψ|ψ⟩ (`energy.hpp:50-55`), so under a CAP
the decaying WP orbital keeps contributing its per-particle mean and inflates
E_total. Use the **real-space** norm (`wp_real_space_stats.norm_check`), not the
momentum-space Parseval constant.

The classical half has no WP orbital, so no norm correction applies — but with the
CAP now ON, its E_absorbed is likewise a retained-excitation lower bound. The
CAP-free control at v = 3.0 measures the gap.

## 7. Validation gates

**Hard (abort):**
- GS: `∫n dV = 100.000`, r_s = 4.183 ± 0.001, SCF converged, no NaN.
- WP t=0 smoke gates (already in `wp/run.cpp`): σ_pz² vs 1/(2σ²), T1−T2 = 3/(4σ²),
  ⟨p_z⟩ = k₀, orthogonalisation-removed weight < 3 % (`LJ_ORTHO_TOL_PC`),
  real-space width consistency.
- Ledger closure every step: `E_SS + E_PS + E_PP = energy_hartree` and
  `E_SB + E_PB = energy_external` (WP); `E_SS = energy_hartree` and
  `E_SB + E_PS = energy_external` (classical). See
  `.claude/rules/decomposed-interaction-energies.md`.
- `deposit_stopping()` must report `complete = True` — it once returned a
  plausible S from an 86-of-3623-step run.

**Informational (report, never abort):**
- E_GS(L_z=105) vs E_GS(L_z=85) = 207.18323 Ha. Adding vacuum should move it very
  little, but a finer/larger box legitimately shifts it, so this is not a gate.
- t=0 CAP loss (predicted 0.020 % / 0.16 %), reproduced by the vacuum controls.
- Velocity drift of the classical projectile — it is a light Ehrenfest particle
  and is *supposed* to decelerate (`.claude/rules/light-projectile-stopping.md`).

## 8. Deliverables

1. `shared/configs/slab_n100_L35x35x105.hpp`
2. `scripts/sigma56_sv/{gs,wp,classical,vac}/run.cpp`
3. `shared/bin/run-s56-{gs,wp,cl,vac,notebooks}.slurm` + `submit-sigma56-sv.sh`
4. `hypotheses/sigma56_sv/` — `s56_stopping.py`, `build_run_notebooks.py`,
   per-run notebooks (with the mandatory density-matrix GIF,
   `.claude/rules/notebook-density-gif.md`), synthesis
5. **The combined S(v) figure**, matching the design of
   `hypotheses/classical_highdensity_sv/dyn_direct/S_of_v_v2_timeavg_sigmar.png`

## 9. Known blocker

`S_of_v_v2_timeavg_sigmar.png` **and its plotting script are not on this
machine** — `hypotheses/classical_highdensity_sv/dyn_direct/` holds only
`S_of_v_direct.csv` and two notebook builders, and
`scripts/classical_highdensity_sv/{dyn,dyn_direct}/results/` are empty. The user
must transfer them before the final figure can match the reference design. This
blocks deliverable 5 only; every run and every other deliverable proceeds.

## 10. Explicitly out of scope

- Re-running σ = 0.5/2/3 at L_z = 105 (user chose the heterogeneous figure with a
  caption caveat).
- Classical twins at σ = 2 and 3.
- v = 4.0 and 4.5 (aliasing-free at these widths and therefore recoverable, but
  the grid is held at 4 points for point-for-point comparability).
- Reproducing the σ = 0.5 classical benchmark — its raw data was lost with the
  `/local/data/public` machine.
