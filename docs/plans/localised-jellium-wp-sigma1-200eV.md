# DRAFT PLAN: main WP projectile on localised jellium — σ=1, 200 eV (normal electron)

**Status: DRAFT — DO NOT EXECUTE YET.** One physics decision is open (launch
standoff, below). Drafted 2026-07-08 at user request ("have this as a draft plan;
let's do this later"). Reverts to the **standard normal-electron wavepacket** run
(NOT the effective-mass fork).

## Why this exists

The effective-mass-fork run (`effmass_12h`, m=2.506 mₑ, 2026-07-08, executed,
notebook `hypotheses/muon_mass_fork/effmass_12h_quantum_run.ipynb`) is **not
trusted** by the user. Decision: **demote the mass-fork run to a LATER campaign
phase** — the mass parameter must be "figured out correctly to fit the needs"
first. For now, revert to the main WP run we have been using (a normal electron,
m=1) with a fresh parameterisation.

## Locked spec (user, 2026-07-08)

| Parameter | Value |
|---|---|
| Projectile | **normal electron wavepacket, m=1** (no fork) |
| σ_WP | **1.0 Bohr** |
| Energy | **200 eV** → k₀ = 3.834, v = 3.834 a.u. |
| CAP η | **−1.0** |
| Total time | **3 × traversal** |
| Checkpoints | **every 1 traversal → 3 checkpoints** (at τ, 2τ, 3τ) |
| Max pre-impact spread | **4%** |
| Grid dx | **0.40 Bohr** |
| Slab thickness | **25 Bohr** (half_width 12.5) |
| Cell | **50 × 50 × L_z**, L_z = 80 (verified) |
| Density (assumed) | r_s = 5.665, **N = 82** (the established main-run SlabN82 density; CONFIRM) |

## Computed sizing (verified 2026-07-08)

- k₀ = 3.834, v = 3.834 a.u.; σ_ρ0 = σ_WP/√2 = 0.707; σ_p = 1/(√2σ_WP) = 0.707.
- **Aliasing @ dx=0.4:** k₀+3σ_p = 5.96 < k_max = π/0.4 = 7.85 ✓ (margin 1.9).
- **dt = 0.05** (E_cut = ½k_max² = 30.8; H·dt = 1.54, below the 2.2 ETRS cliff;
  validated stable at dx=0.4 in the effmass smoke test). **Smoke-test before launch.**
- **τ (one traversal, launch→far slab face) = (D_launch+25)/v ≈ 6.78 a.u.**
- Total = 3τ ≈ 20.3 a.u. → **N_STEPS = 407**.
- **CKPT_EVERY = round(τ/dt) = 136** → checkpoints at steps 136 / 272 / 407
  (t = 6.8 / 13.6 / 20.3 a.u.) = **3 checkpoints** ✓ (matches the run.cpp
  checkpoint logic: intermediate saves at k·CKPT_EVERY < N_STEPS + a final save at N_STEPS).
- Runtime ≈ 35–45 min (dx=0.4, ~50×50×80 grid, N=82 ~ few hundred states → ~5–7 s/step).

## L_z = 80: VERIFIED sufficient ✓

Launch z₀ ≈ −13.5 (1 Bohr before the −12.5 face). WP t=0 extent [−16.3, −10.7].
+z CAP inner face ~ +24 → **11.5 Bohr clean exit runway** past the +12.5 far face.
WP exits the slab at t=6.8, reaches the CAP only at t=9.8 → clean crossing; no t=0
CAP overlap. L_z=90 adds runway (14.5 Bohr) that is not needed here. **Use 80.**

## OPEN DECISION — launch standoff vs the 4% spread cap

A σ_WP=1 packet spreads fast (σ_p=0.71): it reaches 4% spread in only **1.1 Bohr**
of travel, but a clean vacuum standoff needs ~4σ_ρ0 ≈ 2.8 Bohr (24% spread there).
The two cannot both be met for σ=1 at 200 eV. Options put to the user (undecided):

- **A — keep σ=1, soft entry:** launch 1.0 Bohr, spread 3.3%, **~6% of the WP
  inside the slab at t=0**. Honours the exact spec; WP not cleanly vacuum-incident.
- **B — widen to σ=1.5:** launch 2.5 Bohr (2.3σ), ~4% spread, ~1% t=0 overlap.
  Clean incidence; σ ≠ 1.
- **C — keep σ=1, launch 2.8 Bohr (4σ):** ~0% overlap but ~24% spread (violates
  the 4% cap).

**Resolve this before building** (it sets σ_WP and the launch z₀).

## Machinery plan (when executed)

1. **Base on the main standard WP run**, e.g.
   `scripts/qsp_phase4/wp/run.cpp` — it writes the FULL `momentum_distribution.csv`
   (1D n(k)), so the run-notebook gets native n(k)/stopping panels (richer than the
   effmass moments schema). NOT the effmass_12h run.cpp.
2. **New config** `shared/configs/slab_n82_L50x50x80.hpp` (50×50×80, dx grid via
   run, N=82, r_s=5.665, half_width 12.5, EDGE 0). New **GS** at dx=0.40,
   50×50×80 → `shared_gs/slab_n82_L50x50x80_dx0p40`.
3. **Add checkpoint/resume** (the validated `EM_CKPT_EVERY` / `EM_RESUME` +
   native `start_step` mechanism from `effmass_12h/quantum/run.cpp`) to the WP
   run.cpp. Set CKPT_EVERY = 136 (τ/dt) for the 3-checkpoint cadence.
4. CAP η = −1.0, retuned for L_z=80 (mid |z|≈32, inner face ~24).
5. **Smoke-test dt=0.05** (stability + resume bit-faithfulness) before the full run.
6. Launch (~40 min), build the run-notebook (`run-notebook` skill; standard-schema
   run → full battery incl. momentum n(k)).

## Related

- Effective-mass fork run → **later campaign phase** (mass parameter TBD). Its
  spec/plan: `docs/plans/muon-effmass-12h-run.md`; notebook:
  `hypotheses/muon_mass_fork/effmass_12h_quantum_run.ipynb` (energy-method stopping
  panel confounded by CAP absorption — see that notebook's caveat).
- `run-notebook` builder was fixed (detect_type now recognises the `wp_momentum_stats`
  mass-fork schema) — unrelated to this run but recorded.
