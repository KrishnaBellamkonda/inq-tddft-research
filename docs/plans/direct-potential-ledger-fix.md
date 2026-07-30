# Plan: direct-potential projectile — fix the diagnostic ledger + kink-free replica

Status: in progress (2026-07-30)
Owner task: quantum-stopping-power / classical_highdensity_sv

## Problem (diagnosed, evidence in conversation + data)

The classical projectile is a moving external **potential**, but the code routes it
through a **charge → periodic Poisson** solve. Two artifacts follow, both confirmed
in the v-sweep data:

1. **Exit kink.** `gaussian_density` uses absolute distance (no wrap); as the center
   crosses the z-wall (±42.5) the charge clips off the grid (`norm_proj` 1→0 at
   z≈42.5) and the Poisson potential of the clipped, net-charged density lurches →
   kink in `energy_external`/`energy_total` and in the pairwise terms.
2. **e_ps drift (not a gauge).** The projectile is the only net charge; replicated by
   the x,y-periodic Poisson kernel it becomes a **charged sheet** whose potential
   grows linearly with distance. Measured: recorded `e_ps` is linear in z, slope
   −0.35 Ha/Bohr vs sheet prediction −0.51 (R²=0.89), reaching −16 at z=30 instead
   of →0. The **localised** erf/r field on the same density gives `e_ps` positive and
   →0 (7.4 in-slab → 3.1 at z=30). Density depletes at the projectile (repulsion) →
   applied sign is correct; do NOT reverse it.

## What already exists (prior work)

- `inqkit/jellium/gaussian_potential.hpp` — direct erf/r field (no charge/Poisson).
- `inqkit/dynamics/moving_gaussian_projectile_potential.hpp` — direct perturbation.
- `inqkit/dynamics/projectile_force_direct_z` — direct HF force.
- `scripts/classical_highdensity_sv/pilot_direct/run.cpp` — uses the direct
  perturbation + direct force. **energy_total is clean at the wall** (verified).

## Remaining bug

`pilot_direct/run.cpp` still computes the **diagnostic ledger** from the clipping
charge: `compute_coulomb(density, nproj, phiplus)` and
`Uprojbg = -∫nproj·phiplus`. So `e_pp`, `e_ps`, `e_pb`, `energy_proj_bg_ideal`
STILL kink at z≈42 and `e_ps` still drifts to −21 (verified). This is exactly the
symptom the user flagged.

## Fix

1. `interaction_energies.hpp`: add `compute_coulomb_direct(n_slab, v_proj, nplus,
   phiplus, sigma_pot)` — projectile terms from the DIRECT erf/r potential
   `v_proj = gaussian_potential(...)`:
   - `e_pp = 1/(2·σ_pot·√π)` (analytic constant self-energy; no kink)
   - `e_ps = +∫ n_slab·v_proj`  (positive; →0)
   - `e_pb = −∫ nplus·v_proj`   (attraction; →0)
   - `e_ss = ½∫ n_slab·φ_slab`, `e_sb = −∫ n_slab·φ₊`  (slab; Poisson, unchanged)
   Closure preserved: classical E_external(INQ)=e_sb+e_ps, E_hartree(INQ)=e_ss.
2. `scripts/classical_highdensity_sv/dyn_direct/run.cpp` — clone `pilot_direct`
   but use `compute_coulomb_direct` for the ledger and `Uprojbg = e_pb`.
3. Build via `inq-run` from `dyn_direct/`.

## Validation (code-test + simulation-validation)

- **T1 (host, pure):** analytic self-energy `1/(2σ√π)` == numeric Gaussian
  self-Coulomb (radial quadrature). Guards the hardcoded constant.
- **T2 (physical regression):** replica `interactions.csv` `e_ps` is **positive and
  monotonically →0** as projectile recedes; cross-check against the python localised
  ∫n·erf/r (7.4 in-slab → ~3 at z=30). `e_pp` ≈ 0.80 Ha constant.
- **T3 (no-kink):** `energy_total`, `e_ps`, `e_pb`, `energy_proj_bg_ideal` have NO
  curvature spike at z≈42.5 (contrast old v4p5).
- **T4 (differential):** replica vs old v4p5 agree in-slab (trajectory, KE-loss, S)
  — the fix must not move the physics where the old code was valid.
- **T5 (conservation):** (E_elec + KE_proj + U_proj_bg) flat over the whole run.

## Replica run (fastest previous projectile, v=4.5)

Env = exact old v4p5: LX=LY=35 LZ=85 HALF=12.5 N=100 EDGE_W=1.0 PER=2 SPACING=0.5
SIGMA=0.5 MASS=1 DELTA=0.1 DT=0.04 LAUNCH_Z=−24 K0=4.5 N_STEPS=1074 SAVE_EVERY=4
GS=shared_gs/slab_n100_L35x35x85_dx0p5_per2. Detached (setsid). Then a run notebook
with the density GIF at top; expect NO kink.
