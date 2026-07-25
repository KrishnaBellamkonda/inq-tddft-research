# Phase 3 pilot — Run A (perturbation Ehrenfest, v=2, r_s=4.18 slab)

Notebook: `pilot_run_notebook.ipynb` (executed, density GIF embedded).
Run: `scripts/classical_highdensity_sv/pilot/` (1600 steps, dt=0.04, analytic force).

## Result

**1. Central aim — ACHIEVED.**
- Projectile TRANSITED: v 2.0 → 1.40, proj_z −30 → **+70.9** (28 Bohr past the +42.5
  far face) — did NOT stop inside. So v=2 is above the transit floor.
- E_electronic **plateaus perfectly flat after exit** (last-15% ΔE_total std = 0.0000 eV)
  — **no oscillation**. Energy conserved (E_elec+KE_proj+U_proj_bg drift 0.87 eV).
- ⇒ the z-open (`periodicity(2)`) + CAP-free design WORKS. The historical
  non-plateau/oscillation is defeated at a real velocity.

**2. ⚠️ Methodology problem the pilot EXPOSED — E_absorbed = ΔE_total is gauge-broken.**
- ΔE_total = **+445.7 eV** → `S = ΔE/L = 17.8 eV/Bohr` — **UNPHYSICAL** (~20× too high).
- Cause = charged-cell G=0 gauge (the −1 projectile makes the cell net-charged), exactly
  per `reference_charged_cell_hartree_convention`. The projectile's own interaction
  bookkeeping swings by hundreds of eV as it moves and dominates ΔE_total:
  ΔU_proj_bg = −418.8 eV, ΔE_PS = +418.8, ΔE_PB = −418.8 eV, while the projectile only
  lost **−ΔKE_proj = 27.8 eV**.
- ⇒ **Locked design decision 8 ("E_absorbed = ΔE_total trivially clean") is WRONG in
  practice.** The E_absorbed/L headline (Definition 2) cannot be read from raw E_total.

**3. Gauge-clean classical S (the physical answer).**
- Projectile **KE loss across the slab** (z∈[−12.5,+12.5]): 55.16 → 32.04 eV = **23.1 eV**
  over 25 Bohr → **S = 0.93 eV/Bohr** (in-slab −dKE/ds slope = 0.98; mean v = 1.80).
- Sits between Lindhard-point (~0.57) and bulk σ=0.5 (~0.94) at r_s=4.18 — physically sensible.
- Valid classically; **NOT** available for the WP (`feedback_quantum_stopping_not_from_projectile_ke`).

## RESOLVED (2026-07-22): use the GS baseline — Definition 2 works, sweep can proceed

The "gauge-broken E_absorbed" above was a **BASELINE ERROR**, not a broken definition.
The charged-cell G=0 gauge only contaminates E_total **while the −1 projectile is IN the
box**. Both clean endpoints are NEUTRAL-cell:
- GS (no projectile) E_GS = 207.1832 Ha
- post-exit plateau (projectile clipped away at z=70.9) = 208.1715 Ha

**E_absorbed = E_total(plateau) − E_GS = 0.988 Ha = 26.9 eV → S = 1.08 eV/Bohr** (gauge-clean),
matching −ΔKE_proj (27.8 eV, within the 0.86 eV drift) and the in-slab KE-loss S (0.93).
The earlier +445 eV came from the WRONG baseline E_total(0), where the projectile sits at
z=−30 in a charged cell (U_proj_bg/E_PS ±419 eV).

**This E_absorbed method is WP-transferable** (the WP's E_total(plateau)−E_GS is likewise
neutral-cell after the packet exits/absorbs), so the classical↔quantum shared metric is
intact. **Definition 2 = [E_total(plateau) − E_GS] / L_slab, read after full exit.**

⇒ **SWEEP UNBLOCKED.** Run the 6 velocities; extract S via the GS-baseline plateau;
carry the pairwise `interactions.csv` ledger for the still-TBD Definition-1 formula.
Refinement to weigh (non-blocking): dx=0.4 for the fast points (perturbation vs native
∇n differ 1.5%→0.07% there); the pilot dx=0.5 is fine for the floor velocity.

---

## RUN B — native Ehrenfest, real ghost-UPF ion — DISQUALIFIED for this geometry

(Added 2026-07-22 by the Phase-3 pilot agent; notebook
`hypotheses/classical_highdensity_sv/pilot_native/pilot_native_pilot.ipynb`, executed,
density GIF embedded. Contrast run to Run A: real mass-1 ghost-UPF ion moved by INQ's OWN
native Ehrenfest, same slab, own GS with the ghost present, E_GS=217.3 Ha.)

- **Does native Ehrenfest MOVE a z_valence=0 ghost ion? YES** — it is advanced by its
  local HF force. (Settles the open question raised in the campaign resolved-decisions.)
- **But the trajectory is UNPHYSICAL.** Launched z=−30, vz=2 toward the slab, the ghost
  ion decelerates to a stop within ~4 Bohr — **still in vacuum, ~14 Bohr from the slab** —
  then **REVERSES and oscillates** between z≈−58 and −26. It **NEVER reaches the slab**
  (max z=−26.2; near face −12.5). E_total is **not conserved** (swings ~7 Ha ≈ 180 eV).
- **Cause — ghost-UPF long-range +1/r tail** (`reference_ghost_upf_tail_aliasing`): the
  real ghost UPF carries a bare, UNSCREENED +1/r Coulomb tail. In the z-open/xy-periodic
  box (launch z=−30 is only ~12 Bohr from the −42.5 edge) this produces a large spurious
  vacuum force (slab-across-box + periodic images). RUN A avoids it entirely: its drag is
  from (electrons − background) only, so the projectile's own long-range self-field never
  enters the force.

**Native-vs-perturbation verdict:** the perturbation (analytic-force) projectile is
**validated and required** for the production S(v) sweep; the native ghost-UPF ion is
**disqualified** for this geometry. This retires the native-force faithfulness question.
(Note this is a SEPARATE issue from the gauge flaw in item 2 above — the perturbation
Run A trajectory + KE-loss S are correct; only the raw-ΔE_total *metric* is gauge-broken.)
