# Rule: Light projectiles decelerate — size runs and extract S accordingly

Apply to: any electronic-stopping-power / projectile run with a LIGHT projectile
(classical electron of mass m_e, or an electron wavepacket) under FREE EHRENFEST
dynamics — `ResearchProject/systems/**/scripts/**/run.cpp`, run `analyse.py`,
orchestrators (`orchestrate.py`/`finalize.py`), `docs/campaigns/**`, and any S(v)
extraction or pilot/validation gate. Always on.

## The one rule

**A light free-Ehrenfest projectile has tiny kinetic energy (KE = ½·m·v²) and
DECELERATES strongly — it stops within a few Bohr, depositing all its KE, BEFORE a
steady-state wake (5·2π/ω_p) can ever form.** Design the run and the analysis for a
decelerating projectile, not a constant-velocity one.

1. **Do NOT size runs by the wake/traversal criterion** `max(1.5·L_z/v, 5·2π/ω_p)/dt`
   for a light projectile — it yields 10k+ steps the projectile never survives (it
   stopped ~1500 steps in, then sits). Size to capture the **initial-drag window +
   deceleration sweep** (~30–45 a.u.; `ceil(max(30, 100·v)/dt)`).
2. **Extract S(v₀) as the INITIAL drag** = `−d(KE_ion)/ds` over the early
   **near-constant-velocity window** (v ≥ ~0.85·v₀; widen to 0.70/0.50 if sparse),
   from the per-step track. A full-run regression is WRONG — it averages S over the
   whole decelerating velocity range (v₀ down to ~0), not S *at* v₀.
3. **Never abort a pilot/validation gate on the (by-design) velocity drift.** A light
   Ehrenfest projectile is SUPPOSED to decelerate (the user locked free Ehrenfest).
   Gate on a *clean initial-drag slope existing* (finite S, ≥~30 early-window points),
   not on `v-drift < X%`.
4. The decelerating projectile is a feature: one run sweeps a velocity range, so the
   initial KE-loss slope is the friction force at v₀.

## Why

This exact mistake aborted the first overnight cylindrical-jellium S(v) run
(2026-06-28): the pilot gate failed on an 85% v-drift (vz 0.300→0.045 over ~6 Bohr,
real physics), and a 9.5-hour finalizer then polled for production that never ran.
The cause was treating a fast-stopping light electron as a constant-velocity heavy
ion. Heavy projectiles (proton, m≈1836) barely decelerate → constant-velocity S(v)
→ the wake criterion DOES apply; the distinction is the projectile's KE (mass×v²)
versus the stopping power.

## How to apply

- Before sizing a light-projectile run: confirm the velocity stays near v₀ only for
  the early window; size N_STEPS to capture it, not 5 plasma periods.
- In S-extraction: window on `vz ≥ 0.85·v₀`, regress `ke_ion_ha` (or ΔE_system) vs
  path; report S(v₀) with its mean v and point count.
- In gates: replace any `v-drift` abort with a "clean initial-drag S" check.
- Reference implementation:
  `ResearchProject/systems/cylindrical_jellium/scripts/annular_sv/orchestrate.py`
  (`n_steps_for`, `extract_S`, `pilot_gate`). See
  [[reference_light_projectile_deceleration]] and [[reference_cutoff_aliasing_guard]].
