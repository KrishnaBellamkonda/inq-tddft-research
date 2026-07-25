# ADR 0008 — Localised jellium via a static background perturbation

Date: 2026-06-21
Status: accepted (validated — T0 + T1 + T3.4 pass)

## Context

The existing jellium in this repo is **delocalised**: `extra_electrons(N)` with
no ions, where INQ's periodic Poisson solver drops the G=0 Fourier component —
mathematically identical to a uniform positive background spread over the *entire*
cell. This cannot be confined to a region by any parameter change.

We need a **localised** jellium: a finite positive background `n₊(r)` (slab,
sphere, or box) into which a projectile is fired from surrounding vacuum, for
finite-target scattering / stopping-power studies and Lang–Kohn / cluster
benchmarks. The background must be present during the SCF (so electrons localise)
and during real-time propagation (so the target persists while the projectile
flies). `inq/` is immutable.

Two mechanisms were genuinely available:

- **A. Static custom perturbation.** A new `inqkit` class implements INQ's
  perturbation duck-type; its `.potential(t,v)` adds `v_bg = −poisson(n₊)` to the
  KS potential. Both `ground_state::calculate` and `real_time::propagate` accept
  the perturbation, so one object covers SCF and RT.
- **B. Smeared positive pseudo-ions.** Represent `n₊` as many positive Gaussian
  pseudo-charges and let INQ build `vion`.

## Decision

Adopt **Option A**. The localised background is the `inqkit` pair
`jellium/localised_background.hpp` (builds `n₊`) + `jellium/background_
perturbation.hpp` (caches `φ=poisson(n₊)`, adds `−φ` to the KS potential via an
explicit `gpu::run` loop so it works for both the real `inq` potential and the
complex `inq-study` potential). Charge neutrality is the caller's job: set
`extra_electrons(N)` and `n₀ = N/V_inside` so `∫n₊ = N` exactly, which makes the
dropped G=0 cancel exactly (`v_es = poisson(n_elec − n₊)`). This is mathematically
the GPAW jellium recipe. No `inq/` or `inq-study/` edit.

## Consequences

- A true flat-top interior at `n₀` is achievable (Option B's sum-of-Gaussians
  cannot), which the "interior within a few % of n₀" benchmark requires.
- `∫v_bg·n` (electron–background attraction) is folded into INQ's reported energy
  automatically; only `E_self` (background self-energy) must be added by hand,
  because INQ never sees `n₊` as a charge.
- The same object composes with the sin² CAP via `perturbations::sum` for Phase 5.
- **Validated 2026-06-21:** T0 (∫n₊=N slab+sphere; well attractive) passes; the
  slab GS binds electrons inside the slab with interior density flat to 2.0% of
  n₀ (T1) and kinetic energy/electron matching the HEG value at r_s≈4; a 2 au
  static run conserves total energy to 2×10⁻⁸ Ha (T3.4), confirming the
  perturbation is correctly static/Hermitian.
- Reversal cost is high (all run machinery + configs build on this API), which is
  why this is recorded as an ADR.
