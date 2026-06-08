# Rule: jellium base-run specification

Apply to: `ResearchProject/systems/jellium/`

## Standard simulation duration

The **N=138 closed-shell base run at L=30 Bohr cubic, WP_EKIN = 5 eV,
SPACING = 0.85 Bohr, N_STEPS = 990 at dt = 0.020 a.u.** is the canonical
base configuration for the jellium WP–jellium scattering project, set
2026-05-04. Any new variation (different WP energy, different shell
N, different envelope σ, tilted launch, etc.) **must reproduce the same
physical propagation time** of `dt_au × N_STEPS = 19.8 a.u. ≈ 0.479 fs`
unless there is an explicit physics reason to deviate (e.g. a slower
WP that needs a longer flight time, justified in the run-specific
config header).

This duration is chosen to give a single-pass trajectory of
~12 Bohr at k₀ ≈ 0.606 Bohr⁻¹ (5 eV) inside an L = 30 Bohr cubic
box — about 3 Bohr clear of the periodic boundary at L/2 = 15 Bohr
once a 3 σ_r WP envelope is accounted for.

## Why this specific duration

- The WP energy must lie in the **same scale as the lowest electron-hole
  excitations** of the bath for the two systems to couple inelastically.
  At L=30, the kinetic-only ΔE for the |G|²=6 → |G|²=8 transition is
  ≈ 1.19 eV — comparable to the 5 eV WP. This is the regime where
  single-particle excitations dominate the WP slowdown signal.
- A 100 eV WP (the previous run) lives ~80× above the e-h gap and
  couples mostly to plasmonic/collective channels, giving only a few
  per-cent retardation that is easily lost in numerics.
- The flight time of 19.8 a.u. is enough to see several inelastic
  scattering events at this energy scale while still finishing in
  ~25-40 min wall on a single A30 GPU.

## Required Cfg fields for new variants

Every Cfg derived from `Base_N138_L30` must explicitly declare:

- `L_BOHR = 30.0` — *if* the variant changes box length, justify in
  comment (e.g. higher-r_s study).
- `N_ELECTRONS` — must be a closed-shell magic number for the chosen
  L_BOHR (see `docs/sources/free-electron-gas-magic-numbers.md`). Default
  138 unless deliberately exploring open-shell physics.
- `EXTRA_STATES ≥ 20` — covers the next ≥ 2 unoccupied shells so
  that `gamma_transitions` postprocess can build a non-empty occ → unocc
  set (the threshold filter at occ < 0.01 needs states *clearly* above
  Fermi smearing). Below 20 is permitted only if the variant is purely
  about the WP and not about bath e-h dynamics.
- `WP_EKIN_EV` — the only field that should change between variants in
  a "WP-energy sweep". Companion variants live in
  `shared/configs/base_n138_L30_E<NN>.hpp`.
- `SPACING_BOHR` — must satisfy the Nyquist condition `dx ≤ π/(k₀ + 3σ_k)`
  where `σ_k = 1/σ_r` is the WP momentum-space width (Heisenberg).
  Document the calculation in the Cfg header.
- `N_STEPS` — must be retuned per WP energy so the WP traverses
  `~L/2 - 3 σ_r` Bohr in the trajectory. Document the calculation.

## Verification expectations

Every base-run-derived simulation must pass the standard verification
script (`scripts/verify_smoke_outputs.py`) and produce:

- `cod_x ≈ cod_y ≈ 0` (Cartesian frame is centred on the box).
- `cod_z` increasing monotonically with slope close to k₀ at early times.
- Energy drift `< 1 mHa` over the trajectory.
- `density_l2(0) = 0` by construction; subsequent values populated.
- WP momentum peak at `|k| ≈ k₀` at t = 0.

## File-placement reminders

- Configs: `shared/configs/base_n138_L30_E<NN>.hpp` (`Base_N138_L30_E<NN>`
  struct).
- GS saves: `save_gs/gs_L<L>_cubic_N<N>_dx<dx>/` — a fresh GS is needed
  for each new (L, N, SPACING) triple. The same checkpoint can be
  reused across WP-energy variants that share L/N/dx.
- Run dirs: `run_base_n138_L<L>_E<NN>/` for the canonical line; topic-
  specific deviations (tilted, σ-sweep, etc.) get their own
  `run_<topic>/` dirs as before.

## When to deviate

- Higher-r_s study: change L_BOHR (and so density), keep N=138.
  Justify in the Cfg header (e.g. "preserve r_s near sodium").
- Open-shell study: pick a non-magic N near 138 (e.g. N=140); fully
  document that the GS density is no longer flat-by-symmetry and that
  any "drag" signal mixes bath response with broken-symmetry artefacts.
- WP-energy sweep: vary `WP_EKIN_EV` only; SPACING and N_STEPS retuned
  per the formulas above. Companion runs must be in the same
  `run_base_n138_L30_E*/` family for direct comparison.
