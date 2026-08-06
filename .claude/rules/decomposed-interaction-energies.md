# Rule: Every run emits the decomposed pairwise interaction energies

Apply to: EVERY TDDFT run definition with a projectile, in EVERY system —
`ResearchProject/systems/**/scripts/**/run.cpp` (jellium, localised_jellium,
cylindrical_jellium, vacuum, coronene, and any system added later). Always on.
(User decision, 2026-08-01: "I want to ensure that all future runs, does not
matter which system, has these decomposed energies.")

## The one rule

**Every run writes `results/raw/observables/interactions.csv` carrying the
pairwise electrostatic decomposition into the three charge groups —
P (projectile), S (system/bath electrons), B (neutralising background) —
alongside the INQ scalars.**

```
E_SS = ½∫n_S·φ_S     bath-bath
E_PP = ½∫n_P·φ_P     projectile SELF-Hartree      <- the quantum residual
E_PS =  ∫n_S·φ_P     projectile-bath              <- the interaction that stops it
E_SB = -∫n_S·φ₊      bath-background
E_PB = -∫n_P·φ₊      projectile-background
E_BB = ½∫n₊·φ₊       background self (constant, compute once)
```

Header: `inq-stack/include/inqkit/jellium/interaction_energies.hpp`
(`compute_coulomb` for a classical projectile, `compute_coulomb_wp` for a
wavepacket, `orbital_density_field`, `background_self_energy`). Needs
`inqkit/jellium/projectile_background_energy.hpp` for `gaussian_density`.

## Why

INQ's own scalars **cannot** answer the question this decomposition answers,
because the two projectile representations sit in DIFFERENT ledger terms:

| | classical | wavepacket |
|---|---|---|
| projectile enters as | external potential | occupied KS orbital |
| `energy_external` | non-zero | **identically 0** |
| `energy_hartree` | = E_SS | = E_SS + E_PS + E_PP |

So a raw `energy_kinetic` or `energy_hartree` comparison between a classical and
a wavepacket run compares a NET quantity against a GROSS one. The pairwise terms
are representation-independent and ARE comparable.

**E_PP is the projectile self-Hartree — the uncancelled self-interaction of a
wavepacket in LDA.** It has no classical counterpart and is the leading candidate
for the unexplained factor ~2.2 in the bulk-jellium classical/WP stopping ratio
(see `docs/handovers/bulk-jellium-ks-stopping.md`). Without this decomposition it
is not measurable at all; with it, it is a column in a CSV.

## Closure gates (mandatory — they make the terms trustworthy)

The terms must sum back to the INQ scalars. Gate on this, do not assume it:

- classical: `E_SS == energy_hartree`  and  `E_SB + E_PS == energy_external`
- wavepacket: `E_SS + E_PS + E_PP == energy_hartree`  and
  `E_SB + E_PB == energy_external`

`compute_coulomb_wp` returns `e_hartree_check` / `e_external_check` for exactly
this. Verified 2026-08-01 against a completed bulk run: offline E_SS matched
INQ `energy_hartree` to **1.4e-17 Ha**.

## BULK vs SLAB — the background terms differ, and this trips people

- **Slab / localised systems** have an explicit background perturbation, so
  `φ₊ = poisson(n₊)` is non-trivial and E_SB, E_PB, E_BB are all meaningful.
  Build `n₊` via `bg_pert.background_density(basis)`.
- **Bulk jellium** has a UNIFORM background, so `poisson(n₊)` is pure G=0 — which
  INQ drops. **`φ₊` is IDENTICALLY ZERO and E_SB = E_PB = E_BB = 0.** Pass a
  zero-filled field as `phiplus`; still write the columns (as zeros) so the schema
  is identical across systems. In bulk the physics is entirely in E_SS, E_PP, E_PS.

**Gauge caveat:** absolute E_SB / E_PB / E_BB carry the charged-cell G=0 gauge.
Only closure sums and WP-minus-classical DIFFERENCES are gauge-clean. Never quote
an absolute E_PB across systems.

## How to apply

- **New run.cpp:** clone the wiring from a reference —
  `systems/jellium/scripts/bulk_ks_stopping_sigma3/wp/run.cpp` (WP, bulk,
  zero φ₊) and `.../bulk_ks_stopping_sigma3/classical/run.cpp` (classical, bulk,
  `n_P` rebuilt from the ion position each step at `sigma_pot = sigma_WP/√2`);
  or `systems/localised_jellium/scripts/localised_jellium_dynamics/{phase5_wp,proj_dyn}/run.cpp`
  for the slab case with a real background.
- **Cadence:** write it EVERY callback. Two Poisson solves per step is negligible
  against the propagator's per-orbital FFTs — do not thin it to the VTI cadence.
- **Classical `n_P`:** the projectile is not in `n`, so its charge cloud must be
  rebuilt at the CURRENT ion position each step via
  `gaussian_density(basis, center, sigma_pot)`. `sigma_pot = sigma_WP/√2`
  (`.claude/rules/sigma-wp-convention.md`) — using `sigma_WP` directly is wrong by
  √2 and will silently mis-scale E_PP and E_PS.
- **Analysis:** `interactions.csv` sits beside `observables.csv`; concatenate
  segment-suffixed files on resume like every other observable.

## Retrofit status (2026-08-01)

Wired: `bulk_ks_stopping_sigma3/{wp,classical}` (re-running, jobs 32512952/3).
NOT yet wired: the other 10 bulk_ks_stopping runs — their existing results lack
these columns and would need a re-run to gain them. Wire on next touch.
