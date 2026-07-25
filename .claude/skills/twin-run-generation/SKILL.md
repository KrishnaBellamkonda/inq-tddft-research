---
name: twin-run-generation
description: Produce a matched classical(perturbation)+wavepacket TWIN PAIR of localised-jellium runs — identical in every physical parameter except the projectile representation, with the FULL energy decomposition switched on in both. Guarantees parity and observable-completeness via check_twin.py, writing twin_manifest.json (the contract the twin-run-analysis skill consumes). Composes tddft-simulations for the actual build/dispatch.
---

# Twin-Run Generation

Produces a **twin pair**: two localised-jellium runs that are byte-for-byte
identical in configuration EXCEPT the projectile representation —

- `classical/` — projectile = **Gaussian-charge perturbation** (`v_proj = +poisson(n_proj)`
  added to the KS potential; no ghost UPF, no r_cut aliasing). Preferred over the
  ghost-UPF classical because it explains more of the energy (see the handover).
- `wp/` — projectile = **wavepacket** (a real KS electron of width σ_WP, k0).

The pair is the sole input to `twin-run-analysis`. This skill's job is to make the
pair *correctly* and *provably*: same GS, same cell/N/spacing/σ_WP/launch_z/dt, and
— critically — the **full energy decomposition switched on in both runs**, so the
analysis never silently loses a term.

## Why this skill exists

Today, whether a run carries `energy_external` / `energy_proj_bg` depends on which
`ObservableSelection` flags someone flipped. Existing pairs are inconsistent (the
`h0_base_difference` runs emit only total/kinetic/hartree/xc — no external, no
proj_bg). This skill makes "the decomposition is complete and the twins match" a
checked guarantee, not a hope.

## The twin-run contract (what this skill must deliver)

```
<pair_name>/
  wp/         raw/observables/observables.csv   run_summary.txt
  classical/  raw/observables/observables.csv   run_summary.txt
  twin_manifest.json          # written by check_twin.py; valid=true is the gate
```

**Required `observables.csv` columns (both twins):**
`step, time_au, energy_total, energy_kinetic, energy_hartree, energy_xc, energy_external`
(plus `energy_nonlocal`; classical also `energy_ion`).

**Pairwise decomposition — both twins emit `interactions.csv`** (`step,time_au,e_ss,
e_pp,e_ps,e_sb,e_pb,e_bb,…`, Hartree) via `inqkit/jellium/interaction_energies.hpp`
(`compute_coulomb` classical / `compute_coulomb_wp`+`orbital_density_field` WP). These
close exactly to INQ's `E_hartree`/`E_external` and let the analysis skill attribute
every classical-vs-WP difference physically (`reference_twin_pairwise_decomposition`).

For a **dynamic** (Rung-2)
run the classical twin must additionally emit `energy_proj_bg_ideal` per step
(moving projectile → `U_proj_bg` is no longer constant); for a **static** (Rung-1)
run `U_proj_bg` in `run_summary.txt` suffices.

**Parity fields that MUST match** (`check_twin.py` asserts): `periodicity, Lz,
spacing, N, sigma_WP, launch_z, gs_dir`. The ONLY permitted difference is the
projectile block.

**σ convention (do not get this wrong):** the run is labelled by σ_WP; the
classical Gaussian **charge** std is `σ_pot = σ_WP/√2` (so the classical charge std
equals the WP density std). Both twins are reported at the same σ_WP. See
`reference_sigma_matching_convention`.

## Workflow

1. **Pick the config** — geometry/N/spacing/periodicity/Lz, σ_WP, launch_z (→ r
   from the slab face), dt, n_steps, and the **shared GS checkpoint** (both twins
   MUST load the same GS; produce it once if absent).
2. **Build the pair** via `tddft-simulations` (build-once binary + dispatch). Use
   the existing reference run.cpp:
   - classical: `…/localised_jellium_dynamics/proj_perturbation/run.cpp`
     (Gaussian-charge perturbation; writes `U_proj_bg` + full external ledger).
   - wavepacket: the matching `…/phase5_wp`-style WP run.cpp.
   Enable the **full `ObservableSelection`** (energy_external/nonlocal/ion on; for a
   dynamic classical run also `set_proj_bg` per step / the `energy_proj_bg_ideal`
   column). Pin GPU per the session rule (`CUDA_VISIBLE_DEVICES=0` when required).
3. **Gate with `check_twin.py`** (skill-local, stdlib-only):
   ```bash
   /local/data/public/skcb2/tddft/venv/bin/python3 \
     .claude/skills/twin-run-generation/check_twin.py <pair_dir>
   #   or: --wp DIR --classical DIR --pair-dir DIR
   ```
   It verifies completion, parity, that the projectile actually differs, the full
   decomposition columns, and `U_proj_bg` availability — then writes
   `twin_manifest.json`. **A non-zero exit means the pair is NOT ready for
   analysis.** Fix the config and re-dispatch; never hand a failing pair downstream.
4. **Hand off** the `<pair_name>/` dir (with `valid=true` manifest) to
   `twin-run-analysis`.

## Guard rails

- **Same GS, always.** Different ground states silently break every difference.
  `check_twin.py` asserts `gs_dir` parity.
- **Do not tune the two projectiles apart** beyond the representation itself
  (σ_WP, launch_z, k0 must match). To match Gaussian broadening (a σ study), change
  σ_WP in BOTH twins together, not one.
- **Full decomposition or abort.** If either twin is missing a required energy
  column, the pair is invalid — re-run with the flags on rather than patching in
  analysis.
- **Boundary/cadence** rules for the WP run still apply (4σ/1σ launch-stop; VTI
  cadence) per `feedback_jellium_boundary_rule`.

## Files (skill-local, shippable)

| File | Role |
|---|---|
| `check_twin.py` | parity + observable-completeness gate; writes `twin_manifest.json` |

Composes: `tddft-simulations` (build/dispatch/GPU/Gmail), `simulation-validation`
(pre-run tiers). Consumed by: `twin-run-analysis`.
Reference run.cpp: `…/localised_jellium_dynamics/proj_perturbation/run.cpp`.
