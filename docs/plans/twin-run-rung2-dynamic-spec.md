# Spec: twin-run Rung 2 — dynamic + representation-aware

**Parent:** `docs/plans/twin-run-energy-decomposition-skills.md` (Rung 1 built/validated).
**Status:** SPEC (not built). **Depends on:** Rung-1 engine (`twin_decompose.py`),
`check_twin.py`. **No `inq/` edits** — everything is inqkit (wrapper) + run.cpp.
**Refs:** `docs/notes/energy-decomposition-skill.md`,
`docs/notes/gaussian-pertubation-for-classical-simul`,
`.claude/rules/light-projectile-stopping.md`, `.claude/rules/checkpoint-dont-block.md`,
memories `reference_ghost_upf_tail_aliasing`, `reference_inq_propagator_mask_absorber`,
`reference_charged_cell_hartree_convention`.

## Scope

Rung 2 adds two orthogonal things to the validated static engine:

1. **Dynamics** — the projectile moves under its own natural forces (no artificial
   driving); both twins evolve from identical initial conditions and the divergence
   *is* the quantum effect. Adds a projectile trajectory and a projectile-KE store.
2. **Representation-awareness** — the classical twin may be EITHER
   - **`perturbation`** — Gaussian charge → Poisson potential (no r_cut, no
     aliasing); the clean, fully-explained case (residual = WP self-Hartree = 20.81 eV).
   - **`pseudopotential`** — ghost UPF (z_valence=0, erf/r local); the older,
     *mostly*-understood case whose residual is 7.4 eV, with a documented ~14.3 eV
     gap = ghost-UPF tail aliasing (`reference_ghost_upf_tail_aliasing`). Included so
     the skill can analyse both and *explain the difference between them*.

## The two classical representations

| | `perturbation` (Gaussian charge) | `pseudopotential` (ghost UPF) |
|---|---|---|
| Enters KS as | external perturbation `v_proj=+poisson(n_proj)` | a ghost ion's local pseudopotential |
| Motion mechanism | **inqkit `Projectile` Ehrenfest** (custom integrator) | **INQ-native ion Ehrenfest** (ghost ion feels forces, moves) |
| Projectile KE store | `energy_proj_ke = ½·m·V²` (inqkit writes) | `energy_ion_kinetic` (INQ native column, already exists) |
| Projectile position | inqkit `Projectile.R(t)` | INQ ion position |
| `U_proj_bg` | ideal `∫n_proj·v_bg`, per step | **ideal term ONLY** — the impl term is r_cut/grid-pathological (`reference_ghost_upf_tail_aliasing`) |
| `U_proj_bg` sign in residual | `R = d(E_H+E_ext) − U_proj_bg` | `R = d(E_H+E_ext) + U_proj_bg` — ADD (INQ omits the z_valence=0 background-comp. term; `reference_ghost_u_proj_bg_sign`) |
| Aliasing | none | ghost tail aliasing → residual artefact |
| Residual at rest (VERIFIED via G-static) | 20.81 eV (= self-Hartree, clean) | **8.85 eV** (~12 eV below perturbation = ghost aliasing; SIE NOT clean) |
| Role | preferred / clean | reference / mostly-understood |

The energy stores and the residual FORMULA are identical across representations
(`d(E_H+E_ext) − U_proj_bg`); only the *expected value* and the *interpretation of
the remainder* differ. That is what makes one representation-aware engine sufficient.

## Dynamic model

### `perturbation` case — inqkit `Projectile` (Ehrenfest for a rigid Gaussian charge)

New `inqkit/dynamics/projectile.hpp`:
```
struct Projectile {
    double mass, charge;         // a.u.; twin-matched to the WP (m_e=1 unless muon fork)
    vec3   R, V;                 // position, velocity (state)
    void advance(F, dt);         // velocity-Verlet
};
```
- **Force (Hellmann–Feynman, self-force-free):**
  `F = −∇_R [ ∫ n_e·v_proj(r−R) dr  +  U_proj_bg(R) ]`.
  The projectile self-Hartree is R-independent for a rigid symmetric Gaussian → its
  gradient integrated against n_proj is zero by symmetry → **no spurious self-force**
  (a moving ghost UPF would not be so clean; this is a merit of the perturbation).
  Compute from the electron Hartree field + background field (both already built in
  `run.cpp`); a symmetric finite-difference / analytic gradient of the Gaussian.
- **Coupling (wrapper-only, mask-absorber pattern):** the moving perturbation reads
  `R(t)` from the `Projectile` at each substep (`potential(time,…)`); the step
  callback computes `F(t)` and calls `advance`. Order: perturbation uses R at
  substep start; callback updates after the step (standard Ehrenfest ordering).
  Precedent: `reference_inq_propagator_mask_absorber` (callback mutation). NO `inq/` edit.
- **Extend** `gaussian_projectile_perturbation` to take a `Projectile const&` (or a
  `center(time)` callback) instead of a fixed center.

### `pseudopotential` case — INQ-native ghost-ion Ehrenfest

- Add the ghost UPF as an ion with initial velocity `v0`; enable INQ ion (Ehrenfest)
  dynamics. INQ integrates the ion, reports `energy_ion_kinetic` and the ion
  position — **no custom integrator needed**.
- Still must: enable the full `ObservableSelection`; write `energy_proj_bg_ideal`
  per step via `set_proj_bg()` recomputed from the current ion position (the ideal
  term, never impl); pass the cutoff/aliasing guard for the ghost UPF (mandatory,
  `reference_cutoff_aliasing_guard`).

## Unified per-step contract (dynamic additions)

Both classical runs (either representation) add to `observables.csv`:

| column | meaning |
|---|---|
| `energy_proj_ke` | projectile KE ½mV² (perturbation: inqkit; pseudopotential: mirror of `energy_ion_kinetic`) |
| `energy_proj_bg_ideal` | per-step `U_proj_bg` (ideal), replaces the static run_summary constant |
| `proj_z`, `proj_vz` | projectile position / velocity along the beam axis |

The WP run adds:

| column | meaning |
|---|---|
| `wp_centroid_z` | WP centre-of-density along z (inqkit per-step diagnostic on the WP orbital density; or post-hoc via `inqview` `center_of_density`) |
| `wp_sigma_z` | WP spread σ_z(t) (detects broadening — a quantum effect) |

`check_twin.py` gains a `--dynamic` mode asserting these columns are present.

## The conserved-total ledger (resolves the E_kinetic asymmetry)

The WP's *motional* KE lives INSIDE `E_kinetic` (⟨p²⟩/2m of the WP orbital); the
classical's projectile KE is a SEPARATE store. So the like-for-like conserved
quantities are:

- **classical:** `E_conserved = E_total_INQ(electrons) + energy_proj_ke + U_proj_bg`
- **WP:** `E_conserved = E_total_INQ(electrons incl. WP)`

And the localisation comparison must be motional-matched:
`dKin_localisation = dKin − energy_proj_ke(classical)`  (at k0=0 this is just dKin,
recovering Rung 1). The engine computes and reports both.

**Stopping power — CLASSICAL projectile ONLY.** `−d(energy_proj_ke)/ds` over the
initial near-constant-velocity window = the drag on the **classical** projectile;
extract it with the `stopping-power-extraction` skill. Size the run by the
initial-drag window, NOT 5 plasma periods (`light-projectile-stopping.md`).

**Do NOT apply this to the WP (quantum) case.** The WP orbital is not identifiable
as "the projectile", so its kinetic energy is not projectile KE and `−d(E_proj_KE)/ds`
is invalid for the quantum run (`feedback_quantum_stopping_not_from_projectile_ke`).
The **total quantum stopping power** is why the localised-jellium system was chosen:
it is extracted from the **total electronic energy deposited** into the system over
distance, not from any projectile-KE track. The engine therefore computes an
`E_proj_KE` stopping curve for the classical twin only, and reports total-energy
deposition (ΔE_electronic vs projectile path) as the quantum-side stopping proxy.

## Engine extensions (`twin_decompose.py`)

- Read `representation ∈ {perturbation, pseudopotential}` from `twin_manifest.json`.
- Per-step `U_proj_bg` from the `energy_proj_bg_ideal` column (fallback to summary
  constant for static). Already stubbed in `_u_proj_bg_series`.
- Track `energy_proj_ke`; compute `E_conserved(t)` per run; report its drift
  (conservation gate).
- Centroid overlay: `proj_z` vs `wp_centroid_z` → separation Δz(t); residual(t) vs
  Δz(t) and vs `wp_sigma_z(t)`.
- **Representation-aware expected residual:**
  `perturbation` → free-space/open-z self-Hartree (unexplained = open-z gauge, ~0.9 eV);
  `pseudopotential` → self-Hartree − ghost-aliasing; flag the ~14 eV gap as the
  KNOWN ghost-UPF tail-aliasing artefact and cite the note (do NOT report it as
  missing physics).
- `motional-matched` localisation (subtract classical `energy_proj_ke`).

## Interpretation rules added to `twin-run-analysis/SKILL.md`

- **Representation matters.** State which classical twin. For `pseudopotential`,
  the residual is ~7.4 eV and the ~14 eV shortfall vs the 21.7 eV self-Hartree is
  the ghost-UPF tail aliasing — a numerical artefact, resolved by switching to
  `perturbation` (residual → 20.81). This *comparison between representations* is a
  first-class output.
- **Trajectory divergence is a quantum effect.** WP centroid vs classical `proj_z`
  diverging (WP reflects/tunnels/spreads differently) is signal, not noise.
- **Residual drift = WP spreading.** A falling residual with rising `wp_sigma_z`
  = self-Hartree dropping as the WP broadens (wider σ → smaller 1/(σ√2π)).
- **Narrate first few steps explicitly, then the trend** (butterfly cascade), per
  the note — combine residual(t), Δz(t), σ_z(t), and the deceleration curves.

## Validation gates (the no-known-answer regime)

1. **Energy conservation** — `E_conserved(t)` constant within tol in EACH run
   (drift → integrator/force bug). Primary correctness gate for the new dynamics.
2. **t=0 collapse to Rung-1 golden** — first step must reproduce the static numbers
   *per representation*: `perturbation` → 20.81/−16.47/4.34; `pseudopotential` → 7.4
   residual (with the aliasing note).
3. **Force sanity** — `F(t=0)` sign: projectile pulled toward the slab; zero-field →
   constant V; constant-field unit test → uniform acceleration a=F/m (analytic).
4. **Cutoff/aliasing guard** — mandatory pre-launch for the ghost UPF.
5. **Light-projectile sizing** — initial-drag window, gate on a clean drag slope
   existing, never abort on v-drift.
6. **Checkpointing** — interior RT checkpoints every ~200 steps; WARN-not-block on
   projected overrun (`checkpoint-dont-block.md`).

## Smoke / unit tests (skill-local + inqkit two-tier)

- **`Projectile` unit (inqkit test):** constant field → `a=F/m` to machine tol;
  zero field → constant V; symmetric field → self-force = 0; Verlet conserves energy
  on a harmonic well.
- **Ghost-ion parity:** INQ `energy_ion_kinetic` == ½mV² from the reported ion V.
- **Engine dynamic tests (extend `test_twin_decompose.py`):** synthetic dynamic
  fixture with a known linear trajectory → correct `proj_z`, Δz, `E_conserved`
  constant; representation flag routes the expected residual correctly.
- **t=0 twin collapse** for each representation against the golden numbers.

## Test twin pairs to generate

| Pair | classical | purpose | data |
|---|---|---|---|
| **P-static** | perturbation | Rung-1 anchor | EXISTS (golden pair) |
| **G-static** | pseudopotential | reproduce ~8.85 eV residual (+U_proj_bg) + document ~12 eV aliasing | **NEW** — no full-ledger ghost run on disk; cheap static run, GPU0 |
| **P-dyn** | perturbation (moving `Projectile`) | first dynamic quantum-effect narrative | NEW — short (initial-drag window), GPU0 |
| **G-dyn** | pseudopotential (ghost-ion Ehrenfest) | dynamic ghost reference | NEW — short, GPU0 |

Start with **G-static** (validates representation-awareness against a second
known answer, no new dynamics code), then **P-dyn** (validates the `Projectile`
integrator against energy conservation + t=0 collapse), then **G-dyn**.

## Build order (complexity ladder)

1. Engine + `check_twin.py`: add `representation` + per-step `U_proj_bg`/`E_proj_ke`
   handling + dynamic synthetic tests. (No GPU.)
2. **G-static** run (ghost UPF + full ledger) → engine reproduces 7.4 eV + aliasing
   note. Second known-answer anchor. (Cheap GPU.)
3. `inqkit/dynamics/projectile.hpp` + moving perturbation + `Projectile` unit tests.
4. **P-dyn** run → conservation + t=0 gates → first dynamic narrative.
5. Ghost-ion Ehrenfest run.cpp → **G-dyn**.
6. Update both SKILL.md (representation + dynamics rules); dynamic notebook section
   (trajectory + deceleration + residual(t) overlays).

## Files to add/change (all wrapper/skill — NO `inq/`)

| File | Change |
|---|---|
| `inqkit/dynamics/projectile.hpp` | NEW — `Projectile` + Verlet + HF force |
| `inqkit/jellium/gaussian_projectile_perturbation.hpp` | time-dependent center from `Projectile` |
| `inqkit/io/observables_writer.hpp` | `energy_proj_ke, proj_z, proj_vz, wp_centroid_z, wp_sigma_z` |
| `…/proj_perturbation/run.cpp` (+ a ghost variant) | dynamic loop: force→advance; per-step ideal U_proj_bg; k0/v0; checkpoints |
| `.claude/skills/twin-run-analysis/twin_decompose.py` | representation-aware, dynamic |
| `.claude/skills/twin-run-generation/check_twin.py` | `--dynamic` column checks + `representation` in manifest |
| both `SKILL.md` | representation + dynamics interpretation rules |
| `inq-stack/tests/include/inqkit/dynamics/` | `Projectile` Catch2 unit test |
