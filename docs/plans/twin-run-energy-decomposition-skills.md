# Design: twin-run energy-decomposition skills

**Status:** Rung 1 (static, known-answer) BUILT + VALIDATED 2026-07-14. Rung 2
(dynamic) specified, not built.
**Source of intent:** `docs/notes/energy-decomposition-skill.md`,
`docs/notes/gaussian-pertubation-for-classical-simul`.
**Handover:** `docs/handovers/localised-jellium-energy-book-keeping.md`.

## Goal

Two composable skills that turn a matched **classical(perturbation) + wavepacket**
localised-jellium pair into an interpreted, timestep-resolved energy decomposition
— isolating and physically explaining the **quantum effects** of treating the
projectile as a wavepacket rather than a classical point charge.

## Architecture — two skills + a shared contract

```
twin-run-generation ──produces──▶ [TWIN-RUN CONTRACT] ──consumed by──▶ twin-run-analysis
```

- **`twin-run-generation`** — make the pair correctly and *provably*: same GS,
  identical config, full energy decomposition switched on in both; gate with
  `check_twin.py` → `twin_manifest.json`. Composes `tddft-simulations` for the
  actual build/dispatch.
- **`twin-run-analysis`** — deterministic engine (`twin_decompose.py`) does the
  arithmetic + findings table; the agent writes the physical narrative from the
  interpretation rules in the SKILL.md. Skill-local Python (pandas+numpy).

### The twin-run contract (interface both skills bind to)

```
<pair_name>/
  wp/         raw/observables/observables.csv   run_summary.txt
  classical/  raw/observables/observables.csv   run_summary.txt
  twin_manifest.json          # valid=true is the gate
```

- **Required energy columns (both):** `energy_total/kinetic/hartree/xc/external`.
  Dynamic classical twin additionally: `energy_proj_bg_ideal` per step.
- **Parity fields (asserted equal):** `periodicity, Lz, spacing, N, sigma_WP,
  launch_z, gs_dir`. Only the projectile block may differ.
- **σ convention:** labelled by σ_WP; classical charge std `σ_pot = σ_WP/√2`.

## The physics (encoded in twin-run-analysis SKILL.md)

`E_total = E_kin + E_H + E_xc + E_ext (+ nonlocal + ion)`; the classical projectile
is an external potential → same stores as the WP run, different physics inside each.

- **`dKin`** = WP localisation kinetic `3/(4σ²)` (+ `k0²/2`). *(σ=0.5 → 81.6 eV)*
- **Gauge caveat:** individual `dHartree`/`dExt` are charged-cell-gauge-dependent;
  only `d(E_H+E_ext)` is physical.
- **`residual R = d(E_H+E_ext) − U_proj_bg`** = WP self-Hartree `E_H[WP–WP]`;
  r-independent at rest. *(σ=0.5 → 20.81 eV; free-space ref 21.71, open-z ~0.9 lower)*
- **`dXC`** = WP-alone xc, r-independent. *(σ=0.5 → −16.5 eV)*
- **`SIE = R + dXC`** = LDA one-electron self-interaction error (Perdew–Zunger) —
  the irreducible residue. *(σ=0.5 → 4.34 eV)*

## Complexity ladder

- **Rung 1 — static / known-answer (BUILT + VALIDATED).** At-rest pairs; constant
  decomposition; engine reproduces the golden numbers exactly. Regression:
  `tests/test_twin_decompose.py` (synthetic fixture + on-disk golden pair).
  Golden pair: `…/localised_jellium_dynamics/proj_perturbation/{results/proj_pert_dx0p5,
  stress_scratch/s0p5_r12_lz120_p2/results/wp}` → dKin 81.74, dXC −16.47,
  residual 20.81, SIE 4.34 eV.
- **Rung 2 — dynamic + representation-aware (FULL SPEC: `twin-run-rung2-dynamic-spec.md`).**
  Projectile moves under its own natural forces (no driving); both twins evolve from
  identical initial conditions → divergence *is* the quantum effect. Supports TWO
  classical representations: `perturbation` (Gaussian charge, inqkit `Projectile`
  Ehrenfest, clean residual 20.81 eV) and `pseudopotential` (ghost UPF, INQ-native
  ion Ehrenfest, residual 7.4 eV with the known ~14 eV ghost-aliasing gap). Adds an
  `E_proj_KE` store (→ classical-vs-quantum stopping power for free), per-step
  `U_proj_bg`, and WP centroid/σ(t) tracking. Validation without a known answer:
  energy conservation + t=0 collapse to the Rung-1 golden numbers. See the spec for
  the full build ladder and the four test twin pairs (P/G × static/dyn).

## Deliverables (all skill-local, shippable)

| File | Role | Status |
|---|---|---|
| `.claude/skills/twin-run-analysis/twin_decompose.py` | deterministic engine + CLI | done |
| `.claude/skills/twin-run-analysis/twin_notebook_builder.py` | executed analysis `.ipynb` | done |
| `.claude/skills/twin-run-analysis/tests/test_twin_decompose.py` | known-answer tests (7, pass) | done |
| `.claude/skills/twin-run-analysis/SKILL.md` | interpretation rules + workflow | done |
| `.claude/skills/twin-run-generation/check_twin.py` | parity + observable gate + manifest | done |
| `.claude/skills/twin-run-generation/SKILL.md` | contract + workflow + guard rails | done |

## Validation

- `twin_decompose.py`: 7/7 tests pass; CLI reproduces golden numbers; drift ≈ 0
  confirms the static case.
- `check_twin.py`: PASS on golden pair; FAIL correctly on a parity-broken /
  identical-projectile / missing-U_proj_bg pair.
- `twin_notebook_builder.py`: executes end-to-end on the golden pair, 0 errors.

## Open / next

- Rung 2: add per-step `U_proj_bg` to the classical run.cpp; centroid overlay in
  the engine/notebook; generate a first dynamic pair (GPU0) and narrate.
- Consider a wider-σ pair (negligible WP broadening) as a second static datapoint.
