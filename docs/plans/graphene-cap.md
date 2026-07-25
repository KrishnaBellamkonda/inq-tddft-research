# Plan: Graphene + CAP wave-packet/classical scattering (feasibility replica)

Status: **DECISIONS LOCKED** via `/grill-with-docs` 2026-06-18.
Spec: `docs/prompts/absorbing_boundary/graphene_with_cap.md`.
Paper: Yao & Schleife, "Wave-packet electron dynamics ... graphene"
(`ResearchProject/literature/tddft-quantum-projectile/resources/wave-packet-electron-dynamics-on-graphene.pdf`).
Glossary: "Graphene CAP" section of `CONTEXT.md`.
Handover: `docs/handovers/graphene-cap.md`.

This is a **feasibility / methodology replica**, NOT the paper's converged
numbers. Every deviation from Yao & Schleife is listed in §"Deviations" and must
be carried into the notebook + emails as such.

---

## Locked decisions (grill 2026-06-18)

### Scope
- Run the **full pipeline autonomously overnight on the one free GPU (GPU 1;
  GPU 0 is busy with the user's other run)**, heavily instrumented so each stage
  produces observables + plots the user can trust.
- W-tuning sweep is **SKIPPED** — the user already has the data and locked the
  CAP parameters (below).

### System / cell
- **Graphene 4×4 supercell = 32 C atoms**, in-plane 18.6 × 18.6 Bohr
  (|a| = 2.46 Å = 4.651 Bohr; lattice **fixed at experimental, flat**, no relax).
- **z-cell = 60 Bohr** (z ∈ [−30, +30]); graphene sheet at z = 0.
- Cell ≈ 18.6 × 18.6 × 60 Bohr (tall slab). 128 valence e (~64 KS orbitals) + WP.
- **Functional: LDA / ALDA**; **norm-conserving ONCV carbon** pseudopotential;
  **plane-wave cutoff 50 Ha** ⇒ real-space dx ≈ 0.31 Bohr (grid ≈ 59×59×191).
- **Ehrenfest dynamics** — carbon ions move (forces every step).
- Γ-point only (large supercell).

### Absorbing boundary (CAP)
- **Two-sided sin² CAP**, total width **L = 20 Bohr** (10 Bohr each z-end,
  wrapping the periodic boundary at z = ±30). Free region = 40 Bohr (|z| < 20).
- **W = −0.5 Ha** (depth), pre-tuned by the user — no W-sweep this campaign.
- Requires the **complexified `inq-study` engine** (the in-built
  `perturbations::absorbing` only functions against `inq-study`; two-sided via
  `perturbations::sum` of two `absorbing` slabs, as in the vacuum CAP work).
- ⚠ **All CAP results PROVISIONAL until the inq-study engine regression
  (Task #7) passes** — standing governance caveat; do not block on it, but state
  it in every deliverable.

### Wave packet
- **Gaussian, d = 1.1 Å ⇒ σ_r = d/√2 = 1.47 Bohr** (paper-faithful; ψ ∝
  exp(−(r−b)²/2d²), d = √2·σ_r = 2.08 Bohr).
- **E = 100 eV** ⇒ k₀ = √(2E) = 2.711 a.u., along +z.
- **Launch position rule:** 5σ_r inside the near free-region edge, i.e.
  z₀ = −20 + 5σ_r = −20 + 7.35 = **−12.65 Bohr** (≈ 12.65 Bohr / 8.6σ from the
  sheet ⇒ orthonormal to graphene GS, satisfies guard-rail #1). Moves toward +z.

### Trajectories (both perpendicular, paper A/O — NOT the spec's grazing)
- **Standard = centroid (O):** WP/projectile aimed through a carbon **atom**.
- **Channeling (A):** aimed through a hexagon **hollow** site.
- (The spec's literal grazing-along-x was rejected in the grill: incompatible
  with periodic graphene; would need a finite ribbon/flake — deferred.)

### Classical ensemble
- **3 classical trajectories per trajectory-type** (spec said 5; paper >100).
- Classical projectile = **Gaussian-smeared pseudopotential, radial width σ =
  1.47 Bohr** (match the WP — apples-to-apples; overrides the spec's 0.5 Bohr).
  Generate `electron_gaussian_sigma1p47.upf`.
- Sampling: `std::mt19937` Gaussian draws in position **and** momentum, same
  mean/σ as the WP (paper §ensemble); average the ensemble for the classical
  answer.

### Observables (per run)
Baseline = inherited validated sets:
- **WP runs → coronene** required set (RT density {total,system,wp}, WP momentum
  distribution, WP↔GS overlap every step, LEED screens, GS orbital densities,
  initial WP density+WF) **∪ universal core** (energies total/kinetic/hartree/xc,
  density_l2 & Δdensity_l2, GS eigenvalues/occupations, GS density VTI,
  run_summary). Manifest via `minimum_observable_set.hpp`.
- **Classical runs → jellium-classical** set (`electron_track` every step,
  state energies, occupations, density VTIs {total,system}, overlap) ∪ core.

CAP/diffraction extras (ALL added, per grill):
- **Survival fraction ε(t)** = inner-region norm / N₀, + final "electrons
  remaining in vacuum after the WP reaches the CAP" (the W-criterion quantity).
- **Planar-integrated Δn(z,t)** = ∫∫ Δn dx dy (the paper's Fig. 1 trace).
- **Absorbed-fraction(t) + time-to-absorb** (inner-region norm < threshold).
- **Probability current density on the LEED screens** (j·n flux).

Whole-system fields (per grill clarification — NO meaningful "sum of orbitals"):
- **Total density n(r,t)** [real] full cadence (already in set).
- **Total current density j(r,t)** [vector field] — add an `inqkit` grid-current
  VTI writer if none exists (inq-stack only).
- **WP projectile orbital complex WF** full cadence (~1–2 GB/run total).
- **No per-orbital dump** (rejected on storage; would be the only route to
  band-projection analysis).

LEED screens: **8 total = 4 per side @ z = ±4, ±8, ±12, ±16 Bohr** (all inside
the free region; forward = transmission +z, backward = reflection −z), normal = z.

### Numerics
- **Propagator: ETRS** (mandatory with a CAP — Crank–Nicolson renormalises the
  WP each step and silently undoes absorption; established in the vacuum work).
- dt = 0.02 a.u.; n_steps ≈ 800–1000 (propagate until inner-region norm
  plateaus; travel −12.65 → +20 Bohr to far CAP ≈ 16 a.u. + absorption tail).
  Confirm at smoke. WRITE_EVERY ≈ 13 (~60 frames).
- GPU only (GPU 1). NVML/`nvidia-smi` is broken on this host — use the
  `cudaMemGetInfo` probe (`systems/vacuum/gpu_probe`), never nvidia-smi.

### GS validation (Step 1)
- **Standard:** SCF convergence + energy/atom vs an LDA reference + DOS /
  Γ-supercell eigenvalue spectrum showing semimetallic character (gap ≈ 0).
  (A true Dirac-cone band structure needs a separate primitive-cell k-path —
  NOT done this campaign.)

### Deliverables / notifications
- **Per-stage emails** to chiddukanna@gmail.com: (1) setup visualisation up
  front (non-blocking — autonomous, does NOT wait for approval), (2) GS
  validated, (3) no-CAP baseline, (4) standard/centroid ensemble+WP,
  (5) channeling ensemble+WP, (6) final summary. Subject family `[graphene-cap]`.
- **One master study notebook** at
  `ResearchProject/systems/graphene/hypotheses/<sweep>/` (ADR 0007 +
  notebook-making skill), a section per step, rebuilt as stages complete.

---

## Campaign run list (one GPU, ~10–14 h compute)
1. **GS** (1 run, ~min–1 h): graphene 4×4, LDA, 50 Ha, Γ. Validate.
2. **No-CAP baselines** (2 WP runs, one per trajectory) — comparison.
3. **Standard/centroid:** 3 classical + 1 WP (with CAP).
4. **Channeling:** 3 classical + 1 WP (with CAP).
= 1 GS + 2 baselines + 8 trajectory runs = **10 graphene TDDFT runs** + GS.

Per-run ≈ 1.0 h (4×4, 60 Bohr, Ehrenfest); campaign ≈ 10–14 h ± 2×.

---

## Folder layout (ADR 0007, new system `graphene/`)
```
ResearchProject/systems/graphene/
  shared_gs/                 converged graphene GS (reused by all runs)
  shared/                    config headers (Cfg structs, geometry)
  scripts/<sweep>/           run.cpp (build-once, inq-study), dispatch.py, analyse.py
  <sweep>/<run_name>/        per-run outputs (logs gitignored)
  hypotheses/<sweep>/        master notebook, combined CSVs, build_*.py, figs, tests/
shared/pseudopotentials/electron_gaussian_sigma1p47.upf   (generate)
```

## Deviations from Yao & Schleife (carry into notebook + emails)
| Quantity | Paper | This replica | Why |
|---|---|---|---|
| Supercell | 112 C | **24 C (3×2 rect, nx=3)** | one-GPU overnight + nx÷3 folds K→Γ (semimetal). nx=4/32-atom gave a spurious ~2 eV gap; corrected 2026-06-19. |
| In-plane WP images | large cell, minimal | x=13.95 Bohr → notable WP self-image overlap | nx locked to mult. of 3 by physics; nx=6 (ample) too slow. Documented finite-size artifact. |
| z-vacuum | 100 a₀ | 60 Bohr | cost (grid × step-count) |
| Ensemble | >100 | 3 / trajectory | overnight budget |
| Classical width | 0.5 Bohr (spec) | 1.47 Bohr | match WP for fair comparison |
| Trajectories | A/O perpendicular | A/O perpendicular | faithful (spec grazing rejected) |
| Lattice | relaxed | fixed experimental | save relaxation time |
| CAP engine | (paper code) | inq-study complexified | PROVISIONAL until Task #7 |

## Open implementation items (resolve during build, inq-stack/inq-study only)
- Verify carbon ONCV pseudo available (reuse coronene's); generate σ=1.47 UPF.
- Confirm/add grid current-density `j(r,t)` VTI writer in `inqkit` (none → add).
- Confirm LEED `plane_screen` supports the 8-screen layout + flux accumulation.
- Smoke one short run to fix n_steps / WRITE_EVERY and per-run wall time, then
  re-estimate the campaign before committing the full dispatch.
- Build run binary against **inq-study** (`INQ_SOURCE=…/inq-study inq-run`).
