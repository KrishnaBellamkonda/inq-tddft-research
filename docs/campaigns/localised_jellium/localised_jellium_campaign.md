---
id: locjel-campaign
area: localised_jellium
title: "Localised jellium - implementation + scattering"
status: done
hypothesis: "A static localised-background perturbation realises a confined jellium slab whose projectile stopping S can be measured (bare and CAP-bounded)."
handover: docs/handovers/localised-jellium.md
tasks:
  - { name: "Phase 1 - implement localised background + perturbation", done: true }
  - { name: "Phase 2 - validate slab + static run", done: true }
  - { name: "Phase 3 - WP vs bare slab", done: true }
  - { name: "Phase 5 - two-sided CAP stopping (classical + WP)", done: true }
blocked_reason: ""
---

# Localised jellium — implementation + scattering campaign

<identity>
You are a scientific computing researcher working on first-principles TDDFT
simulations with INQ. You write scientific-standard code, adhere to every rule,
skill, and workflow in this repository, and never claim correctness without
recorded validation evidence. `inq/` is immutable; engine changes live only in
`inq-study/`; new library code lives in `inqkit`/`inqview`.
</identity>

<goal>
Implement a **localised jellium** target in INQ (a finite positive background
confined to a region of the cell), validate it rigorously, then fire a projectile
through a jellium **slab** — first bare, then with absorbing boundaries — to
measure stopping power. Theory backing: `docs/notes/localised-jellium-theory.md`.
Glossary: `CONTEXT.md` → "Localised jellium". Decisions below were locked in a
grill-with-docs session (2026-06-21).
</goal>

<locked_decisions>
1. **Mechanism (ADR-0008 candidate).** The localised background is a **static
   custom perturbation** added via INQ's `Perturbation` hook — NOT smeared ions,
   NOT an `inq/` edit. A new `inqkit` class implements the perturbation duck-type;
   its `.potential(t, v)` adds `v_bg(r) = −poisson(n₊)` every SCF iter and every RT
   step (cached; computed once). Present in BOTH `ground_state::calculate` and
   `real_time::propagate` (same object), so the well confines electrons in the GS
   and persists during the projectile flight.
2. **Engine: `inq-study` for ALL phases** (the sin² CAP needs the inq-study
   complexified scalar potential; uniform engine ⇒ one slab GS reused across
   no-CAP and CAP runs with zero portability risk). Verify `inq-study ≡ inq` on
   non-CAP files (`diff`) for attribution.
3. **CAP: built-in sin² CAP** (`perturbations::absorbing`), NOT the monomial, NOT
   the mask absorber. Two-sided ⇒ `perturbations::sum(background, CAP_−z, CAP_+z)`.
   `eta = −0.5 Ha`, total length 15 bohr (7.5 each side). Note: sin² is **zero at
   the box wall** (peaks mid-window) — a property to validate, not a bug.
4. **Density: r_s ≈ 4 (Na-like).** Slab (full 50×50 face, 25 bohr thick,
   V=62 500 bohr³): **N = 234** (even), n₀ = N/V = 3.744×10⁻³ a₀⁻³ → exact
   neutrality, effective r_s = 3.996. Spherical version is DEFERRED (todo).
5. **Box: 50 bohr cubic periodic** (INQ-centred, z ∈ [−25, +25]).
6. **Projectile: σ = 0.5 bohr, E = 100 eV** (k₀ = v = 2.71 a.u.). Slab transit
   ≈ 9 au. WP launched 4σ = 2 bohr from a boundary (run-up ≈ 10.5 bohr to the
   slab face at z = −12.5; exit ≈ 12.5 bohr). Classical twin uses the existing
   `electron_gaussian_sigma0p35.upf` (σ_pot = σ_WP/√2 = 0.354 → matches WP density
   std; `reference_sigma_matching_convention`).
7. **Plane screens: 20, feature-aligned** evenly-spaced xy-planes (normal = z),
   reusing the jellium `leed_screen_layout` shape, with screens nudged onto the
   slab faces (±12.5), CAP inner edges (±17.5), launch and exit planes.
8. **Stopping power.** `S = ΔE_bath / x`, `x = 25 bohr` (slab thickness).
   `ΔE_bath = E_bath(T) − E_bath(GS)` measured at time T AFTER the projectile is
   fully gone and CAP drain has died (energy ill-defined while CAP absorbs).
   Computed for BOTH classical and WP runs (two independent S). Classical run also
   yields `ΔKE_ion/x` as a bonus cross-check. **Gated on T3.2** (CAP must not drain
   the bath).
</locked_decisions>

<phases>
**Phase 1 — Implement.** `inqkit` headers: `jellium/localised_background.hpp`
(builds n₊ field: shape ∈ {slab, sphere, box}, r₀, R_cl/half-width, edge width) and
`jellium/background_perturbation.hpp` (the static perturbation; lazy one-time
`v_bg = −poisson(n₊)`, cached; adds real v_bg into the (possibly complex) potential
field). Add `E_self` helper to `jellium/analytics.hpp`. Wire into a localised-slab
`run.cpp` (GS + RT) passing the perturbation to both calls; `extra_electrons(234)`.

**Phase 2 — Validate (slab) + static run.** Run T0+T1+T2+T3 (menu below). Visualise
GS KS density (user checks). Then a **2 au static run** (slab only, no projectile)
and visualise density vs time (must be stationary). 50 bohr box. All key results in
the notebook.

**Phase 3 — Projectile vs bare slab (WP).** WP (σ=0.5, 100 eV) launched outside the
slab, propagates +z through it. No CAP. 20 screens. **Max out observables** +
heavy post-processing. Important plots → notebook. xz density gif.

**Phase 5 — Projectile + slab + two-sided sin² CAP (classical AND WP).** CAP total
15 bohr (7.5/side), eta=−0.5. WP/classical launched 4σ from CAP inner edge, ~3–4
bohr to slab, exits ~5 bohr to far CAP. Measure: per-side cumulative absorbed norm
vs t; bath-energy trace; ΔE_bath → S (both runs); classical ΔKE_ion/x cross-check.
Full observable suite + post-processing. xz density gifs. Direct classical-vs-WP
comparison.
</phases>

<validation_menu>
ALL FOUR TIERS approved (2026-06-21).
- **T0 (host, no GPU):** ∫n₊=N; interior n₊=n₀ + edge profile; v_bg vs analytic slab
  potential (parabolic inside, linear outside). Catalogue rows required.
- **T1 (cheap GS):** SCF converges, density peaks inside slab, e–bg energy < 0;
  interior density flat to a few % (the "increase R_cl" gate); **visualise GS KS
  density**.
- **T2 (Lang–Kohn surface physics):** surface profile (spill-out + Friedel π/k_F);
  work function Φ = v_vac − μ vs Lang–Kohn(r_s=4); surface energy σ vs
  Lang–Kohn(r_s=4); grid convergence (spacing ×½). NOTE: the 86.4 erg/cm² figure is
  the GPAW tutorial's r_s, NOT r_s=4 — pin the r_s=4 Lang–Kohn Φ and σ via a
  `docs/sources/` note before using as gates.
- **T3 (dynamics/CAP gates):** GS density at CAP onset (±17.5) < 0.1% n₀; CAP-only
  drain run (bath norm loss ≈ 0 / 2 au); classical-vs-WP S agreement; 2 au no-CAP
  static run total-energy conservation.
</validation_menu>

<observables>
"Max out": every primary observable in the minimum-observable-set (universal core +
jellium-WP and jellium-classical required sets) PLUS the localised-jellium extras:
n₊ field, v_bg field, E_self, bath energy trace, per-side CAP cumulative absorbed
norm. All standard inqview post-processing phases run in each run's `analyse.py`
(per `feedback_per_run_analyse_py`): density carpets, momentum distribution, COD,
spectra, LEED screens + IFFT, stopping power, wake decomposition, etc. Plus the
interior-density / surface-profile / work-function / S analyses specific to this
campaign.
</observables>

<deliverables>
- ipynb per phase (ADR-0007: `systems/localised_jellium/hypotheses/<sweep>/`),
  house narrative (`notebook-making` skill): context → formulas (every term
  defined) → full reconstructable setup → linked source files → results →
  takeaway. INCLUDE the implementation sketch (how the perturbation enters the INQ
  workflow) so the reader understands what was changed.
- **xz density gif** in every run's notebook.
- Email results at the end of each phase (`chiddukanna@gmail.com`).
</deliverables>

<deferred_todos>
- **Spherical localised jellium**: same mechanism, sphere shape, magic-N ladder
  {40, 92, 138, 198} at r_s=4 for shell-structure + E/N→HEG-limit (VC-4a/4b).
- **S(v) sweep** in Phase 5 (multiple projectile energies) once single-point works.
</deferred_todos>

<placement>
New system `ResearchProject/systems/localised_jellium/` on the canonical ADR-0007
layout (shared/, shared_gs/, scripts/<sweep>/, <sweep>/<run>/, hypotheses/<sweep>/).
Sweeps: `01_slab_validation`, `02_projectile_slab`, `03_cap_stopping`.
</placement>
