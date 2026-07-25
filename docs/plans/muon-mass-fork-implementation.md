# Plan: per-state mass fork (inq-study) + validation suite

Status: IN PROGRESS (started 2026-07-06). Companion:
`docs/campaigns/muon_projectile/inq_study_engine_notes.md` (design + call chain),
`docs/notes/muon_plan_understanding.md` (user's understanding + feedback).

Goal: give any Kohn–Sham orbital a tunable mass via a per-state inverse-mass array,
so INQ can simulate (a) a muon projectile in electron jellium, (b) an all-muon
jellium, (c) mass-tuned band structure. `inq/` is immutable; all edits in `inq-study`.

## Part 1 — Implementation (convert the PROPOSAL annotations to code)

Design = Strategy A (per-state inverse mass on `electrons`, read by `ks_hamiltonian`,
applied in `operations/laplacian`). Alignment fact that makes it safe: **per-state
mass only matters in RT on the full state set (muon-in-electron-jellium), where the
Hamiltonian's factor array aligns with `phi`'s partition. In GS the mass is uniform
(all-1 electron, or global muon), so any block/subset indexing is harmless.**

| Step | File | Concrete change |
|---|---|---|
| 1 | `systems/electrons.hpp` | member `gpu::array<double,2> inverse_mass_`; accessors; `reextent`+`fill(1.0)` at both alloc sites. |
| 2 | `operations/laplacian.hpp` | add `_states` variants (distinct names — NO overload ambiguity): `laplacian_states`, `laplacian_add_states`, `laplacian_expectation_value_states`, taking `gpu::array<double,1> const & factors`; kernel `factors[ist]*(…)`. |
| 3 | `hamiltonian/ks_hamiltonian.hpp` | member `gpu::array<double,1> kinetic_factor_`; ctor gains defaulted `inverse_mass` arg → build `kinetic_factor_ = -0.5*inverse_mass` (or all -0.5); call the `_states` variants at :202/:235/:245. |
| 4 | `real_time/propagate.hpp` | pass `electrons.inverse_mass()[0]` into ctor (RT opt-in). |
| 5 | `ground_state/calculator.hpp` | pass `electrons.inverse_mass()[0]` into ctor (GS opt-in; all-muon/bands). |
| — | `ground_state/initial_guess.hpp` | NO CHANGE (overlap-only Hamiltonian). |
| 6 | consumer `run.cpp` | set `inverse_mass()`: WP slot only (projectile) or all entries (muon jellium). |

Design decision: use **distinct `_states` function names** rather than overloads, so
every existing caller of scalar `laplacian*` is untouched (protects the Tier-0
bit-for-bit invariant).

## Part 2 — Validation suite (expected values fixed UP FRONT)

Analytic facts used as oracles (atomic units, ħ=1):
- Plane-wave kinetic eigenvalue: `T e^{ik·r} = (k²/2m) e^{ik·r}`.
- Free Gaussian spreading (density std): `σ_ρ(t)² = σ_ρ(0)² + (t/(2 m σ_ρ(0)))²`
  → spreading rate ∝ **1/m**.  (σ_WP = √2·σ_ρ; sigma-wp-convention.)
- Free WP group velocity: `v = k₀/m` → centroid slope ∝ **1/m**.
- Particle in a box: `E_n = n²π²/(2 m L²)` ∝ **1/m**.
- One-component HEG rescaling: at fixed r_s, `E ∝ m` (eff. units), `S ∝ m²` (real).

### Tier 0 — INERT WHEN OFF (the foundation; cheap)
| T0.1 | With `inverse_mass ≡ 1`, the existing electron-jellium **GS energy is bit-for-bit** identical to a pre-fork build. | MUST |
| T0.2 | inq-study `ctest` suite passes unchanged. | MUST |
| T0.3 | An existing electron RT run: energy trace identical to pre-fork. | MUST |

### Tier 1 — KERNEL UNIT TESTS (Catch2, cheap; `_engine` where a live set is needed)
| T1.1 | Plane wave at known k: `laplacian_states` with factor −0.5 → eigenvalue `k²/2`; with −0.5/207 → `k²/414`. Ratio = 1/207. | analytic |
| T1.2 | `laplacian_expectation_value_states`: ⟨T⟩ per state = `k²/2m`; muon slot and electron slots in the SAME set give the right per-state values. | analytic |
| T1.3 | **Electrons-unchanged:** electron slots' ⟨T⟩ identical whether or not a muon slot is present (no leakage). | analytic |
| T1.4 | **Wrong-slot guard:** only the flagged slot has non-unit mass; all others report m=1. | analytic |
| T1.5 | **Ledger consistency:** the factor used by apply (:202/:235) equals the one used by expectation (:245) — assert same array object. | structural |

### Tier 2 — VACUUM FREE PROPAGATION (user's σ(t) test + analytic overlays; GPU)
| T2.1 | Stationary Gaussian (k₀=0), plot σ_ρ(t) for m∈{1,207}; overlay analytic `σ_ρ(t)`. Muon spreads ~207× slower. | <5% vs analytic |
| T2.2 | **Group velocity:** moving WP (k₀≠0), centroid slope = k₀/m (electron vs muon). | <2% |
| T2.3 | **KE conservation:** ⟨T⟩ = k₀²/2m constant, drift <0.1%. | analytic |
| T2.4 | **Norm conservation:** ‖ψ‖ drift <1e-4 (Hermiticity of the scaled operator). | unitarity |
| T2.5 | **Time reversal:** propagate +N steps then −N; recover initial WP (‖Δ‖ small). | reversibility |

### Tier 3 — MUON vs ELECTRON JELLIUM (user's density comparison; GPU)
| T3.1 | **Localised-jellium slab GS density profile**, electron vs muon bath: heavier mass → LESS edge spill-out (lower kinetic pressure). *(Note: a UNIFORM jellium GS density is mass-independent — the SLAB's spill-out is the discriminating observable.)* | qualitative + monotone |
| T3.2 | **Rescaling gate:** muon jellium at fixed dimensionless r_s reproduces an electron run's energetics scaled by the known power of m (`E ∝ m` eff. units). | quantitative |

### Tier 4 — ADDITIONAL ROBUSTNESS CHECKS (suggested — see chat)
| T4.1 | **Mass-dial continuity:** sweep m∈{1,2,5,10,50,207}; spreading rate vs 1/m is linear through origin (R²≈1). | linearity |
| T4.2 | **GPU vs CPU agreement** (6+ sig figs) for a mixed-mass case — catches per-state indexing bugs. | 1e-6 |
| T4.3 | **MPI partition correctness:** 1 rank vs N state-parallel ranks give identical results; muon (last state) keeps its mass whichever rank owns it. | 1e-9 |
| T4.4 | **Particle-in-a-box** (mass m in a finite well): `E_n ∝ 1/m` — analytic quantum check WITH a potential. | <1% |
| T4.5 | **RT energy drift, mixed masses:** muon WP in electron jellium, total-E drift <0.1% over the run (end-to-end ledger test). | 0.1% |
| T4.6 | **Free-particle dispersion E(k):** empty-lattice eigenvalues scale as `k²/2m` (band-structure tie-in). | <1% |

Tier 0/1/(4.2–4.4) are cheap and run first. Tier 2/3/(4.1,4.5,4.6) are GPU sims →
`simulation-validation` menu, user-approved before launch.

## Part 2b — Physics-sim configs (draft for approval; NOT launched)

Both are GPU runs (validation-gates: GPU default). Both need the GPU build first.

### SIM 1 — vacuum σ(t) spreading  (covers T2.1 σ(t), T2.2 v_group, T2.3 KE cons.)
- **System:** empty cubic box, NO jellium background, NO occupied electrons — a
  single injected Gaussian WP as the only orbital. Free Ehrenfest.
- **Method:** inject WP (σ_WP, centre = box centre, k₀); set `inverse_mass` for the
  WP slot; `real_time::propagate` (ETRS — no exact exchange); record
  `wp_real_space_stats` (σ_ρ(t)), `wp_momentum_stats`, centroid, ⟨T⟩ each step.
- **Two masses, same everything else:** m = 1 (electron) and m = 206.77 (muon).
- **Parameters (proposed):** σ_WP = 0.5 Bohr (σ_ρ = σ_WP/√2 ≈ 0.354);
  L = 48 Bohr; spacing 0.4 Bohr; dt = 0.02 a.u.
  - **Panel A (spreading):** k₀ = 0, run to t ≈ 12 a.u. Electron σ_ρ grows to
    ~8 Bohr (well inside L/2 = 24); muon stays ≈ 0.354 (rigid). Overlay analytic
    `σ_ρ(t)² = σ_ρ(0)² + (t/2mσ_ρ0)²`. Fit slope of σ_ρ² vs t² → extract m; **check
    m to <5%** and the muon/electron spreading-rate ratio = 206.77.
  - **Panel B (muon spreads):** muon-only, k₀=0, run to t ≈ 120 a.u. so the muon
    reaches σ_ρ ~ 0.6; confirm it follows the SAME law at its own rate.
  - **Panel C (group velocity, T2.2):** k₀ = 0.5 Bohr⁻¹ along z; centroid slope =
    k₀/m (electron 0.5, muon 0.0024 Bohr/a.u.); **check <2%**.
- **Checks/adheres to:** analytic dispersion law (family 1) + KE & norm
  conservation (family 2). Box-size guard: WP tail must stay < boundary (reuse the
  4σ rule); abort on NaN.

### SIM 2 — muon vs electron jellium density  (covers T3.1 spill-out, T3.2 rescale)
- **System:** the existing localised-jellium slab (`ResearchProject/systems/
  localised_jellium`), GS only. Two runs: global mass m = 1 vs m = 206.77 (ALL bath
  particles muon — set every `inverse_mass` entry). GS path now carries mass
  (calculator.hpp).
- **Method:** converge GS for each mass; extract the density z-profile (VTI →
  inqview z-profile); compare edge spill-out. Also log E_GS, E_kinetic.
- **Check (T3.1):** heavier bath → LOWER kinetic pressure → density hugs the
  background more tightly → **less edge spill-out** (monotone, qualitative).
  *(A uniform jellium density is mass-independent — the slab edge is the
  discriminating observable, per the plan's T3.1 note.)*
- **Check (T3.2 rescaling gate):** hold the DIMENSIONLESS r_s fixed (rescale the
  box/spacing by 1/m) and confirm E_GS(muon) = m · E_GS(electron) in effective
  units — the HEG scaling law → a quantitative correctness gate on the GS-mass path.
- **Caveat to surface:** at the SAME physical density the muon slab is strongly
  correlated (r_s × 207) → LDA is stretched; the spill-out trend is still
  qualitatively valid but absolute energetics are approximate (Q4 scenario 2).

Pre-launch: `simulation-validation` tier menu + user approval; boundary/NaN
guards; provenance in the handover (run IDs, σ_ρ-fit m, E_GS values).

## Part 3 — Provenance
- catalogue rows in `docs/validation/test-catalogue.md`;
- formula-validation subagent for the spreading-law + plane-wave oracles;
- handover `docs/handovers/muon-mass-fork.md`.
