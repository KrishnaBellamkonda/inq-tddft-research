# Plan — Localised jellium implementation + scattering campaign

Authoritative companion to `docs/prompts/localised_jellium/localised_jellium_campaign.md`.
Theory: `docs/notes/localised-jellium-theory.md`. Glossary: `CONTEXT.md` →
"Localised jellium". Locked in grill-with-docs 2026-06-21.

## Status (2026-06-21)
- Design: **LOCKED** (8 decisions, 4 validation tiers — see prompt).
- Code: **NOT STARTED**.
- GPU: **both cards busy** with user's own `run`/`orted` jobs (GPU0 PID 2069269,
  GPU1 PID 2071679). GS launch is gated on a free card.

## Engine / immutability
- All runs on `inq-study` (sin² CAP needs its complexified scalar potential).
- **No `inq/` edits. No `inq-study/` edits** — the background is a wrapper-only
  `inqkit` perturbation; the sin² CAP already exists in inq-study.
- Verify `diff -q` inq vs inq-study on non-CAP files for attribution.

## Phase 1 — Implementation (NO GPU)
Files (new unless noted):
1. `inq-stack/include/inqkit/jellium/localised_background.hpp` — build `n₊` field on
   a given basis. Params: shape{slab,sphere,box}, centre r₀, R_cl/half-width, edge
   width w, n₀. Fill via `gpu::run` + `point_op.rvector_cartesian`. Slab: |z−z₀| <
   half-width; sphere: |r−r₀| < R_cl. Edge: sharp Θ now, erf/Fermi option later.
2. `inq-stack/include/inqkit/jellium/background_perturbation.hpp` — perturbation
   duck-type (copy `perturbations::none` skeleton + `has_potential()→true`).
   `mutable std::optional<field> v_bg_`. `.potential(t,v)`: if unset, build n₊ on
   `v.basis()`, `v_bg_ = −poisson::solve(n₊)`, cache; then `increment(v, v_bg_)`.
   Must compile against real (inq) AND complex (inq-study) potential field — add a
   real field into the (complex) potential's real part.
3. `inq-stack/include/inqkit/jellium/analytics.hpp` — add `e_self_sphere = 0.6 N²/R`
   and `e_self_slab` (per-area) helpers.
4. `ResearchProject/systems/localised_jellium/` skeleton (ADR-0007). Slab GS +
   RT `run.cpp` built on the jellium `run_template` (pass `pert` to GS + RT;
   `extra_electrons(234)`; write n₊, v_bg, E_self to run_summary/observables).
   Config header carries r_s, shape, half-width(=12.5), edge w, box, projectile.

## Phase 1 tests — T0 (host-only, NO GPU)
`inq-stack/tests/include/inqkit/jellium/test_localised_background.cpp` (Catch2,
pure host): ∫n₊=N (within quadrature tol); interior n₊=n₀; v_bg vs analytic slab
potential (parabolic inside, linear outside, value match). Add rows to
`docs/validation/test-catalogue.md`. **These run before any GPU work.**

## Phase 2 — Slab validation + static run (GPU)
- GS of the 234-e slab → `shared_gs/`. T1: SCF convergence, density peaks in slab,
  e–bg < 0, interior flat to few % (gate), **GS KS density viz**.
- T2: surface profile (Friedel π/k_F), Φ, σ, grid ×½ convergence. Pin r_s=4
  Lang–Kohn Φ/σ via `docs/sources/lang-kohn-1970.md` first.
- T3.4: 2 au static run (no projectile, no CAP) → total energy conserved; density
  stationary in time (viz). Notebook `hypotheses/01_slab_validation/`.

## Phase 3 — WP vs bare slab (GPU)
- Reuse slab GS. WP inject (σ=0.5, 100 eV) at z=−23; propagate +z. No CAP.
- ~20 au / dt=0.02 (~1000 steps). 20 screens. Full observable suite + analyse.py
  all phases. xz gif. Notebook `hypotheses/02_projectile_slab/`.

## Phase 5 — Projectile + slab + two-sided sin² CAP (GPU, two runs)
- `sum(background, CAP_−z, CAP_+z)`; CAP each side: mid_pos = ±21.25/50 fractional,
  width = 7.5/50 fractional, amplitude = −0.5 Ha.
- T3.1 (GS density at ±17.5 < 0.1% n₀) + T3.2 (CAP-only drain ≈ 0) BEFORE the
  stopping runs.
- Run A: WP projectile. Run B: classical projectile (`electron_gaussian_sigma0p35`
  ion, v=2.71). Same geometry, launched concurrently when 2 GPUs free.
- Measure per-side cumulative absorbed norm, bath-energy trace; ΔE_bath/x (both),
  ΔKE_ion/x (classical). T3.3: the two S agree. Notebook
  `hypotheses/03_cap_stopping/` with direct classical-vs-WP comparison + gifs.

## Execution order (autonomous)
1. [now, no GPU] Phase 1 headers → T0 host tests pass → build slab run.cpp (CPU
   compile check).
2. [gated on 1 free GPU] Phase 2 GS + T1/T2/T3.4 + static run.
3. [gated] Phase 3 WP run.
4. [gated on 2 free GPUs] Phase 5 classical + WP concurrent.
Each phase: analyse.py → notebook → email → handover update.

## Open / to-pin
- r_s=4 Lang–Kohn Φ and σ values (literature-review; 86.4 erg/cm² is a different
  r_s).
- ADR-0008 (perturbation mechanism) to be written.
- Edge profile (sharp Θ vs erf/Fermi) — start sharp, soften if Gibbs ringing
  appears in v_bg (worksheet Part 8).
