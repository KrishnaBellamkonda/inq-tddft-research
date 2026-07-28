# Plan: validate the CAP energy-normalization hypothesis (vacuum WP suite)

Date: 2026-07-28. Companion to `docs/notes/inq-energy-normalization-error.md`.
Designed by a planning agent. All runs are cheap vacuum WP sims on **GPU 1**
(`CUDA_VISIBLE_DEVICES=1`; GPU 0 runs the jellium GS/replica).

## Hypothesis under test
The CAP "energy shoots up" entirely because INQ reports the per-particle
(norm-normalized) kinetic energy (`energy.hpp:50-55`, `occ_sum` divides by
`<psi|psi>`). The extensive energy `E_ext = E_reported*norm` decays smoothly with
the absorbed norm. The effect is independent of CAP geometry (one/two-sided),
strength eta, width W, and periodic wrap; it depends ONLY on the orbital norm
decaying under a non-Hermitian potential.

Falsifiers: (a) `E_ext` fails to track norm; (b) residual scales with eta/W (→
reflection); (c) a NORM-PRESERVING absorber still shows the rise; (d) rise appears
without any norm loss.

## Key source facts confirmed by the planner
- Vacuum non-interacting WP: `E_reported == energies.csv:kinetic == e_kin_ha`
  (per-particle mean); all other KS terms are 0.
- `wp_momentum_stats.csv` logs `e_kin_ha` (= the norm-divided kinetic) and a norm;
  `wp_real_space_stats.csv` logs the real-space norm. `E_ext = E_reported*norm` is
  computable per-step (energy_diagnostics.py already does this).
- **MaskAbsorber** (`inqkit/absorbers/mask_absorber.hpp`) multiplies psi by sin²
  each step. Under **ETRS** the amplitude removal survives → **norm-losing** mask.
  Under **Crank-Nicolson** (`crank_nicolson.hpp:139,147,162,165` orthonormalizes +
  renormalizes density to num_electrons each step) the removal is undone →
  **norm-preserving** absorber. Toggling the propagator is the decisive test.
- Existing env-driven scripts: `wp_traversal_energy` (one-sided CAP), `cap_sweep`
  (eta/W sweep + `inner_region_norm`), `twosided_cap_vs_mask`, `mfa_sweep` (ETRS
  mask). `wp_traversal_energy/run.cpp` hardcodes ETRS + `perturbations::absorbing`
  (no mask, no propagator knob).

## Global conventions
- Prefix every run `CUDA_VISIBLE_DEVICES=1`.
- Aliasing guard: `k_max=pi/h > k0 + 4*dk`, `dk=1/(sigma0*sqrt2)`. sigma0=3 → dk=0.236;
  h=0.4 → k_max=7.854 → k0 ≤ 6.9. Raise energy only with k0≤6.0 or drop h to 0.3.
- `WP_MOM_EVERY=1` for per-step norm. Ignore `energies.csv:wp_norm` (NaN); use stats CSVs.
- Prediction key for every CAP/mask run: `E_reported` stays pinned/rises, `norm`
  decays, `E_ext=E_reported*norm` decays ∝ norm.

## Execution order (dependency-sequenced)
1. **Phase 0** — baselines: 0a no-CAP control, 0b one-sided CAP (the phenomenon).
2. **Phase 6** — post-process 0b: confirm `energies.csv:kinetic == e_kin_ha` (not
   `e_kin_ha*norm`) → proves INQ prints the /norm quantity. No run.
3. **Phase 3 (DECISIVE)** — mask absorber, matched geometry:
   - 3a mask+ETRS (norm-losing) → should reproduce the CAP artifact.
   - 3b mask+CN (norm-preserving, norm≈1) → MUST show NO rise. If it rises, hypothesis dies.
   - 3c "norm loss without a CAP" complement (mask+ETRS / hard wall).
4. **Phase 2** — partial-absorption ladder: weak CAP so norm ends ~0.5/0.3/0.1;
   `E_ext/E0` vs norm must be the identity line y=x.
5. **Phase 1** — geometry independence: 1a eta∈{-0.3,-0.7,-1,-2,-3.5}, 1b W∈{10,15,20,25},
   1c two-sided. Residual must NOT scale with eta/W (else reflection).
6. **Phase 5** — spectral-width rise-rate: vary sigma0∈{2,3,5} and k0∈{2.711,5.421};
   fractional rise of `E_reported` scales with dk/k0, independent of eta/W.
7. **Phase 4** — numerics: h∈{0.5,0.4,0.3}, dt∈{0.02,0.01,0.005}; artifact must persist.

## Code changes required
- **For 3a/3b in matched geometry** (only change needed): add to
  `wp_traversal_energy/run.cpp` (a) `WP_ABS` switch (`cap`|`mask`) → when `mask`,
  build `inqkit::absorbers::MaskAbsorber(2, z_cap0, CAP_L, wp_idx)` and call
  `absb.apply(electrons)` in the callback for step>0 (mirror `mfa_sweep:180`),
  dropping the `perturbations::absorbing` arg to propagate; (b) `WP_PROP` switch
  (`etrs`|`cn`) appending `.crank_nicolson()`. ~10 lines each, patterns exist in
  `mfa_sweep/run.cpp`.
- Optional: add `e_kin_extensive=kinetic*norm` column + populate `wp_norm`. Not
  required (diagnostics reconstruct E_ext).

## New diagnostics (energy_diagnostics.py)
- `E_ext/E0` vs `norm` identity-line panel (Phase 2).
- Overlay `energies.csv:kinetic`, `e_kin_ha`, `e_kin_ha*norm` (Phase 6).
- `reflection_residual = inner_region_norm/N0` panel (Phase 1, separates reflection).

## The single decisive experiment
**3a vs 3b** (mask+ETRS vs mask+CN): identical spatial clipping, propagator toggles
only whether norm is preserved. Norm-losing → artifact; norm-preserving → no
artifact ⇒ isolates norm decay as the sole cause. Do this early (step 3).
