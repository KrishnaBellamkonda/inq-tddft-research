# E_proj_bg — dual-route validation (energy book-keeping campaign, B1)

Date: 2026-07-11. Status: **PASSED** (advisor-approved protocol; formula-validation
subagent waived for this post-hoc analysis scalar — advisor ruling 2, logged in
`docs/handovers/localised-jellium-energy-book-keeping.md`).

## Quantity

E_proj_bg(z_p) = ∫ ρ_ghost(r) φ₊(r) d³r — the projectile↔positive-background Coulomb
term deliberately omitted by the classical insertion runs
(`campaign_autorun/classical/run.cpp`, `ghost_background_term_omitted = true`).
ρ_ghost = −1·Gaussian(σ_pot = σ_WP/√2 = 0.354, entered independently in each route —
the √2 trap is the known failure mode). φ₊ = potential of the slab background
(n₀ = 1.312e-3, half-width 12.5) in the stated convention.

## Convention (advisor condition iii)

Primary: **3D-periodic mean-zero** (INQ drops G=0; mirrors
`inqkit::jellium::background_perturbation.hpp` → `solvers::poisson::solve(n₊)`).
For p2 (2D-periodic + open z) the per-G kernel solver in
`hypotheses/campaign_autorun_study/b1_decomposition.py::poisson_p2` is used
(G_xy ≠ 0: (2π/G)e^{−G|z−z′|}; G_xy = 0: open sheet −2π|z−z′|, gauge documented).
Only same-convention combinations are differenced (charged-cell Hartree caveat).

## Routes (independent implementations, no shared φ code)

- Route A: closed-form piecewise-parabolic periodic mean-zero φ₊ + 1D Gaussian
  marginal quadrature (pure formula, no grids).
- Route B: 3D FFT Poisson on a grid (48×48×1200; G=0 dropped), 3D Gaussian charge
  normalised on-grid.

## Results (Lz = 120, the h0 ledger box)

| r (Bohr from face) | A (eV) | B (eV) | \|A−B\| |
|---:|---:|---:|---:|
| 0 (centre) | −79.52 | −79.52 | 0.00 |
| 4 | −34.79 | −34.59 | 0.20 |
| 12 | −5.26 | −5.09 | 0.17 |
| 20 | +18.30 | +18.43 | 0.13 |
| 28 | +35.87 | +35.96 | 0.09 |
| 36 | +47.46 | +47.51 | 0.05 |
| 40 | +51.01 | +51.05 | 0.04 |

Max |A−B| = 0.20 eV on an 80 eV scale (0.25% of scale; the % vs local value peaks at a
zero crossing, which is a percentage artefact). Limiting cases: point-charge check at
z_p = −50 (σ_pot 0.354 vs 0.05): 48.966 vs 48.972 eV ✓; slab-centre closed-form ✓.

Re-validation in the qsp_phase3 geometry (Lz = 90, for E_pb(t) along the classical
track): |A−B| ≤ 0.23 eV at z_p ∈ {−23.75, −18.0, −14.7} — PASSED.

## Integration check (known-case gate, advisor ruling 4)

The exact four-term t=0 decomposition (E_wb + E_selfH + E_bgw − E_ghb, all terms in the
run's own p2 Poisson convention, ghost from the PARSED UPF truncated at its 50-Bohr mesh
end with lateral images) reproduces the measured d(H+E) of the h0_p2 ledger at
r = {4, 12, 28, 40} to within ±4 eV on 40–170 eV terms (≈2%); the ghost-model ablation is
decisive (no-images +96 eV off; untruncated −200…−510 eV off). See the campaign notebook
§B1 and `b1_decomposition.py`.

## Catalogue

Test rows added to `docs/validation/test-catalogue.md` (E_proj_bg dual-route; B1
decomposition known-case gate).
