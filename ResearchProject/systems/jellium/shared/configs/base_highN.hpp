// ============================================================================
// shared/configs/base_highN.hpp  (jellium, closed-shell high-density)
//
// Base_HighN inherits the L=60 cubic cell, WP injection (200 eV, σ=0.53 Å,
// +z launch, centred at origin), real-time grid, and screen layout from
// `Base`, and overrides only:
//
//   * N_ELECTRONS = 514  — closest closed-shell magic number above
//     4 × 128 = 512. Fully fills the |G|² ≤ 16 shell (cumulative count of
//     paired (n_x,n_y,n_z) integer triples). Density 2.380e-3 e/bohr³ is
//     4.015× the Base density (5.926e-4); r_s ≈ 4.64 bohr, sodium/
//     potassium regime. See `docs/sources/free-electron-gas-magic-numbers.md`.
//   * EXTRA_STATES = 8  — covers the next |G|²=17 shell (48 spatial
//     orbitals) without over-populating; absorbs Fermi smearing without
//     introducing partial-shell artefacts.
//   * SCF_TOL_HA = 1.0e-6 — tightened from 1e-4 to rule out
//     under-convergence as a confounding factor for the orbital-symmetry
//     inspection on the closed-shell case (user instruction 2026-05-03).
//
// All other Cfg fields (cell, WP, dt, N_STEPS, screens) are inherited
// verbatim from Base. The legacy partial-shell `Base` (N=128) is kept for
// reproducibility of the previous run.
// ============================================================================
#pragma once

#include "base.hpp"

namespace jellium::config {

struct Base_HighN : Base {
    static constexpr int    N_ELECTRONS  = 514;
    static constexpr int    EXTRA_STATES = 8;
    static constexpr double SCF_TOL_HA   = 1.0e-6;
};

}  // namespace jellium::config
