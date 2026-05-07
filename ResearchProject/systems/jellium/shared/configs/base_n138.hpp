// ============================================================================
// shared/configs/base_n138.hpp  (jellium, closed-shell, low-energy WP)
//
// Base_N138 inherits the L=60 cubic cell, WP centring at origin, and screen
// layout from `Base`. Overrides:
//
//   * N_ELECTRONS = 138  — closed-shell magic number at |G|² ≤ 6 (cumulative
//     paired-(n_x,n_y,n_z) count). Density 6.39e-4 e/bohr³, r_s ≈ 7.21 bohr —
//     near the legacy N=128 r_s of 7.38, but on a fully-filled shell so the
//     GS density is uniform by symmetry. Source:
//     `docs/sources/free-electron-gas-magic-numbers.md`.
//   * EXTRA_STATES = 8  — covers the next |G|²=8 shell (12 spatial orbitals)
//     without over-populating; absorbs Fermi smearing.
//   * SCF_TOL_HA = 1.0e-6 — tightened SCF.
//   * WP_EKIN_EV = 100 (vs Base 200) → k₀ = 2.711 Bohr⁻¹ (vs 3.834).
//   * SPACING_BOHR = 0.55 (vs Base 0.50). Coarsening is bounded by the
//     Nyquist condition that the WP momenta up to k₀ + 3σ_k must be
//     resolved (3σ envelope of |φ̃_wp(k)|²):
//
//         σ_k = 1/σ_r = 1/1.0015  ≈ 0.999  Bohr⁻¹
//         k_max = k₀ + 3σ_k       = 2.711 + 2.996 = 5.707 Bohr⁻¹
//         Required dx ≤ π / k_max = 0.5505 Bohr
//         Chosen dx  = 0.55          (k_Nyquist = π/0.55 = 5.712 Bohr⁻¹)
//
//     This is exactly at the threshold; if the WP develops energetic
//     scattered components beyond 3σ, raise the cutoff or revert to 0.50.
//
//   * N_STEPS = 320 unchanged — at k₀ = 2.711 the WP travels
//     2.711·0.02·320 = 17.3 Bohr in the trajectory, well clear of the
//     periodic boundary at L/2 = 30 Bohr. Plenty of margin.
// ============================================================================
#pragma once

#include "base.hpp"

namespace jellium::config {

struct Base_N138 : Base {
    static constexpr int    N_ELECTRONS    = 138;
    static constexpr int    EXTRA_STATES   = 8;
    static constexpr double SCF_TOL_HA     = 1.0e-6;
    static constexpr double SPACING_BOHR   = 0.55;

    // WP at 100 eV → k₀ = √(2 · 100 / 27.21138625) ≈ 2.7110 Bohr⁻¹
    static constexpr double WP_EKIN_EV      = 100.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;
};

}  // namespace jellium::config
