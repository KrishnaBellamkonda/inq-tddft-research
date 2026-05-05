// ============================================================================
// shared/configs/base_n138_L50_E1p5.hpp
//
// Variant of the project base aimed at a regime where:
//   * inelastic e-h coupling can be cleanly observed without the
//     periodic-box revival contamination of the L=30 base run;
//   * the WP momentum spread is much smaller than k_0 so the WP is a
//     well-defined moving Gaussian and not a barely-moving stationary
//     blob (sigma_r enlarged from 1.0 to 5.0 Bohr).
//
// Inherits from Base. Overrides:
//
//   * L_BOHR = 50 (vs project base 30). Revival timescale ~ (m sigma L^2)/(2 pi)
//     scales with L^2, so L=50 buys ~2.8x longer until the kinematic
//     revival becomes visible — enough to interpret the WP slowdown as
//     genuinely inelastic.
//
//   * N_ELECTRONS = 138 (closed shell |G|^2 <= 6, same shell as the
//     project base). Density 138/L^3 = 1.104e-3 e/Bohr^3, r_s ~= 5.97
//     Bohr. Lithium-like regime, between the L=30 sodium-like and the
//     legacy L=60 dilute case.
//
//   * EXTRA_STATES = 20 (project standard).
//
//   * WP_EKIN_EV = 1.5 (vs base 5). Lowest e-h gap at L=50 is
//     0.5 * (2 pi/50)^2 * 2 = 0.01579 Ha = 0.430 eV, so 1.5 eV WP is
//     ~3.5x the gap — well in the resonant single-particle excitation
//     regime, but not so close to the gap that the WP can only excite
//     one transition.
//
//   * WP_SIGMA_BOHR = 5.0 (vs base 1.0015). sigma_k = 1/sigma_r = 0.20
//     Bohr^-1. With k_0 = 0.332 Bohr^-1 (1.5 eV), the WP momentum is
//     resolved at k_0/sigma_k = 1.66 - moderately well-defined moving
//     Gaussian. At sigma=1 the WP was 50% backward-going by k content;
//     at sigma=5 only ~5% is.
//
//   * SPACING_BOHR = 1.0 (vs base 0.85). Nyquist:
//        sigma_k = 1/sigma_r        = 0.200 Bohr^-1
//        k_max   = k_0 + 3 sigma_k = 0.332 + 0.600 = 0.932 Bohr^-1
//        required dx <= pi/k_max   = 3.37 Bohr
//        chosen dx = 1.0           (k_Nyquist = 3.14 Bohr^-1, ~3.4x margin)
//
//   * N_STEPS = 1500 at dt=0.020 -> t_final = 30.0 a.u. = 0.726 fs
//     (vs project standard 19.80 a.u.). The slower WP needs longer
//     flight: travel = k_0 dt N_STEPS = 0.332 * 0.020 * 1500 = 9.96 Bohr.
//     Combined with the 3 sigma_r = 15 Bohr WP envelope, the leading
//     edge reaches z ~= 25 Bohr = L/2. This deviation from the project
//     standard is permitted per .claude/rules/jellium-base-run-spec.md
//     "When to deviate" -> "WP-energy sweep ... documented in the
//     run-specific config header".
//
//   * WP centre = (0,0,0) (box centre, INQ centred-Cartesian).
//   * SCF_TOL_HA = 1e-6 (project standard).
//
// New observable enabled in this run via the run_template:
//   - inqkit::observables::OccupationsWriter dumps f_i(t) at every
//     5*WRITE_EVERY = 10 steps, so the energy_balance postprocess can
//     compute occupation-weighted Delta-E_bath without joining with a
//     separate GS occupations CSV.
// ============================================================================
#pragma once

#include "base.hpp"

namespace jellium::config {

struct Base_N138_L50_E1p5 : Base {
    static constexpr double L_BOHR  = 50.0;
    static constexpr double LX_BOHR = L_BOHR;
    static constexpr double LY_BOHR = L_BOHR;
    static constexpr double LZ_BOHR = L_BOHR;

    static constexpr int    N_ELECTRONS    = 138;
    static constexpr int    EXTRA_STATES   = 20;
    static constexpr double SCF_TOL_HA     = 1.0e-6;
    static constexpr double SPACING_BOHR   = 1.0;

    // WP at 1.5 eV → k₀ = √(2 · 1.5 / 27.21138625) ≈ 0.33212 Bohr⁻¹
    static constexpr double WP_EKIN_EV      = 1.5;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // Wider WP — addresses the σ_k ~ k₀ pathology at base-run σ.
    static constexpr double WP_SIGMA_ANG  = 5.0 / ANG_TO_BOHR;     // ~2.65 Å
    static constexpr double WP_SIGMA_BOHR = 5.0;

    static constexpr double WP_CX_BOHR = 0.0;
    static constexpr double WP_CY_BOHR = 0.0;
    static constexpr double WP_CZ_BOHR = 0.0;

    static constexpr double DT_AU             = 0.020;
    static constexpr int    N_STEPS           = 1500;
    static constexpr int    WRITE_EVERY       = 2;
    static constexpr int    SCREEN_SNAP_EVERY = 6;

    static constexpr int    N_SCREENS = 20;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = DT_AU * N_STEPS / FS_TO_AU;
    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
};

}  // namespace jellium::config
