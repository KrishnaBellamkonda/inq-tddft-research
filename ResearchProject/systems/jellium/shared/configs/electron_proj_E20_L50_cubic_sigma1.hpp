// ============================================================================
// electron_proj_E20_L50_cubic_sigma1.hpp — Run-5 σ=1 follow-up (2026-05-19).
//
// Same bath as Run-5 (L=50, N=162, dx=0.40) but with the σ=1 concentrated
// WP. Pairs with run_classical_n162_L50_E20 (already done, 4.82h wall).
// Self-spread capped at 475 steps per the σ=1 convention.
//
// Goal: test whether σ=1 gives TDDFT-classical agreement at v=1.21 the
// way it did at v=2.71 (1% agreement, σ=1 task #23). If yes, the
// σ_w ≪ r_s rule is universal across the velocity range.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_E20_L50_sigma1 : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;
    static constexpr double DT_AU        = 0.020;

    static constexpr double WP_EKIN_EV      = 20.0;
    static constexpr double WP_SIGMA_BOHR   = 1.0;        // override base 5.0
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 1.212
    static constexpr double WP_KX = 0.0, WP_KY = 0.0, WP_KZ = WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =  0.0;
    static constexpr double WP_CY_BOHR =  0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -21

    // Self-spread cap: σ=1 density 3σ tail hits the box face at
    // t ≈ 9.5 a.u. regardless of v.
    static constexpr int N_STEPS     = 475;
    static constexpr int WRITE_EVERY = boundary::write_every_for(N_STEPS);

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E20_L50_sigma1_WP_dx0p40 : Common_E20_L50_sigma1 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

}  // namespace jellium::config
