// ============================================================================
// electron_proj_E25_L50_cubic_sigma1.hpp — Run-1 σ=1 follow-up (2026-05-19).
//
// Same bath as Run-1 (L=50, N=162, dx=0.40) but σ=1 WP. Pairs with
// run_classical_n162_L50_E25 (done as Run-9, 4.84h wall, S=0.88 eV/Bohr).
//
// Goal: σ=1 vs classical agreement at v=1.356 — strengthens the
// σ_w ≪ r_s campaign rule, complements σ=1 task #23 (v=2.71).
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_E25_L50_sigma1 : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;
    static constexpr double DT_AU        = 0.020;

    static constexpr double WP_EKIN_EV      = 25.0;
    static constexpr double WP_SIGMA_BOHR   = 1.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 1.356
    static constexpr double WP_KX = 0.0, WP_KY = 0.0, WP_KZ = WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =  0.0;
    static constexpr double WP_CY_BOHR =  0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -21

    static constexpr int N_STEPS     = 475;
    static constexpr int WRITE_EVERY = boundary::write_every_for(N_STEPS);

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E25_L50_sigma1_WP_dx0p40 : Common_E25_L50_sigma1 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

}  // namespace jellium::config
