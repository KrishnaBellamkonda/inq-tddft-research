// ============================================================================
// shared/configs/highdens_n162_L30_E50_sigma1.hpp — σ=1 high-density run at
// E=50 eV (L=30, r_s≈3.41). Part of the density-dependence σ=1 sweep
// alongside E={100, 200, 300} at the same density.
//
// SELF-SPREAD BOUND (σ=1, L=30, v=1.916):
//   -11 + 1.916·t + 2.121·√(1+t²) = 15  →  t* ≈ 6.4 a.u.  →  N_STEPS = 320.
//
// Nyquist: k₀=1.916, σ_k=0.707. k_max=4.04. k_Nyq=7.85. Ratio 0.51 — OK.
//
// GS reuse: gs_L30_cubic_N162_dx0p40.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_HighDens_N162_L30_E50_sigma1 : Base_N162_L50_E1p5 {
    static constexpr double L_BOHR  = 30.0;
    static constexpr double LX_BOHR = L_BOHR;
    static constexpr double LY_BOHR = L_BOHR;
    static constexpr double LZ_BOHR = L_BOHR;
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;

    static constexpr double DT_AU = 0.020;

    static constexpr double WP_EKIN_EV      = 50.0;
    static constexpr double WP_SIGMA_BOHR   = 1.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV); // 1.916
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =  0.0;
    static constexpr double WP_CY_BOHR =  0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -11

    static constexpr int N_STEPS     = 320;
    static constexpr int WRITE_EVERY = 2;              // 160 frames

    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;
    static constexpr double PROJ_LAUNCH_X = 0.0;
    static constexpr double PROJ_LAUNCH_Y = 0.0;
    static constexpr double PROJ_LAUNCH_Z = WP_CZ_BOHR;
    static constexpr double PROJ_VEL_X    = 0.0;
    static constexpr double PROJ_VEL_Y    = 0.0;
    static constexpr double PROJ_VEL_Z    = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct HighDens_N162_L30_E50_sigma1_WP_dx0p40 : Common_HighDens_N162_L30_E50_sigma1 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct HighDens_N162_L30_E50_sigma1_Classical_dx0p40 : Common_HighDens_N162_L30_E50_sigma1 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
