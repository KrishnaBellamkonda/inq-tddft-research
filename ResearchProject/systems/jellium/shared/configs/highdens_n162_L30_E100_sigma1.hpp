// ============================================================================
// shared/configs/highdens_n162_L30_E100_sigma1.hpp — σ=1 companion to the
// highdens (L=30, r_s=3.41) WP run. Same density as
// highdens_n162_L30_E100.hpp, but with σ=1 Bohr instead of σ=0.5, so the
// density-effect can be isolated from the σ-effect when compared to the
// standard-density σ=1 run (run_wp_n162_L50_E100_sigma1).
//
// SELF-SPREAD BOUND at σ=1 (same physics as the L=50 σ=1 case): the WP
// density σ_density(t) = (1/√2)·√(1 + t²) grows in vacuum independent of v.
// The 3σ_density tail reaches the far box face (+L/2=+15) at the time t*
// satisfying
//
//     launch_z + v·t* + (3/√2)·√(1+t²) = +15
//     -11 + 2.711·t* + 2.121·√(1+t²)  = 15  →  t* ≈ 5.35 a.u.
//
// N_STEPS = 268 at dt=0.020.
//
// Nyquist: σ_p = 1/(σ√2) = 0.707 Bohr⁻¹; k_max_WP = k₀ + 3σ_p = 2.71 + 2.12
// = 4.83. k_Nyquist at dx=0.40 = 7.85. Ratio 0.61 — comfortable.
//
// GS reuse: same checkpoints/gs_L30_cubic_N162_dx0p40 as the σ=0.5 highdens.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_HighDens_N162_L30_E100_sigma1 : Base_N162_L50_E1p5 {
    static constexpr double L_BOHR  = 30.0;
    static constexpr double LX_BOHR = L_BOHR;
    static constexpr double LY_BOHR = L_BOHR;
    static constexpr double LZ_BOHR = L_BOHR;
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;

    static constexpr double DT_AU = 0.020;

    static constexpr double WP_EKIN_EV      = 100.0;
    static constexpr double WP_SIGMA_BOHR   = 1.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV); // 2.7111
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =  0.0;
    static constexpr double WP_CY_BOHR =  0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -11

    // Self-spread-capped (see header note). t* ≈ 5.35 a.u. → N_STEPS = 268.
    static constexpr int N_STEPS     = 268;
    static constexpr int WRITE_EVERY = 2;              // 134 frames

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

struct HighDens_N162_L30_E100_sigma1_WP_dx0p40 : Common_HighDens_N162_L30_E100_sigma1 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

}  // namespace jellium::config
