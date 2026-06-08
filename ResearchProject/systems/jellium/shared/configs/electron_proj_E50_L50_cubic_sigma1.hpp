// ============================================================================
// shared/configs/electron_proj_E50_L50_cubic_sigma1.hpp  — σ=1 WP + classical
// companion at E=50 eV on the standard L=50 bath (r_s=5.69).
//
// Fills the E=50 gap in the σ=1 energy sweep (existing: E={20,25,100,200,300}).
//
// SELF-SPREAD BOUND (σ=1, L=50, v=1.916 Bohr/a.u.):
// σ_density(t) = (1/√2)·√(1+t²). The 3σ_density tail reaches the far box
// face (+25) when:
//   launch_z + v·t + (3/√2)·√(1+t²) = +25
//   -21 + 1.916·t + 2.121·√(1+t²) = 25
//   → t* ≈ 11.4 a.u.  →  N_STEPS = 570.
//
// Nyquist: k₀ = 1.916, σ_k = 1/(σ√2) = 0.707. k_max = 1.916 + 3×0.707
// = 4.04. k_Nyquist at dx=0.40 = 7.85. Ratio 0.51 — comfortable.
//
// GS reuse: same gs_L50_cubic_N162_dx0p40 as all other L=50 runs.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_E50_L50_sigma1 : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;

    static constexpr double DT_AU        = 0.020;

    static constexpr double WP_EKIN_EV      = 50.0;
    static constexpr double WP_SIGMA_BOHR   = 1.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV); // 1.916
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -21

    // Self-spread-capped (see header note). t* ≈ 11.4 a.u. → N_STEPS = 570.
    static constexpr int N_STEPS     = 570;
    static constexpr int WRITE_EVERY = 2;              // 285 frames

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

struct Electron_Proj_E50_L50_sigma1_WP_dx0p40 : Common_E50_L50_sigma1 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Electron_Proj_E50_L50_sigma1_Classical_dx0p40 : Common_E50_L50_sigma1 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
