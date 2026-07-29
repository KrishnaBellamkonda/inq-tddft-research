// ============================================================================
// shared/configs/electron_proj_E20_L50_cubic.hpp — Run-5: 20 eV pair (2026-05-18)
//
// WP+Classical pair at E=20 eV (DNA-radiolysis energy scale, Knudsen
// regime). v = k0 = 1.212 Bohr/atu. Standard 4σ/1σ rule:
//   launch_z = -5, stop_z = +20, traversal = 25, N_STEPS ≈ 1032.
//
// Reuses gs_L50_cubic_N162_dx0p40 (no new GS).
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_E20_L50_cubic : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;
    static constexpr double DT_AU        = 0.020;

    static constexpr double WP_EKIN_EV      = 20.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV); // 1.212
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =  0.0;
    static constexpr double WP_CY_BOHR =  0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -5

    static constexpr int N_STEPS =
        boundary::n_steps_for(WP_SIGMA_BOHR, L_BOHR, WP_K0, DT_AU);
    static constexpr int WRITE_EVERY = boundary::write_every_for(N_STEPS);

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

struct Electron_Proj_E20_L50_cubic_WP_dx0p40 : Common_E20_L50_cubic {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Electron_Proj_E20_L50_cubic_Classical_dx0p40 : Common_E20_L50_cubic {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
