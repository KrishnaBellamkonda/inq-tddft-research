// ============================================================================
// shared/configs/electron_proj_backlog_L50_cubic.hpp — backlog re-runs added
// 2026-05-20.
//
// Re-runs of the legacy E={50, 300, 600} eV WP cases with the new
// WPMomentumStats / WPRealSpaceStats / energy_balance observables that
// the legacy versions lacked. Also adds a fresh classical companion at
// E=600 (legacy didn't have one).
//
// All four runs share the L=50 cubic N=162 dx=0.40 bath (reuses the
// existing gs_L50_cubic_N162_dx0p40 checkpoint — no new GS needed).
//
// dt=0.020, σ=5. Standard boundary_rule (launch_z = -L/2 + 4σ = -5).
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

template <int ENERGY_EV>
struct Common_Backlog_L50_cubic : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;
    static constexpr double DT_AU        = 0.020;

    static constexpr double WP_EKIN_EV      = double(ENERGY_EV);
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);
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

    // Classical companion
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

#define DECLARE_BACKLOG_CFG(EV)                                                  \
struct Backlog_L50_cubic_WP_E##EV : Common_Backlog_L50_cubic<EV> {              \
    static constexpr bool WP_ENABLED = true;                                     \
    static constexpr int  N_IONS     = 0;                                        \
};                                                                               \
struct Backlog_L50_cubic_Classical_E##EV : Common_Backlog_L50_cubic<EV> {       \
    static constexpr bool WP_ENABLED = false;                                    \
    static constexpr int  N_IONS     = 1;                                        \
};

DECLARE_BACKLOG_CFG(50)
DECLARE_BACKLOG_CFG(300)
DECLARE_BACKLOG_CFG(600)

#undef DECLARE_BACKLOG_CFG

}  // namespace jellium::config
