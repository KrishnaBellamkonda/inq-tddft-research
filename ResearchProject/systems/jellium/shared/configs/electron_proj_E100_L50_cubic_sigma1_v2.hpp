// ============================================================================
// shared/configs/electron_proj_E100_L50_cubic_sigma1_v2.hpp — v2 rerun with:
//   * dt = 0.01 a.u. (halved from 0.02 for better temporal resolution)
//   * WP orbital wavefunction + density saving enabled
//   * All-orbital dump at final IFW step for secondary electron analysis
//
// Self-spread bound (same physics, but dt=0.01 → 2× N_STEPS):
//   t* ≈ 9.5 a.u. → N_STEPS = 950 (was 475 at dt=0.02)
//
// Nyquist: unchanged (dx=0.40, k_max=4.83 → dx_max=0.650 > 0.40 ✓)
// Memory: unchanged (same grid, same GS)
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_E100_L50_sigma1_v2 : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;

    static constexpr double DT_AU        = 0.010;   // v2: halved from 0.020

    static constexpr double WP_EKIN_EV      = 100.0;
    static constexpr double WP_SIGMA_BOHR   = 1.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -21

    // Self-spread-capped: t* ≈ 9.5 a.u. at dt=0.01 → 950 steps
    static constexpr int N_STEPS     = 950;
    static constexpr int WRITE_EVERY = 3;              // ~317 density frames

    // WP wavefunction saving cadence (heavier I/O — less frequent)
    static constexpr int WF_WRITE_EVERY = 10;          // ~95 wavefunction frames

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

struct Electron_Proj_E100_L50_sigma1_v2_WP : Common_E100_L50_sigma1_v2 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Electron_Proj_E100_L50_sigma1_v2_Classical : Common_E100_L50_sigma1_v2 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
    // Boundary-rule N_STEPS: traversal=45 / (v=2.711 * dt=0.01) → 1660
    static constexpr int N_STEPS     = 1660;
    static constexpr int WRITE_EVERY = 6;              // ~277 density frames
};

}  // namespace jellium::config
