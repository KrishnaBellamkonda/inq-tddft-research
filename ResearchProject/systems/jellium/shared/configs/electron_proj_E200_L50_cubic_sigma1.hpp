// ============================================================================
// shared/configs/electron_proj_E200_L50_cubic_sigma1.hpp — σ=1 WP at E=200 eV
// for the 2026-05-21 meeting follow-up. Companion to
// electron_proj_E100_L50_cubic_sigma1.hpp; only WP_EKIN_EV (100 → 200) and
// N_STEPS (475 → 385) change.
//
// SELF-SPREAD BOUND at σ=1 (same physics as the E=100 σ=1 case): the WP
// density σ_density(t) = (1/√2)·√(1 + t²) grows in vacuum independent of v.
// The 3σ tail of the WP reaches the far box face (+25 Bohr) at the time t*
// satisfying
//
//     launch_z + v·t* + 3·σ_density(t*) = +25
//     −21 + 3.840·t* + (3/√2)·√(1 + t²) = 25
//
// Solving: t* ≈ 7.70 a.u. (vs 9.50 at E=100). N_STEPS = 385 at dt=0.020.
//
// Nyquist: σ_p = 1/(σ√2) = 0.707 Bohr⁻¹; k_max_WP = k₀ + 3σ_p = 3.84 + 2.12
// = 5.96 Bohr⁻¹. k_Nyquist at dx=0.40 = 7.85 Bohr⁻¹. Ratio 0.76 — comfortable.
//
// GS reuse: same checkpoints/gs_L50_cubic_N162_dx0p40 as the E=100 σ=1 run.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_E200_L50_sigma1 : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;

    static constexpr double DT_AU        = 0.020;

    static constexpr double WP_EKIN_EV      = 200.0;
    static constexpr double WP_SIGMA_BOHR   = 1.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV); // 3.840
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -21

    // Self-spread-capped (see header note). t* ≈ 7.70 a.u. → N_STEPS = 385.
    static constexpr int N_STEPS     = 385;
    static constexpr int WRITE_EVERY = 2;              // 193 frames

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

struct Electron_Proj_E200_L50_sigma1_WP_dx0p40 : Common_E200_L50_sigma1 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

}  // namespace jellium::config
