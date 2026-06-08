// ============================================================================
// shared/configs/electron_proj_E25_L50_cubic.hpp  — stretch run, WP only.
//
// Per plan §8.4: this Cfg is intended for an UNMATCHED WP run probing the
// quantum-scattering regime (κ = 1.47) where Born scattering theory breaks
// down. No matched classical sibling is budgeted (classical 25 eV would
// cost ~6.4 h alone, exceeding the 9 h GPU-0 budget).
//
// Only the WP struct is declared; the Classical struct is omitted on
// purpose. Any future run.cpp using this header should use
// `Electron_Proj_E25_L50_cubic_WP_dx0p40`.
//
// At E=25 eV: KE_Ha = 0.919, v = sqrt(2 * 0.919) = 1.356 bohr/atu.
// v/v_F = 4.02 (low end of the Bragg-peak band).
// κ = 2/v = 1.475 (Bohr classical regime).
// k₀σ = 6.78 (lower end of the classical-packet limit; structure may begin
//   to appear, see plan §6.2 σ_z diagnostic).
//
// NYQUIST: k_max = 1.356 + 0.60 = 1.956. k_Nyquist = 7.854. Ratio 0.25 —
//   very clean.
//
// PROPAGATION (updated 2026-05-17 to the universal boundary_rule):
//   launch_z = -L/2 + 4σ = -5      (was -10 in the pre-rule draft)
//   stop_z   = +L/2 -  σ = +20
//   traversal = L - 5σ = 25 Bohr
//   t_total  = 25 / k0(25 eV) = 25 / 1.35553 = 18.443 a.u.
//   N_STEPS  = boundary::n_steps_for(σ=5, L=50, k0, dt=0.02) = 923
//              (plan §"Exact next steps" quotes 922 with v=1.356 rounded;
//              constexpr at full k0 precision gives 923. One extra dt
//              step ⇒ 0.027 Bohr past stop_z, negligible.)
//   WRITE_EVERY = boundary::write_every_for(923) = 3
//              ⇒ 307 frames, on target.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_E25_L50_cubic : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;       // log only

    static constexpr double DT_AU       = 0.020;

    static constexpr double WP_EKIN_EV      = 25.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 1.356
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                     // -5

    static constexpr int    N_STEPS     =
        boundary::n_steps_for(WP_SIGMA_BOHR, L_BOHR, WP_K0, DT_AU);    // 922
    static constexpr int    WRITE_EVERY =
        boundary::write_every_for(N_STEPS);                            // 3

    // Classical-electron projectile companion (Run-9, added 2026-05-18).
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;
    static constexpr double PROJ_LAUNCH_X = 0.0;
    static constexpr double PROJ_LAUNCH_Y = 0.0;
    static constexpr double PROJ_LAUNCH_Z = WP_CZ_BOHR;                    // -5
    static constexpr double PROJ_VEL_X    = 0.0;
    static constexpr double PROJ_VEL_Y    = 0.0;
    static constexpr double PROJ_VEL_Z    = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E25_L50_cubic_WP_dx0p40 : Common_E25_L50_cubic {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

// Run-9: classical companion to Run-1 at E=25 eV
struct Electron_Proj_E25_L50_cubic_Classical_dx0p40 : Common_E25_L50_cubic {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
