// ============================================================================
// shared/configs/electron_proj_E300_L50_cubic.hpp  — mid-Bethe pair.
//
// Sits between the E=100 and E=600 pairs and the E=1500 anchor. Fills the
// regime-classification velocity sweep (docs/plans/jellium-regime-constrained-
// simulations.md §8.2/8.3) at the κ ≈ 0.43, v/v_F ≈ 13.9 point — squarely in
// the Bethe regime (mid-Bethe), still in the classical-projectile limit
// (κ < 1) and the classical-packet limit (k₀σ = 23.5).
//
// Two Cfgs sharing one cubic 50 × 50 × 50 Bohr periodic jellium bath
// (N=162, dx=0.40 Bohr, dt=0.020 a.u., N_STEPS=374):
//
//   Electron_Proj_E300_L50_cubic_WP_dx0p40         — Gaussian wave packet
//   Electron_Proj_E300_L50_cubic_Classical_dx0p40  — classical electron
//                                                    (custom UPF + mass override)
//
// Launch (0, 0, -10) Bohr (IDENTICAL across the E sweep).
// At E=300 eV: KE_Ha = 11.03, v = sqrt(2 * 11.03) = 4.6960 bohr/atu.
//
// NYQUIST: k_max = k₀ + 3 σ_k = 4.696 + 0.60 = 5.296. k_Nyquist = π/0.40
//   = 7.854. Ratio 0.67 — clean.
//
// PROPAGATION: t_total = 35 / 4.696 = 7.453 a.u. N_STEPS = ceil(7.453 /
//   0.020) = 374.
//
// EXPECTED PHYSICS (plan §4.3 Bethe-Lindhard prediction):
//   S = 0.107 eV/Bohr. After ~10 % box deficit (plan §6.1, q_min·v/ω_p
//   = 4.7), expect 0.094 – 0.100 eV/Bohr measured.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"

namespace jellium::config {

struct Common_E300_L50_cubic : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;       // log only

    static constexpr double DT_AU       = 0.020;
    static constexpr int    N_STEPS     = 374;          // total t = 7.48 a.u.
    static constexpr int    WRITE_EVERY = 2;            // 187 density frames

    static constexpr double WP_EKIN_EV      = 300.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 4.696
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR = -10.0;

    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;

    static constexpr double PROJ_LAUNCH_X =   0.0;
    static constexpr double PROJ_LAUNCH_Y =   0.0;
    static constexpr double PROJ_LAUNCH_Z = -10.0;

    static constexpr double PROJ_VEL_X =   0.0;
    static constexpr double PROJ_VEL_Y =   0.0;
    static constexpr double PROJ_VEL_Z = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E300_L50_cubic_WP_dx0p40 : Common_E300_L50_cubic {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Electron_Proj_E300_L50_cubic_Classical_dx0p40 : Common_E300_L50_cubic {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
