// ============================================================================
// shared/configs/electron_proj_E50_L50_cubic.hpp  — Bragg-peak-onset pair.
//
// This is the *quantum-dominated* point in the regime sweep. At κ ≈ 1.04 and
// v/v_F ≈ 5.69, both Born and Bohr scattering pictures partially apply, and
// the host response is in the Bragg-peak crossover region where the full
// Lindhard dielectric matters. The largest WP-vs-classical Δ(E) of the
// sweep is expected here (plan §6.3, §9.4 Test 4).
//
// Two Cfgs sharing one cubic 50 × 50 × 50 Bohr periodic jellium bath
// (N=162, dx=0.40 Bohr, dt=0.020 a.u., N_STEPS=913):
//
//   Electron_Proj_E50_L50_cubic_WP_dx0p40         — Gaussian wave packet
//   Electron_Proj_E50_L50_cubic_Classical_dx0p40  — classical electron
//
// Launch (0, 0, -10) Bohr (IDENTICAL across the sweep).
// At E=50 eV: KE_Ha = 1.838, v = sqrt(2 * 1.838) = 1.917 bohr/atu.
//
// NYQUIST: k_max = k₀ + 3 σ_k = 1.917 + 0.60 = 2.517. k_Nyquist = π/0.40
//   = 7.854. Ratio 0.32 — very clean.
//
// PROPAGATION: t_total = 35 / 1.917 = 18.26 a.u. N_STEPS = ceil(18.26 /
//   0.020) = 913.
//
// EXPECTED PHYSICS (plan §4.3 + §4.5):
//   Bethe-Lindhard prediction S = 0.426 eV/Bohr. With Bloch correction at
//   κ = 1.04: S = 0.407 eV/Bohr. After ≤5 % box deficit (plan §6.1,
//   q_min·v/ω_p = 1.92), expect 0.32 – 0.41 eV/Bohr measured.
//
// k₀σ = 1.917 × 5 = 9.6 — still in the classical-packet limit by a thin
// margin (plan §1, Q2). Note: the packet's 3σ envelope ≈ 15 Bohr is
// comparable to the launch-offset |z₀| = 10, so the leading edge of the
// packet is at z = +5 at t=0. Plan §8.3 Sim 4 cost note: use
// t_start ≈ 1.6 a.u. (Δz = 3 Bohr) for windowed fit.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"

namespace jellium::config {

struct Common_E50_L50_cubic : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;       // log only

    static constexpr double DT_AU       = 0.020;
    static constexpr int    N_STEPS     = 913;          // total t = 18.26 a.u.
    static constexpr int    WRITE_EVERY = 2;            // 456 density frames

    static constexpr double WP_EKIN_EV      = 50.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 1.917
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

struct Electron_Proj_E50_L50_cubic_WP_dx0p40 : Common_E50_L50_cubic {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Electron_Proj_E50_L50_cubic_Classical_dx0p40 : Common_E50_L50_cubic {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
