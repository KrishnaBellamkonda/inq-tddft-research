// ============================================================================
// shared/configs/electron_proj_E600_L50_cubic.hpp  — higher-energy companion
//   in the meeting case-study sweep. Sits at v ≈ 6.640 a.u. between the
//   v_F-peak (~0.5–1.0) and the Bethe asymptote (10.5).
//
// Two Cfgs sharing one cubic 50 x 50 x 50 Bohr periodic jellium bath
// (N=162, dx=0.40 Bohr, dt=0.020 a.u., N_STEPS=264):
//
//   Electron_Proj_E600_L50_cubic_WP_dx0p40         — Gaussian wave packet
//   Electron_Proj_E600_L50_cubic_Classical_dx0p40  — classical electron
//                                                    (custom UPF + mass override)
//
// Launch position IDENTICAL to the E=100 / E=1500 pairs: (0, 0, -10) Bohr.
// Velocity is the only thing that differs across the sweep.
//
// At E=600 eV: KE_Ha = 22.05, v = sqrt(2 * 22.05) = 6.6404 bohr/atu.
//   - WP:        k_0 = 6.6404 bohr⁻¹, sigma = 5.0 Bohr.
//   - Classical: v_z = 6.6404 bohr/atu, m = m_e.
//
// NYQUIST CHECK (chosen by the user to avoid aliasing):
//   k_Nyquist = π / 0.40 = 7.854 bohr⁻¹.
//   k_max (WP) = k_0 + 3 sigma_k = 6.6404 + 0.60 = 7.240 bohr⁻¹.
//   Ratio: 7.240 / 7.854 = 0.92 — 8 % Nyquist headroom (vs E=1000 which
//   would be 1.09 / over-Nyquist, and E=1500 which was 1.41 / heavily
//   aliased). E=600 is the highest energy at which the WP is physically
//   clean at dx = 0.40.
//
// PROPAGATION TIME (per user instruction "projectile reaches boundary
// at initial speed"):
//   t_total = (L_box - z_launch) / v_initial
//           = (25 - (-10)) / 6.6404
//           = 35 / 6.6404 = 5.271 a.u.
//   N_STEPS = ceil(5.271 / 0.020) = 264.
//
// After t = 5.271 a.u., projectile centroid reaches z = -10 + 6.6404 *
// 5.271 = +25 Bohr (= +z box face). Classical decelerates ~few % under
// .ehrenfest() so actually reaches z ≈ +24.something.
//
// TIME-STEP CHOICE:
//   dt = 0.020 a.u. matches the E=100 and the L=50 plasmon runs.
//   ETRS stability: dt * v = 0.020 * 6.640 = 0.133 ≪ dx = 0.40 ✓
//   (~3× margin; the E=1500 case at dt=0.005 was 4× finer because
//    v=10.5 there pushed dt*v closer to the dx limit).
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"

namespace jellium::config {

// Common bath/grid; both WP and classical Cfgs share these.
struct Common_E600_L50_cubic : Base_N162_L50_E1p5 {

    // ----- Grid: same dx=0.40 as the E=100 / E=1500 classical → reuse
    //   checkpoints/gs_L50_cubic_N162_dx0p40/ . No new GS needed.
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;       // pi^2/(2 dx^2);
                                                        // log only.

    // ----- Real-time -----------------------------------------------------
    static constexpr double DT_AU       = 0.020;
    // N_STEPS = ceil(35 / v_initial / dt) = ceil(35 / 6.6404 / 0.02) = 264.
    static constexpr int    N_STEPS     = 264;          // total t = 5.28 a.u.
    static constexpr int    WRITE_EVERY = 2;            // 132 density frames

    // ----- Wave-packet projectile (E_kin = 600 eV) ----------------------
    static constexpr double WP_EKIN_EV      = 600.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 6.6404
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // WP_SIGMA_BOHR inherited (5.0 Bohr) from Base_N138_L50_E1p5.

    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR = -10.0;     // same as E=100, E=1500

    // ----- Classical-electron projectile (same KE, m = m_e) -------------
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;

    static constexpr double PROJ_LAUNCH_X =   0.0;
    static constexpr double PROJ_LAUNCH_Y =   0.0;
    static constexpr double PROJ_LAUNCH_Z = -10.0;     // matches WP launch

    static constexpr double PROJ_VEL_X =   0.0;
    static constexpr double PROJ_VEL_Y =   0.0;
    static constexpr double PROJ_VEL_Z = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    // Time-window placeholders.
    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E600_L50_cubic_WP_dx0p40 : Common_E600_L50_cubic {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Electron_Proj_E600_L50_cubic_Classical_dx0p40 : Common_E600_L50_cubic {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
