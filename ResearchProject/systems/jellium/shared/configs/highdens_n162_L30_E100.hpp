// ============================================================================
// shared/configs/highdens_n162_L30_E100.hpp  — Run-6 of the 2026-05-21
// meeting campaign: high-density jellium probe.
//
// Shrinks the box from L=50 to L=30 at fixed N=162, raising the density
// from 1.296e-3 e/Bohr^3 (r_s=5.69, Li-like) to 6.000e-3 e/Bohr^3
// (r_s=3.41, between Li and Al — sodium-to-aluminium regime). Same WP
// energy as Run-3 (100 eV) so the only changed variable is the bath
// density. N=162 remains closed-shell at any cubic L per
// docs/sources/free-electron-gas-magic-numbers.md (filled shells through
// |G|^2 <= 6 in the gamma-only k-point convention).
//
// SIGMA SHRUNK TO 0.5 Bohr (vs 5 Bohr at L=50). Reason: σ=5 in L=30 would
// have its 3σ_density tail of 10.6 Bohr stretching across more than half
// the box at t=0 — geometrically infeasible. σ=0.5 keeps the WP a true
// localised probe in the smaller bath, with 3σ_density tail of 1.06 Bohr
// well clear of the L/2=15 box face. The injector's σ parameter convention
// (psi ~ exp(-r^2/(2σ²))) means density σ_r is σ/√2 = 0.354 Bohr.
//
// Cross-validation point: the σ-sweep family (Run-7) includes σ=0.5 at
// L=50 / E=100 eV, so the Run-6 high-density vs Run-7 σ=0.5 isolates the
// pure density effect.
//
// PROPAGATION (boundary_rule, σ=0.5, L=30, v=2.711, dt=0.02):
//   launch_z   = -L/2 + 4σ = -13.0  Bohr
//   stop_z     = +L/2 -   σ = +14.5 Bohr
//   traversal  = L - 5σ = 27.5 Bohr
//   N_STEPS    = ceil(27.5 / 2.711 / 0.02) = 508
//   WRITE_EVERY = write_every_for(508) ≈ 2 → 254 frames.
//
// NYQUIST: σ_p = 1/(σ√2) = 1.414 Bohr^-1; k_max_WP = k0 + 3σ_p = 2.711 +
// 4.24 = 6.95 Bohr^-1. k_Nyquist at dx=0.40 = π/0.40 = 7.85, so 6.95/7.85
// = 88% — comfortable.
//
// SHARED GS: save_gs/gs_L30_cubic_N162_dx0p40/ (new, ~30 min on one A30).
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_HighDens_N162_L30_E100 : Base_N162_L50_E1p5 {
    // Override geometry: smaller box, same N → r_s = 3.41
    static constexpr double L_BOHR  = 30.0;
    static constexpr double LX_BOHR = L_BOHR;
    static constexpr double LY_BOHR = L_BOHR;
    static constexpr double LZ_BOHR = L_BOHR;
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;       // pi^2/(2 dx^2)

    static constexpr double DT_AU       = 0.020;

    // WP @ 100 eV, σ shrunk for the smaller box.
    static constexpr double WP_EKIN_EV      = 100.0;
    static constexpr double WP_SIGMA_BOHR   = 0.5;       // override base 5.0
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV); // 2.7111
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =  0.0;
    static constexpr double WP_CY_BOHR =  0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -13

    static constexpr int N_STEPS =
        boundary::n_steps_for(WP_SIGMA_BOHR, L_BOHR, WP_K0, DT_AU);   // 508
    static constexpr int WRITE_EVERY =
        boundary::write_every_for(N_STEPS);                           // 2

    // Classical-electron projectile (matched KE, m = m_e).
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;

    static constexpr double PROJ_LAUNCH_X = 0.0;
    static constexpr double PROJ_LAUNCH_Y = 0.0;
    static constexpr double PROJ_LAUNCH_Z = WP_CZ_BOHR;                // -13
    static constexpr double PROJ_VEL_X    = 0.0;
    static constexpr double PROJ_VEL_Y    = 0.0;
    static constexpr double PROJ_VEL_Z    = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct HighDens_N162_L30_E100_WP_dx0p40 : Common_HighDens_N162_L30_E100 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct HighDens_N162_L30_E100_Classical_dx0p40 : Common_HighDens_N162_L30_E100 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
