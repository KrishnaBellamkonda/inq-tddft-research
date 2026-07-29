// ============================================================================
// shared/configs/knudsen_sweep_L50_cubic.hpp — Run-8 velocity sweep, 2026-05-18.
//
// 5 WP+Classical pairs at E ∈ {700, 800, 900, 1000, 1100} eV. Bridges the
// v ≈ 2.7 → v ≈ 10.5 gap on the campaign's stopping-power curve and
// brackets the regime where Knudsen et al. (arXiv 2605.12854) claim
// classical-WP stopping-power convergence. All share:
//
//   bath:      L=50 Bohr cubic, N=162, dx=0.30 Bohr (existing GS)
//   WP:        σ_w = 5.0 Bohr, dt = 0.010 a.u.
//   launch_z:  -L/2 + 4σ = -5 Bohr (boundary_rule standard)
//   stop_z:    +L/2 - σ = +20 Bohr
//   traversal: L - 5σ = 25 Bohr
//
// dt choice: 0.01 a.u. (plan default; the dt-convergence subtest at
// E=1100 with dt=0.005 was deferred since the rollup figure can be
// recomputed if a future subtest indicates noticeable drift). At
// dt=0.01 the ETRS condition dt·v_max < dx → 0.01 · 8.999 = 0.090 ≪
// 0.30 is comfortably satisfied at the highest sweep energy.
//
// Nyquist: σ_p = 1/(σ√2) = 0.1414 Bohr⁻¹; k_max_WP = k0 + 3σ_p. At E=1100,
//          k_max = 9.42 vs k_Nyquist(dx=0.30) = 10.47 → 90 % usage.
//          Comfortable.
//
// Reused GS: checkpoints/gs_L50_cubic_N162_dx0p30/ — already exists.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

template <int ENERGY_EV>
struct Common_KnudsenSweep_L50_cubic : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.30;
    static constexpr double CUTOFF_HA    = 54.83;       // pi^2/(2 dx^2)

    static constexpr double DT_AU       = 0.010;

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
    static constexpr int WRITE_EVERY =
        boundary::write_every_for(N_STEPS);

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

// WP and Classical at each of the 5 sweep energies.
#define DECLARE_KNUDSEN_SWEEP_CFG(EV)                                          \
struct KnudsenSweep_L50_cubic_WP_E##EV : Common_KnudsenSweep_L50_cubic<EV> {   \
    static constexpr bool WP_ENABLED = true;                                   \
    static constexpr int  N_IONS     = 0;                                      \
};                                                                             \
struct KnudsenSweep_L50_cubic_Classical_E##EV : Common_KnudsenSweep_L50_cubic<EV> { \
    static constexpr bool WP_ENABLED = false;                                  \
    static constexpr int  N_IONS     = 1;                                      \
};

DECLARE_KNUDSEN_SWEEP_CFG(700)
DECLARE_KNUDSEN_SWEEP_CFG(800)
DECLARE_KNUDSEN_SWEEP_CFG(900)
DECLARE_KNUDSEN_SWEEP_CFG(1000)
DECLARE_KNUDSEN_SWEEP_CFG(1100)

#undef DECLARE_KNUDSEN_SWEEP_CFG

}  // namespace jellium::config
