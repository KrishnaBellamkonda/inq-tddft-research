// shared/configs/sv_ladder_L50_sigma0p5.hpp
// ----------------------------------------------------------------------------
// S(v) ladder, free-Ehrenfest classical electron in the r_s=5.69 jellium bath.
//
// Decided in the grill-with-docs session 2026-06-12 (see
// docs/plans/overnight-gaussian-classical-jellium.md):
//   - projectile = -1 / sigma=0.5 Bohr erf-smoothed Gaussian charge
//     (electron_gaussian_sigma0p5.upf, V(0)=+1.596 Ha, repulsive),
//   - mass = m_e (NOT a fictitious mass), dynamics = free ehrenfest
//     (the projectile genuinely decelerates),
//   - launch +z on-axis, ONE-traversal path cap,
//   - reuses the validated GS checkpoints/gs_L50_cubic_N162_dx0p40.
//
// Velocity v0 (a.u.), N_STEPS, WRITE_EVERY and the output subdir are supplied
// at RUNTIME via environment variables (PROJ_V0 / SV_N_STEPS / SV_WRITE_EVERY /
// SV_OUT_SUBDIR) so a single build serves the whole ladder. The constants here
// are only fallback defaults.
// ----------------------------------------------------------------------------
#pragma once

#include "electron_proj_E100_L50_cubic.hpp"

namespace jellium::config {

struct SV_Ladder_L50_sigma0p5 : Common_E100_L50_cubic {
    // erf-smoothed sigma=0.5 electron (replaces the near-point ONCV psp).
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron_gaussian_sigma0p5.upf";

    // launch near the -z face with a small margin; travel +z one traversal.
    static constexpr double PROJ_LAUNCH_X = 0.0;
    static constexpr double PROJ_LAUNCH_Y = 0.0;
    static constexpr double PROJ_LAUNCH_Z = -20.0;   // 5 Bohr from the -25 face

    // fallback defaults (overridden at runtime):
    static constexpr double PROJ_VEL_Z_DEFAULT = 1.0;     // a.u.
    static constexpr int    N_STEPS_DEFAULT    = 400;
    static constexpr int    WRITE_EVERY_DEFAULT = 50;
};

// sigma=0.4 sensitivity sibling (only the psp path differs).
struct SV_Ladder_L50_sigma0p4 : SV_Ladder_L50_sigma0p5 {
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron_gaussian_sigma0p4.upf";
};

}  // namespace jellium::config
