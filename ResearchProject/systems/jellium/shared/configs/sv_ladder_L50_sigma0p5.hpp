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

// ----------------------------------------------------------------------------
// sigma=3.0 Bohr large-width probe (2026-06-15 extension; plan §12).
//
// A broad Gaussian (V(0)=0.532 Ry in-file ≈ 0.266 Ha physical) whose form
// factor e^{−q²σ²} strongly suppresses all q≳0.3 Bohr⁻¹ — expected to deviate
// markedly BELOW the point-charge Lindhard reference. Unlike the small-σ runs,
// the trailing 4σ=12 Bohr tail forces the STANDARD boundary rule launch:
//   launch_z = -L/2 + 4σ = -25 + 12 = -13 Bohr  (NOT the -20 of the small-σ
//   ladder, which would push the σ=3 tail 9 Bohr past the -25 face).
// stop_z = +22, traversal 35 Bohr; max centroid reach at v0=3 (const v) = +5.
// Needs its own build (launch_z is compile-time, not env-driven).
// ----------------------------------------------------------------------------
struct SV_Ladder_L50_sigma3p0 : SV_Ladder_L50_sigma0p5 {
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron_gaussian_sigma3p0.upf";
    static constexpr double PROJ_LAUNCH_Z = -13.0;   // -L/2 + 4σ, σ=3, L=50
};

}  // namespace jellium::config
