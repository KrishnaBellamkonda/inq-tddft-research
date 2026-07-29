// ============================================================================
// shared/configs/positive_ion_L50_v0p33.hpp
//
// Positive-ion (proton) charge-conjugate companion to the negative-WP runs.
// Matched-velocity comparison: v_proton = v_F ≈ 0.337 a.u. (same as WP k_0
// in the L=50 closed-shell N=162 run).
//
// Aim: produce the textbook accumulation-wake behind a positive projectile,
// to compare directly against the depletion-anti-wake behind the negative WP.
// See docs/plans/jellium_positive_ion_companion.md for the full motivation.
//
// Differs from Base_N162_L50_E1p5: no WP injection; one H atom added at
// the launch position with impulsive velocity along +z; finer dt and more
// steps to keep total time at 30 a.u. = 0.726 fs.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"

namespace jellium::config {

struct Positive_Ion_L50_v0p33 : Base_N162_L50_E1p5 {
    // ----- Disable WP injection -----------------------------------------
    static constexpr bool        WP_ENABLED       = false;

    // ----- Add one proton at (0, 0, -L/4), velocity (0, 0, +v_F) --------
    static constexpr int         N_IONS           = 1;
    static constexpr const char* ION_SPECIES      = "H";
    static constexpr double      ION_LAUNCH_X     = 0.0;
    static constexpr double      ION_LAUNCH_Y     = 0.0;
    static constexpr double      ION_LAUNCH_Z     = -12.5;        // -L/4 Bohr
    static constexpr double      ION_VELOCITY_X   = 0.0;
    static constexpr double      ION_VELOCITY_Y   = 0.0;
    static constexpr double      ION_VELOCITY_Z   = 0.3320360921808203;  // = v_F at N=162 L=50, matches WP k_0

    // ----- Propagation grid (finer dt for ion stability) ----------------
    static constexpr double      DT_AU            = 0.005;        // 4× finer than WP run
    static constexpr int         N_STEPS          = 6000;         // total_time = 30 a.u.
    static constexpr int         WRITE_EVERY      = 8;            // ~same GIF cadence as WP run

    // ----- Inherited from Base_N162_L50_E1p5: L=50, N_ELECTRONS=162,
    //       SPACING_BOHR=1.0, EXTRA_STATES=20, xc=LDA, T=100 K -----------
};

}  // namespace jellium::config
