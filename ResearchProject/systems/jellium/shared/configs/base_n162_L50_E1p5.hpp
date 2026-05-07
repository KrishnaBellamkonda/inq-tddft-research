// ============================================================================
// shared/configs/base_n162_L50_E1p5.hpp
//
// True closed-shell variant: N = 162 fills the |G|^2 = 6 shell completely
// (24 spatial states added on top of the 57 from |G|^2 <= 5 ⇒ 81 spatial,
// 162 paired electrons; see docs/sources/free-electron-gas-magic-numbers.md).
// At L=50 this gives density 1.296e-3 e/bohr^3, r_s ≈ 5.69 bohr (slightly
// denser than the N=138 case at this box, still lithium-like).
//
// All other parameters inherited from Base_N138_L50_E1p5: L=50, dx=1.0,
// EXTRA_STATES=20, SCF_TOL_HA=1e-6, WP_EKIN_EV=1.5, WP_SIGMA_BOHR=5.0,
// dt=0.020, N_STEPS=1500.
//
// The N=138 run on the same box had partially-filled |G|^2=6 (only 12 of
// 24 spatial states), giving the smeared 0.95-1.0 occupations seen in
// occupations_vs_time.csv. With N=162 every shell is fully filled, so
// f_i should be exactly 2.0 for the lowest 81 states and 0 for everything
// above. This is the configuration intended by the project's
// "uniform-density jellium" rationale.
// ============================================================================
#pragma once

#include "base_n138_L50_E1p5.hpp"

namespace jellium::config {

struct Base_N162_L50_E1p5 : Base_N138_L50_E1p5 {
    static constexpr int N_ELECTRONS = 162;
};

}  // namespace jellium::config
