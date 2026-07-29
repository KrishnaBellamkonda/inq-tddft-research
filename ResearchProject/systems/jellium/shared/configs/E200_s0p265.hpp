// shared/configs/E200_s0p265.hpp — N=38, σ=0.265 Å, E=200 eV, +z (narrow WP)
#pragma once
#include "base.hpp"
namespace jellium::config {
struct E200_s0p265 : Base {
    static constexpr double WP_SIGMA_ANG  = 0.265;
    static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;     // ≈ 0.5008
    // Inherits WP_CZ_BOHR = 0 (box centre) from Base.
    // Same single-pass reasoning as Base: travel ≈ k0*dt*N_STEPS ≈ 24.5
    // bohr in 320 steps; the narrower WP just makes the front sharper.
    static constexpr int    N_STEPS       = 320;
    static constexpr double T2_AU         = DT_AU * N_STEPS;
};
}  // namespace jellium::config
