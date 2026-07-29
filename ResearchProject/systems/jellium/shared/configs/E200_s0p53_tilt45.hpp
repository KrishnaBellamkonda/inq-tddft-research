// shared/configs/E200_s0p53_tilt45.hpp — N=38, σ=0.53 Å, E=200 eV, 45° xz-plane tilt
#pragma once
#include "base.hpp"
namespace jellium::config {
struct E200_s0p53_tilt45 : Base {
    // Same |k| as Base (E=200 eV); 45° tilt in xz-plane.
    // k_x = k_z = k₀ / √2.
    static constexpr double WP_KX = +Base::WP_K0 / 1.41421356237309504880;
    static constexpr double WP_KY = 0.0;
    static constexpr double WP_KZ = +Base::WP_K0 / 1.41421356237309504880;
    static constexpr double WP_KZ_MAGNITUDE = Base::WP_K0;
    static constexpr int    N_STEPS = 350;
    static constexpr double T2_AU   = DT_AU * N_STEPS;
};
}  // namespace jellium::config
