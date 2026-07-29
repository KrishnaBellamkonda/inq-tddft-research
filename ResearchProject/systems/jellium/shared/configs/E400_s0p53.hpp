// shared/configs/E400_s0p53.hpp — N=38, σ=0.53 Å, E=400 eV, +z (high-energy WP)
#pragma once
#include "base.hpp"
namespace jellium::config {
struct E400_s0p53 : Base {
    static constexpr double WP_EKIN_EV = 400.0;
    static constexpr double WP_K0      = k0_from_ev(WP_EKIN_EV);
    static constexpr double WP_KZ      = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;
    static constexpr int    N_STEPS    = 295;
    static constexpr double T2_AU      = DT_AU * N_STEPS;
};
}  // namespace jellium::config
