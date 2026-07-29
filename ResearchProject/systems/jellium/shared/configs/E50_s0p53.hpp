// shared/configs/E50_s0p53.hpp — N=38, σ=0.53 Å, E=50 eV, +z (low-energy WP)
#pragma once
#include "base.hpp"
namespace jellium::config {
struct E50_s0p53 : Base {
    static constexpr double WP_EKIN_EV = 50.0;
    static constexpr double WP_K0      = k0_from_ev(WP_EKIN_EV);
    static constexpr double WP_KZ      = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;
    static constexpr int    N_STEPS    = 834;     // ½ k₀ ⇒ 2× steps for same range
    static constexpr double T2_AU      = DT_AU * N_STEPS;
};
}  // namespace jellium::config
