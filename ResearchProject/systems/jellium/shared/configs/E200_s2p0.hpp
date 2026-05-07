// shared/configs/E200_s2p0.hpp — N=38, σ=2.0 Å, E=200 eV, +z (wide WP)
#pragma once
#include "base.hpp"
namespace jellium::config {
struct E200_s2p0 : Base {
    static constexpr double WP_SIGMA_ANG  = 2.0;
    static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;     // ≈ 3.779
    // Wider WP — larger envelope means we want extra margin from +L/2.
    // Inherits WP_CZ_BOHR = 0 (box centre) from Base. Single-pass:
    // k0 * dt * N_STEPS ≈ 3.83 * 0.020 * 280 = 21.5 bohr; remaining
    // L/2 - 21.5 - 5σ ≈ 30 - 21.5 - 19 = -10.5 bohr ⟹ half the WP envelope
    // wraps even at t = T_final, so reduce N_STEPS so the WP front sits at
    // ~L/2 - 3σ ≈ 18.6 bohr from origin (~242 steps).
    static constexpr int    N_STEPS       = 240;
    static constexpr double T2_AU         = DT_AU * N_STEPS;
};
}  // namespace jellium::config
