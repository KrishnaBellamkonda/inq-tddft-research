// Larger box (Lz = 80 Bohr) + b = 6 Bohr (= 0.5 x base, near approach).
#pragma once
#include "cell_35x35x80.hpp"

namespace coronene::config {

struct b6_35x35x80 : cell_35x35x80 {
    static constexpr double WP_OFFSET_BOHR = 6.0;
    static constexpr double WP_CZ_BOHR     = +WP_OFFSET_BOHR;
    static constexpr int    N_STEPS = compute_n_steps(
        LZ_BOHR, WP_OFFSET_BOHR, WP_SIGMA_BOHR, WP_K0, DT_AU);
};

}  // namespace coronene::config
