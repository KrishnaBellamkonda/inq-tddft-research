// Larger box (Lz = 80 Bohr) + b = 18 Bohr (= 1.5 x base).
#pragma once
#include "cell_35x35x80.hpp"

namespace coronene::config {

struct b18_35x35x80 : cell_35x35x80 {
    static constexpr double WP_OFFSET_BOHR = 18.0;
    static constexpr double WP_CZ_BOHR     = +WP_OFFSET_BOHR;
    static constexpr int    N_STEPS = compute_n_steps(
        LZ_BOHR, WP_OFFSET_BOHR, WP_SIGMA_BOHR, WP_K0, DT_AU);
};

}  // namespace coronene::config
