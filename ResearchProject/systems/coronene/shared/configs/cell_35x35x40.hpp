// Smaller box: Lz = 40 Bohr (= 2/3 x base). b stays at base.
#pragma once
#include "tsubonoya_2014_base.hpp"

namespace coronene::config {

struct cell_35x35x40 : Base {
    static constexpr double LZ_BOHR = 40.0;
    static constexpr int    N_STEPS = compute_n_steps(
        LZ_BOHR, WP_OFFSET_BOHR, WP_SIGMA_BOHR, WP_K0, DT_AU);
};

}  // namespace coronene::config
