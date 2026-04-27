// sigma = 3 x base (wide WP). Everything else inherited from Base.
#pragma once
#include "tsubonoya_2014_base.hpp"

namespace coronene::config {

struct s3 : Base {
    static constexpr double WP_SIGMA_BOHR = Base::WP_SIGMA_BOHR * 3.0;
    static constexpr int    N_STEPS = compute_n_steps(
        LZ_BOHR, WP_OFFSET_BOHR, WP_SIGMA_BOHR, WP_K0, DT_AU);
};

}  // namespace coronene::config
