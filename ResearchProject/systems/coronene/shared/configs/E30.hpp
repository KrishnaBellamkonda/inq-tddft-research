// E = 30 eV (low-energy WP). Everything else inherited from Base.
#pragma once
#include "tsubonoya_2014_base.hpp"

namespace coronene::config {

struct E30 : Base {
    static constexpr double WP_EKIN_EV      = 30.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);
    static constexpr double WP_KZ           = -WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;
    // |k| changes -> recompute end-of-box arrival time.
    static constexpr int    N_STEPS = compute_n_steps(
        LZ_BOHR, WP_OFFSET_BOHR, WP_SIGMA_BOHR, WP_K0, DT_AU);
};

}  // namespace coronene::config
