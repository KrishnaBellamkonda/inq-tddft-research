// shared/configs/E200_s0p53_N40.hpp — N=40 open-shell jellium, σ=0.53 Å, E=200 eV, +z
#pragma once
#include "base.hpp"
namespace jellium::config {
struct E200_s0p53_N40 : Base {
    // Open-shell variant. Originally N=40 at L=40; rescaled to N=135 at
    // L=60 to preserve the open-shell character (one electron above the
    // even-N closed shell at the same r_s). Uses its own checkpoint
    // gs_L60_cubic_N135.
    static constexpr int N_ELECTRONS = 135;
};
}  // namespace jellium::config
