// ============================================================================
// shared/configs/cc_bond_35x35x60.hpp
//
// Off-center WP impact: beam aimed at a Type B (radial spoke) C-C bond
// midpoint instead of the molecular centre. All other parameters identical
// to Base (Tsubonoya 2014 paper parameters) for a clean single-variable
// comparison with run_base.
//
// Bond: C₁ = (1.421, 0.0, 0.0) Å  →  C₂ = (2.842, 0.0, 0.0) Å
//        (shared/geometry/coronene.xyz lines 3–4)
// Midpoint in Å:  (2.1315, 0, 0)
// Midpoint in Bohr: 2.1315 × 1.88973 = 4.028 Bohr
//
// The WP is shifted +4.028 Bohr in x from the molecular centre, with the
// same z-offset (12 Bohr above the molecular plane). The 3σ tail at
// x = 4.028 + 3×1.0015 = 7.033 Bohr is well inside the cell boundary
// at Lx/2 = 17.39 Bohr. N_STEPS is unchanged: the z-travel time is the
// same regardless of the x-offset.
// ============================================================================
#pragma once

#include "tsubonoya_2014_base.hpp"

namespace coronene::config {

struct CC_bond_35x35x60 : Base {
    static constexpr double WP_CX_BOHR = 2.1315 * ANG_TO_BOHR;  // 4.028 Bohr
};

}  // namespace coronene::config
