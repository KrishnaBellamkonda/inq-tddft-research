// ============================================================================
// broadening_35x35x80.hpp — long pre-collision flight for Gaussian wave-packet
// spreading, with the centroid kept INSIDE the box (no boundary crossing →
// no periodic interference) for the whole run.
//
// Purpose: primary deliverable is the PRE-COLLISION Gaussian-broadening plot
// (sigma_z(t) vs the analytic free-spread law), plus momentum redistribution
// before/after collision. Max temporal resolution: WRITE_EVERY = 1.
//
// Launch: z = +30 Bohr, moving -z toward the molecule at z = 0. The far
// boundary is -40 Bohr. Pre-collision flight = 30 Bohr (t ~ 7.8 a.u.); sigma
// grows ~1 -> ~4 Bohr by the time the centroid reaches the molecule.
//
// Centroid-safe stop: end when the centroid reaches z = -(Lz/2 - margin), so
// it NEVER reaches the boundary. This OVERRIDES compute_n_steps() (which would
// run until centroid+sigma cleared the far face, i.e. PAST the boundary).
// NOTE: the 3*sigma envelope still wraps in the late post-collision phase, so
// post-collision analysis must be time-windowed; the pre-collision segment is
// fully interference-free.
// ============================================================================
#pragma once
#include "cell_35x35x80.hpp"

namespace coronene::config {

struct broadening_35x35x80 : cell_35x35x80 {
    // Long pre-collision flight (centre impact: WP_CX = WP_CY = 0 from Base).
    static constexpr double WP_OFFSET_BOHR = 30.0;
    static constexpr double WP_CZ_BOHR     = +WP_OFFSET_BOHR;

    // Max resolution.
    static constexpr int    WRITE_EVERY       = 1;
    static constexpr int    SCREEN_SNAP_EVERY = 10;

    // Centroid-safe N_STEPS: t_end = (offset + Lz/2 - margin) / k0.
    static constexpr double SAFE_MARGIN_BOHR = 5.0;
    static constexpr int    N_STEPS =
        static_cast<int>(((WP_OFFSET_BOHR + 0.5 * LZ_BOHR - SAFE_MARGIN_BOHR)
                          / WP_K0) / DT_AU + 0.5);
};

}  // namespace coronene::config
