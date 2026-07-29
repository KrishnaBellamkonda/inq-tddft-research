// ============================================================================
// systems/localised_jellium/shared/configs/slab_n82_L50x50x120.hpp
//
// Source-of-truth Cfg for the localised-jellium GS parameter-study campaign,
// run-set 2 / H0 (base WP-vs-classical E_total(0) gap, then E_total(r) sweep).
//
// IDENTICAL slab/density/electrons to slab_n82_L50x50x90.hpp — ONLY the z extent
// grows 90 -> 120 Bohr so the WP/classical projectile can be placed up to
// r = 40 Bohr from the near slab face (z = -52.5) inside the box, and so the
// periodicity-2 (open-z) Rozzi cutoff rc = L_z = 120 comfortably exceeds the
// ~67 Bohr slab+WP z-extent.
//
//   * Cell: 50 x 50 x 120 Bohr orthorhombic (INQ-centred, z in [-60,+60]).
//   * Slab: full 50x50 face, 25 Bohr thick (half-width 12.5), centred at z=0.
//     V_inside = 50*50*25 = 62500 Bohr^3 (UNCHANGED by the z-extension).
//   * N = 82 electrons. n0 = 82/62500 = 1.312e-3 a0^-3 -> r_s ~ 5.665 (UNCHANGED).
//   * Edge: sharp Theta (edge_width 0) for H0 (w-sweep in H1 sets the clean w).
//   * Projectile: sigma_WP = 0.5 Bohr. For H0 the WP is STATIONARY (k0 = 0);
//     launch z is passed at run time (LJ_LAUNCH_Z): r=4 -> z=-16.5, r=40 -> z=-52.5
//     (r measured from the near slab face at z = -12.5).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_120 = 27.21138625;

inline constexpr double const_sqrt_120(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 30; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_120(double e_ev) {
	return const_sqrt_120(2.0 * e_ev / HA_TO_EV_120);
}

struct SlabN82_L50x50x120 {
	// Cell (orthorhombic: x,y = 50 preserve in-plane density; z = 120 vacuum)
	static constexpr double LX_BOHR      = 50.0;
	static constexpr double LY_BOHR      = 50.0;
	static constexpr double LZ_BOHR      = 120.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background (full x,y; confined along z = axis 2) — UNCHANGED vs 90-box
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.0;        // sharp Theta (H1 sets clean w)

	// Electronic structure — UNCHANGED vs 90-box (slab volume identical)
	static constexpr int    N_ELECTRONS    = 82;           // even
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 62500
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.312e-3
	static constexpr int    EXTRA_STATES   = 20;
	static constexpr double TEMPERATURE_EV = 0.00862;      // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (sigma_WP = 0.5). H0 uses k0 = 0 (stationary); launch z via env.
	static constexpr double WP_SIGMA_BOHR = 0.5;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_K0         = k0_from_ev_120(WP_EKIN_EV);   // ~2.71 (control only)
};

} // namespace localised_jellium::config
