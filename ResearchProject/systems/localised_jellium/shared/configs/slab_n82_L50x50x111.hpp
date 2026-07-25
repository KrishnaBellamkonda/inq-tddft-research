// ============================================================================
// systems/localised_jellium/shared/configs/slab_n82_L50x50x111.hpp
//
// Cfg for the `wide-wavepacket-lowspread` campaign, ENLARGED-BOX / open-z variant
// (2026-07-03 grilling decisions). Derived from slab_n82_L50x50x101.hpp; differs
// ONLY in the z extent and the boundary/CAP treatment:
//
//   * z extent 101 -> 111 Bohr (z in [-55.5, +55.5]). +10 Bohr gives the fast wide
//     packet more CAP runway so it is fully absorbed BEFORE it can wrap around.
//   * BOUNDARY: periodicity 2 (open-z: periodic in x,y, non-periodic in z) — set in
//     the run.cpp cell (`.periodicity(2)`), an ELECTROSTATICS choice (2D-truncated
//     Coulomb, no z-image interaction; poisson.hpp:190).
//     CAVEAT (user-accepted, to debug later): the localised-jellium GS study found
//     open-z gives a net-charge G=0 monopole self-energy (0.5*rc^2, rc=L_z;
//     poisson.hpp:49) that biases E_total by ~L_z^2*Q^2. It is a CONSTANT while Q is
//     fixed but STEPS once the CAP drains the WP (Q:83->82), so the long-time energy
//     PLATEAU sits below the true deposited energy until corrected. Q(t) is logged
//     (electron_number.csv) so the monopole is recoverable post-hoc.
//   * CAP (in run.cpp, NOT this header): two-sided sin^2, eta = -1.0 Ha (was -0.7),
//     14 Bohr/side (was 10), region [+/-41.5, +/-55.5], inner faces +/-41.5.
//
//   * Cell: 50 x 50 x 111 Bohr orthorhombic, periodicity 2 (open-z).
//   * Slab: full 50x50 face, 25 Bohr thick (half-width 12.5), centred z=0.
//     V_inside = 50*50*25 = 62500 Bohr^3 (UNCHANGED). N=82, n0=1.312e-3,
//     r_s~5.665 (UNCHANGED -- slab interior is box-independent).
//   * Region layout (z): slab [-12.5,12.5] . free [+/-12.5,+/-41.5] .
//     CAP [+/-41.5,+/-55.5] (inner faces +/-41.5).
//   * Projectile: sigma_WP=3.5 Bohr, E=300 eV, k0~4.696, launched +z at z=-26.5
//     (14.0 Bohr = 4 sigma from the -12.5 slab face; 15.0 Bohr to the -41.5 CAP
//     inner face -- WP tail clears the CAP at t=0).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_111 = 27.21138625;

inline constexpr double const_sqrt_111(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 30; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_111(double e_ev) {
	return const_sqrt_111(2.0 * e_ev / HA_TO_EV_111);
}

struct SlabN82_L50x50x111 {
	// Cell (orthorhombic: x,y = 50 preserve in-plane density; z = 111, open-z)
	static constexpr double LX_BOHR      = 50.0;
	static constexpr double LY_BOHR      = 50.0;
	static constexpr double LZ_BOHR      = 111.0;   // 101 -> 111 (+10 Bohr CAP runway)
	static constexpr double SPACING_BOHR = 0.40;

	// Slab background (full x,y; confined along z = axis 2) — UNCHANGED vs 101-box
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;        // erfc-softened (>= dx=0.40)

	// Electronic structure — UNCHANGED vs 101-box (slab volume identical)
	static constexpr int    N_ELECTRONS    = 82;           // even
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 62500
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.312e-3
	static constexpr int    EXTRA_STATES   = 20;
	static constexpr double TEMPERATURE_EV = 0.00862;      // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (wide WP sigma_WP=3.5, operating point E=300 eV) — 4-sigma launch
	static constexpr double WP_SIGMA_BOHR = 3.5;
	static constexpr double WP_EKIN_EV    = 300.0;
	static constexpr double WP_K0         = k0_from_ev_111(WP_EKIN_EV);   // ~4.696
	static constexpr double WP_CZ_BOHR    = -26.5;                        // 4-sigma launch
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
