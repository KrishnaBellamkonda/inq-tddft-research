// ============================================================================
// systems/localised_jellium/shared/configs/slab_n82_L50x50x101.hpp
//
// Source-of-truth Cfg for the `wide-wavepacket-lowspread` campaign (wide,
// near-rigid WP, sigma_WP=3.5, matched classical). IDENTICAL slab/density/
// electrons to slab_n82_L50x50x90.hpp; differs ONLY in:
//   * z extent 90 -> 101 Bohr, so the equidistant launch sits EXACTLY 4*sigma
//     (=14 Bohr at sigma_WP=3.5) from BOTH the slab face and the CAP inner face.
//   * spacing 0.50 -> 0.40 Bohr, set by the cutoff/aliasing guard:
//     k_max=pi/dx must exceed k0 + 4*sigma_p, sigma_p=1/(sqrt2*sigma_WP)=0.202.
//     dx=0.40 -> k_Nyq=7.85, E_cut=839 eV -> >=6 sigma_p margin to E=600 eV
//     (phase-5's dx=0.50 fails at E>=500). Clear of the dx=0.30 WP-init deadlock.
//
//   * Cell: 50 x 50 x 101 Bohr orthorhombic periodic (z in [-50.5, +50.5]).
//   * Slab: full 50x50 face, 25 Bohr thick (half-width 12.5), centred at z=0.
//     V_inside = 50*50*25 = 62500 Bohr^3 (UNCHANGED). N=82, n0=1.312e-3,
//     r_s~5.665 (UNCHANGED -- slab interior is box-independent).
//   * Region layout (z): slab [-12.5,12.5] . free [+/-12.5,+/-40.5] .
//     CAP [+/-40.5,+/-50.5] (inner faces +/-40.5).  CAP = SAME as phase-5:
//     two-sided sin2, eta=-0.7 Ha, 10 Bohr/side.
//   * Projectile (operating point): sigma_WP=3.5 Bohr, E=300 eV, k0~4.696,
//     launched +z at z=-26.5 (equidistant: 14.0 Bohr = 4 sigma from both the
//     -12.5 slab face and the -40.5 CAP inner face). The matched classical UPF is
//     electron_gaussian_wpsigma3p5.upf (sigma_pot=2.475). Phase-1 sweeps E via
//     env (LJ_K0 / LJ_EKIN); k0 and tau are real-time-only -> GS reuse OK.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_101 = 27.21138625;

inline constexpr double const_sqrt_101(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 30; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_101(double e_ev) {
	return const_sqrt_101(2.0 * e_ev / HA_TO_EV_101);
}

struct SlabN82_L50x50x101 {
	// Cell (orthorhombic: x,y = 50 preserve in-plane density; z = 101)
	static constexpr double LX_BOHR      = 50.0;
	static constexpr double LY_BOHR      = 50.0;
	static constexpr double LZ_BOHR      = 101.0;
	static constexpr double SPACING_BOHR = 0.40;   // cutoff-guard refined (see header)

	// Slab background (full x,y; confined along z = axis 2) — UNCHANGED vs 90-box
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;        // erfc-softened (>= dx=0.40) — kills the
	                                                       // sharp-edge Gibbs aliasing that appears
	                                                       // because |z|=12.5 is NOT a grid node at
	                                                       // dx=0.40 (it is at dx=0.50). Per GS-study H1.

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

	// Projectile (wide WP sigma_WP=3.5, operating point E=300 eV) — 4-sigma launch
	static constexpr double WP_SIGMA_BOHR = 3.5;
	static constexpr double WP_EKIN_EV    = 300.0;
	static constexpr double WP_K0         = k0_from_ev_101(WP_EKIN_EV);   // ~4.696
	static constexpr double WP_CZ_BOHR    = -26.5;                        // 4-sigma equidistant launch
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
