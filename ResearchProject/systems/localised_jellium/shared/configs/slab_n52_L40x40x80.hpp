// ============================================================================
// systems/localised_jellium/shared/configs/slab_n52_L40x40x80.hpp
//
// Source-of-truth Cfg for the effmass_sigma1 LEAN re-run (concentrated
// sigma_WP=1 chirped-focus WP + matched classical twin, m=2.10).
// Plan: docs/plans/effmass-sigma1-lean-rerun.md.
//
// Derived from slab_n82_L50x50x101.hpp by the 2026-07-09 user directive:
// keep dx=0.333 (chirp Nyquist requirement), shrink x=y to 40 Bohr and z to
// 80 Bohr, and HOLD THE SLAB DENSITY fixed (N scales with transverse area):
//   * n0 = 1.312e-3 (identical to the sigma=2 run) over 40*40*25 = 40000 Bohr^3
//     -> N = 52.5 -> 52 (even). r_s = 5.679 (0.2% from the 82-electron 5.667).
//   * Grid at dx=0.33333: 120 x 120 x 240 = 3.46 M pts (FFT-friendly).
//     E_cut = (pi/dx)^2/2 = 44.4 Ha; dt=0.04 -> H*dt = 1.78 < 2.2 cliff.
//   * Region layout (z, box [-40,+40]): slab [-12.5,12.5] . free . CAP
//     [+/-25,+/-40] (15 Bohr/side, reflectivity-tuned width, eta=-1.0).
//     Launch z0=-16.5 (4*sigma standoff). Packet z-budget verified:
//     sigma_rho_z(exit)=3.2 -> 3sigma front +22 < 25; back tail -20 > -25.
//   * Aliasing: k0 + 3 sigma_p = 5.693 + 2.121 = 7.81 = 0.83 * k_Nyq(9.42). OK
//     (the vacuum-validated chirp regime; dx=0.40 FAILS this -> under-focus).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_L40 = 27.21138625;

struct SlabN52_L40x40x80 {
	// Cell (orthorhombic; z in [-40, +40])
	static constexpr double LX_BOHR      = 40.0;
	static constexpr double LY_BOHR      = 40.0;
	static constexpr double LZ_BOHR      = 80.0;
	static constexpr double SPACING_BOHR = 0.33333;  // chirp Nyquist requirement

	// Slab background (full x,y; confined along z) — same density/thickness as
	// the sigma=2 run's slab
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;   // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;    // erfc-softened (per GS-study H1)

	// Electronic structure — density-matched: n0 held at the sigma=2 value
	static constexpr int    N_ELECTRONS    = 52;       // 1.312e-3 * 40000 = 52.5 -> even
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 40000
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.300e-3
	static constexpr int    EXTRA_STATES   = 10;       // lean (T=100K needs ~10; 12h-run precedent)
	static constexpr double TEMPERATURE_EV = 0.00862;  // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (concentrated sigma_WP=1, chirped focus; m=2.10, v=2.711)
	static constexpr double WP_SIGMA_BOHR = 1.0;
	static constexpr double WP_K0         = 5.693;     // m*v = 2.10*2.711
	static constexpr double WP_CZ_BOHR    = -16.5;     // 4*sigma standoff from slab face
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
