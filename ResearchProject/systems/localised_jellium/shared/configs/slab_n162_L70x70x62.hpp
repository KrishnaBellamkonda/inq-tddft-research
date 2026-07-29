// ============================================================================
// systems/localised_jellium/shared/configs/slab_n162_L70x70x62.hpp
//
// Source-of-truth Cfg for the GENUINE 162-electron localised-jellium mass-pair
// WP campaign (mass_pair_n162): matched quantum wavepacket runs at projectile
// mass 1 (m1) and mass 2 (m2), sigma_WP=1, E=100 eV, CAP eta=-1.0.
// Plan: docs/plans/mass-pair-n162-sigma1-cap.md (user spec 2026-07-19).
//
// DESIGN (why these numbers):
//   * GENUINE N=162 (not density-matched). 162 is a magic number. The transverse
//     box is sized to hold exactly 162 electrons in the 25-Bohr slab at the
//     canonical localised-jellium density: Lx=Ly=70.4 -> n0 = 162/(70.4^2*25) =
//     1.3054e-3 a0^-3 -> r_s = 5.68 (0.3% from the reference 5.665). Neutral to
//     162 electrons exactly (n0 defined FROM N and the box, so no net cell charge).
//   * dx = 0.40 Bohr: the FINEST spacing allowed by the 24 GB A30 memory ceiling
//     (user directive). Grid 176 x 176 x 156 = 4.83 M pts. E_cut=(pi/dx)^2/2 =
//     30.8 Ha = 839 eV >> 100 eV. Aliasing gate PASSED for both masses
//     (p0+3sigma_p = 4.83 (m1) / 5.96 (m2) < k_Nyq=7.854).
//   * Slab z in [-12.5,12.5] (25 Bohr), erfc edge_width=1.0 (smoothening).
//   * Lz=62.4 (156 pts): 4-sigma geometry. WP launch z0=-16.5 (=slab_face-4sigma).
//     CAP region [+/-21.2, +/-31.2] (10 Bohr/side, eta=-1.0), inner face +/-21.2
//     => WP standoff 4.7 Bohr = 4.7 sigma from the CAP (>= 4 sigma). See plan.
//   * Projectile at E=100 eV: k0 = sqrt(2*m*E_Ha). m1: k0=2.711, v=2.711.
//     m2 (inverse_mass=0.5): k0=3.834, v=1.917. Set per-run in the WP dispatcher.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_N162 = 27.21138625;

struct SlabN162_L70x70x62 {
	// Cell (orthorhombic; z in [-31.2, +31.2])
	static constexpr double LX_BOHR      = 70.4;
	static constexpr double LY_BOHR      = 70.4;
	static constexpr double LZ_BOHR      = 62.4;
	static constexpr double SPACING_BOHR = 0.40;      // A30 memory ceiling (finest allowed)

	// Slab background (full x,y; confined along z)
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;  // 25 Bohr thick, z in [-12.5,12.5]
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;   // erfc smoothening (user spec)

	// Electronic structure — GENUINE 162 electrons; n0 defined FROM N and the box
	static constexpr int    N_ELECTRONS    = 162;
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 123904
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.3075e-3
	static constexpr int    EXTRA_STATES   = 18;      // larger DOS than N=82; pilot-confirmed
	static constexpr double TEMPERATURE_EV = 0.00862; // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (sigma_WP=1, E=100 eV; k0/inverse_mass set per-run in the dispatcher)
	static constexpr double WP_SIGMA_BOHR = 1.0;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_CZ_BOHR    = -16.5;    // 4*sigma standoff from slab face
	// m=1 default kinematics (m=2 overrides via env EM_K0/EM_INV_MASS):
	static constexpr double WP_K0         = 2.711;    // sqrt(2*1*100/27.2114)
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
