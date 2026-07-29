// ============================================================================
// systems/localised_jellium/shared/configs/slab_n120_L60x60x62.hpp
//
// 120-electron localised-jellium slab for the mass-pair WP campaign — the
// memory-relieved successor to slab_n162_L70x70x62 (user directive 2026-07-19:
// "let's try 120 electrons or so"). The genuine N=162 system (99 KS states,
// 176x176x156) saturated the 24 GB A30 during ETRS propagation (min GPU free =
// 0 MB, intermittent step-0 CUDA illegal-access). N=120 -> 78 states / 152x152x156
// grid ~= 59% of the memory (~14 GB used, ~10 GB free) => reliable propagation.
//
// EVERYTHING in the z-direction is IDENTICAL to the n162 config: slab z in
// [-12.5,12.5], erfc edge_width=1.0, Lz=62.4, WP launch -16.5 (4 sigma), CAP
// region +/-[21.2,31.2] (10 Bohr/side, eta=-1.0), dx=0.40. Only the transverse
// box shrinks to hold 120 e at r_s~=5.69.
//   * Lx=Ly=60.8 (152 pts) -> n0=120/(60.8^2*25)=1.2985e-3 -> r_s=5.686.
//   * Grid 152x152x156 = 3.60 M pts. E_cut=(pi/0.4)^2/2=839 eV >> 100 eV.
//   * Aliasing gate PASSED (same sigma_WP/energy as n162; dx unchanged).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_N120 = 27.21138625;

struct SlabN120_L60x60x62 {
	// Cell (orthorhombic; z in [-31.2, +31.2])
	static constexpr double LX_BOHR      = 60.8;
	static constexpr double LY_BOHR      = 60.8;
	static constexpr double LZ_BOHR      = 62.4;
	static constexpr double SPACING_BOHR = 0.40;

	// Slab background (identical to n162)
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;  // 25 Bohr thick, z in [-12.5,12.5]
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;   // erfc smoothening

	// Electronic structure — 120 electrons; n0 defined FROM N and the box
	static constexpr int    N_ELECTRONS    = 120;
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 92416
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.2985e-3
	static constexpr int    EXTRA_STATES   = 18;      // 60 occ + 18 -> 78 states (headroom)
	static constexpr double TEMPERATURE_EV = 0.00862; // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (sigma_WP=1, E=100 eV; identical z-kinematics to n162)
	static constexpr double WP_SIGMA_BOHR = 1.0;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_CZ_BOHR    = -16.5;    // 4*sigma standoff from slab face
	static constexpr double WP_K0         = 2.711;    // m=1: sqrt(2*1*100/27.2114)
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
