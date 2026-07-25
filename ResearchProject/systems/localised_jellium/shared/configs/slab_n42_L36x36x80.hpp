// ============================================================================
// systems/localised_jellium/shared/configs/slab_n42_L36x36x80.hpp
//
// Source-of-truth Cfg for the effective-mass muon-WP run RE-PLANNED for <=12 h
// (N3/A, 2026-07-08). Same r_s density and 25 Bohr slab thickness as the
// slab_n82_L50x50x90 reference, but the TRANSVERSE box is shrunk 50 -> 36 and
// z 90 -> 80 to bring the wall time under 12 h on one GPU. Holding the density
// n0 fixed with an EVEN electron count gives N = 42 (n0 = 42/32400 = 1.296e-3,
// r_s = 5.69) at the coarser dx = 0.40 grid.
//
//   * Cell: 36 x 36 x 80 Bohr orthorhombic periodic (INQ-centred, z in [-40,+40]).
//   * Slab: full 36x36 face, 25 Bohr thick (half-width 12.5), centred at z=0.
//     V_inside = 36*36*25 = 32400 Bohr^3.
//   * N = 42 electrons (even). n0 = 42/32400 = 1.296e-3 a0^-3 -> r_s = 5.69.
//   * dx = 0.40 Bohr (k_max = pi/dx = 7.854; grid 90x90x200 = 1.62 M pts).
//   * Region layout (z): slab [-12.5,12.5] . free . CAP near |z|=32 (mid 0.40*Lz),
//     8 Bohr/side (width 0.10*Lz), inner face ~|z|=28.
//   * Projectile (quantum WP, set at run time via env): sigma_WP = 2.0,
//     k0 = 6.7933 (= k_max - 3 sigma_p), m_eff = 2.506 m_e, v = 2.711 a.u.
//     (= 100 eV electron velocity -> same S(v)), E = 251 eV, launch z = -16.389.
//
// Transverse-box caveat: 36<50 risks periodic wake-wrap at r_s~5.7 (weak
// screening); the momentum-peak S(v0) readout is largely immune. See
// docs/plans/muon-effmass-12h-run.md.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_42 = 27.21138625;

struct SlabN42_L36x36x80 {
	// Cell (orthorhombic: x,y = 36; z = 80)
	static constexpr double LX_BOHR      = 36.0;
	static constexpr double LY_BOHR      = 36.0;
	static constexpr double LZ_BOHR      = 80.0;
	static constexpr double SPACING_BOHR = 0.40;

	// Slab background (full x,y; confined along z = axis 2) — density held vs 50-box
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick (unchanged)
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.0;        // sharp Θ (soften if Gibbs)

	// Electronic structure — N chosen even to hold n0 at the shrunk transverse box
	static constexpr int    N_ELECTRONS    = 42;           // even; 42/32400 = 1.296e-3
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 32400
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.296e-3
	static constexpr int    EXTRA_STATES   = 10;           // 21 occ + 10; lighter than the 61-state run
	static constexpr double TEMPERATURE_EV = 0.00862;      // ≈100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile reference values (set at run time via env; here for provenance)
	static constexpr double WP_SIGMA_BOHR = 2.0;
	static constexpr double WP_K0         = 6.7933;        // k_max - 3 sigma_p
	static constexpr double WP_INV_MASS   = 0.39907;       // 1 / 2.506
	static constexpr double WP_VELOCITY   = 2.711;         // = 100 eV electron
	static constexpr double WP_CZ_BOHR    = -16.389;       // 2.75 sigma_rho0 before slab face
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
