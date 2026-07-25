// ============================================================================
// systems/localised_jellium/shared/configs/slab_n100_L35x35x85.hpp
//
// Source-of-truth Cfg for the HIGH-DENSITY classical S(v) benchmark campaign
// (docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md,
//  id classical-highdensity-sv). Grill-locked 2026-07-21.
//
// Denser localised jellium slab than the r_s=5.68 baseline: transverse box
// shrunk to 35x35 at fixed 25-Bohr thickness with N=100 to raise the density.
//   * n0 = 100 / (35*35*25) = 100/30625 = 3.2653e-3 a0^-3  ->  r_s = 4.183
//     (Na-like moderate-density metal; a real bump from the 82e r_s=5.667).
//   * Grid at dx=0.50: 70 x 70 x 170 = 833k pts. E_cut = (pi/dx)^2/2 = 19.74 Ha
//     = 537 eV single-particle Nyquist. dt=0.04 -> H*dt well below the 2.2 cliff.
//   * Boundary: periodicity(2) — x,y periodic (infinite slab), z OPEN/finite so
//     the moving Gaussian-charge projectile LEAVES the box (no wraparound) and,
//     with NO CAP, energy is conserved -> exact post-exit E_electronic plateau.
//   * Region layout (z, box [-42.5,+42.5]): slab [-12.5,+12.5] . vacuum each side
//     (~30 Bohr) for launch standoff + clean exit + wake room. Lock exact launch_z
//     per velocity at the pilot.
//   * Projectile = mass-1 electron Gaussian CHARGE (perturbation, not ghost UPF),
//     sigma_WP=0.5 -> sigma_pot=0.35355 (charge std == WP density std;
//     sigma-wp-convention). Velocity grid = high-v transit floor + up (pilot-set).
//   * dx=0.50 is PROVISIONAL: sigma_pot=0.354 and momentum transfers ~2v at the
//     fast end can approach k_Nyq=pi/0.5=6.28; cutoff_guard.py is a MANDATORY
//     per-velocity pre-launch gate (drop dx to ~0.35-0.40 for the fastest points
//     if it fails).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_L35 = 27.21138625;

struct SlabN100_L35x35x85 {
	// Cell (orthorhombic; z in [-42.5, +42.5]); periodicity(2) applied in run.cpp
	static constexpr double LX_BOHR      = 35.0;
	static constexpr double LY_BOHR      = 35.0;
	static constexpr double LZ_BOHR      = 85.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background (full x,y; confined along z) — 25 Bohr thick, denser
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;   // 25 Bohr thick (L_slab)
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;    // erfc-softened (GS-study H1)

	// Electronic structure — denser: N=100 over 30625 Bohr^3
	static constexpr int    N_ELECTRONS    = 100;
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 30625
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 3.2653e-3
	static constexpr int    EXTRA_STATES   = 24;       // ~50 occupied + headroom (T=100 K)
	static constexpr double TEMPERATURE_EV = 0.00862;  // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (mass-1 electron Gaussian charge; sigma matched to sigma_WP=0.5)
	static constexpr double WP_SIGMA_BOHR  = 0.5;                 // sigma_WP (label)
	static constexpr double PROJ_SIGMA_POT = 0.5 / 1.4142135624;  // sigma_pot = sigma_WP/sqrt(2) = 0.35355
	static constexpr double PROJ_MASS      = 1.0;                 // classical electron
	static constexpr double PROJ_CHARGE    = -1.0;                // electron
};

} // namespace localised_jellium::config
