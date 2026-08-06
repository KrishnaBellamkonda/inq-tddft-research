// ============================================================================
// systems/localised_jellium/shared/configs/slab_n40_L35x35x85.hpp
//
// LOW-DENSITY partner of slab_n100_L35x35x85.hpp for the wrap-around slab
// KS-orbital stopping study (docs/plans/slab-ks-orbital-stopping-wrap.md).
// User decision 2026-07-31: two densities, same box, N the only variable.
//
// EVERY geometric and numerical parameter is identical to the N=100 config —
// same 35 x 35 x 85 cell, same 25-Bohr slab, same erfc edge, same dx, dt,
// smearing, SCF settings, periodicity(2). ONLY the electron count changes, so
// the density is the single independent variable and any difference between the
// two campaigns is attributable to it and to nothing else.
//
//   * n0 = 40 / (35*35*25) = 40/30625 = 1.3061e-3 a0^-3  ->  r_s = 5.675
//     Density ratio against the N=100 partner: 2.50x.
//   * Derived HEG scales: k_F = v_F = 0.3382, E_F = 1.56 eV,
//     hbar omega_p = 3.49 eV  ->  T_plasmon = 49.0 a.u.
//   * 20 occupied spatial states + 24 extra = 44 states (vs 74 at N=100), so
//     these runs are roughly 0.6x the cost of the N=100 ones.
//
// WHY THIS DENSITY. It is doubly anchored:
//   (a) it is the project's long-standing localised-jellium reference density
//       (slab_n82_L50x50x*.hpp: n0 = 1.312e-3, r_s = 5.665 — 0.2 % away), and
//   (b) it matches the LOW-density point of the bulk KS-orbital stopping study
//       (r_s = 5.702, docs/handovers/bulk-jellium-ks-stopping.md), which is the
//       bulk result these slab runs are being tested against.
// So the two-density slab pair spans the same lever the bulk pair did (2.50x
// here vs 2.92x there) with both endpoints tied to existing measurements.
//
// dx = 0.40 is the PRODUCTION grid here (not 0.50 as in the N=100 header's
// original classical campaign): E_cut = (pi/dx)^2/2 = 30.84 Ha = 839 eV, and at
// sigma_WP = 2 the packet's momentum width is sigma_p = 1/(sqrt2 sigma) = 0.354
// against k_Nyq = pi/0.4 = 7.85, so every moment is exact to <1e-4 %.
//
// BOUNDARY: periodicity(2) — x,y periodic (infinite slab), z open
// ELECTROSTATICALLY. Note that the wavefunction basis is a plain 3-D FFT and is
// periodic in ALL THREE directions regardless, so KS orbitals wrap in z. That is
// the mechanism the wrap-around study relies on; see the plan, section 1.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_L35_N40 = 27.21138625;

struct SlabN40_L35x35x85 {
	// Cell (orthorhombic; z in [-42.5, +42.5]); periodicity(2) applied in run.cpp
	static constexpr double LX_BOHR      = 35.0;
	static constexpr double LY_BOHR      = 35.0;
	static constexpr double LZ_BOHR      = 85.0;
	static constexpr double SPACING_BOHR = 0.40;   // production grid

	// Slab background — IDENTICAL to the N=100 partner
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;   // 25 Bohr thick (L_slab)
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;    // erfc-softened

	// Electronic structure — the ONLY difference: 40 electrons instead of 100
	static constexpr int    N_ELECTRONS    = 40;
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 30625
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.3061e-3
	static constexpr int    EXTRA_STATES   = 24;       // 20 occupied + headroom -> 44 states
	static constexpr double TEMPERATURE_EV = 0.00862;  // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile — sigma_WP = 2 for this study (sigma-wp-convention: sigma always
	// means the WAVEPACKET width; the classical charge std is the derived
	// sigma_pot = sigma_WP/sqrt2 = 1.41421, surfaced only inside the binaries).
	static constexpr double WP_SIGMA_BOHR  = 2.0;
	static constexpr double PROJ_SIGMA_POT = 2.0 / 1.4142135624;  // 1.41421
	static constexpr double PROJ_MASS      = 1.0;
	static constexpr double PROJ_CHARGE    = -1.0;
};

} // namespace localised_jellium::config
