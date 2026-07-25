// ============================================================================
// systems/localised_jellium/shared/configs/slab_n102_L25x25x140_w0p5.hpp
//
// Source-of-truth Cfg for the CAP energy-plateau diagnostic campaign
// (docs/campaigns/localised_jellium/wp_cap_energy_plateau.md). Localised jellium
// SLAB (fills the periodic x,y face; localised in z) with a WP projectile,
// run WITHOUT and WITH a two-sided CAP to measure how much energy the CAP drains
// from the plateauing energy_total.
//
//   * Cell: 25 × 25 × 140 Bohr orthorhombic periodic (INQ-centred, z ∈ [−70,+70]).
//   * Slab: full 25×25 face, 25 Bohr thick (half-width 12.5), centred z=0,
//     faces ±12.5. erfc edge softening w = 0.5 Bohr.
//     V_inside = 25·25·25 = 15625 Bohr³.
//   * N = 102 electrons. n₀ = 102/15625 = 6.528e-3 a₀⁻³ → r_s ≈ 3.32,
//     E_F ≈ 4.5 eV, ħω_p ≈ 7.8 eV.
//   * Grid spacing h = 0.5 Bohr (must pass cutoff_guard.py for σ_WP=1, E=100 eV).
//   * Projectile: σ_WP=1 Bohr, E=100 eV, k₀≈2.711, mass 1 (electron), launched
//     +z at z=−20.5 (8 Bohr from the −12.5 slab face).
//   * CAP (run 2 only): TWO-SIDED sin² absorber, η=−0.7 Ha, 10 Bohr/side at the
//     far ends z∈[±60,±70] (fractional mid ±0.4643, width 0.07143). Sits in
//     vacuum → no bath over-drain. Functional ONLY on inq-study.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_N102 = 27.211386245988;

inline constexpr double const_sqrt_n102(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 40; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_n102(double e_ev) {
	return const_sqrt_n102(2.0 * e_ev / HA_TO_EV_N102);
}

struct SlabN102_L25x25x140_w0p5 {
	// Cell (orthorhombic: x,y = 25 set in-plane density; z = 140 long vacuum)
	static constexpr double LX_BOHR      = 25.0;
	static constexpr double LY_BOHR      = 25.0;
	static constexpr double LZ_BOHR      = 140.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background (full x,y; confined along z = axis 2)
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.5;        // w = 0.5 erfc softening

	// Electronic structure
	static constexpr int    N_ELECTRONS    = 102;          // even
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 15625
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 6.528e-3
	static constexpr int    EXTRA_STATES   = 24;
	static constexpr double TEMPERATURE_EV = 0.00862;      // ≈100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (σ_WP=1, 100 eV, mass 1)
	static constexpr double WP_SIGMA_BOHR = 1.0;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_K0         = k0_from_ev_n102(WP_EKIN_EV);   // ≈2.711
	static constexpr double WP_CZ_BOHR    = -20.5;                          // 8 Bohr from −12.5 face
	static constexpr double WP_KZ         = +WP_K0;

	// Two-sided CAP (run 2 only; fractions of LZ, physical z = frac·LZ from centre)
	static constexpr double CAP_ETA_HA    = -0.7;
	static constexpr double CAP_L_BOHR    = 10.0;                            // per side
	static constexpr double CAP_WIDTH_FRAC= CAP_L_BOHR / LZ_BOHR;            // 0.07143
	static constexpr double CAP_MID_FRAC  = (LZ_BOHR/2.0 - CAP_L_BOHR/2.0) / LZ_BOHR; // 0.4643
};

} // namespace localised_jellium::config
