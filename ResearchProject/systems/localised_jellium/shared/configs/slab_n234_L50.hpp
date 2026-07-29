// ============================================================================
// systems/localised_jellium/shared/configs/slab_n234_L50.hpp
//
// Source-of-truth Cfg for the localised jellium SLAB campaign (grill-locked
// 2026-06-21). A finite positive jellium background confined to a slab,
// injected as a static perturbation (inqkit::jellium::localised_background_
// perturbation); see docs/campaigns/localised_jellium/localised_jellium_campaign.md
// and docs/notes/localised-jellium-theory.md.
//
//   * Cell: 50 Bohr cubic periodic (INQ-centred, z ∈ [−25,+25]).
//   * Slab: full 50×50 face, 25 Bohr thick (half-width 12.5), centred at z=0.
//     V_inside = 50·50·25 = 62500 Bohr³.
//   * N = 234 electrons (even). n₀ = N/V = 3.7440e-3 a₀⁻³ → r_s = 3.996 (Na).
//     Exact neutrality (∫n₊ = N) makes the G=0 cancellation exact.
//   * Spacing 0.50 Bohr: Nyquist π/0.5 = 6.28 > WP k₀ = 2.71; resolves the
//     surface Friedel wavelength π/k_F ≈ 6.5 Bohr; ≥0.40 avoids the WP-init
//     deadlock. dx tightened for T2 grid-convergence only.
//   * Projectile (Phase 3/5): σ=0.5 Bohr, E=100 eV, k₀=2.71, launched +z at
//     z=−23 (4σ from the −z wall), run-up ≈10.5 Bohr to the slab face.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double ANG_TO_BOHR = 1.8897259886;
inline constexpr double HA_TO_EV    = 27.21138625;

inline constexpr double const_sqrt(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 30; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev(double e_ev) {
	return const_sqrt(2.0 * e_ev / HA_TO_EV);
}

struct SlabN234L50 {
	// Cell
	static constexpr double L_BOHR       = 50.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background (full x,y; confined along z = axis 2)
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.0;        // sharp Θ (soften if Gibbs)

	// Electronic structure
	static constexpr int    N_ELECTRONS    = 234;
	static constexpr double V_INSIDE_BOHR3 = L_BOHR * L_BOHR * (2.0 * SLAB_HALF_WIDTH);
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3; // 3.744e-3
	static constexpr int    EXTRA_STATES   = 20;
	static constexpr double TEMPERATURE_EV = 0.00862;      // ≈100 K  TODO: Can be decreases even further
	static constexpr double SCF_TOL_HA     = 1.0e-4; // TODO: Can be decreased to 1.0e-6
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (Phase 3/5; unused by GS)
	static constexpr double WP_SIGMA_BOHR = 0.5;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_K0         = k0_from_ev(WP_EKIN_EV);   // ≈2.71
	static constexpr double WP_CZ_BOHR    = -23.0;                    // 4σ from −z wall
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
