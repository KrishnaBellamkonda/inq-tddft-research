// ============================================================================
// systems/localised_jellium/shared/configs/slab_n82_L50x50x70.hpp
//
// Source-of-truth Cfg for the `quantum-stopping-power` campaign, Phase 1
// (docs/campaigns/jellium_wp_stopping/quantum-stopping-power.md, locked 2026-06-24).
// A localised positive jellium slab at the n162-reference density (r_s≈5.67),
// in a NON-cubic 50×50×70 Bohr box (z extended to give CAP + far-launch room).
//
//   * Cell: 50 × 50 × 70 Bohr orthorhombic periodic (INQ-centred, z ∈ [−35,+35]).
//   * Slab: full 50×50 face, 25 Bohr thick (half-width 12.5), centred at z=0.
//     V_inside = 50·50·25 = 62500 Bohr³ (UNCHANGED by the z-extension).
//   * N = 82 electrons (even). n₀ = 82/62500 = 1.312e-3 a₀⁻³ → r_s ≈ 5.665,
//     matching the long-standing n162-in-50³ jellium (r_s=5.69) for S(v)
//     comparability. Exact neutrality (∫n₊ = N) → exact G=0 cancellation.
//   * Region layout (z): slab [−12.5,12.5] · free [±12.5,±25] · CAP [±25,±35].
//   * Projectile (Phase 1 SIE / Phase 2): σ=0.5 Bohr, E=100 eV, k₀≈2.71,
//     launched +z at z=−32 (far; ≈19.5 Bohr from the slab face, 3 Bohr from the
//     −35 wall). Phase 1 runs CAP OFF (only the t=0 energy is needed).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_82 = 27.21138625;

inline constexpr double const_sqrt_82(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 30; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_82(double e_ev) {
	return const_sqrt_82(2.0 * e_ev / HA_TO_EV_82);
}

struct SlabN82_L50x50x70 {
	// Cell (orthorhombic: x,y = 50 preserve in-plane density; z = 70 vacuum)
	static constexpr double LX_BOHR      = 50.0;
	static constexpr double LY_BOHR      = 50.0;
	static constexpr double LZ_BOHR      = 70.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background (full x,y; confined along z = axis 2)
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.0;        // sharp Θ (soften if Gibbs)

	// Electronic structure
	static constexpr int    N_ELECTRONS    = 82;           // even
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 62500
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.312e-3
	static constexpr int    EXTRA_STATES   = 20;
	static constexpr double TEMPERATURE_EV = 0.00862;      // ≈100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (Phase 1 SIE diagnostic; CAP off)
	static constexpr double WP_SIGMA_BOHR = 0.5;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_K0         = k0_from_ev_82(WP_EKIN_EV);   // ≈2.71
	static constexpr double WP_CZ_BOHR    = -32.0;                       // far launch
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
