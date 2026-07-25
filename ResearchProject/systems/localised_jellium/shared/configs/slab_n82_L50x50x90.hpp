// ============================================================================
// systems/localised_jellium/shared/configs/slab_n82_L50x50x90.hpp
//
// Source-of-truth Cfg for the `quantum-stopping-power` campaign, Phase 3
// (big-box / long-time σ=0.5 WP+classical pair, energy-method stopping).
// Identical slab/density/electrons to slab_n82_L50x50x70.hpp — ONLY the z
// extent grows 70 → 90 Bohr (more vacuum) so that:
//   * the WP launches equidistant (11.25 Bohr) from the slab face and the CAP,
//     killing the t=0 absorption seen at the 70-box launch, and
//   * collective (plasmon) flux reaches the CAP later (≈66 a.u.), keeping the
//     retained-energy "system" cleaner over the τ=100 a.u. run.
//
//   * Cell: 50 × 50 × 90 Bohr orthorhombic periodic (INQ-centred, z ∈ [−45,+45]).
//   * Slab: full 50×50 face, 25 Bohr thick (half-width 12.5), centred at z=0.
//     V_inside = 50·50·25 = 62500 Bohr³ (UNCHANGED by the z-extension).
//   * N = 82 electrons. n₀ = 82/62500 = 1.312e-3 a₀⁻³ → r_s ≈ 5.665 (UNCHANGED).
//   * Region layout (z): slab [−12.5,12.5] · free [±12.5,±35] · CAP [±35,±45].
//   * CAP: TWO-SIDED sin² absorber (stock engine, cap_lo+cap_hi at mid=±40/90), η=−0.7 Ha,
//     10 Bohr/side, region [±35,±45] (inner faces ±35), each bump peaks at ±40. The
//     benchmarked "known devil" (~1.3% reflection); a seam-centred variant was reverted.
//   * Projectile: σ=0.5 Bohr, E=100 eV, k₀≈2.71, launched +z at z=−23.75
//     (equidistant: 11.25 Bohr from both the −12.5 slab face and the −35 CAP face).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_90 = 27.21138625;

inline constexpr double const_sqrt_90(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 30; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_90(double e_ev) {
	return const_sqrt_90(2.0 * e_ev / HA_TO_EV_90);
}

struct SlabN82_L50x50x90 {
	// Cell (orthorhombic: x,y = 50 preserve in-plane density; z = 90 vacuum)
	static constexpr double LX_BOHR      = 50.0;
	static constexpr double LY_BOHR      = 50.0;
	static constexpr double LZ_BOHR      = 90.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background (full x,y; confined along z = axis 2) — UNCHANGED vs 70-box
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.0;        // sharp Θ (soften if Gibbs)

	// Electronic structure — UNCHANGED vs 70-box (slab volume identical)
	static constexpr int    N_ELECTRONS    = 82;           // even
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 62500
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.312e-3
	static constexpr int    EXTRA_STATES   = 20;
	static constexpr double TEMPERATURE_EV = 0.00862;      // ≈100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (σ=0.5, 100 eV) — equidistant launch in the bigger box
	static constexpr double WP_SIGMA_BOHR = 0.5;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_K0         = k0_from_ev_90(WP_EKIN_EV);   // ≈2.71
	static constexpr double WP_CZ_BOHR    = -23.75;                      // equidistant launch
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
