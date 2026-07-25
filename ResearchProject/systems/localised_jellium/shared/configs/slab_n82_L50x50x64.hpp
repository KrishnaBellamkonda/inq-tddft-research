// ============================================================================
// systems/localised_jellium/shared/configs/slab_n82_L50x50x64.hpp
//
// Source-of-truth Cfg for the `sigma1_massonly` run (2026-07-09): a concentrated
// σ=1 WP projectile, mass-only (accept dispersion), 100 eV, through the r_s≈5.67
// slab, sized to the user's 4σ standoffs so it fits the 1–2 h budget.
//
// Identical slab/density/electrons to slab_n82_L50x50x90.hpp — ONLY the z extent
// SHRINKS 90 → 64 Bohr (the 4σ_WP standoffs on both sides need far less vacuum)
// and the grid is dx=0.40 (not the legacy 0.50). Region layout (z, box [−32,+32]):
//   * Slab  [−12.5, +12.5]  (25 Bohr, half-width 12.5, centred at 0)
//   * Launch z₀ = −16.5      (4σ_WP=4 Bohr before the −12.5 near face)
//   * CAP   [±20.5, ±32]     sin² band, inner edge ±20.5 (4σ_WP=4 Bohr behind
//                            launch), peak ±26.25, outer edge = box edge ±32.
//   * V_inside = 50·50·25 = 62500 Bohr³ (UNCHANGED) → N=82, n₀=1.312e-3, r_s≈5.67.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_64 = 27.21138625;

inline constexpr double const_sqrt_64(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 30; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_64(double e_ev) {
	return const_sqrt_64(2.0 * e_ev / HA_TO_EV_64);
}

struct SlabN82_L50x50x64 {
	// Cell (orthorhombic: x,y = 50 preserve in-plane density; z = 64, box [−32,+32])
	static constexpr double LX_BOHR      = 50.0;
	static constexpr double LY_BOHR      = 50.0;
	static constexpr double LZ_BOHR      = 64.0;
	static constexpr double SPACING_BOHR = 0.40;

	// Slab background (full x,y; confined along z = axis 2) — UNCHANGED vs 90-box
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.0;        // sharp Θ (soften if Gibbs)

	// Electronic structure — UNCHANGED vs 90-box (slab volume identical)
	static constexpr int    N_ELECTRONS    = 82;           // even
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 62500
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.312e-3
	static constexpr int    EXTRA_STATES   = 20;
	static constexpr double TEMPERATURE_EV = 0.00862;      // ≈100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (σ=1, 100 eV) — mass-tuned externally (m_eff≈3.45 via inverse_mass);
	// WP_K0 here is the m=1 electron value, NOT the fork k₀=√(2·E·m). The RT run
	// sets its own EM_K0 (=5.0356 for m=3.45); these fields are provenance only.
	static constexpr double WP_SIGMA_BOHR = 1.0;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_K0         = k0_from_ev_64(WP_EKIN_EV);   // ≈2.71 (m=1 ref)
	static constexpr double WP_CZ_BOHR    = -16.5;                       // 4σ launch
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
