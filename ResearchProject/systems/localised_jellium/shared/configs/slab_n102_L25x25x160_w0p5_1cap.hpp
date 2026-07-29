// ============================================================================
// systems/localised_jellium/shared/configs/slab_n102_L25x25x160_w0p5_1cap.hpp
//
// REPLICA (2026-07-28, user-specified) of SlabN102_L25x25x140_w0p5 with:
//   * LZ 140 -> 160 Bohr (longer box; z ∈ [−80,+80]).
//   * CAP: ONE-SIDED, 20 Bohr TOTAL at the +z end (z ∈ [60,80]) — same total CAP
//     length as the old two-sided 10+10, but on one side only. WP launches on the
//     −z side and traverses the slab into it.
//   * N_STEPS default 8000 (t=160 a.u.) — set in run.cpp, not here.
// Everything else IDENTICAL to the 140 config: LX=LY=25, slab half-width 12.5,
// N=102 -> r_s=3.32 (SAME density; LZ does not affect density), WP σ=1/E=100 eV,
// launch z=−20.5, grid h=0.5, 24 extra states, T≈100 K.
//
// The "energy fix" (INQ reports per-particle kinetic; corrected extensive energy =
// E_total − e_kin_ha·(1−norm)) is applied in the run's ANALYSIS from the per-step
// e_kin_ha (wp_momentum_stats) + norm — no engine change. See
// [[reference_inq_reports_normalized_energy]].
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

struct SlabN102_L25x25x160_w0p5_1cap {
	// Cell (LX=LY=25 as before; LZ 140 -> 160)
	static constexpr double LX_BOHR      = 25.0;
	static constexpr double LY_BOHR      = 25.0;
	static constexpr double LZ_BOHR      = 160.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background (UNCHANGED -> same r_s = 3.32)
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.5;

	// Electronic structure (UNCHANGED). V_inside independent of LZ.
	static constexpr int    N_ELECTRONS    = 102;
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 15625
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 6.528e-3
	static constexpr int    EXTRA_STATES   = 24;
	static constexpr double TEMPERATURE_EV = 0.00862;
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (UNCHANGED)
	static constexpr double WP_SIGMA_BOHR = 1.0;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_K0         = k0_from_ev_n102(WP_EKIN_EV);   // ≈2.711
	static constexpr double WP_CZ_BOHR    = -20.5;
	static constexpr double WP_KZ         = +WP_K0;

	// ONE-SIDED CAP: single 20-Bohr band at the +z end, z ∈ [60,80] of [−80,80].
	// mid_z = 70 -> mid_frac = 70/160 = 0.4375; width = 20/160 = 0.125.
	static constexpr double CAP_ETA_HA    = -0.7;
	static constexpr double CAP_L_BOHR    = 20.0;                             // total (one-sided)
	static constexpr double CAP_WIDTH_FRAC= CAP_L_BOHR / LZ_BOHR;             // 0.125
	static constexpr double CAP_MID_FRAC  = (LZ_BOHR/2.0 - CAP_L_BOHR/2.0) / LZ_BOHR; // 0.4375 (+z only)
};

} // namespace localised_jellium::config
