// ============================================================================
// systems/localised_jellium/shared/configs/slab_n92_L35x35x120_w0p5.hpp
//
// extkin_plateau_E100 (2026-07-29, user-interviewed design — decision log in
// docs/plans/norm-corrected-stopping-power.md "Run design (2026-07-29)"):
//   * Box 35×35×120 Bohr (z = traversal axis), h=0.5 → 70×70×240 grid.
//   * Slab at r_s = 4.0 EXACTLY (n0 = 3/(4π·64)); thickness DERIVED from
//     N=92: t = N/(n0·A) = 20.13 Bohr (half-width 10.067, faces ±10.07).
//     N=92 = closest even to the user's ≥20-Bohr thickness directive.
//   * WP σ=1.5 / E=100 eV (k0=2.711), launch z=−17.5 (5σ from the slab face,
//     18σ from the −z CAP inner edge). Compact-projectile choice: free
//     dispersion ×1.58 at slab entry accepted by the user.
//   * TWO-SIDED CAP: 15 Bohr each end (inner edges ±45), η=−1.0 (replica
//     η=−0.7 rescaled 20→15 Bohr; predicted survival ~8e-4).
//   * dt=0.04, N_STEPS=1500 (t=60 a.u.) — first CAP run at this dt (H·dt=0.79).
//   * The norm-division fix runs IN-RUN via OrbitalKineticStats (all states,
//     every step) — see [[reference_inq_reports_normalized_energy]].
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_N92 = 27.211386245988;
inline constexpr double PI_N92       = 3.14159265358979323846;

inline constexpr double const_sqrt_n92(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 40; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_n92(double e_ev) {
	return const_sqrt_n92(2.0 * e_ev / HA_TO_EV_N92);
}

struct SlabN92_L35x35x120_w0p5 {
	// Cell (z in [-60, +60])
	static constexpr double LX_BOHR      = 35.0;
	static constexpr double LY_BOHR      = 35.0;
	static constexpr double LZ_BOHR      = 120.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background: r_s = 4.0 exact; half-width derived from N
	static constexpr int    N_ELECTRONS    = 92;
	static constexpr double RS_BOHR        = 4.0;
	static constexpr double N0             = 3.0 / (4.0 * PI_N92 * RS_BOHR * RS_BOHR * RS_BOHR); // 3.7301e-3
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  =
		double(N_ELECTRONS) / (2.0 * N0 * LX_BOHR * LY_BOHR);   // 10.067 → 20.13 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.5;

	// Electronic structure
	static constexpr int    EXTRA_STATES   = 16;      // user: >10; 35% of the 46 occupied
	static constexpr double TEMPERATURE_EV = 0.00862; // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (compact σ=1.5; expansion accepted — decision log)
	static constexpr double WP_SIGMA_BOHR = 1.5;
	static constexpr double WP_EKIN_EV    = 100.0;
	static constexpr double WP_K0         = k0_from_ev_n92(WP_EKIN_EV);  // ≈2.711
	static constexpr double WP_CZ_BOHR    = -17.5;    // 5σ from slab face (−10.07)
	static constexpr double WP_KZ         = +WP_K0;

	// TWO-SIDED CAP: 15-Bohr bands at both ends, z ∈ ±[45,60] of [−60,60].
	// mid_z = ±52.5 → mid_frac = 52.5/120 = 0.4375; width = 15/120 = 0.125.
	static constexpr double CAP_ETA_HA     = -1.0;
	static constexpr double CAP_L_BOHR     = 15.0;                       // per side
	static constexpr double CAP_WIDTH_FRAC = CAP_L_BOHR / LZ_BOHR;       // 0.125
	static constexpr double CAP_MID_FRAC   = (LZ_BOHR/2.0 - CAP_L_BOHR/2.0) / LZ_BOHR; // 0.4375

	// Dynamics defaults (overridable via env in run.cpp)
	static constexpr double DT_AU_DEFAULT   = 0.04;
	static constexpr int    N_STEPS_DEFAULT = 1500;   // t = 60 a.u.
};

} // namespace localised_jellium::config
