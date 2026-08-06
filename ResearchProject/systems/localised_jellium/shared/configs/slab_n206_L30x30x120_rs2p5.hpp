// ============================================================================
// systems/localised_jellium/shared/configs/slab_n206_L30x30x120_rs2p5.hpp
//
// Source-of-truth Cfg for the NAZAROV-GROSS MASS LADDER campaign.
// Plan: docs/plans/nazarov-gross-slab-mass-ladder.md
//
// WHY A NEW, DENSER SLAB. The campaign needs a projectile of mass ~1 to CROSS
// the slab at a velocity BELOW the Bragg peak (~1.6 v_F). Deposit/KE at
// v = 1.4 v_F, sigma_WP = 4, 15 Bohr slab (linear-response sizing):
//     r_s = 5.665 (the existing slab)  ->  1.63   projectile STOPS
//     r_s = 4.0                        ->  0.86   marginal
//     r_s = 2.5 (this config)          ->  0.30   crosses, ~16% v loss
// KE ~ v_F^2 ~ 1/r_s^2 rises steeply with density while the form-factor-
// suppressed deposit stays nearly flat (4.3 -> 5.0 eV over that whole range),
// so density is the only lever that buys traversal. The existing r_s = 5.665
// slab CANNOT run this campaign.
//
//   * Cell: 30 x 30 x 120 Bohr orthorhombic periodic (INQ-centred, z in [-60,+60]).
//   * Slab: full 30x30 face, 15 Bohr thick (half-width 7.5), centred at z = 0.
//     V_inside = 30*30*15 = 13500 Bohr^3.
//   * N = 206 electrons. n0 = 206/13500 = 1.52593e-2 a0^-3 -> r_s = 2.5011.
//     v_F = k_F = 0.76733, E_F = 8.011 eV, omega_p = 0.43790 (T_p = 14.35 a.u.).
//     Bragg peak ~ 1.6 v_F = 1.228 a.u.
//   * Region layout (z): slab [-7.5,7.5] . free [+/-7.5,+/-45] . CAP [+/-45,+/-60].
//   * CAP: TWO-SIDED sin^2 absorber, eta = -1.0 Ha, 15 Bohr/side, region
//     [+/-45,+/-60], each bump peaking at +/-52.5. eta = -1.0 is the VALIDATED
//     value from the muon-mass-fork CAP study (2026-07-11): complete WP
//     absorption, 0.13% reflection, zero bath drain (N_min = 51.999 of 52).
//     Do NOT weaken below -1.0: eta = -0.4 gave 4.1% edge reflection.
//   * Projectile: WIDE packet sigma_WP = 4.0 Bohr, launched +z at z = -25.0
//     (20 Bohr = 5 sigma clear of the CAP inner face; 17.5 Bohr run-up to the
//     slab face). The classical twin uses sigma_pot = sigma_WP/sqrt(2) =
//     2.828427 (.claude/rules/sigma-wp-convention.md).
//   * VELOCITY, not energy, is the sweep invariant: every mass runs at the SAME
//     v0 = 1.40 v_F = 1.074266 a.u. (0.875 of the Bragg peak). Each run's
//     k0 = MASS * V0_AU is therefore computed in run.cpp, not here.
//
// GRID / TIMESTEP (both calibrated in the plan, sections 4.1 and 4.2):
//   * Aliasing needs pi/h >= M*v0 + 3/(2*sigma_WP). At sigma_WP = 4 the
//     coarsest allowed h is 2.15 (M=1) / 0.87 (M=3) -- all far above the
//     h = 0.50 used here, so the BATH sets the grid, not the projectile.
//   * dt <= 0.08 * M * h^2  (= 0.02*M at h = 0.5). Calibrated exactly on two
//     working runs: p3 (h=0.5, M=1, dt=0.02) and sigma1_masspair (h=0.5, M=2,
//     dt=0.04). => M = 3 / 1.2 / 1.0 / 0.5 run at dt = 0.06 / 0.024 / 0.02 /
//     0.01. INFERENCE from two points: confirm with the Phase-2 drift gate.
//     NOTE this means cost ~ 1/M: light masses are MORE expensive, not less.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_206 = 27.21138625;

struct SlabN206_L30x30x120_rs2p5 {
	// Cell (orthorhombic: x,y = 30 preserve in-plane density; z = 120 vacuum+CAP)
	static constexpr double LX_BOHR      = 30.0;
	static constexpr double LY_BOHR      = 30.0;
	static constexpr double LZ_BOHR      = 120.0;
	static constexpr double SPACING_BOHR = 0.50;   // -> 60 x 60 x 240 = 864k points
	// periodicity(2): x,y periodic (infinite slab), z OPEN. The projectile cannot
	// wrap, so the own-wake re-entry that capped every window in the channeling
	// study at t ~ 22 a.u. cannot occur here; the CAPs terminate z instead.
	static constexpr int    PERIODICITY  = 2;

	// Slab background (full x,y; confined along z = axis 2)
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 7.5;   // 15 Bohr thick (L_slab)
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;   // erfc-softened, >= dx = 0.50

	// Electronic structure
	static constexpr int    N_ELECTRONS    = 206;     // even; 103 occupied
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 13500
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.52593e-2
	static constexpr int    EXTRA_STATES   = 20;      // 103 + 20 (+1 WP) = 124 states
	static constexpr double TEMPERATURE_EV = 0.00862; // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Derived electron-gas scales (recorded so runs/analysis never re-derive them).
	// Self-checked against closed forms by shared/configs/tests/ — see that test
	// for the tolerances; every value below is exact to the digits shown.
	static constexpr double RS_BOHR    = 2.5010708;  // (3/(4 pi N0))^(1/3)
	static constexpr double KF_AU      = 0.7673347;  // = v_F
	static constexpr double EF_EV      = 8.01107;
	static constexpr double OMEGA_P_AU = 0.4378967;  // sqrt(3/r_s^3), T_p = 14.35 a.u.
	static constexpr double V_BRAGG_AU = 1.2277355;  // ~1.6 v_F, the Bragg-peak proxy

	// CAP: two-sided sin^2, region [+/-45,+/-60] -> fractional mid/width for
	// perturbations::absorbing(eta, mid_frac, width_frac), matching the
	// slab_n82_L50x50x90 convention (mid = region centre, width = FULL extent).
	static constexpr double CAP_ETA_HA    = -1.0;
	static constexpr double CAP_MID_FRAC  = 52.5 / LZ_BOHR;  // 0.4375
	static constexpr double CAP_WIDTH_FRAC = 15.0 / LZ_BOHR; // 0.125
	static constexpr double CAP_INNER_FACE_BOHR = 45.0;

	// Projectile — velocity is the sweep invariant; k0 = MASS * V0_AU per run
	static constexpr double WP_SIGMA_BOHR   = 4.0;
	static constexpr double WP_SIGMA_POT    = 2.8284271;  // sigma_WP/sqrt(2), classical twin
	static constexpr double V0_AU           = 1.0742685;  // 1.40 v_F, 0.875 of the Bragg peak
	static constexpr double WP_CZ_BOHR      = -25.0;      // 5 sigma clear of the CAP face
	static constexpr double ANALYSIS_Z_END  = +30.0;      // stop before the far CAP: 55 Bohr path
	// KE = 0.5 * M * V0_AU^2:  M=3 -> 47.11 eV, 1.2 -> 18.84, 1.0 -> 15.70, 0.5 -> 7.85 eV
};

} // namespace localised_jellium::config
