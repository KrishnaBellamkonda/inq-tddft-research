// shared/configs/channeling_tube_rs3.hpp
// ----------------------------------------------------------------------------
// LOCKED geometry + physics table for the CHANNELING TWIN study:
// a matched classical / wavepacket pair shot on-axis down the hollow bore of a
// periodic r_s = 3 annular jellium tube.
//
// Plan: docs/plans/cylindrical-channeling-ks-stopping.md   (design locked by the
// user; this header is the single source of truth every binary compiles against).
//
// AIM. Show that a quantum electron WAVEPACKET, used as a CHANNELING projectile,
// reproduces the classical projectile's stopping power — thereby validating a
// KS-orbital definition of S against the established classical delta-E/ds one.
//
// WHY A TUBE AND WHY THIS FAST. In bulk, the wavepacket's kinetic energy was
// contaminated by an INTERACTION-driven growth of var(p): the momentum-spread
// term rose while the drift term 1/2 <p>^2 stayed flat. var(p) is CONSERVED
// under free evolution, so that growth is interaction, not dispersion.
// Channeling suppresses the interaction (the packet flies through vacuum and
// couples to the wall only via the smooth image force), so dT/dt collapses onto
// the drift channel — which IS the classical projectile kinetic energy. The high
// velocity (v/v_F = 3.0) is not optional: it gets the packet across before
// dispersion delocalises it, and keeps the fractional energy loss small so the
// projectile stays at ~constant velocity.
//
// EVERY NUMBER BELOW IS DERIVED, NOT ASSUMED. N = 160 is the electron count that
// makes n0 = N/V_annulus land on r_s = 3.000000 for this exact geometry; the run
// binaries then set n0 = N/V so integral n_+ = N EXACTLY (the G=0 cancellation
// requirement). Verified numerically 2026-08-01.
//
//   V_annulus = pi (R_out^2 - R_in^2) L_z = 18095.573685 Bohr^3
//   n0        = 160 / V_annulus           = 8.8419413e-3 a0^-3
//   r_s       = (3/(4 pi n0))^(1/3)       = 3.000000
//   omega_p   = sqrt(4 pi n0)             = 0.333333 a.u.  = 9.0705 eV
//   T_plasmon = 2 pi / omega_p            = 18.850 a.u.
//   k_F = v_F = (3 pi^2 n0)^(1/3)         = 0.639719 a.u.
//   E_F                                   = 0.204620 Ha   = 5.5680 eV
//   k_TF      = sqrt(4 k_F / pi)          = 0.902505 -> screening length 1.108 Bohr
//
// PROJECTILE (identical in both twins except for the representation):
//   E_proj    = 50 eV = 1.83746611 Ha  ->  v = k0 = 1.91701127 a.u.
//   v / v_F   = 2.997          (fast, well clear of the Lindhard peak)
//   lambda_p  = 2 pi v / omega_p = 36.14 Bohr  (fits inside L_z = 60; do NOT
//               shorten L_z below ~48 or the wake self-overlaps)
//   sigma_WP  = 4.0 Bohr  <- THE LABEL. Both twins are reported at this sigma.
//   sigma_pot = sigma_WP/sqrt2 = 2.82842712 Bohr  <- the classical CHARGE std,
//               an internal quantity (.claude/rules/sigma-wp-convention.md). It
//               is the value that makes the classical charge std equal the WP
//               DENSITY std, i.e. makes the two projectiles source the SAME
//               potential. Never label a run or an axis by sigma_pot.
//   k0 sigma_WP = 7.67   (drift dominates the momentum spread by 10.8 sigma_p)
//
// DISPERSION BUDGET (sigma_d(t) = sqrt(sigma_WP^2/2 + t^2/(2 sigma_WP^2))):
//   t =  0 a.u.  sigma_d = 2.83   2 sigma_d =  5.66   (clears the R_in=10 bore)
//   t = 13       sigma_d = 3.64   2 sigma_d =  7.29   (end of the fit window)
//   t = 23.3     sigma_d = 5.00   2 sigma_d = 10.00   <- 2 sigma_d REACHES the wall
//   t = 30       sigma_d = 6.01   2 sigma_d = 12.02   (end of the run)
//   6 sigma_d = L_xy = 40 at t = 34.1 a.u. > 30, so the packet NEVER overlaps its
//   own transverse periodic image during the run — a limitation the slab study
//   had to accept and this geometry does not.
//   => the KS-stopping fit window is t in [transient, ~13 a.u.]; f_bore(t) from
//      inqkit::observables::radial_occupancy is the MEASURED version of this
//      budget and is what the analysis actually gates on.
//
// GRID / ALIASING. dx = 0.5 -> k_Nyq = 6.2832, E_cut = 537 eV. The WP momentum
// std is sigma_p = 1/(sqrt2 sigma_WP) = 0.1768, so k0 + 3 sigma_p = 2.45 << k_Nyq
// and the aliased tail is 0.00 %. cutoff_guard.py: PASS for both twins
// (verified 2026-08-01).
//
// TRAJECTORY. Launch on-axis at z = -28 (2 Bohr inside the -z face) moving +z.
// Over 1500 steps of dt = 0.02 (t = 30 a.u.) the projectile covers 57.5 Bohr and
// ends at z = +29.5, i.e. ONE traversal of the periodic tube with no wrap. The
// launch point is nevertheless within 0.7 sigma_pot of the face, so BOTH twins
// must use the minimum-image kernels (charge, force, occupancy, centroid) — the
// wavepacket wraps exactly on the FFT basis and the classical charge must too.
// ----------------------------------------------------------------------------
#pragma once

#include <cmath>

namespace cylindrical_jellium::config {

struct ChannelingTubeRs3 {

	// ---- cell -------------------------------------------------------------
	static constexpr double LX_BOHR  = 40.0;
	static constexpr double LY_BOHR  = 40.0;
	static constexpr double LZ_BOHR  = 60.0;
	static constexpr double SPACING_BOHR = 0.5;      // 80 x 80 x 120 grid
	static constexpr int    PERIODICITY  = 3;        // fully periodic: the tube is
	                                                 // infinite along z and the
	                                                 // transverse cell is a lattice

	// ---- tube -------------------------------------------------------------
	static constexpr double R_IN_BOHR   = 10.0;      // hollow bore radius
	static constexpr double R_OUT_BOHR  = 14.0;      // wall outer radius (4 Bohr wall
	                                                 // = 3.6 screening lengths)
	static constexpr double EDGE_W_BOHR = 0.5;       // erfc softening at BOTH radial
	                                                 // edges (~1 grid cell)
	static constexpr int    TUBE_AXIS   = 2;         // z
	static constexpr double CENTER_X = 0.0, CENTER_Y = 0.0, CENTER_Z = 0.0;

	// ---- electrons --------------------------------------------------------
	static constexpr int    N_ELECTRONS   = 160;     // bath; -> r_s = 3.000000
	static constexpr int    EXTRA_STATES  = 24;      // 80 occupied + 24 -> 104 states
	                                                 // (the WP occupies the LAST one)
	static constexpr double TEMPERATURE_EV = 0.00862;
	static constexpr double SCF_TOL_HA     = 1.0e-7;
	static constexpr int    SCF_MAX_STEPS  = 400;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.3;

	// ---- projectile (shared by both twins) --------------------------------
	static constexpr double PROJ_ENERGY_EV = 50.0;
	static constexpr double PROJ_V0        = 1.91701127;   // = k0 = sqrt(2 E), m = 1
	static constexpr double SIGMA_WP_BOHR  = 4.0;          // THE LABEL
	static constexpr double LAUNCH_Z_BOHR  = -28.0;
	static constexpr double PROJ_MASS      = 1.0;          // electron mass, a.u.
	static constexpr double PROJ_CHARGE    = -1.0;

	// ---- propagation ------------------------------------------------------
	static constexpr double DT_AU      = 0.02;
	static constexpr int    N_STEPS    = 1500;             // t = 30 a.u.
	static constexpr int    SAVE_EVERY = 5;                // -> 300 density frames
	static constexpr int    WF_EVERY   = 75;               // -> 20 wavefunction dumps
	static constexpr double FD_DELTA   = 0.05;             // classical force step

	// ---- derived (constexpr so a binary can print them and a test can pin them)
	static double v_annulus()  { return M_PI * (R_OUT_BOHR*R_OUT_BOHR - R_IN_BOHR*R_IN_BOHR) * LZ_BOHR; }
	static double n0(int n_elec = N_ELECTRONS) { return double(n_elec) / v_annulus(); }
	static double r_s(int n_elec = N_ELECTRONS) { return std::cbrt(3.0 / (4.0*M_PI*n0(n_elec))); }
	static double omega_p(int n_elec = N_ELECTRONS) { return std::sqrt(4.0*M_PI*n0(n_elec)); }
	static double k_fermi(int n_elec = N_ELECTRONS) { return std::cbrt(3.0*M_PI*M_PI*n0(n_elec)); }

	// sigma_pot = sigma_WP / sqrt2 — the classical charge std that makes the two
	// projectiles source the same potential. INTERNAL ONLY (sigma-wp-convention).
	static double sigma_pot() { return SIGMA_WP_BOHR / std::sqrt(2.0); }

	// Free-particle density spread of the wavepacket at time t (Bohr).
	static double sigma_d(double t_au) {
		return std::sqrt(SIGMA_WP_BOHR*SIGMA_WP_BOHR/2.0
		                 + t_au*t_au/(2.0*SIGMA_WP_BOHR*SIGMA_WP_BOHR));
	}
};

}  // namespace cylindrical_jellium::config
