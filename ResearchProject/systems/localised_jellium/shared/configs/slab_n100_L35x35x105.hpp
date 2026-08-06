// ============================================================================
// systems/localised_jellium/shared/configs/slab_n100_L35x35x105.hpp
//
// Source-of-truth Cfg for the sigma_WP = 5 and 6 Bohr twin campaign
// (docs/plans/sigma56-sv-twin.md, sweep `sigma56_sv`). User-locked 2026-08-02.
//
// IDENTICAL to SlabN100_L35x35x85 except LZ_BOHR 85 -> 105 and the production
// spacing 0.50 -> 0.40. The extra 20 Bohr is PURE VACUUM added along z:
//   * slab thickness (25 Bohr), transverse box (35 x 35) and N = 100 are
//     unchanged, so V_inside = 30625 Bohr^3, n0 = 3.2653e-3 and r_s = 4.183 are
//     BIT-IDENTICAL to the 85-Bohr config. The bath is the same physical system.
//   * Only E_GS shifts (a different box is a different calculation), which is
//     why the GS run reports it as INFO and does not gate on the 85-Bohr value.
//
// WHY THE BOX GREW. A wavepacket of width sigma_WP has a t=0 density std of
// sigma_WP/sqrt(2) = 4.243 Bohr at sigma = 6. In the 85-Bohr box the launch point
// z = -24 sits only 6 Bohr from the CAP inner edge (-30) and 11.5 Bohr from the
// slab face (-12.5), so a sigma = 6 packet cannot be placed without either 7.9 %
// of it starting inside the absorber or ~10 % of it starting inside the slab:
// 3 sigma_d of clearance on each side needs 25.5 Bohr and the box offers 17.5.
// At LZ = 105 the CAP inner edge moves to -40, the gap becomes 27.5 Bohr, and
// launch z = -27.5 leaves 12.5 Bohr to the CAP (2.95 sigma_d, 0.16 % of the
// packet inside) and 15.0 Bohr to the slab (3.54 sigma_d, 0.020 % inside) — at
// or below the 0.23 % t=0 CAP loss already accepted for the sigma = 3 campaign.
//
// Grid at dx = 0.40: 88 x 88 x 263 (LZ/dx = 262.5). E_cut = (pi/dx)^2/2 =
// 30.84 Ha = 839 eV single-particle Nyquist; sigma_p = 1/(sqrt2 sigma) = 0.141
// (sigma=5) / 0.118 (sigma=6) against k_Nyq = 7.85, so momentum aliasing is
// ~1e-14 % at every velocity in the grid — the constraint that cut v = 4.0/4.5
// from the sigma = 0.5 campaign does not bite here.
//
// Boundary: periodicity(2) — x,y periodic (infinite slab), z open for the
// Poisson solve. NOTE this does NOT make z open for the ORBITALS (the kinetic
// operator is a plain 3-D FFT), which is why both halves of this campaign carry
// absorbing bands; see docs/handovers/wavepacket-highdensity-sv-twin.md.
//
// Projectile sigma is NOT pinned here — it is an env knob (LJ_SIGMA) because the
// campaign sweeps it. WP_SIGMA_BOHR below is only the default. sigma_pot is
// always derived as sigma_WP/sqrt(2) inside the run binaries
// (.claude/rules/sigma-wp-convention.md).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

struct SlabN100_L35x35x105 {
	// Cell (orthorhombic; z in [-52.5, +52.5]); periodicity(2) applied in run.cpp
	static constexpr double LX_BOHR      = 35.0;
	static constexpr double LY_BOHR      = 35.0;
	static constexpr double LZ_BOHR      = 105.0;
	static constexpr double SPACING_BOHR = 0.40;   // production grid

	// Slab background — UNCHANGED from the 85-Bohr config (25 Bohr thick)
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;   // 25 Bohr thick (L_slab)
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 1.0;    // erfc-softened (GS-study H1)

	// Electronic structure — UNCHANGED: N=100 over the same 30625 Bohr^3 interior
	static constexpr int    N_ELECTRONS    = 100;
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 30625
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 3.2653e-3 -> r_s = 4.183
	static constexpr int    EXTRA_STATES   = 24;       // ~50 occupied + headroom (T=100 K)
	static constexpr double TEMPERATURE_EV = 0.00862;  // ~100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Absorbing bands — perturbations::absorbing takes FRACTIONAL cell coords,
	// NOT Bohr (it compares point_op.rvector()[2] in [-0.5, 0.5)). Passing Bohr
	// would place the CAP at z ~ 0.4 Bohr, straight through the slab centre.
	static constexpr double CAP_L_BOHR     = 12.5;                        // per z face
	static constexpr double CAP_ETA_HA     = -1.0;                        // < 0 absorbs
	static constexpr double CAP_WIDTH_FRAC = CAP_L_BOHR / LZ_BOHR;        // 0.119047619048
	static constexpr double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;  // 0.440476190476 (= 46.25 Bohr)
	static constexpr double CAP_Z_INNER    = LZ_BOHR / 2.0 - CAP_L_BOHR;  // +/- 40.0 Bohr

	// Projectile. sigma is swept (LJ_SIGMA); this is only the campaign default.
	// The classical twin's Gaussian width is sigma_pot = sigma_WP/sqrt(2), derived
	// in the run binaries, and BOTH halves are labelled by sigma_WP.
	static constexpr double WP_SIGMA_BOHR  = 6.0;      // sigma_WP (label), default
	static constexpr double PROJ_MASS      = 1.0;      // classical electron
	static constexpr double PROJ_CHARGE    = -1.0;     // electron
	static constexpr double LAUNCH_Z_BOHR  = -27.5;    // common to both halves and both sigma
};

} // namespace localised_jellium::config
