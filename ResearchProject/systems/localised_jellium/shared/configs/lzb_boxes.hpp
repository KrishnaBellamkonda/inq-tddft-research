// ============================================================================
// systems/localised_jellium/shared/configs/lzb_boxes.hpp
//
// The FOUR box presets of the slab->bulk L_slab sweep (`lz_bulk_sweep`).
// Plan: docs/plans/jellium-slab-extend-Lz.md (user-locked 2026-08-05).
//
// The sweep varies the slab THICKNESS at fixed density (n0 = 3.2653e-3,
// r_s = 4.183 — bit-identical to every prior slab campaign) and extrapolates
// the deposit stopping power in 1/L_slab. Geometry is PER-sigma-FAMILY so each
// family's new runs match its existing L_slab = 25 anchor exactly:
//
//   sigma = 0.5 family — 85-box layout: launch standoff 11.5 Bohr from the
//     face, face->CAP gap 17.5, so L_z = L_slab + 60. Anchor: the completed
//     wp_highdensity_sv sigma = 0.5 runs (launch z = -24 in the 85-Bohr box).
//   sigma = 5 family — 105-box layout: standoff 15, gap 27.5, L_z = L_slab
//     + 80. Anchor: the completed sigma56_sv sigma = 5 twins (launch -27.5).
//
// Holding the standoff fixed WITHIN a family makes the packet's arrival width
// a function of (sigma, v) only — L_slab drops out — which is the user's
// comparability requirement. The residual in-transit <sigma_d> drift with L is
// analytic and small at sigma = 5 (<= ~10 %).
//
// Unlike SlabN100_L35x35x105 these presets are selected at RUNTIME (env
// LZB_CFG), so ONE binary per run type serves all four boxes — the pattern the
// classical/vac sigma56 binaries already use for their env-driven geometry.
// Everything not listed per box is shared and unchanged from the sigma56
// config: dx = 0.40, L_xy = 35, erfc edge 1.0, CAP 12.5 Bohr per face at
// eta = -1 Ha, T = 100 K, LDA.
//
// EXTRA_STATES scales with N_e (same ~48 % headroom as 100 e- -> 74 states):
// 60 e- -> 30 occ + 15 = 45 states; 140 e- -> 70 occ + 34 = 104 states.
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cstdlib>
#include <iostream>
#include <string>

namespace localised_jellium::config {

// ---- shared across every box of the sweep (matches SlabN100_L35x35x105) ----
struct LzbShared {
	static constexpr double LX_BOHR        = 35.0;
	static constexpr double LY_BOHR        = 35.0;
	static constexpr double SPACING_BOHR   = 0.40;
	static constexpr double EDGE_WIDTH_BOHR= 1.0;
	static constexpr double CAP_L_BOHR     = 12.5;   // per z face
	static constexpr double CAP_ETA_HA     = -1.0;   // < 0 absorbs
	static constexpr double TEMPERATURE_EV = 0.00862;
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;
	static constexpr double PROJ_MASS      = 1.0;
	static constexpr double PROJ_CHARGE    = -1.0;
	static constexpr int    SLAB_AXIS      = 2;
	static constexpr double SLAB_CENTER    = 0.0;
};

struct LzbBox {
	std::string name;          // preset key = run-name prefix
	double      LZ_BOHR;       // box length along z
	double      SLAB_HALF;     // L_slab / 2
	int         N_ELECTRONS;   // = n0 * 35 * 35 * L_slab (n0 held fixed)
	int         EXTRA_STATES;
	double      LAUNCH_Z;      // = -(SLAB_HALF + family standoff)
	double      SIGMA_DEFAULT; // the family's sigma_WP
	std::string GS_TAG;        // shared_gs/<GS_TAG> checkpoint directory

	double n0()        const { return double(N_ELECTRONS) / (LzbShared::LX_BOHR * LzbShared::LY_BOHR * 2.0 * SLAB_HALF); }
	double l_slab()    const { return 2.0 * SLAB_HALF; }
	double cap_z_inner() const { return LZ_BOHR / 2.0 - LzbShared::CAP_L_BOHR; }
	double standoff()  const { return -SLAB_HALF - LAUNCH_Z; }
};

// Preset table. n0 check: 60/(1225*15) = 140/(1225*35) = 100/(1225*25)
// = 3.2653e-3 exactly — the gs binary gates on r_s = 4.183 from this.
inline LzbBox lzb_box(std::string const& key) {
	if (key == "s0p5_L15") return {key,  75.0,  7.5,  60,  15, -19.0, 0.5, "slab_n60_L35x35x75_dx0p4_per2"};
	if (key == "s0p5_L35") return {key,  95.0, 17.5, 140,  34, -29.0, 0.5, "slab_n140_L35x35x95_dx0p4_per2"};
	if (key == "s5p0_L15") return {key,  95.0,  7.5,  60,  15, -22.5, 5.0, "slab_n60_L35x35x95_dx0p4_per2"};
	if (key == "s5p0_L35") return {key, 115.0, 17.5, 140,  34, -32.5, 5.0, "slab_n140_L35x35x115_dx0p4_per2"};
	std::cerr << "FATAL: unknown LZB_CFG '" << key << "' (want s0p5_L15, s0p5_L35, s5p0_L15, s5p0_L35)\n";
	std::exit(2);
}

inline LzbBox lzb_box_from_env() {
	const char* v = std::getenv("LZB_CFG");
	if (!v || std::string(v).empty()) {
		std::cerr << "FATAL: LZB_CFG is not set — this binary serves four boxes and "
		             "refuses to guess which one you meant.\n";
		std::exit(2);
	}
	return lzb_box(v);
}

} // namespace localised_jellium::config
