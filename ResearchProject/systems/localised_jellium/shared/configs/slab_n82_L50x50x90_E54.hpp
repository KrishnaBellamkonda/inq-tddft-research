// ============================================================================
// systems/localised_jellium/shared/configs/slab_n82_L50x50x90_E54.hpp
//
// Phase-4 energy-matched S(v) point: σ_WP = 0.5 projectile at v = 2.0 a.u.
// (E = ½ m v² = 54.42 eV), an ON-GRID classical S(v) energy
// (jellium/hypotheses/06_sigma_convergence → sv_convergence_energy.png).
//
// IDENTICAL to slab_n82_L50x50x90.hpp (same 90-box, slab, density, electrons,
// GS, CAP, equidistant launch z = −23.75) — ONLY the projectile energy changes
// 100 eV → 54.42 eV, so v drops 2.711 → 2.0. The GS in
// shared_gs/slab_n82_L50x50x90 is REUSED unchanged (energy is a real-time-only
// quantity). σ is labelled as σ_WP (UPF stays electron_gaussian_sigma0p35.upf,
// σ_pot = σ_WP/√2).
//
// All quantities atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace localised_jellium::config {

inline constexpr double HA_TO_EV_E54 = 27.21138625;

inline constexpr double const_sqrt_E54(double x) {
	if (x <= 0.0) return 0.0;
	double g = x > 1.0 ? x : 1.0;
	for (int i = 0; i < 30; ++i) g = 0.5 * (g + x / g);
	return g;
}
inline constexpr double k0_from_ev_E54(double e_ev) {
	return const_sqrt_E54(2.0 * e_ev / HA_TO_EV_E54);
}

struct SlabN82_L50x50x90_E54 {
	// Cell (orthorhombic: x,y = 50 preserve in-plane density; z = 90 vacuum)
	static constexpr double LX_BOHR      = 50.0;
	static constexpr double LY_BOHR      = 50.0;
	static constexpr double LZ_BOHR      = 90.0;
	static constexpr double SPACING_BOHR = 0.50;

	// Slab background (full x,y; confined along z = axis 2) — UNCHANGED
	static constexpr int    SLAB_AXIS        = 2;
	static constexpr double SLAB_HALF_WIDTH  = 12.5;       // 25 Bohr thick
	static constexpr double SLAB_CENTER_BOHR = 0.0;
	static constexpr double EDGE_WIDTH_BOHR  = 0.0;

	// Electronic structure — UNCHANGED (slab volume identical; GS reused)
	static constexpr int    N_ELECTRONS    = 82;
	static constexpr double V_INSIDE_BOHR3 = LX_BOHR * LY_BOHR * (2.0 * SLAB_HALF_WIDTH); // 62500
	static constexpr double N0             = double(N_ELECTRONS) / V_INSIDE_BOHR3;        // 1.312e-3
	static constexpr int    EXTRA_STATES   = 20;
	static constexpr double TEMPERATURE_EV = 0.00862;      // ≈100 K
	static constexpr double SCF_TOL_HA     = 1.0e-4;
	static constexpr int    SCF_MAX_STEPS  = 300;
	static constexpr int    SCF_MIX_NDIM   = 8;
	static constexpr double SCF_MIX_ALPHA  = 0.1;

	// Projectile (σ_WP = 0.5, v = 2.0 a.u. → 54.42 eV) — equidistant launch
	static constexpr double WP_SIGMA_BOHR = 0.5;
	static constexpr double WP_EKIN_EV    = 54.422772491976;   // ½·2.0²·HA_TO_EV → v=2.0
	static constexpr double WP_K0         = k0_from_ev_E54(WP_EKIN_EV);   // = 2.0
	static constexpr double WP_CZ_BOHR    = -23.75;                       // equidistant launch
	static constexpr double WP_KZ         = +WP_K0;
};

} // namespace localised_jellium::config
