// ============================================================================
// shared/configs/tsubonoya_2014_base.hpp
//
// Tsubonoya, Hu & Watanabe, Phys. Rev. B 90, 035416 (2014) — coronene LEED
// reference parameters. Source-of-truth for every coronene replication run
// under ResearchProject/systems/coronene/run_*/.
//
// Cfg pattern: every variant (E30, E800, s0p33, ...) is a struct that inherits
// from Base and shadows only the constants that differ. The propagation
// driver is a single function template `run_propagation<Cfg>(...)` so each
// run.cpp picks its variant struct, instantiates the template, and is done.
//
// Convention: INQ uses centred orthorhombic cells [-L/2, +L/2]. The molecule
// (shared/geometry/coronene.xyz) sits at z = 0; the wave packet is launched
// from z = +b along -z.
//
// All quantities in atomic units / Bohr / eV.
// ============================================================================
#pragma once

#include <cmath>

namespace coronene::config {

// ---- Unit conversions (CODATA 2018) ---------------------------------------
inline constexpr double ANG_TO_BOHR = 1.8897259886;
inline constexpr double HA_TO_EV    = 27.21138625;
inline constexpr double FS_TO_AU    = 41.341374575751;

// Compile-time sqrt (Newton iteration). std::sqrt is not constexpr under
// CUDA, so we ship our own. 30 iterations is well past double-precision
// convergence for x in any plausible scientific range.
inline constexpr double const_sqrt(double x) {
    if (x <= 0.0) return 0.0;
    double g = x > 1.0 ? x : 1.0;
    for (int i = 0; i < 30; ++i) {
        g = 0.5 * (g + x / g);
    }
    return g;
}

// Helper: |k| from kinetic energy in eV.
inline constexpr double k0_from_ev(double ekin_ev) {
    // |k| = sqrt(2 * E_Ha) = sqrt(2 * E_eV / HA_TO_EV)
    return const_sqrt(2.0 * ekin_ev / HA_TO_EV);
}

// Total propagation time = time for (centroid + 1*sigma) of the WP to reach
// the far end of the box at speed |k|. The WP starts at z = +offset moving
// in -z; the far end is at z = -Lz/2. So the trailing edge (z = +offset+sigma)
// has cleared the box when z_cen(t) = -Lz/2 - sigma, i.e. after travelling
//
//     (offset + sigma) - (-Lz/2) = offset + sigma + Lz/2
//
// at speed |k|. N_steps is the closest multiple of dt to that end time.
inline constexpr int compute_n_steps(double lz, double offset,
                                     double sigma, double k0, double dt) {
    const double end_time_au = (offset + sigma + 0.5 * lz) / k0;
    int n = static_cast<int>(end_time_au / dt + 0.5);
    return n > 0 ? n : 1;
}

// ---- Base config: the Tsubonoya 2014 paper parameters --------------------
struct Base {
    // Cell (Bohr)
    static constexpr double LX_BOHR = 18.4 * ANG_TO_BOHR;   // 34.7710
    static constexpr double LY_BOHR = 18.4 * ANG_TO_BOHR;   // 34.7710
    static constexpr double LZ_BOHR = 31.7 * ANG_TO_BOHR;   // 59.9043

    // DFT / SCF
    static constexpr int    EXTRA_STATES   = 8;
    static constexpr double CUTOFF_HA      = 40.0;
    static constexpr double SCF_TOL_HA     = 1.0e-6;
    static constexpr int    SCF_MAX_STEPS  = 1000;
    static constexpr int    SCF_MIX_NDIM   = 8;
    static constexpr double SCF_MIX_ALPHA  = 0.1;

    // Wave packet (paper Eq. 1)
    static constexpr double WP_SIGMA_BOHR  = 0.53 * ANG_TO_BOHR;   // 1.0015
    static constexpr double WP_OFFSET_BOHR = 6.35 * ANG_TO_BOHR;   // 12.000
    static constexpr double WP_EKIN_EV     = 200.0;
    static constexpr double WP_K0          = k0_from_ev(WP_EKIN_EV);

    // WP centre and direction
    static constexpr double WP_CX_BOHR = 0.0;
    static constexpr double WP_CY_BOHR = 0.0;
    static constexpr double WP_CZ_BOHR = +WP_OFFSET_BOHR;
    static constexpr double WP_KX = 0.0;
    static constexpr double WP_KY = 0.0;
    static constexpr double WP_KZ = -WP_K0;
    // Magnitude alias for log lines that don't want the sign
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // Real-time propagation (dt fixed at the paper value; N_STEPS is derived
    // from physics — see compute_n_steps above).
    static constexpr double DT_AU             = 0.020;
    static constexpr int    N_STEPS           = compute_n_steps(
        LZ_BOHR, WP_OFFSET_BOHR, WP_SIGMA_BOHR, WP_K0, DT_AU);
    static constexpr int    WRITE_EVERY       = 10;
    static constexpr int    SCREEN_SNAP_EVERY = 30;

    // LEED accumulator. Per-screen physics-derived windows are computed in
    // run_template.hpp from compute_screen_window(); the paper window is the
    // single global window kept for paper-figure comparison.
    static constexpr int    N_SCREENS = 20;
    static constexpr double T1_FS = 0.077;
    static constexpr double T2_FS = 0.25;
    static constexpr double T1_AU = T1_FS * FS_TO_AU;
    static constexpr double T2_AU = T2_FS * FS_TO_AU;

    // How many sigmas of the WP envelope are considered as the "extent"
    // for screen-window timing. n=1 reproduces the Phase-3 logic
    // (centroid +/- 1*sigma); for these runs we use n=2 so the full WP
    // envelope has cleared the screen before integration begins (forward)
    // and the rebound trailing edge has more headroom (back).
    // compute_screen_window scales sigma by this factor when computing
    // t_start (forward + back) and t_end (back).
    static constexpr double WP_ENVELOPE_SIGMAS = 2.0;
};

}  // namespace coronene::config
