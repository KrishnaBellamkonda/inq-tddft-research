// ============================================================================
// inqkit/config/tsubonoya_2014_coronene.hpp
//
// Compile-time configuration for an exact replica of the coronene (C24H12)
// electron-wavepacket scattering simulation reported in:
//
//   K. Tsubonoya, C. Hu, K. Watanabe,
//   "Time-dependent density-functional theory simulation of electron
//    wave-packet scattering with nanoflakes",
//   Phys. Rev. B 90, 035416 (2014).
//
// All parameters are in atomic units. Constants taken verbatim from the
// paper (Eqs. 1, 5, and §III), converted to a.u. as noted.
//
// Coordinate convention: INQ uses [-L/2, +L/2] for orthorhombic cells
// (see inq/src/systems/cell.hpp:212). The flake therefore sits at z=0,
// matching coronene_centred.xyz (Qball-parity geometry).
// ============================================================================
#pragma once

#include <cmath>

namespace inqkit::config::tsubonoya_2014 {

// Unit conversions (CODATA 2018)
inline constexpr double ANG_TO_BOHR = 1.8897259886;
inline constexpr double HA_TO_EV    = 27.21138625;
inline constexpr double FS_TO_AU    = 41.341374575751;  // 1 fs = 41.34... a.u.

// ---- Cell ------------------------------------------------------------------
// Paper: 18.4 x 18.4 x 31.7 A^3
inline constexpr double LX_BOHR = 18.4 * ANG_TO_BOHR;   // 34.7710
inline constexpr double LY_BOHR = 18.4 * ANG_TO_BOHR;   // 34.7710
inline constexpr double LZ_BOHR = 31.7 * ANG_TO_BOHR;   // 59.9043

// ---- DFT setup -------------------------------------------------------------
// Paper: ALDA, Troullier-Martins NC PSPs, spin-unpolarised, fixed ions.
// We use INQ's default norm-conserving pseudopotentials as a substitute.
inline constexpr int    EXTRA_STATES   = 8;
inline constexpr double CUTOFF_HA      = 54.0;       // ~54 Ha (~108 Ry, paper grid spacing 0.16 A => Δr ~ 0.30 Bohr)
inline constexpr double SCF_TOL_HA     = 1.0e-6;
inline constexpr int    SCF_MAX_STEPS  = 1000;
inline constexpr int    SCF_MIX_NDIM   = 8;
inline constexpr double SCF_MIX_ALPHA  = 0.1;

// ---- Wave-packet (Eq. 1) ---------------------------------------------------
// Paper: d = 0.53 A, b at 6.35 A above the flake, perpendicular to the
//        flake plane (i.e. along z). E_kin = 200 eV. Direction toward flake.
inline constexpr double WP_SIGMA_BOHR  = 0.53 * ANG_TO_BOHR;   // ~1.0015
inline constexpr double WP_OFFSET_BOHR = 6.35 * ANG_TO_BOHR;   // ~12.000
inline constexpr double WP_EKIN_EV     = 200.0;
inline constexpr double WP_EKIN_HA     = WP_EKIN_EV / HA_TO_EV; // ~7.3499 Ha
inline const     double WP_K0          = std::sqrt(2.0 * WP_EKIN_HA); // ~3.834 Bohr^-1

// WP centre: directly above the flake at z = 0 (centred-cell convention).
// Direction of motion: -z (toward flake at z=0 then onward to z<0).
inline constexpr double WP_CX_BOHR = 0.0;
inline constexpr double WP_CY_BOHR = 0.0;
inline constexpr double WP_CZ_BOHR = +WP_OFFSET_BOHR;

inline constexpr double WP_KX = 0.0;
inline constexpr double WP_KY = 0.0;
inline const     double WP_KZ = -WP_K0;

// ---- Real-time propagation (Eq. 5) -----------------------------------------
// Paper: Δt = 4.84e-4 fs = 0.020 a.u. for the LEED window.
//        10000 steps total (LEED phase).
inline constexpr double DT_AU                = 0.020;            // ~4.84e-4 fs
// 600 steps = 12.0 a.u. = 0.290 fs. WP traverses from z=+12 Bohr to
// z=-29.95 Bohr (cell -Lz/2) at v=|k|=3.834 a.u., which takes ~10.94 a.u.;
// 600 steps is the WP transit time + small margin.
inline constexpr int    N_STEPS              = 600;
inline constexpr int    WRITE_EVERY          = 10;               // density frame cadence
inline constexpr int    SCREEN_SNAP_EVERY    = 30;               // per-step screen snapshot cadence
inline constexpr int    N_SCREENS            = 20;

// LEED time-integration window (paper Eq. 5):
// t1 = 0.077 fs (WP arrival at flake), t2 = 0.25 fs (WP at box boundary).
inline constexpr double T1_FS = 0.077;
inline constexpr double T2_FS = 0.25;
inline constexpr double T1_AU = T1_FS * FS_TO_AU;   // ~3.183
inline constexpr double T2_AU = T2_FS * FS_TO_AU;   // ~10.335

// ---- Screen positions ------------------------------------------------------
// 20 z-positions spread across [-Lz/2, +Lz/2]. Boundaries close to ±Lz/2
// minus a small margin; 18 interior positions with a small deterministic
// jitter so multiple screens never coincide on a grid plane.
//
// Note: cell origin is at -Lz/2; +Lz/2 = +29.95 Bohr in our centred frame.
// Layout convention parallels jellium-wp-rt/run_01_base/run.cpp but
// rescaled to the centred cell.

}  // namespace inqkit::config::tsubonoya_2014
