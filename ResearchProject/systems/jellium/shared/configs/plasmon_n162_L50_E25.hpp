// ============================================================================
// shared/configs/plasmon_n162_L50_E25.hpp
//
// Run E — 25 eV spectroscopy companion to run_plasmon_n162_L50_E15 (15 eV)
// and run_plasmon_n162_L50_E3p4_varyv (3.4 eV). Same box / GS / propagator
// length; only WP_EKIN_EV (-> 25.0) and WP_SIGMA_BOHR (-> 3, matching E3.4)
// change. Purpose: an INDEPENDENT, frequency-resolved loss function L(q,w)
// at a third projectile velocity so the medium-property assumption can be
// tested directly (M8: "how different are the loss functions?"). L(q,w) is a
// property of the jellium, not the projectile; if linear response holds the
// E15/E3.4/E25 loss functions should overlay.
//
// WHY a long, COARSE-grid run (not the fine-grid E25 transit run):
//   - Frequency resolution dE = 2*pi/T_sim is set ONLY by total propagation
//     time. T_sim = 2000 a.u. -> dE = 0.086 eV (resolves omega_p = 3.47 eV
//     ~40x over). The existing E25 transit run is T = 11 a.u. -> dE = 15.5 eV
//     (< 1 bin below the plasmon) and is UNUSABLE for spectroscopy.
//   - dt is deliberately NOT decreased: at fixed N_STEPS, a smaller dt SHORTENS
//     T_sim and WORSENS dE. Resolution comes from long T_sim, so we keep the
//     proven stable dt = 0.020 (inherited) and N_STEPS = 100000 -> T = 2000 a.u.
//   - The slow/wide WP (sigma=3) lets dx=1.0 (50^3 grid) satisfy Nyquist, so
//     100000 steps cost ~6.8 h wall (matched to E15/E3.4: wall_time_s=24571).
//     A fine-grid (dx=0.4, 125^3) run of equal duration would be ~16x costlier.
//
// Nyquist (dx = 1.0 inherited):
//   k_0 = k0_from_ev(25) = sqrt(2*25/27.2114) = 1.3555 Bohr^-1
//   sigma_k = 1/sigma = 1/3 = 0.3333 ;  k_max = k_0 + 3*sigma_k = 2.355 Bohr^-1
//   Nyquist k at dx=1.0 is pi = 3.1416 Bohr^-1 > 2.355  -> OK
//   backward-going content at k=0 sits at -k_0/sigma_k = -4.07 sigma -> negligible
//
// Cadence: WRITE_EVERY = 10 (vs 200 for E15/E3.4) for HIGH-cadence scalar
//   observables (momentum stats / energies / occupations at 0.2 a.u. spacing,
//   10000 samples) to support the momentum/scattering analysis. This does NOT
//   affect spectral resolution (dE fixed by T_sim) and only RAISES the field
//   Nyquist ceiling to pi/(dt*WRITE_EVERY) = pi/0.2 = 15.7 Ha = 427 eV (well
//   above the e-h continuum). Field-VTI disk ~ 80 GB (10000 frames, 50^3) on
//   2.7 TB free. The E15/E3.4 loss functions (WRITE_EVERY=200) are recovered
//   identically over 0-21 eV, so the M8 comparison remains apples-to-apples.
//
// Predicted kinematic vs plasmon discriminator at v(25 eV) = 1.3555:
//   omega_kin(m=1) = v*q_1 = 1.3555 * 0.1257 = 0.1704 a.u. = 4.64 eV
//   omega_BG(q_1)  = 3.59 eV (unchanged medium plasmon)
//   -> the kinematic and plasmon peaks are well separated (1.05 eV = 12x dE).
// ============================================================================
#pragma once

#include "plasmon_n162_L50_E15.hpp"

namespace jellium::config {

struct Plasmon_N162_L50_E25 : Plasmon_N162_L50_E15 {
    static constexpr double WP_EKIN_EV      = 25.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // ~ 1.3555
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // sigma=3 matches E3.4 (clean WP-shape-controlled comparison); Nyquist OK at dx=1.0.
    static constexpr double WP_SIGMA_BOHR = 3.0;
    static constexpr double WP_SIGMA_ANG  = 3.0 / ANG_TO_BOHR;

    // T_sim = 2000 a.u. (dE = 0.086 eV) — matched to E15/E3.4 for FFT resolution.
    static constexpr int    N_STEPS         = 100000;  // T_sim = 2000 a.u. ~ 48.4 fs
    // High-cadence scalars for momentum/scattering (20x finer than E15/E3.4).
    static constexpr int    WRITE_EVERY     = 10;      // ~10000 frames, ~80 GB
};

}  // namespace jellium::config
