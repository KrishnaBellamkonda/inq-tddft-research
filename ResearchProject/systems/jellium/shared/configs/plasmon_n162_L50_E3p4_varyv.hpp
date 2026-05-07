// ============================================================================
// shared/configs/plasmon_n162_L50_E3p4_varyv.hpp
//
// Run D — variable-velocity follow-up to run_plasmon_n162_L50_E15. Same box
// (L=50, N=162 closed shell), but the WP velocity is OFF the m=1 plasmon
// resonance so that the kinematic ("wrap" / rigid-WP-translation) and bath-
// plasmon predictions for the n_q_m(t) FFT peaks are no longer degenerate.
//
// Discriminator math (units: a.u. except where stated):
//   - Kinematic peak in m-th channel:  omega_kin(m) = m * v * q_1 = m * v * 2*pi/L
//   - Bohm-Gross plasmon peak:         omega_BG(q_m) = sqrt(omega_p^2 + (3/5)*v_F^2*q_m^2 + q_m^4/4)
//   These two coincide at m=1 only when v = v_res^{m=1} = omega_BG(q_1)/q_1.
//   Choosing v != v_res^{m=1} breaks the degeneracy.
//
// At L=50, N=162: omega_p = 0.1276, k_F = v_F = 0.337, q_1 = 0.1257.
//   omega_BG(q_1) = 0.1320 a.u. = 3.59 eV   v_res^{m=1} = 1.050  E_WP = 15.0 eV
//   omega_BG(q_2) = 0.1469 a.u. = 4.00 eV
//
// This Cfg picks v = 0.5 a.u. (E_WP = 3.40 eV — same energy region as the
// plasmon, well off both v_res^{m=1} = 1.05 and v_res^{m=2} = 0.585):
//   - omega_kin(m=1) at v=0.5: 0.5 * 0.1257 = 0.0628 a.u. = 1.71 eV
//   - omega_BG(q_1) (unchanged):                            3.59 eV
//   - separation 1.88 eV = 22x FFT resolution at T_sim=2000 a.u. (dE=0.086 eV)
//
//   - omega_kin(m=2) at v=0.5: 1.0  * 0.1257 = 0.1257 a.u. = 3.42 eV
//   - omega_BG(q_2) (unchanged):                            4.00 eV
//   - separation 0.58 eV = 6.7x FFT resolution
//
// The verdict: if Run D's m=1 channel peak is at 1.71 eV the previous Run B
// peak was a kinematic / rigid-WP / wrap-period artefact. If Run D's m=1
// channel peak is at 3.59 eV (or both peaks present but plasmon-dominant)
// the bath plasmon detection in Run B is real.
//
// WP shape: sigma_r = 3 (slightly tighter than the sigma=5 of previous runs).
//   sigma_k = 1/sigma = 0.333 -> backward-going Gaussian content at -k_0/sigma_k = -1.5
//   = ~6.7 % backward weight. Acceptable for a discriminator run.
//   k_max = k_0 + 3*sigma_k = 0.5 + 1.0 = 1.5 < pi (Nyquist at dx=1.0).
//
// All other parameters inherited from Plasmon_N162_L50_E15 (same propagator
// length T_sim = 2000 a.u. for matched FFT resolution).
// ============================================================================
#pragma once

#include "plasmon_n162_L50_E15.hpp"

namespace jellium::config {

struct Plasmon_N162_L50_E3p4_VaryV : Plasmon_N162_L50_E15 {
    static constexpr double WP_EKIN_EV      = 3.40;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // ≈ 0.500
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // tighter envelope per user request — sigma=3 instead of 5
    static constexpr double WP_SIGMA_BOHR = 3.0;
    static constexpr double WP_SIGMA_ANG  = 3.0 / ANG_TO_BOHR;

    // Reuse N_STEPS = 100000, WRITE_EVERY = 200 from Plasmon_N162_L50_E15
    // (T_sim = 2000 a.u., 500 density frames). For the smoke test the
    // override below sets N_STEPS = 100 (toggle for full run).
    static constexpr int    N_STEPS         = 100000;  // T_sim = 2000 a.u. ≈ 48.4 fs
    static constexpr int    WRITE_EVERY     = 200;     // ~500 density frames
};

}  // namespace jellium::config
