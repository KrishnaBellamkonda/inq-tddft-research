// ============================================================================
// shared/configs/plasmon_n162_L50_E15.hpp
//
// Plasmon-detection variant of Base_N162_L50_E1p5. The only changes are
// (a) WP_EKIN_EV = 15.0 (vs 1.5) so the WP velocity matches the m=1 axial
// plasmon resonance v_res = omega(q1)/q1 ≈ 1.05 a.u., and (b) a longer
// propagation (T_sim = 2000 a.u. ≈ 48.4 fs) for FFT resolution
// dE = 2π/T_sim · Ha2eV ≈ 0.086 eV, comfortably finer than the m=1 vs m=2
// plasmon separation (≈ 0.41 eV).
//
// Predicted m=1 plasmon (Bohm-Gross dispersion ω² = ω_p² + (3/5)v_F²q² + q⁴/4):
//   q_1   = 2π/L = 0.1257 Bohr⁻¹
//   ω_p   = √(4π·n) = 0.1276 Ha,    n = 162/50³ = 1.296e-3 e/Bohr³
//   ω(q1) = 0.1320 Ha = 3.59 eV,     period T_p = 2π/ω = 47.6 a.u. = 1.15 fs
//   v_res = ω(q1)/q1 = 1.050 a.u.,   E_WP = ½v² = 0.5513 Ha = 15.0 eV  → k₀
//
// Numerical-stability sanity:
//   k_0 = 1.05;  σ_r = 5;  σ_k = 1/σ_r = 0.20
//   k_max = k_0 + 3·σ_k = 1.65 Bohr⁻¹    (Nyquist π/dx at dx=1.0 → 3.14 Bohr⁻¹, OK)
//   λ_dB = 2π/k_0 = 5.98 Bohr  ⇒ ~6 grid points per de Broglie wavelength at dx=1
//
// See `docs/plans/jellium_plasmon_detection.md` for the full validation
// table, Run A / Run B / Run C definitions, and verdict criteria.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"

namespace jellium::config {

struct Plasmon_N162_L50_E15 : Base_N162_L50_E1p5 {
    // ----- Higher-energy WP — match m=1 plasmon resonance ---------------
    static constexpr double WP_EKIN_EV      = 15.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // ≈ 1.0498
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // ----- Longer propagation for plasmon FFT resolution ---------------
    // Default: 100 000 steps at dt=0.020 → T_sim = 2000 a.u. = 48.4 fs,
    // FFT resolution ΔE = 2π/T_sim·Ha2eV ≈ 0.086 eV, well below the
    // m=1↔m=2 plasmon separation 0.41 eV.
    static constexpr int    N_STEPS         = 100000;  // T_sim = 2000 a.u. ≈ 48.4 fs
    static constexpr int    WRITE_EVERY     = 200;     // ~500 density frames
};

}  // namespace jellium::config
