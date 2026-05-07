// ============================================================================
// shared/configs/base_n138_L30_E5.hpp  (jellium, closed-shell, low-energy WP)
//
// Base_N138_L30_E5 inherits from Base and overrides the cell, grid spacing,
// WP kinetic energy, and propagation length to set up the canonical
// e-h-coupling regime defined in `.claude/rules/jellium-base-run-spec.md`.
// This is the new project base configuration as of 2026-05-04.
//
// === Why these numbers ===
//
//   * L_BOHR = 30 (vs Base 60). Half the box. The WP travels half the
//     distance in the trajectory; the kinetic-only e-h transition gap
//     |G|^2=6 -> |G|^2=8 is dE = 1/2 * (2pi/L)^2 * 2 = 0.04386 Ha = 1.19 eV
//     at L=30 (vs 0.30 eV at L=60), placing the lowest unoccupied shell
//     in the same scale as the 5 eV WP.
//
//   * N_ELECTRONS = 138 (closed shell |G|^2 <= 6). Same shell as the
//     L=60 base. Density at L=30 is 138/27000 = 5.11e-3 e/bohr^3
//     (r_s = 3.61 bohr, sodium-like); compare 6.39e-4 e/bohr^3 at L=60
//     (r_s = 7.21). The closed-shell symmetry is preserved; only the
//     density changes.
//
//   * EXTRA_STATES = 20 (vs Base 4). With kT ~ 0.0086 eV at T=100 K,
//     Fermi smearing only penetrates 1-2 states; 20 extras cover the
//     |G|^2 = 8 (12-fold) and |G|^2 = 9 (24-fold) shells comprehensively
//     so the gamma_transitions postprocess builds a non-empty occ -> unocc
//     set, and the no-WP KS-energy diagnostics see the bath redistribute
//     across many low-lying transitions.
//
//   * WP_EKIN_EV = 5. Resonant with the lowest e-h transitions. k_0 =
//     sqrt(2 * 5 / 27.211) = 0.6063 Bohr^-1 (vs 2.711 at 100 eV; 4.47x slower).
//
//   * SPACING_BOHR = 0.85. Nyquist:
//        sigma_k = 1 / sigma_r   = 1 / 1.0015     = 0.999  Bohr^-1
//        k_max   = k_0 + 3 sigma_k = 0.606 + 2.996 = 3.602  Bohr^-1
//        Required dx <= pi / k_max = 0.872 Bohr
//        Chosen dx  = 0.85          (k_Nyquist = pi/0.85 = 3.696 Bohr^-1;
//                                    2.6 % above k_max)
//
//   * N_STEPS = 990 at dt = 0.020 a.u. -> 19.80 a.u. (0.479 fs) of
//     propagation. WP travels k_0 * dt * N_STEPS = 0.6063 * 0.020 * 990
//     = 12.0 Bohr from the origin, leaving ~3 Bohr clear of L/2 = 15
//     even with the 3 sigma envelope. THIS DURATION IS THE PROJECT
//     STANDARD (see .claude/rules/jellium-base-run-spec.md).
//
//   * WP centre = (0,0,0) (box centre in INQ centred-Cartesian frame).
//   * SCF_TOL_HA = 1e-6 (tightened, project-wide).
// ============================================================================
#pragma once

#include "base.hpp"

namespace jellium::config {

struct Base_N138_L30_E5 : Base {
    static constexpr double L_BOHR  = 30.0;
    static constexpr double LX_BOHR = L_BOHR;
    static constexpr double LY_BOHR = L_BOHR;
    static constexpr double LZ_BOHR = L_BOHR;

    static constexpr int    N_ELECTRONS    = 138;
    static constexpr int    EXTRA_STATES   = 20;
    static constexpr double SCF_TOL_HA     = 1.0e-6;
    static constexpr double SPACING_BOHR   = 0.85;

    // WP at 5 eV → k₀ = √(2 · 5 / 27.21138625) ≈ 0.60633 Bohr⁻¹
    static constexpr double WP_EKIN_EV      = 5.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // WP centre = box centre = (0,0,0) (inherited from Base, made explicit)
    static constexpr double WP_CX_BOHR = 0.0;
    static constexpr double WP_CY_BOHR = 0.0;
    static constexpr double WP_CZ_BOHR = 0.0;

    // Real-time. 990 steps × 0.020 a.u. = 19.80 a.u. — the project
    // standard duration; see .claude/rules/jellium-base-run-spec.md.
    static constexpr double DT_AU             = 0.020;
    static constexpr int    N_STEPS           = 990;
    static constexpr int    WRITE_EVERY       = 2;
    static constexpr int    SCREEN_SNAP_EVERY = 6;

    static constexpr int    N_SCREENS = 20;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = DT_AU * N_STEPS / FS_TO_AU;
    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
};

}  // namespace jellium::config
