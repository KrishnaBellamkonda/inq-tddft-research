// ============================================================================
// shared/configs/electron_proj_E100_L50_cubic_sigma1.hpp  — concentrated-WP
// companion to electron_proj_E100_L50_cubic_v2 (user request 2026-05-17).
//
// SAME bath (L=50 cubic, N=162, dx=0.40) and SAME WP energy (100 eV), but
// the WP is much more spatially concentrated (σ=1 Bohr instead of σ=5 Bohr).
//
// The concentrated WP probes a fundamentally different limit of the
// projectile-bath coupling: with σ_density(0) = 1/√2 = 0.707 Bohr the WP is
// smaller than the bath's inter-electron spacing (r_s = 5.69 Bohr), so it
// behaves more like a point projectile in real space (smaller k-space
// uncertainty, σ_p = 1/(σ√2) = 0.707 Bohr⁻¹). Contrast σ=5: density
// σ = 3.54 Bohr, σ_p = 0.141 Bohr⁻¹ — five times larger spatial and
// five times smaller momentum spread.
//
// SELF-SPREAD CONSTRAINT (the binding bound on N_STEPS for σ=1):
// Analytic free-particle σ_density(t) = (σ/√2) · √(1 + (t/σ²)²). For σ=1,
// σ_density(9.5 a.u.) ≈ 6.75 Bohr. With centroid at -21 + 2.7111·9.5 =
// +4.75 Bohr, the 3σ_density(t) tail reaches +4.75 + 20.3 ≈ +25.0 — the
// far box face. Propagation beyond t ≈ 9.5 a.u. produces wrap-around
// contamination on the density second moment (verified empirically in
// the Run-2b Python toy bug-hunt).
//
// We therefore cap N_STEPS at 475 (t_total = 9.50 a.u.), distinct from
// the boundary_rule.n_steps_for value of 830. The corresponding IFW end
// is set by self-spread, NOT the centroid criterion:
//
//   t_IFW_self = 9.5 a.u.    (where 3σ_density(t) reaches the far face)
//
// vs. the boundary-rule centroid-criterion value
//
//   t_IFW_centroid = (ifw_end_z - launch_z) / v = (+25 - 3·1 - (-21)) / 2.711
//                  = 43 / 2.711 = 15.86 a.u.
//
// The binding bound is whichever is smaller — at σ=1 the self-spread wins
// at t ≈ 9.5 a.u., well before the centroid criterion fires.
//
// Nyquist: σ_p = 0.707 ⇒ k_max_WP = k0 + 3σ_p = 2.71 + 2.12 = 4.83 Bohr⁻¹.
// k_Nyquist at dx=0.40 = 7.85 Bohr⁻¹. Ratio 0.61 — comfortable.
//
// GS reuse: same checkpoints/gs_L50_cubic_N162_dx0p40 as Run-3 v2; no new
// GS needed.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_E100_L50_sigma1 : Base_N162_L50_E1p5 {
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;

    static constexpr double DT_AU        = 0.020;

    static constexpr double WP_EKIN_EV      = 100.0;
    static constexpr double WP_SIGMA_BOHR   = 1.0;      // override 5.0 from base
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV); // 2.7111
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -21

    // SELF-SPREAD-CAPPED N_STEPS (see header note above). The boundary_rule
    // ifw_end criterion would allow 793 steps; the actual safe window is
    // 475 steps because σ_density(t) blows up much faster at σ=1.
    static constexpr int N_STEPS     = 475;            // t_total = 9.5 a.u.
    static constexpr int WRITE_EVERY = 2;              // 237 frames

    // Classical-electron companion (matched KE, m = m_e).
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;
    static constexpr double PROJ_LAUNCH_X = 0.0;
    static constexpr double PROJ_LAUNCH_Y = 0.0;
    static constexpr double PROJ_LAUNCH_Z = WP_CZ_BOHR;                // -21
    static constexpr double PROJ_VEL_X    = 0.0;
    static constexpr double PROJ_VEL_Y    = 0.0;
    static constexpr double PROJ_VEL_Z    = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E100_L50_sigma1_WP_dx0p40 : Common_E100_L50_sigma1 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Electron_Proj_E100_L50_sigma1_Classical_dx0p40 : Common_E100_L50_sigma1 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
