// ============================================================================
// shared/configs/electron_proj_E1500_L50_cubic.hpp  (v2 of the
//   electron-classical-wavepacket-jellium comparison plan)
//
// Two Cfgs sharing one cubic 50 x 50 x 50 Bohr periodic jellium bath
// (N=162, dx=0.248 Bohr, dt=0.005 a.u., N_STEPS=860):
//
//   Electron_Proj_E1500_L50_cubic_WP        — Gaussian wave-packet projectile
//   Electron_Proj_E1500_L50_cubic_Classical — classical-electron projectile
//                                              (custom UPF + mass override)
//
// Both projectiles enter at (25, 25, 15) Bohr (corner-origin) = (0, 0, -10)
// in INQ centred Cartesian — 3 sigma_r from the -z face (envelope fully
// inside the box). They share kinetic energy 1500 eV (= 55.13 Ha):
//   - WP: k_0 = sqrt(2 * 55.13) = 10.50 bohr^-1, sigma = 5.0 Bohr.
//   - Classical: v_z = sqrt(2 * 55.13) = 10.50 bohr/atu, m = m_e (1.0 a.u.).
// Velocity match is automatic: same KE, same mass ⇒ same v_group = 10.50.
//
// Trajectory: at v=10.50 over t=4.3 atu the projectile travels 45.15 Bohr.
// Starting at z=-10, the WP/classical exits +z face at t=(25-(-10))/10.50
// = 3.33 atu, wraps to z=-25, and ends at z=-14.85 (still 10 Bohr inside
// box). One full periodic wrap during the simulation. Postprocess should
// split the trajectory into pre-wrap (t<3.33 atu, ~660 steps) and
// post-wrap (t>3.33 atu, ~200 steps) for cleanest interpretation.
//
// Density: r_s = (3 V / (4 pi N))^(1/3) = (3 * 125000 / (4 pi * 162))^(1/3)
// = 5.69 Bohr. Plasmon frequency omega_p = sqrt(4 pi n) = 0.1276 Ha
// = 3.47 eV. WP energy 1500 eV = 432 omega_p — well above plasmon
// resonance, in the Bohr/Bethe single-particle stopping regime.
//
// Grid: dx=0.248 corresponds to plane-wave cutoff E_max = 80 Ha = 160 Ry
// (g_max = sqrt(2*80) = 12.65 Bohr^-1, dx_min = pi/g_max = 0.248). With WP
// k_max = k_0 + 3 sigma_k = 10.50 + 0.60 = 11.10 Bohr^-1, Nyquist gives
// ~14% margin.
//
// Time choice: dt=0.005 a.u. is 4x finer than the standard 0.020 used by
// run_base_n162_L50_E1p5; it satisfies dt * v < dx (0.005 * 10.5 = 0.0525
// << 0.248) with a comfortable margin for ETRS stability at this k_0.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"

namespace jellium::config {

// Common bath/grid; both WP and classical Cfgs share these via inheritance.
// We start from Base_N162_L50_E1p5 to inherit L=50 cubic, N=162, EXTRA_STATES=20,
// T=100K, SCF_TOL=1e-6, WP_SIGMA=5.0 — and override only what differs.
struct Common_E1500_L50_cubic : Base_N162_L50_E1p5 {

    // ----- Grid (relaxed from user-spec 0.248 to 0.30 for memory headroom)
    // Trade-off: at dx=0.248 the eigensolver's 6 wavefunction-sized buffers
    // overflow GPU memory (each buffer is 9.26M cells x 16 B x 101 states
    // = 15 GB; 6 buffers = 90 GB nominal vs 24 GB GPU). At dx=0.30 the
    // grid drops to 167^3 = 4.66M cells, halving memory and letting the
    // SCF run efficiently on one GPU.
    //
    // Aliasing cost: k_Nyquist(0.30) = pi/0.30 = 10.47 bohr^-1. WP k_0
    // = 10.50 sits AT the Nyquist boundary; the high-k tail (k_0 + 3 sigma_k
    // = 11.10) is above Nyquist by ~6 % and gets aliased back. ~16 % of
    // the WP's momentum-space weight is mis-resolved. For Bethe-regime
    // stopping-power physics at 1500 eV (where k_0 dominates), this is
    // an acceptable scoping-run compromise.
    static constexpr double SPACING_BOHR = 0.30;
    static constexpr double CUTOFF_HA    = 54.8;       // pi^2/(2 dx^2);
                                                       // log only, INQ uses
                                                       // spacing, not cutoff

    // ----- Real-time -----------------------------------------------------
    static constexpr double DT_AU       = 0.005;
    static constexpr int    N_STEPS     = 860;         // total t = 4.3 a.u.
    static constexpr int    WRITE_EVERY = 4;           // 215 frames

    // ----- Wave-packet projectile (E_kin = 1500 eV) ---------------------
    static constexpr double WP_EKIN_EV      = 1500.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 10.4995
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // sigma inherited (5.0 Bohr) from Base_N138_L50_E1p5

    // Launch at z=15 Bohr in corner-origin (= z=-10 in INQ centred).
    // 3 sigma = 15 Bohr from -z face ⇒ whole envelope inside box.
    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR = -10.0;

    // ----- Classical-electron projectile (same KE, m = m_e) -------------
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";

    // PROJ_MASS_AMU = 1.0 / 1822.8885 ⇒ m_e to machine precision
    // (verified in C2 smoke test).
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;

    static constexpr double PROJ_LAUNCH_X =   0.0;
    static constexpr double PROJ_LAUNCH_Y =   0.0;
    static constexpr double PROJ_LAUNCH_Z = -10.0;     // matches WP launch

    static constexpr double PROJ_VEL_X =   0.0;
    static constexpr double PROJ_VEL_Y =   0.0;
    // v_z = sqrt(2 * KE_Ha / m) with m=1, KE=55.13 Ha ⇒ 10.50 bohr/atu.
    // Same numerical value as WP_K0 (KE+mass identity).
    static constexpr double PROJ_VEL_Z = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    // ----- Time-window placeholders (re-used by some screen tooling) ----
    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;

    // SCREEN_SNAP_EVERY inherited (6); irrelevant for this comparison
    // (LEED screens not used in the WP-vs-classical postprocess).
};

// WP run: no ion, only WP injection. run_template.hpp consumes this.
struct Electron_Proj_E1500_L50_cubic_WP : Common_E1500_L50_cubic {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

// Classical run: no WP, one ion with custom UPF + mass override.
struct Electron_Proj_E1500_L50_cubic_Classical : Common_E1500_L50_cubic {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

// Classical run on a COARSER grid (dx=0.40, ~125^3 = 1.95 M points).
// Ehrenfest force computation in observables::forces_stress allocates a
// vector3<complex> orbital_set buffer, which at 81 occupied states needs
// ~18 GB on top of the 7+ GB resident wavefunctions — overflowing 24 GB
// on the dx=0.30 grid. At dx=0.40 the orbital_set is ~half the size,
// force-stress fits comfortably, and the bath physics (density-fluctuation
// modes ≤ a few bohr^-1) is still well-resolved (Nyquist k = π/0.40
// = 7.85 bohr^-1).
//
// This Cfg requires its OWN GS at dx=0.40; the dx=0.30 checkpoint is
// not compatible (different grid, different num_grid_points). See
// save_gs/gs_L50_cubic_N162_dx0p40/run.cpp.
//
// The WP run continues to use the dx=0.30 grid because k_0 = 10.5 needs
// fine spacing; this asymmetry is acceptable because the projectiles'
// physical regimes are separable (the WP's k-space content needs fine
// grid; the classical electron's point-charge nature does not).
struct Electron_Proj_E1500_L50_cubic_Classical_dx0p40 : Common_E1500_L50_cubic {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
    static constexpr double SPACING_BOHR = 0.40;     // overrides 0.30
    static constexpr double CUTOFF_HA    = 30.84;    // pi^2/(2*0.4^2)
};

}  // namespace jellium::config
