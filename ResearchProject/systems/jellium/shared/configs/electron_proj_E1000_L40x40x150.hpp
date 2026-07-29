// ============================================================================
// shared/configs/electron_proj_E1000_L40x40x150.hpp
//
// Two Cfgs sharing one elongated 40 x 40 x 150 Bohr orthorhombic-periodic
// jellium bath (N=162, dx=0.30 Bohr, dt=0.005 a.u., N_STEPS=2800):
//
//   Electron_Proj_E1000_L40x40x150_WP        — Gaussian wave-packet projectile
//   Electron_Proj_E1000_L40x40x150_Classical — classical-electron projectile
//                                              (custom UPF + mass override)
//
// Both projectiles enter at (20, 20, 25) Bohr in the user's corner-origin
// convention, which converts to (0, 0, -50) in INQ's centred-Cartesian
// coordinates. They share kinetic energy 1000 eV (= 36.75 Ha):
//   - WP: k_0 = sqrt(2 * 36.75) = 8.5732 bohr^-1, σ = 5.0 Bohr.
//   - Classical: v_z = sqrt(2 * 36.75) = 8.5732 bohr/atu, m = m_e (1.0 a.u.).
// Velocity match is automatic: same KE, same mass ⇒ same v_group = 8.5732.
//
// Open issue (flagged in docs/plans/the-objective-in-this-dapper-moon.md):
// at 40x40x150 with N=162, r_s = (3 V / (4 pi N))^(1/3) = 7.07 Bohr — much
// more dilute than the L=50 N=162 base run (r_s = 5.69). Direct comparisons
// of plasmon ω_p, k_F, etc. with prior journal entries will need re-derivation
// at this density.
//
// Grid choice: dx=0.30 is dictated by Nyquist for k_0=8.57:
//   k_max = k_0 + 3 σ_k = 8.57 + 3·(1/5) = 9.17 Bohr^-1
//   dx_max = π / k_max = 0.343 Bohr
// dx=0.30 gives k_Nyquist = π/0.30 = 10.47 Bohr^-1, ~14 % margin above k_max.
//
// Time choice: WP/projectile travels v · t = 8.57 · 14 = 120 Bohr in 14 a.u.
// ⇒ z_final ≈ -50 + 120 = +70 Bohr, leaving 5 Bohr to the +z face. dt=0.005
// (4× tighter than the standard 0.020) is needed because the WP/classical
// projectile traverses many grid points per step at v=8.57; rule of thumb
// dt · v < dx for stability ⇒ dt < 0.30/8.57 ≈ 0.035; we use 0.005 for safety.
// 14 a.u. / 0.005 = 2800 steps.
// ============================================================================
#pragma once

#include "base.hpp"

namespace jellium::config {

// Base for both Cfgs — all bath, grid, WP/projectile parameters live here.
struct Common_E1000_L40x40x150 : Base {

    // ----- Cell (orthorhombic, fully periodic) ---------------------------
    static constexpr double LX_BOHR = 40.0;
    static constexpr double LY_BOHR = 40.0;
    static constexpr double LZ_BOHR = 150.0;
    static constexpr double L_BOHR  = LX_BOHR;     // legacy alias for log
                                                   // strings; the actual cell
                                                   // construction in run.cpp
                                                   // uses LX/LY/LZ directly.

    // ----- Bath -----------------------------------------------------------
    static constexpr int    N_ELECTRONS    = 162;          // r_s ≈ 7.07 in this volume
    static constexpr int    EXTRA_STATES   = 20;
    static constexpr double SPACING_BOHR   = 0.30;
    static constexpr double SCF_TOL_HA     = 1.0e-6;
    // (TEMPERATURE_EV, SCF_MAX_STEPS, SCF_MIX_NDIM, SCF_MIX_ALPHA inherited
    //  from Base — same defaults the L=50 N=162 GS used.)

    // ----- Real-time -----------------------------------------------------
    static constexpr double DT_AU       = 0.005;
    static constexpr int    N_STEPS     = 2800;            // total t = 14.0 a.u.
    static constexpr int    WRITE_EVERY = 8;               // 350 frames

    // ----- Wave-packet projectile (E_kin = 1000 eV) ---------------------
    static constexpr double WP_EKIN_EV      = 1000.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 8.5732
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_SIGMA_BOHR = 5.0;
    static constexpr double WP_SIGMA_ANG  = WP_SIGMA_BOHR / ANG_TO_BOHR;

    // Launch position. User spec (20, 20, 25) is in corner-origin coords;
    // INQ uses centred Cartesian, so subtract L/2 from each component:
    //   (20, 20, 25) - (20, 20, 75) = (0, 0, -50)
    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR = -50.0;            // 25 bohr from -z face

    // ----- Classical-electron projectile (same KE, m = m_e) -------------
    // The species is a copy of anti-proton-ONCV-1.2.upf (-1 Coulomb tail,
    // Z_valence=0) with the ion mass overridden host-side to electron mass.
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "shared/pseudopotentials/electron-ONCV-1.2.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";  // INQ needs an
                                                             // element symbol; the UPF
                                                             // and mass override pin
                                                             // down the actual physics.

    // PROJ_MASS_AMU = 1.0 / 1822.8885 ⇒ INQ stores 5.485799e-4 amu and
    // returns mass() = 1822.8885 * 5.485799e-4 = 1.0000 atomic units = m_e
    // (verified to machine precision in smoke_C2).
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;

    static constexpr double PROJ_LAUNCH_X =   0.0;
    static constexpr double PROJ_LAUNCH_Y =   0.0;
    static constexpr double PROJ_LAUNCH_Z = -50.0;         // matches WP launch

    // v_z = sqrt(2 * KE / m) with m=1, KE=36.75 Ha ⇒ 8.5732 bohr/atu.
    // Same numerical value as WP_K0 by KE+mass identity.
    static constexpr double PROJ_VEL_X =   0.0;
    static constexpr double PROJ_VEL_Y =   0.0;
    static constexpr double PROJ_VEL_Z = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    // ----- Time-window placeholders (re-used by some screen tooling) ----
    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

// WP run: no ion, only the WP injection.
struct Electron_Proj_E1000_L40x40x150_WP : Common_E1000_L40x40x150 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

// Classical run: no WP, one ion with custom UPF + mass override.
struct Electron_Proj_E1000_L40x40x150_Classical : Common_E1000_L40x40x150 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

}  // namespace jellium::config
