// ============================================================================
// shared/configs/bulk_ks_stopping_L46x46x80_sigma3.hpp
//
// SIGMA-SWEEP variant of the L46x46x80 bulk KS-stopping twin pair:
// sigma_WP = 3.0 Bohr instead of 2.0. Everything else -- box, bath, grid,
// projectile ENERGY, dt -- is held identical, so the packet width is the only
// variable. Plan: docs/plans/bulk-jellium-ks-stopping.md (user decision,
// 2026-07-31).
//
// WHY. The r_s = 5.702 and r_s = 3.987 pairs both show S_classical / S_WP ~ 6x,
// and the density lever moved it only 13 %. Width is the remaining candidate:
// a wavepacket presents a SMEARED charge whose coupling to the bath is weaker
// than a compact classical projectile's. This sweep tests that axis.
//
// ---------------------------------------------------------------------------
// COUNTER-INTUITIVE, AND THE REASON THE SWEEP IS NOT MONOTONE IN sigma:
// A NARROWER PACKET SPREADS FASTER. sigma_d(t) = sqrt(sigma^2/2 + t^2/(2 sigma^2)),
// so small sigma has a LARGER t^2 coefficient. sigma = 1 and sigma = 2 cross at
// t = 2.0 a.u.; beyond that sigma = 1 is the BROADER projectile (13.05 vs 6.67
// Bohr at t = 18.4). The width that minimises the ARRIVAL width is
// sigma_opt = sqrt(t) ~ 3-4 for this flight. So on IN-FLIGHT width the ordering
// is  sigma=1 (widest) > sigma=2 > sigma=3 (narrowest)  -- NOT the t=0 ordering.
// Any interpretation of this sweep must use in-flight sigma_d, not sigma(0).
//
// ---------------------------------------------------------------------------
// LAUNCH POSITION follows the BOUNDARY RULE z0 = -L_z/2 + 4 sigma (user
// decision, 2026-07-31, overriding an earlier plan to pin z0 = -32 for all
// sigma). Each packet therefore keeps a constant 4-sigma clearance from the
// launch face: z0 = -28.0 for sigma = 3.0. Path lengths consequently DIFFER
// between sigmas, so cross-sigma comparison must be done on a COMMON TIME
// window, not on raw path.
//
// CLEAN WINDOW for this variant: T_IFW (longitudinal, dispersion-aware)
// = 19.48, T_TRANSVERSE (periodic images) = 31.26; FIT_T1 = min = 19.48.
// The LONGITUDINAL +z face binds.
//
// THE COMMON CROSS-SIGMA WINDOW IS [4.0, 9.37] a.u. -- set by sigma = 1 in the
// 40x40x80 box, the fastest spreader in the smallest transverse cell. All
// sigma=2 and sigma=3 runs have valid data there. Verified 2026-07-31 by
// re-fitting the completed sigma=2 runs on it: S2 shifts only 5-8 %, systematics
// roughly triple but stay ~10 % of the value. USE [4, 9.37] FOR ANY CROSS-SIGMA
// COMPARISON; this header's FIT_T1 is the widest window valid for THIS run alone.
//
// ---------------------------------------------------------------------------
// LOCALISATION ENERGY T1 - T2 = 3/(4 sigma^2) = 2.27 eV at t = 0, i.e. 2 %
// of the 100 eV drift energy (sigma = 2 gave 5.10 eV = 5.1 %). The S1-vs-S2
// offset established on the density sweep scales with this, so expect the
// T1 channel to be LESS contaminated here.
//
// GRID. dx = 0.4 unchanged. k_max = k0 + 3 sigma_p = 3.418 against
// k_Nyquist = 7.854 (2.30x margin). VERIFIED numerically on this exact
// grid (2026-07-31): <p_z>, T1 and T1-T2 reproduce analytic values to MACHINE PRECISION.
//
// PSEUDOPOTENTIAL. electron_gaussian_wpsigma3p0.upf, already present in the repository (sigma_charge = 3/sqrt2 = 2.12132).
// ============================================================================
#pragma once

#include "base.hpp"
#include "boundary_rule.hpp"
#include "bulk_ks_stopping_L46x46x80.hpp"

namespace jellium::config {

struct Common_Bulk_KS_L46x46x80_sigma3 : Common_Bulk_KS_Stopping_L46x46x80 {
    static constexpr double WP_SIGMA_BOHR = 3.0;

    // Re-derived from the new sigma (the base's values are sigma = 2).
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, LZ_BOHR);              // -28.0
    static constexpr int N_STEPS =
        boundary::n_steps_for(WP_SIGMA_BOHR, LZ_BOHR, WP_K0, DT_AU);  // 600
    static constexpr double T_IFW_AU = boundary::ifw_end_t_dispersive(
        WP_SIGMA_BOHR, LZ_BOHR, WP_K0, WP_CZ_BOHR);              // 19.48
    static constexpr double T_TRANSVERSE_AU =
        boundary::transverse_clean_t(WP_SIGMA_BOHR, LX_BOHR);    // 31.26
    static constexpr double FIT_T0_AU = 4.0;
    static constexpr double FIT_T1_AU = T_IFW_AU < T_TRANSVERSE_AU
                                          ? T_IFW_AU : T_TRANSVERSE_AU;
    // The window to use when comparing ACROSS sigma (see header).
    static constexpr double COMMON_FIT_T1_AU = 9.37;

    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;

    static constexpr double PROJ_LAUNCH_Z = WP_CZ_BOHR;
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/"
        "systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma3p0.upf";
};

struct Bulk_KS_L46x46x80_sigma3_WP : Common_Bulk_KS_L46x46x80_sigma3 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Bulk_KS_L46x46x80_sigma3_Classical : Common_Bulk_KS_L46x46x80_sigma3 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

// ---- Compile-time guards --------------------------------------------------
static_assert(Common_Bulk_KS_L46x46x80_sigma3::N_STEPS == 600,
              "sigma3 L46x46x80: N_STEPS changed - re-derive the header numbers");
static_assert(Common_Bulk_KS_L46x46x80_sigma3::WP_CZ_BOHR == -28.0,
              "sigma3 L46x46x80: launch z should be -L_z/2 + 4 sigma");
// NOTE: which constraint binds DIFFERS across the sweep, so this guard asserts
// only that the window is non-empty and long enough to fit a slope - it does
// NOT assert which of the two limits is smaller (that flips with sigma).
static_assert(Common_Bulk_KS_L46x46x80_sigma3::FIT_T1_AU
                  - Common_Bulk_KS_L46x46x80_sigma3::FIT_T0_AU > 5.0,
              "sigma3 L46x46x80: fit window shorter than 5 a.u. - too short for a slope");
// The cross-sigma common window must actually be inside this run's valid window.
static_assert(Common_Bulk_KS_L46x46x80_sigma3::COMMON_FIT_T1_AU
                  <= Common_Bulk_KS_L46x46x80_sigma3::FIT_T1_AU + 1e-9,
              "sigma3 L46x46x80: the common cross-sigma window is NOT clean for this "
              "run - it extends past this run's own interference limit");
// Nyquist guard at the new (larger) momentum width.
static_assert(Common_Bulk_KS_L46x46x80_sigma3::WP_K0
                  + 3.0 / (2.0 * Common_Bulk_KS_L46x46x80_sigma3::WP_SIGMA_BOHR)
                  < 3.14159265358979
                        / Common_Bulk_KS_L46x46x80_sigma3::SPACING_BOHR,
              "sigma3 L46x46x80: grid too coarse for the WP bandwidth");

}  // namespace jellium::config
