// ============================================================================
// shared/configs/bulk_ks_stopping_L46x46x80.hpp
//
// Bulk-jellium KS-orbital stopping power: a matched classical + wavepacket TWIN
// PAIR of 100 eV electron projectiles in a fully periodic orthorhombic jellium
// bath. Plan: docs/plans/bulk-jellium-ks-stopping.md (user design lock,
// 2026-07-30).
//
//   Bulk_KS_Stopping_L46x46x80_WP         — Gaussian wave-packet projectile
//   Bulk_KS_Stopping_L46x46x80_Classical  — classical electron projectile
//                                           (Gaussian UPF + m_e mass override)
//
// PURPOSE. Extract S = -dT/ds from FOUR KS-orbital-dependent definitions:
//   T1 = <p^2>/2m   (full orbital kinetic energy; WPMomentumStats e_kin_ha)
//   T2 = <p>^2/2m   (drift only; discards localisation + scattering energy)
//   s3 = density centroid of the WP KS orbital  (WPRealSpaceStats, circular)
//   s4 = integral of <p_z> dt                   (post-processed from pz_mean)
// plus the classical reference S from the Ehrenfest ion trajectory.
//
// T1 - T2 = (3/2) sum_d sigma_pd^2 = 3/(8 sigma^2) = 2.55 eV at t=0 for
// sigma = 2. Its CHANGE along the trajectory is the momentum-broadening
// (angular-scattering) contribution to apparent stopping — the physics contrast
// this run is built to measure.
//
// NOTE on s3 vs s4: the WP run has NO ions, so the KS Hamiltonian is purely
// local (kinetic + Hartree + ALDA) and Ehrenfest gives d<z>/dt = <p_z>/m
// EXACTLY. There is no CAP here either (contrast the qsp5 runs, where CAP
// non-unitarity broke the identity at t ~ 5). s3 and s4 must therefore agree to
// numerical precision: their comparison is a VALIDATION CHECK, not a second
// independent physics channel.
//
// ---------------------------------------------------------------------------
// GEOMETRY DERIVATION (why 46 x 46 x 80 and not the originally-proposed
// 35 x 35 x 110 — see the plan section 3 for the full working)
//
// A free Gaussian spreads as sigma_d(t) = sqrt(sigma^2/2 + t^2/(2 sigma^2)).
// At sigma = 2 Bohr the packet grows 1.41 -> 7.2 Bohr over the flight, so TWO
// independent constraints bound the interference-free window:
//
//   LONGITUDINAL  z0 + v t + 3 sigma_d(t) < L_z/2      ->  t = 18.97 a.u.
//   TRANSVERSE    L_xy >= 6 sigma_d(t)                 ->  t = 21.31 a.u. at L_xy=46
//
// (both computed by boundary_rule.hpp::ifw_end_t_dispersive / transverse_clean_t
// and pinned there by static_assert). L_xy = 46 puts the transverse limit safely
// PAST the longitudinal one, so L_z is the binding constraint and widening the
// box further would buy nothing. L_xy = 35 would have bound first, at 16.0 a.u.
//
// The minimum achievable arrival width is sigma_d,min = sqrt(t) at sigma = sqrt(t),
// so shortening L_z (110 -> 80) is the only lever with unbounded returns.
//
// WARNING — the STATIC rule is wrong here. boundary_rule.hpp::t_ifw_au (the
// legacy +L/2 - 3*sigma_launch form) claims 24.34 a.u. for this geometry: a 22 %
// over-estimate, because it ignores spreading entirely. Use the dispersive
// helper for any run whose flight is long compared with sigma^2.
//
// ---------------------------------------------------------------------------
// BATH. N = 218 electrons in 46 x 46 x 80 = 169,280 Bohr^3 gives
// n = 1.2878e-3 e/Bohr^3, r_s = 5.702 Bohr — essentially the r_s = 5.69 of the
// legacy L=50 / N=162 cubic runs, so this run is directly comparable with that
// body of work. N = 218 sits in a 0.24 eV free-electron shell-closure gap for
// THIS box (enumerated numerically, 2026-07-30), far above the 100 K smearing,
// so the SCF should give exactly 109 doubly-occupied states and nothing above.
//
// hbar omega_p = 3.462 eV, v_F = 0.3366, E_F = 0.05664 Ha. v/v_F = 8.05 — well
// above the Fermi velocity, plasmon excitation kinematically allowed.
//
// KNOWN LIMITATION (stated in every downstream notebook): 2 pi / omega_p =
// 49.4 a.u. exceeds the whole 25.8 a.u. run. The bath cannot complete one plasma
// oscillation, so the extracted S is an INITIAL-DRAG stopping power, not a
// converged steady-state S(v). This is forced: a light 100 eV electron crosses
// this box in 26 a.u. and the wake criterion (5 plasma periods) is geometrically
// unreachable. See .claude/rules/light-projectile-stopping.md.
//
// ---------------------------------------------------------------------------
// GRID. dx = 0.40 (user decision 2026-07-30) -> 115 x 115 x 200 = 2.65 M points,
// cutoff 30.84 Ha, matching every prior jellium run. Nyquist k = pi/0.40 = 7.85
// against k_max = k0 + 3 sigma_k = 2.711 + 3*(1/(2*2)) = 3.46 — 2.3x headroom.
//
// SIGMA CONVENTION (.claude/rules/sigma-wp-convention.md). sigma = 2.0 Bohr is
// the WAVEPACKET sigma: psi ~ exp(-r^2 / 2 sigma^2), so the density std is
// sigma/sqrt2 = 1.414. The classical twin uses the Gaussian UPF generated with
// generate_gaussian_psp(sigma_wp=2.0), whose CHARGE std is sigma/sqrt2 = 1.414 —
// so both projectiles present the IDENTICAL charge cloud exp(-r^2/sigma^2) to
// the bath. Both are labelled sigma = 2.
// ============================================================================
#pragma once

#include "base.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

// Base for both Cfgs — bath, grid, projectile and timing all live here so the
// twin pair cannot drift apart.
struct Common_Bulk_KS_Stopping_L46x46x80 : Base {

    // ----- Cell (orthorhombic, fully periodic in x, y AND z) --------------
    static constexpr double LX_BOHR = 46.0;
    static constexpr double LY_BOHR = 46.0;
    static constexpr double LZ_BOHR = 80.0;
    static constexpr double L_BOHR  = LZ_BOHR;   // legacy alias used in log strings
                                                 // ONLY; run.cpp builds the cell
                                                 // from LX/LY/LZ directly.

    // ----- Bath -----------------------------------------------------------
    static constexpr int    N_ELECTRONS  = 218;   // closed shell, gap 0.24 eV
    static constexpr int    EXTRA_STATES = 20;    // degenerate buffer + WP slot
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84; // pi^2/(2 dx^2), for the log
    static constexpr double SCF_TOL_HA   = 1.0e-6;

    // ----- Projectile (shared by both representations) --------------------
    static constexpr double WP_EKIN_EV    = 100.0;
    static constexpr double WP_SIGMA_BOHR = 2.0;                    // WAVEPACKET sigma
    static constexpr double WP_K0         = k0_from_ev(WP_EKIN_EV); // 2.7111
    static constexpr double WP_KX = 0.0;
    static constexpr double WP_KY = 0.0;
    static constexpr double WP_KZ = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR = 0.0;
    static constexpr double WP_CY_BOHR = 0.0;
    // Standard boundary rule: -L_z/2 + 4 sigma = -32.0
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, LZ_BOHR);

    // ----- Real time ------------------------------------------------------
    static constexpr double DT_AU = 0.040;
    // n_steps_for(sigma=2, L=80, v=2.7111, dt=0.04) = 646: centroid travels the
    // 70-Bohr traversal from -32 to stop_z = +38.
    static constexpr int N_STEPS =
        boundary::n_steps_for(WP_SIGMA_BOHR, LZ_BOHR, WP_K0, DT_AU);

    // Cadence (user decision 2026-07-30): every 2 steps = 323 density frames.
    // NOT boundary::write_every_for(646) (= 2 as it happens) — pinned explicitly
    // because the user chose the frame count, not the campaign default.
    // Propagation-callback stride. The ENERGY series (observables.csv) is
    // written from that callback, so it follows THIS, not WRITE_EVERY.
    // Kept fine (2) because the CLASSICAL stopping power is fitted to
    // energy_total(x): thinning it to the VTI cadence would cut the fit
    // window from ~68 points to ~17. Storage is saved on the VTI writes,
    // which are gated separately on WRITE_EVERY inside the callback.
    static constexpr int OBS_EVERY = 2;

    static constexpr int WRITE_EVERY    = 8;    // density_total / _wp / _delta
    static constexpr int WF_WRITE_EVERY = 32;    // complex wavefunction (108 frames)
    static constexpr int STATS_EVERY    = 1;    // momentum + real-space moments:
                                                // THE stopping data, every step

    // Interior RT checkpoint cadence (.claude/rules/checkpoint-dont-block.md).
    static constexpr int CKPT_EVERY = 200;

    // ----- Analysis window (consumed by analyse.py and the notebooks) -----
    // Dispersion-aware interference-free window; see the header comment.
    static constexpr double T_IFW_AU = boundary::ifw_end_t_dispersive(
        WP_SIGMA_BOHR, LZ_BOHR, WP_K0, WP_CZ_BOHR);            // 18.97
    static constexpr double T_TRANSVERSE_AU =
        boundary::transverse_clean_t(WP_SIGMA_BOHR, LX_BOHR);  // 21.31
    // Lower edge drops the launch/orthogonalisation transient.
    static constexpr double FIT_T0_AU = 4.0;
    static constexpr double FIT_T1_AU = T_IFW_AU < T_TRANSVERSE_AU
                                          ? T_IFW_AU : T_TRANSVERSE_AU;

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;

    // ----- Classical projectile -------------------------------------------
    // Gaussian electron UPF whose CHARGE std is WP_SIGMA_BOHR/sqrt2 = 1.414, so
    // the classical cloud is identical to the WP density at t = 0. Generated by
    // inqview.io.gaussian_psp.generate_gaussian_psp(sigma_wp=2.0) on 2026-07-30
    // and validated against the analytic erf form to 5e-11.
    static constexpr const char* PROJ_PSEUDO_PATH =
        "/rds/user/skcb2/hpc-work/tddft/inq-tddft-research/ResearchProject/"
        "systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma2p0.upf";
    static constexpr const char* PROJ_SPECIES_SYMBOL = "H";
    static constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;   // = m_e in a.u.
    static constexpr double PROJ_LAUNCH_X = 0.0;
    static constexpr double PROJ_LAUNCH_Y = 0.0;
    static constexpr double PROJ_LAUNCH_Z = WP_CZ_BOHR;        // identical launch
    static constexpr double PROJ_VEL_X    = 0.0;
    static constexpr double PROJ_VEL_Y    = 0.0;
    static constexpr double PROJ_VEL_Z    = WP_K0;             // v = k0 since m = m_e
};

struct Bulk_KS_Stopping_L46x46x80_WP : Common_Bulk_KS_Stopping_L46x46x80 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Bulk_KS_Stopping_L46x46x80_Classical : Common_Bulk_KS_Stopping_L46x46x80 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
    // Same N_STEPS as the WP twin: parity of the time axis is what makes the
    // pair comparable, and the classical projectile follows the same nominal
    // trajectory (it decelerates under free Ehrenfest, so it arrives late —
    // that deceleration IS the measurement).
};

// ---- Compile-time guards on the derived geometry ---------------------------
static_assert(Common_Bulk_KS_Stopping_L46x46x80::N_STEPS == 646,
              "bulk KS stopping: N_STEPS should be 646 at dt=0.04");
static_assert(Common_Bulk_KS_Stopping_L46x46x80::WP_CZ_BOHR == -32.0,
              "bulk KS stopping: launch z should be -32 (= -L_z/2 + 4 sigma)");
// L_z, not L_xy, must be the binding constraint — the whole point of L_xy = 46.
static_assert(Common_Bulk_KS_Stopping_L46x46x80::T_IFW_AU
                  < Common_Bulk_KS_Stopping_L46x46x80::T_TRANSVERSE_AU,
              "bulk KS stopping: transverse box too narrow — images bind before "
              "the +z face is reached; widen LX/LY or shorten the run");
// The fit window must be long enough for a meaningful slope.
static_assert(Common_Bulk_KS_Stopping_L46x46x80::FIT_T1_AU
                  - Common_Bulk_KS_Stopping_L46x46x80::FIT_T0_AU > 10.0,
              "bulk KS stopping: fit window shorter than 10 a.u.");
// Nyquist guard (.claude cutoff-aliasing rule): k_max = k0 + 3/(2 sigma).
static_assert(Common_Bulk_KS_Stopping_L46x46x80::WP_K0
                  + 3.0 / (2.0 * Common_Bulk_KS_Stopping_L46x46x80::WP_SIGMA_BOHR)
                  < 3.14159265358979 / Common_Bulk_KS_Stopping_L46x46x80::SPACING_BOHR,
              "bulk KS stopping: grid too coarse for the WP bandwidth");

}  // namespace jellium::config
