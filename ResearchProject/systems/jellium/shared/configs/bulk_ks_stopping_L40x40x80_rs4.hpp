// ============================================================================
// shared/configs/bulk_ks_stopping_L40x40x80_rs4.hpp
//
// DENSITY REPLICA of the bulk KS-stopping twin pair. Identical projectile,
// identical propagation, ~2.9x DENSER bath. Plan:
// docs/plans/bulk-jellium-ks-stopping.md (user decision, 2026-07-31).
//
//   Bulk_KS_Stopping_L40x40x80_rs4_WP         — Gaussian wave-packet projectile
//   Bulk_KS_Stopping_L40x40x80_rs4_Classical  — classical Gaussian-electron
//                                               projectile (UPF + m_e mass)
//
// ---------------------------------------------------------------------------
// WHY THIS RUN EXISTS
//
// The r_s = 5.702 pair (bulk_ks_stopping_L46x46x80) gave
//     S_classical = 0.377 eV/Bohr   vs   S_WP(drift) = 0.057 eV/Bohr
// a factor of 6.6. The standing hypothesis (user, 2026-07-31) is that the WP
// INJECTION evacuates bath density around the packet, so the wavepacket travels
// through a locally rarefied medium and under-reports the drag; the slab runs
// reportedly show a much smaller discrepancy. See the handover section
// "OPEN HYPOTHESIS: initialisation density-clearing".
//
// DENSITY IS THE LEVER. If clearing drives the gap, raising n shortens the
// screening length and the WP/classical ratio should shrink markedly. If the
// ratio is roughly density-independent, the cause lies in the KE definitions or
// in free-packet dispersion instead. THE PROJECTILE IS THEREFORE HELD FIXED —
// sigma, energy, launch z, dt, N_STEPS and L_z are all bit-identical to the
// r_s = 5.702 pair, so n is the only variable that moves.
//
// ---------------------------------------------------------------------------
// GEOMETRY. L_z = 80 is UNCHANGED (same flight, same launch, same step count).
// L_xy drops 46 -> 40 purely to keep the state count affordable at the higher
// density; see the cost note below.
//
// CONSEQUENCE — THE TRANSVERSE CONSTRAINT NOW BINDS. At L_xy = 46 the ordering
// was T_IFW (18.97) < T_TRANSVERSE (21.31), so L_z bound the clean window. At
// L_xy = 40 it INVERTS: T_TRANSVERSE = 18.43 < T_IFW = 18.97, so the periodic
// images in x/y are what end the trustworthy window. FIT_T1 takes the min of
// the two either way, so the analysis is correct; but the window is 2.8 %
// shorter than the r_s = 5.702 pair's (18.43 vs 18.97 a.u.) and any comparison
// of the two pairs must either accept that or re-fit the old run on [4, 18.43].
// The static_assert below pins the INVERSION deliberately, so that a future
// edit widening L_xy trips it and forces this comment to be re-read.
//
// ---------------------------------------------------------------------------
// BATH. N = 482 electrons in 40 x 40 x 80 = 128,000 Bohr^3 gives
// n = 3.7656e-3 e/Bohr^3, r_s = 3.9874 — a 2.92x density increase over the
// r_s = 5.702 pair. N = 482 is a genuine closed shell FOR THIS BOX with a
// 0.252 eV gap (enumerated numerically 2026-07-31), slightly BETTER than the
// 0.244 eV of the original run and ~29x the 100 K smearing.
//
// hbar omega_p = 5.919 eV (was 3.462), v_F = 0.4813 (was 0.3366),
// E_F = 0.11583 Ha. v/v_F = 5.63 — still comfortably above the Fermi velocity,
// so plasmon excitation remains kinematically allowed.
//
// IMPROVED, THOUGH STILL NOT STEADY STATE. 2 pi / omega_p = 28.88 a.u. against
// the 25.84 a.u. run: 0.89 plasma periods, up from 0.52 at r_s = 5.702. The
// bath now has time for nearly a full collective oscillation. S is STILL an
// initial-drag stopping power, not a converged S(v) — see
// .claude/rules/light-projectile-stopping.md — but materially less far off.
//
// ---------------------------------------------------------------------------
// GRID. dx = 0.50 (user decision 2026-07-31, "if you can get away with it")
// -> 80 x 80 x 160 = 1.024 M points, cutoff 19.74 Ha, k_Nyquist = 6.283.
//
// VERIFIED, not assumed (2026-07-31): the WP momentum moments were computed on
// this exact grid and reproduce the analytic values to MACHINE PRECISION —
// <p_z> = 2.711063 (0 error), T1 = 105.1021 eV (0 error), T1 - T2 = 5.1021 eV
// (0 error). The packet's spectral content reaches k0 + 3 sigma_k = 3.461,
// which is 1.82x inside Nyquist; the bath's k_F = 0.4813 is 13x inside; and the
// classical Gaussian UPF's form factor exp(-k^2 sigma_pot^2 / 2) is ~1e-17 at
// the Nyquist edge. Nothing in this run has spectral weight near the cutoff.
//
// DEVIATION TO DECLARE. Every other jellium run in this project used dx = 0.40.
// The WP/classical RATIO extracted here is internally consistent (both halves
// share this grid), but a cross-pair comparison of ratios carries a second-order
// grid difference. State this wherever the two pairs are compared.
//
// ---------------------------------------------------------------------------
// COST. 241 occupied + 20 extra = 261 states (+1 WP slot) on a 1.024 M grid is
// ~4.3 GB per orbital copy, ~13 GB for the three copies ETRS propagation needs
// — comfortable on an 80 GB A100. Projected ~1.6x the r_s = 5.702 pair:
// WP ~2.4 h, classical ~3.2 h, both well inside the 24 h wall (so the classical
// half stays non-resumable per the user decision of 2026-07-30).
//
// Rejected alternative for the record: r_s = 3 in the ORIGINAL 46 x 46 x 80 box
// needs N = 1514, i.e. 778 states and ~99 GB — it does not fit on an A100 at
// all. That is why the density target is r_s = 4, not r_s = 3.
//
// SIGMA CONVENTION (.claude/rules/sigma-wp-convention.md). sigma = 2.0 Bohr is
// the WAVEPACKET sigma; the classical UPF's CHARGE std is sigma/sqrt2 = 1.414,
// so both projectiles present the identical charge cloud. Both are labelled
// sigma = 2. Same UPF file as the r_s = 5.702 pair — the projectile is the
// control variable and must not be regenerated.
// ============================================================================
#pragma once

#include "base.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

struct Common_Bulk_KS_Stopping_L40x40x80_rs4 : Base {

    // ----- Cell (orthorhombic, fully periodic in x, y AND z) --------------
    static constexpr double LX_BOHR = 40.0;
    static constexpr double LY_BOHR = 40.0;
    static constexpr double LZ_BOHR = 80.0;   // UNCHANGED from the r_s=5.702 pair
    static constexpr double L_BOHR  = LZ_BOHR;   // legacy alias, log strings only

    // ----- Bath -----------------------------------------------------------
    static constexpr int    N_ELECTRONS  = 482;   // closed shell, gap 0.252 eV
    static constexpr int    EXTRA_STATES = 20;    // degenerate buffer + WP slot
    static constexpr double SPACING_BOHR = 0.50;
    static constexpr double CUTOFF_HA    = 19.74; // pi^2/(2 dx^2), for the log
    static constexpr double SCF_TOL_HA   = 1.0e-6;

    // ----- Projectile — IDENTICAL to the r_s=5.702 pair (control variable) -
    static constexpr double WP_EKIN_EV    = 100.0;
    static constexpr double WP_SIGMA_BOHR = 2.0;                    // WAVEPACKET sigma
    static constexpr double WP_K0         = k0_from_ev(WP_EKIN_EV); // 2.7111
    static constexpr double WP_KX = 0.0;
    static constexpr double WP_KY = 0.0;
    static constexpr double WP_KZ = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    static constexpr double WP_CX_BOHR = 0.0;
    static constexpr double WP_CY_BOHR = 0.0;
    // Standard boundary rule: -L_z/2 + 4 sigma = -32.0 (L_z unchanged)
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, LZ_BOHR);

    // ----- Real time — IDENTICAL to the r_s=5.702 pair ---------------------
    static constexpr double DT_AU = 0.040;
    static constexpr int N_STEPS =
        boundary::n_steps_for(WP_SIGMA_BOHR, LZ_BOHR, WP_K0, DT_AU);  // 646

    // Propagation-callback stride. The ENERGY series (observables.csv) is
    // written from that callback, so it follows THIS, not WRITE_EVERY.
    // Kept fine (2) because the CLASSICAL stopping power is fitted to
    // energy_total(x): thinning it to the VTI cadence would cut the fit
    // window from ~68 points to ~17. Storage is saved on the VTI writes,
    // which are gated separately on WRITE_EVERY inside the callback.
    static constexpr int OBS_EVERY = 2;

    static constexpr int WRITE_EVERY    = 8;    // density_total / _wp / _delta
    static constexpr int WF_WRITE_EVERY = 32;    // complex wavefunction
    static constexpr int STATS_EVERY    = 1;    // THE stopping data, every step

    static constexpr int CKPT_EVERY = 200;

    // ----- Analysis window -------------------------------------------------
    static constexpr double T_IFW_AU = boundary::ifw_end_t_dispersive(
        WP_SIGMA_BOHR, LZ_BOHR, WP_K0, WP_CZ_BOHR);            // 18.97
    static constexpr double T_TRANSVERSE_AU =
        boundary::transverse_clean_t(WP_SIGMA_BOHR, LX_BOHR);  // 18.43  <-- BINDS
    static constexpr double FIT_T0_AU = 4.0;
    static constexpr double FIT_T1_AU = T_IFW_AU < T_TRANSVERSE_AU
                                          ? T_IFW_AU : T_TRANSVERSE_AU;

    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;

    // ----- Classical projectile — same UPF, same launch, same velocity -----
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

struct Bulk_KS_Stopping_L40x40x80_rs4_WP
    : Common_Bulk_KS_Stopping_L40x40x80_rs4 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Bulk_KS_Stopping_L40x40x80_rs4_Classical
    : Common_Bulk_KS_Stopping_L40x40x80_rs4 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

// ---- Compile-time guards on the derived geometry ---------------------------
static_assert(Common_Bulk_KS_Stopping_L40x40x80_rs4::N_STEPS == 646,
              "rs4 pair: N_STEPS must stay 646 — parity with the r_s=5.702 pair "
              "is what makes the density comparison controlled");
static_assert(Common_Bulk_KS_Stopping_L40x40x80_rs4::WP_CZ_BOHR == -32.0,
              "rs4 pair: launch z should be -32 (= -L_z/2 + 4 sigma)");
// DELIBERATE INVERSION vs the L=46 pair: here the TRANSVERSE images bind first.
// If this assert ever fires, L_xy has been widened past ~42 and the header's
// window discussion is stale — re-read it before proceeding.
static_assert(Common_Bulk_KS_Stopping_L40x40x80_rs4::T_TRANSVERSE_AU
                  < Common_Bulk_KS_Stopping_L40x40x80_rs4::T_IFW_AU,
              "rs4 pair: transverse limit no longer binds — L_xy was widened, so "
              "the header's 'transverse constraint now binds' note is stale");
// The fit window must still be long enough for a meaningful slope.
static_assert(Common_Bulk_KS_Stopping_L40x40x80_rs4::FIT_T1_AU
                  - Common_Bulk_KS_Stopping_L40x40x80_rs4::FIT_T0_AU > 10.0,
              "rs4 pair: fit window shorter than 10 a.u.");
// Nyquist guard at the COARSER dx = 0.50: k_max = k0 + 3/(2 sigma) = 3.461
// against pi/0.5 = 6.283. Verified numerically to machine precision, see header.
static_assert(Common_Bulk_KS_Stopping_L40x40x80_rs4::WP_K0
                  + 3.0 / (2.0 * Common_Bulk_KS_Stopping_L40x40x80_rs4::WP_SIGMA_BOHR)
                  < 3.14159265358979
                        / Common_Bulk_KS_Stopping_L40x40x80_rs4::SPACING_BOHR,
              "rs4 pair: grid too coarse for the WP bandwidth");
// The whole point of the run: the bath must actually be denser.
static_assert(Common_Bulk_KS_Stopping_L40x40x80_rs4::N_ELECTRONS
                  / (Common_Bulk_KS_Stopping_L40x40x80_rs4::LX_BOHR
                     * Common_Bulk_KS_Stopping_L40x40x80_rs4::LY_BOHR
                     * Common_Bulk_KS_Stopping_L40x40x80_rs4::LZ_BOHR)
                  > 2.0 * 1.287807e-3,
              "rs4 pair: bath is not at least 2x denser than the r_s=5.702 pair — "
              "the density lever is the entire purpose of this run");

}  // namespace jellium::config
