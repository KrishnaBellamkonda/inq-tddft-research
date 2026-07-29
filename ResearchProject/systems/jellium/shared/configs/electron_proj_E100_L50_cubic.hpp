// ============================================================================
// shared/configs/electron_proj_E100_L50_cubic.hpp  — low-energy companion to
//   shared/configs/electron_proj_E1500_L50_cubic.hpp, intended for the
//   meeting-prep classical-vs-WP case study (`docs/plans/over-the-last-two-
//   sharded-muffin.md`).
//
// Two Cfgs sharing one cubic 50 x 50 x 50 Bohr periodic jellium bath
// (N=162, dx=0.40 Bohr, dt=0.020 a.u., N_STEPS=200):
//
//   Electron_Proj_E100_L50_cubic_WP_dx0p40         — Gaussian wave packet
//   Electron_Proj_E100_L50_cubic_Classical_dx0p40  — classical electron
//                                                    (custom UPF + mass override)
//
// Both projectiles enter at (0, 0, -10) Bohr (INQ centred Cartesian) —
// IDENTICAL launch position to the E=1500 pair, satisfying the user's
// "same initial position" requirement. They share kinetic energy
// 100 eV (= 3.6749 Ha):
//   - WP:        k_0 = sqrt(2 * 3.6749) = 2.7111 bohr^-1, sigma = 5.0 Bohr.
//   - Classical: v_z = sqrt(2 * 3.6749) = 2.7111 bohr/atu, m = m_e (1.0 a.u.).
// Velocity match is automatic: same KE, same mass ⇒ same v_group = 2.7111.
//
// PHYSICS RATIONALE (for the meeting):
// At v=2.71 a.u., the projectile sits **between** the v_F Lindhard peak
// (v_F ≈ 0.55 at r_s=5.69) and the high-v Bethe asymptote (v=10.5 in
// the E=1500 sibling pair). It is the regime where collective (plasmon)
// and single-particle (e-h) channels can both contribute, and where the
// classical electron's lack of momentum-space spread (a delta function
// in k) departs most from the WP's finite-width (sigma_k=0.20 Bohr^-1)
// distribution. Matching the bath, geometry, dt window, and launch
// position of the E=1500 reference isolates the projectile-kind effect
// at fixed bath.
//
// PRIOR DATA POINT:
// run_base_n138 (E=100 eV, but L=60 cubic, N=138, dx=0.55) ran successfully
// at dt=0.020, N_STEPS=320. The new pair places E=100 eV in the canonical
// L=50 N=162 family (same bath as the plasmon runs and the E=1500
// classical), so it is *not* a re-run — it is the missing matched-pair
// data point in the L=50 family.
//
// GRID / NYQUIST CHECK:
// dx = 0.40 ⇒ k_Nyquist = π/0.40 = 7.854 Bohr^-1.
// WP k_max = k_0 + 3 sigma_k = 2.711 + 0.600 = 3.311 Bohr^-1.
// 3.311 / 7.854 = 0.42 ⇒ very comfortable Nyquist headroom, NO aliasing
// (contrast E=1500 WP where k_0+3sigma_k=11.10 was 6% over Nyquist).
//
// TIME-STEP CHOICE:
// dt = 0.020 a.u. (vs E=1500's 0.005). The ETRS stability constraint
// dt * v < dx gives 0.020 * 2.711 = 0.0542 ≪ 0.40 — comfortable. Matches
// the dt used by run_base_n138 (E=100 eV) and the L=50 base/plasmon
// runs, so we keep the established stability margin.
//
// PROPAGATION TIME:
// N_STEPS = 646 ⇒ t_total = 12.92 a.u. ≈ 0.313 fs. Travel at initial v:
// 2.711 * 12.92 = 35.03 Bohr ⇒ projectile centroid reaches z = -10 +
// 35.03 = +25 Bohr (the +z box face) exactly at t_end. This is the
// "maximum-traversal" propagation length, set by user instruction
// (2026-05-13 dapper-moon plan addendum): t_total = (L_box - z_launch) /
// v_initial, taking the projectile from the launch position to the
// opposite box face at its initial speed. The WP centroid (no
// deceleration) reaches z = +25; the classical projectile under
// .ehrenfest() decelerates ~1% so reaches z ≈ +24.65.
//
// After t_end the WP envelope leading edge has crossed into the periodic
// image, but at that point the propagation is over — no wrap
// contamination during the trajectory itself.
// ============================================================================
#pragma once

#include "base_n162_L50_E1p5.hpp"
#include "boundary_rule.hpp"

namespace jellium::config {

// Common bath/grid; both WP and classical Cfgs share these via inheritance.
// Reuses Base_N162_L50_E1p5 for L=50 cubic, N=162, EXTRA_STATES=20, T=100K,
// SCF_TOL=1e-6, WP_SIGMA=5.0 — overrides only what differs at E=100 eV.
struct Common_E100_L50_cubic : Base_N162_L50_E1p5 {

    // ----- Grid: dx=0.40 to reuse the existing classical GS at
    //   checkpoints/gs_L50_cubic_N162_dx0p40/ (no new GS required for
    //   either run in this pair). Nyquist headroom is very comfortable
    //   at E=100 eV.
    static constexpr double SPACING_BOHR = 0.40;
    static constexpr double CUTOFF_HA    = 30.84;       // pi^2/(2 dx^2);
                                                        // log only.

    // ----- Real-time -----------------------------------------------------
    static constexpr double DT_AU       = 0.020;
    // N_STEPS = ceil((L_box - z_launch) / v_initial / dt) = 35 / 2.7111 /
    // 0.02 ≈ 646. Sets t_total so the projectile centroid reaches the
    // +z box face (z = +25 Bohr) at its initial velocity.
    static constexpr int    N_STEPS     = 646;          // total t = 12.92 a.u.
    static constexpr int    WRITE_EVERY = 2;            // 323 density frames

    // ----- Wave-packet projectile (E_kin = 100 eV) ----------------------
    static constexpr double WP_EKIN_EV      = 100.0;
    static constexpr double WP_K0           = k0_from_ev(WP_EKIN_EV);  // 2.7111
    static constexpr double WP_KX           = 0.0;
    static constexpr double WP_KY           = 0.0;
    static constexpr double WP_KZ           = +WP_K0;
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // WP_SIGMA_BOHR inherited (5.0 Bohr) from Base_N138_L50_E1p5.

    // Launch IDENTICAL to E=1500 pair: z=-10 Bohr in INQ centred
    // Cartesian (= 15 Bohr from the -z face at z=-25, 3 sigma_r). This
    // satisfies the user's "same initial position" requirement.
    static constexpr double WP_CX_BOHR =   0.0;
    static constexpr double WP_CY_BOHR =   0.0;
    static constexpr double WP_CZ_BOHR = -10.0;

    // ----- Classical-electron projectile (same KE, m = m_e) -------------
    // Pseudopotential and mass override IDENTICAL to E=1500 pair.
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
    // v_z = sqrt(2 * KE_Ha / m) with m=1, KE=3.6749 Ha ⇒ 2.7111 bohr/atu.
    // Same numerical value as WP_K0 (KE+mass identity).
    static constexpr double PROJ_VEL_Z = const_sqrt(2.0 * WP_EKIN_EV / HA_TO_EV);

    // ----- Time-window placeholders (used by some screen tooling) -------
    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;

    // SCREEN_SNAP_EVERY inherited (6); irrelevant for this comparison.
};

// WP run: no ion, only WP injection.
struct Electron_Proj_E100_L50_cubic_WP_dx0p40 : Common_E100_L50_cubic {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

// Classical run: no WP, one ion with custom UPF + mass override.
struct Electron_Proj_E100_L50_cubic_Classical_dx0p40 : Common_E100_L50_cubic {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

// ============================================================================
// v2 Cfgs (added 2026-05-17 for the meeting campaign Run-3)
//
// Same bath, same energy, same projectile species — but the launch + stop
// + N_STEPS are pinned to the universal boundary_rule (launch=-5, stop=+20,
// traversal=25 Bohr, N_STEPS=462). The v2 runs also enable the new
// WPMomentumStats / WPRealSpaceStats observables (handled in run.cpp).
//
// The v1 pair (Common_E100_L50_cubic above) used launch=-10 and a 35-Bohr
// traversal that violated the 1σ stop rule. Kept for legacy comparability.
// ============================================================================

struct Common_E100_L50_cubic_v2 : Common_E100_L50_cubic {
    // Override launch + propagation length to boundary_rule.
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                     // -5
    static constexpr double PROJ_LAUNCH_Z = WP_CZ_BOHR;                // -5

    static constexpr int N_STEPS =
        boundary::n_steps_for(WP_SIGMA_BOHR, L_BOHR, WP_K0, DT_AU);    // 462
    static constexpr int WRITE_EVERY = boundary::write_every_for(N_STEPS);

    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E100_L50_cubic_WP_dx0p40_v2 : Common_E100_L50_cubic_v2 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Electron_Proj_E100_L50_cubic_Classical_dx0p40_v2 : Common_E100_L50_cubic_v2 {
    static constexpr bool WP_ENABLED = false;
    static constexpr int  N_IONS     = 1;
};

// ─── Run-4: extra-states basis-completeness test (2026-05-18) ────────────────
// Tests whether the unaccounted-overlap / "missing electrons" puzzle in the
// WP runs is a basis-completeness artefact. Same Cfg as v2 but with
// extra_states = 40 and 80 (vs the default 20). New GSes required.
// Three data points total: extra ∈ {20 (= v2), 40, 80}.

struct Common_E100_L50_cubic_x40 : Common_E100_L50_cubic_v2 {
    static constexpr int EXTRA_STATES = 40;
};

struct Electron_Proj_E100_L50_cubic_WP_dx0p40_x40 : Common_E100_L50_cubic_x40 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

struct Common_E100_L50_cubic_x80 : Common_E100_L50_cubic_v2 {
    static constexpr int EXTRA_STATES = 80;
};

struct Electron_Proj_E100_L50_cubic_WP_dx0p40_x80 : Common_E100_L50_cubic_x80 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

// ─── Run-7: σ-sweep at E=100 eV (added 2026-05-18) ─────────────────────────
// Same bath as v2, vary the WP wavefunction σ. We already have:
//   σ=5  via Common_E100_L50_cubic_v2 (Run-3 WP_v2)
//   σ=1  via electron_proj_E100_L50_cubic_sigma1.hpp (σ=1 task)
//   σ=0.5 at L=30 (Run-6 high-density — different bath)
// Adding: σ=0.5 at L=50 (same bath as v2 — pure σ effect, separates from
// Run-6's density effect), σ=3, σ=8. σ=0.25 deferred (would need dx=0.30
// which has the known WP-init CUDA memory issue).
//
// SELF-SPREAD CAPPING is required at small σ: σ_density(t) = (σ/√2) ·
// √(1 + (t/σ²)²). For σ=0.5, at t=4.8 a.u. the 3σ density tail hits the
// far face (L/2 - launch_z=20+25=45, σ_density≈4.5 at t=4.8 → 13.5 →
// the launch wall, with launch_z=-L/2+4σ=-23). N_STEPS capped at 240
// (t_total=4.8 a.u. — still ~10 Bohr of centroid travel at v=2.71).
// For σ=3 and σ=8 the standard boundary_rule is fine.

struct Common_E100_L50_sigma3 : Common_E100_L50_cubic_v2 {
    static constexpr double WP_SIGMA_BOHR = 3.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -13
    static constexpr int N_STEPS =
        boundary::n_steps_for(WP_SIGMA_BOHR, L_BOHR, WP_K0, DT_AU);   // 461 ish
    static constexpr int WRITE_EVERY = boundary::write_every_for(N_STEPS);
    // wf-saving cadence for the 2026-05-31 σ-sweep rerun (momentum + loss fn +
    // bath-only density via total_excluding_orbital). = WRITE_EVERY so the WP
    // density is available every density frame for free bath subtraction.
    static constexpr int WF_WRITE_EVERY = WRITE_EVERY;
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E100_L50_sigma3_WP_dx0p40 : Common_E100_L50_sigma3 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

// σ=0.5 at L=50 (same bath as Run-3 — pure σ effect at standard density).
// Self-spread-capped at 240 steps.
struct Common_E100_L50_sigma0p5 : Common_E100_L50_cubic_v2 {
    static constexpr double WP_SIGMA_BOHR = 0.5;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);                    // -23
    static constexpr int N_STEPS     = 240;        // self-spread cap
    static constexpr int WRITE_EVERY = boundary::write_every_for(N_STEPS);
    static constexpr int WF_WRITE_EVERY = WRITE_EVERY;   // σ-sweep rerun (see σ3)
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E100_L50_sigma0p5_WP_dx0p40 : Common_E100_L50_sigma0p5 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

// σ=8 with the relaxed 3σ/1σ boundary rule (σ=8 in L=50 leaves no room
// for the standard 4σ/1σ rule — see boundary_rule.hpp).
struct Common_E100_L50_sigma8 : Common_E100_L50_cubic_v2 {
    static constexpr double WP_SIGMA_BOHR = 8.0;
    static constexpr double WP_CZ_BOHR =
        boundary::launch_z_relaxed(WP_SIGMA_BOHR, L_BOHR);            // -1
    static constexpr int N_STEPS =
        boundary::n_steps_for_relaxed(WP_SIGMA_BOHR, L_BOHR, WP_K0, DT_AU);
    static constexpr int WRITE_EVERY = boundary::write_every_for(N_STEPS);
    static constexpr int WF_WRITE_EVERY = WRITE_EVERY;   // σ-sweep rerun (see σ3)
    static constexpr double T2_AU = DT_AU * N_STEPS;
    static constexpr double T2_FS = T2_AU / FS_TO_AU;
};

struct Electron_Proj_E100_L50_sigma8_WP_dx0p40 : Common_E100_L50_sigma8 {
    static constexpr bool WP_ENABLED = true;
    static constexpr int  N_IONS     = 0;
};

}  // namespace jellium::config
