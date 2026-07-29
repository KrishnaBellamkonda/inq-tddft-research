// ============================================================================
// shared/configs/boundary_rule.hpp  (jellium)
//
// Universal launch/stop/traversal/cadence helpers for every new jellium WP
// run started 2026-05-17 onwards (campaign for the 2026-05-21 Emilio
// meeting and beyond).
//
// Rule sources:
//   - Plan: docs/plans/jellium-meeting-2026-05-21.md §"Universal rules"
//   - Design journal: docs/journals/researchproject/
//                     2026-05-17_jellium_meeting_design.md §"Key resolved
//                     ambiguities" (#1, #3, #10)
//   - Memory: feedback_jellium_boundary_rule.md
//
// Two geometry rules:
//
//   1. STANDARD (used by 99 % of runs):
//        launch_z(σ, L)   = -L/2 + 4σ
//        stop_z (σ, L)    = +L/2 -  σ
//        traversal        = L - 5σ
//        t_IFW end is when centroid reaches +L/2 - 4σ
//
//   2. RELAXED (only when σ is so large that the standard rule is
//      geometrically infeasible in the box — e.g. σ=8 in L=50):
//        launch_z_relaxed(σ, L) = -L/2 + 3σ
//        stop_z_relaxed  (σ, L) = +L/2 -  σ
//        traversal_relaxed      = L - 4σ
//        t_IFW end is when centroid reaches +L/2 - 3σ (same formula
//        as standard — IFW is a Gaussian-tail-at-far-face criterion,
//        independent of launch convention).
//
// IFW physics note: the IFW end is when the WP centroid is 3σ from the
// far box face — at that point the WP density at the face is
// exp(-9/2) ≈ 0.011 of peak, the leading-edge tail starts overlapping
// with its +z periodic image, and the run can no longer be considered
// interference-free. This is independent of where the WP was launched
// from, so both rules share +L/2 - 3σ as the IFW-end criterion.
//
// The relaxed rule must be opted-in deliberately per Cfg, with a comment
// explaining why the standard rule is infeasible. See task #31 in the
// 2026-05-17 task list for the σ→∞ open question.
//
// Cadence rule:
//   target_frames = 300 (default, comparable gif lengths + disk usage
//                        across the whole campaign).
//   WRITE_EVERY = max(1, round(N_STEPS / target_frames))
//
// All quantities atomic units / Bohr / eV. No external citation.
// ============================================================================
#pragma once

#include <algorithm>
#include <cmath>

namespace jellium::config::boundary {

// ----- Geometry: standard 4σ / 1σ rule --------------------------------------

// Launch z-coordinate (WP centroid at t=0) under the standard rule.
inline constexpr double launch_z(double sigma_r_bohr, double L_bohr) {
    return -0.5 * L_bohr + 4.0 * sigma_r_bohr;
}

// Stop z-coordinate (WP centroid at t=t_end) under the standard rule.
// Both rules share the same stop: +L/2 - σ.
inline constexpr double stop_z(double sigma_r_bohr, double L_bohr) {
    return 0.5 * L_bohr - sigma_r_bohr;
}

// Traversal length L - 5σ under the standard rule.
inline constexpr double traversal(double sigma_r_bohr, double L_bohr) {
    return L_bohr - 5.0 * sigma_r_bohr;
}

// Centroid position at which the interference-free window ends.
// The IFW criterion is "Gaussian 3σ tail reaches the far box face" —
// a physics criterion, independent of launch convention. So this is the
// same formula for both standard and relaxed rules: +L/2 - 3σ.
inline constexpr double ifw_end_z(double sigma_r_bohr, double L_bohr) {
    return 0.5 * L_bohr - 3.0 * sigma_r_bohr;
}

// ----- Geometry: relaxed 3σ / 1σ rule (large σ only) ------------------------

inline constexpr double launch_z_relaxed(double sigma_r_bohr, double L_bohr) {
    return -0.5 * L_bohr + 3.0 * sigma_r_bohr;
}

inline constexpr double stop_z_relaxed(double sigma_r_bohr, double L_bohr) {
    return 0.5 * L_bohr - sigma_r_bohr;
}

inline constexpr double traversal_relaxed(double sigma_r_bohr,
                                          double L_bohr) {
    return L_bohr - 4.0 * sigma_r_bohr;
}

// Relaxed-rule IFW end is the same formula (Gaussian-tail physics):
// +L/2 - 3σ. Provided as an alias for symmetry; identical to ifw_end_z.
inline constexpr double ifw_end_z_relaxed(double sigma_r_bohr,
                                          double L_bohr) {
    return ifw_end_z(sigma_r_bohr, L_bohr);
}

// ----- Time domain: standard rule -------------------------------------------

// Total propagation time at constant initial velocity.
inline constexpr double t_total_au(double sigma_r_bohr, double L_bohr,
                                   double v_bohr_per_au) {
    return traversal(sigma_r_bohr, L_bohr) / v_bohr_per_au;
}

// t_IFW under standard rule: travel from launch_z to ifw_end_z at v.
inline constexpr double t_ifw_au(double sigma_r_bohr, double L_bohr,
                                 double v_bohr_per_au) {
    return (ifw_end_z(sigma_r_bohr, L_bohr)
            - launch_z(sigma_r_bohr, L_bohr)) / v_bohr_per_au;
}

// N_STEPS = ceil((L - 5σ) / (v · dt)).
inline constexpr int n_steps_for(double sigma_r_bohr, double L_bohr,
                                 double v_bohr_per_au, double dt_au) {
    // No std::ceil in constexpr pre-C++23; inline an integer-friendly
    // version: ceil(a/b) = (a + b - 1) / b when a, b > 0. Here we work
    // with doubles, so use the floor + 1 trick after casting.
    const double t_total = t_total_au(sigma_r_bohr, L_bohr, v_bohr_per_au);
    const double n_double = t_total / dt_au;
    const int n_floor = static_cast<int>(n_double);
    return n_double > static_cast<double>(n_floor) ? n_floor + 1 : n_floor;
}

// ----- Time domain: relaxed rule --------------------------------------------

inline constexpr double t_total_relaxed_au(double sigma_r_bohr,
                                           double L_bohr,
                                           double v_bohr_per_au) {
    return traversal_relaxed(sigma_r_bohr, L_bohr) / v_bohr_per_au;
}

inline constexpr double t_ifw_relaxed_au(double sigma_r_bohr,
                                         double L_bohr,
                                         double v_bohr_per_au) {
    return (ifw_end_z_relaxed(sigma_r_bohr, L_bohr)
            - launch_z_relaxed(sigma_r_bohr, L_bohr)) / v_bohr_per_au;
}

inline constexpr int n_steps_for_relaxed(double sigma_r_bohr,
                                         double L_bohr,
                                         double v_bohr_per_au,
                                         double dt_au) {
    const double t_total =
        t_total_relaxed_au(sigma_r_bohr, L_bohr, v_bohr_per_au);
    const double n_double = t_total / dt_au;
    const int n_floor = static_cast<int>(n_double);
    return n_double > static_cast<double>(n_floor) ? n_floor + 1 : n_floor;
}

// ----- Cadence: write-every for ~300-frame target ---------------------------

// Target N_FRAMES per run across the campaign.
inline constexpr int DEFAULT_TARGET_FRAMES = 300;

// Returns max(1, round(n_steps / target_frames)).
inline constexpr int write_every_for(int n_steps,
                                     int target_frames
                                       = DEFAULT_TARGET_FRAMES) {
    if (n_steps <= 0 || target_frames <= 0) return 1;
    // Integer round-half-up: (n_steps + target_frames/2) / target_frames.
    const int rounded = (n_steps + target_frames / 2) / target_frames;
    return rounded < 1 ? 1 : rounded;
}

// ----- Static smoke tests: verify the rule constants at compile time -------
// The standard rule at σ=5 Bohr, L=50 Bohr, v=2.711 Bohr/atu, dt=0.02:
//   launch_z = -5  ; stop_z = +20 ; traversal = 25
//   t_total  = 9.222 a.u.  ; t_ifw = 3.689 a.u.  ; N_STEPS = 462
// (User design lock 2026-05-17 quoted N_STEPS=461 — this is a single-step
// rounding difference from the integer-ceil convention; the helper
// answer rounds up to 462 which is the safer floor for boundary reach.)
//
// Relaxed at σ=8, L=50, v=2.711:
//   launch_z_relaxed = -1 ; stop_z = +17 ; traversal_relaxed = 18
//   t_total  = 6.640 a.u.  ; t_ifw = 5.532 a.u.  ; N_STEPS_relaxed = 332
//
// Cadence (defaults 300): write_every_for(461)=2 ; write_every_for(922)=3 ;
//   write_every_for(348)=1 ; write_every_for(1290)=4 (was 2 with WE=2
//   convention from the legacy 25-eV Cfg).
//
// These are documented expected values from the plan run-inventory table;
// not run-as-tests here (Cfg headers are header-only). The Tutorial smoke
// test that exercises these via static_assert lives in
// Tutorial/wp-boundary-rule-test/ (Infra-10 follow-up, not yet created).
// ----------------------------------------------------------------------------

static_assert(launch_z(5.0, 50.0) == -5.0,
              "launch_z(5, 50) should be -5");
static_assert(stop_z(5.0, 50.0) == 20.0,
              "stop_z(5, 50) should be +20");
static_assert(traversal(5.0, 50.0) == 25.0,
              "traversal(5, 50) should be 25");
static_assert(ifw_end_z(5.0, 50.0) == 10.0,
              "ifw_end_z(5, 50) should be +10");

static_assert(launch_z_relaxed(8.0, 50.0) == -1.0,
              "launch_z_relaxed(8, 50) should be -1");
static_assert(stop_z_relaxed(8.0, 50.0) == 17.0,
              "stop_z_relaxed(8, 50) should be +17");
static_assert(traversal_relaxed(8.0, 50.0) == 18.0,
              "traversal_relaxed(8, 50) should be 18");
// IFW physics is rule-independent: σ=8, L=50 ⇒ +L/2 - 3σ = +25 - 24 = +1.
// At σ=8 launched from -1 (relaxed) the IFW window is only 2 Bohr / ~10 %
// of the 18-Bohr traversal — most of the σ=8 run is post-IFW. This is a
// known consequence of using σ=8 in L=50; the run is still informative
// for the σ→∞ professor's-claim question (TODO #31).
static_assert(ifw_end_z_relaxed(8.0, 50.0) == 1.0,
              "ifw_end_z_relaxed(8, 50) should be +1");

static_assert(write_every_for(461) == 2,
              "write_every_for(461, 300) should be 2");
static_assert(write_every_for(922) == 3,
              "write_every_for(922, 300) should be 3");
static_assert(write_every_for(348) == 1,
              "write_every_for(348, 300) should be 1");
static_assert(write_every_for(1) == 1,
              "write_every_for(1, 300) should clamp to 1");

}  // namespace jellium::config::boundary
