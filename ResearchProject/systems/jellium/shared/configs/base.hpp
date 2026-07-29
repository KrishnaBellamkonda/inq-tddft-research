// ============================================================================
// shared/configs/base.hpp  (jellium)
//
// Source-of-truth Cfg for every jellium run. Updated 2026-05-01 for the
// L=60 / centred-launch redesign (plan: .claude/plans/
// in-this-task-we-lively-meerkat.md):
//
//   * Cubic-periodic cell, L = 60 bohr (≈ 31.75 Å).
//   * N = 128 electrons (r_s = (3/(4π n))^{1/3} ≈ 7.38 a₀ preserved by
//     keeping n = N/V = 128/216000 = 5.926e-4 e/bohr³, matching the
//     legacy L=40, N=38 density to 0.3 %). Closed-shell behaviour is
//     established by the SCF; EXTRA_STATES provides degenerate buffer.
//   * Wave packet: σ = 0.53 Å (1.0015 bohr), E = 200 eV,
//     k₀ = 3.834 bohr⁻¹, launched at +z from the box centre
//     (30, 30, 30) bohr. Single-pass to z ≈ 55 bohr leaves a 5-bohr
//     margin from the periodic boundary.
//   * Real-time: dt = 0.020 a.u., N_STEPS = 320 (single-pass:
//     k₀·dt·N_STEPS ≈ 24.5 bohr travel from z=30).
//
// Cfg pattern: each variant is a `struct : Base { … overrides … };` so each
// override is statically known. The propagation driver
// `jellium::run_template::run_propagation<Cfg>(...)` instantiates one
// template per variant.
//
// All quantities in atomic units / Bohr / eV. No published reference paper
// is on file for this exact configuration; treat as internal-reference.
// ============================================================================
#pragma once

#include <cmath>

namespace jellium::config {

// ---- Unit conversions (CODATA 2018) ---------------------------------------
inline constexpr double ANG_TO_BOHR = 1.8897259886;
inline constexpr double HA_TO_EV    = 27.21138625;
inline constexpr double FS_TO_AU    = 41.341374575751;

// Compile-time sqrt (Newton iteration). std::sqrt is not constexpr under
// CUDA, so we ship our own. 30 iterations is past double-precision
// convergence for any plausible scientific range.
inline constexpr double const_sqrt(double x) {
    if (x <= 0.0) return 0.0;
    double g = x > 1.0 ? x : 1.0;
    for (int i = 0; i < 30; ++i) {
        g = 0.5 * (g + x / g);
    }
    return g;
}

// |k| from kinetic energy in eV. |k| = sqrt(2 * E_Ha) = sqrt(2 * E_eV / HA_TO_EV).
inline constexpr double k0_from_ev(double ekin_ev) {
    return const_sqrt(2.0 * ekin_ev / HA_TO_EV);
}

// ---- Base config: closed-shell N=38 jellium, 200 eV, σ=0.53 Å, +z --------
struct Base {
    // Cell (cubic, periodic). LX_BOHR/LY_BOHR/LZ_BOHR aliases preserved so
    // the run_template body matches coronene's.
    static constexpr double L_BOHR  = 60.0;
    static constexpr double LX_BOHR = L_BOHR;
    static constexpr double LY_BOHR = L_BOHR;
    static constexpr double LZ_BOHR = L_BOHR;

    // Electronic structure. N_ELECTRONS scales with V to preserve r_s
    // (38 e at L=40 → 128 e at L=60). EXTRA_STATES bumped to absorb any
    // near-degeneracy at the Fermi level under the larger box.
    static constexpr int    N_ELECTRONS    = 128;
    static constexpr int    EXTRA_STATES   = 4;
    static constexpr double SPACING_BOHR   = 0.50;
    static constexpr double TEMPERATURE_EV = 0.00862;     // ≈ 100 K
    static constexpr double SCF_TOL_HA     = 1.0e-4;
    static constexpr int    SCF_MAX_STEPS  = 300;
    static constexpr int    SCF_MIX_NDIM   = 8;
    static constexpr double SCF_MIX_ALPHA  = 0.1;

    // Stub for log compatibility with coronene template (jellium uses
    // spacing, not cutoff).
    static constexpr double CUTOFF_HA = 0.0;

    // Wave packet (Eq.-1 Gaussian envelope, paper-style normalisation)
    static constexpr double WP_SIGMA_ANG  = 0.53;
    static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;   // ≈ 1.0015
    static constexpr double WP_EKIN_EV    = 200.0;
    static constexpr double WP_K0         = k0_from_ev(WP_EKIN_EV);

    // WP centre and direction. INQ uses Cartesian coordinates
    // r ∈ [-L/2, +L/2] (cell centred at origin), so the box centre
    // is (0, 0, 0). The WP has L/2 ≈ 30 bohr of free space ahead of
    // it in +z before periodic loop-back (single-pass; see N_STEPS).
    static constexpr double WP_CX_BOHR = 0.0;
    static constexpr double WP_CY_BOHR = 0.0;
    static constexpr double WP_CZ_BOHR = 0.0;

    static constexpr double WP_KX = 0.0;
    static constexpr double WP_KY = 0.0;
    static constexpr double WP_KZ = +WP_K0;            // +z launch (jellium convention)
    static constexpr double WP_KZ_MAGNITUDE = WP_K0;

    // Real-time. N_STEPS sized for single-pass propagation:
    // travel = k0 * dt * N_STEPS ≈ 3.83 * 0.020 * 320 = 24.5 bohr,
    // so the WP front reaches z ≈ 30 + 24.5 = 54.5 bohr (leaves ≈5.5 bohr
    // from the periodic boundary at z=L=60).
    static constexpr double DT_AU             = 0.020;
    static constexpr int    N_STEPS           = 320;
    static constexpr int    WRITE_EVERY       = 2;
    static constexpr int    SCREEN_SNAP_EVERY = 3;

    // Paper window — placeholder for jellium (no reference paper). Set to
    // [0, total_time_au] so the paper accumulator equals the full-time
    // accumulator.
    static constexpr int    N_SCREENS = 20;
    static constexpr double T1_FS = 0.0;
    static constexpr double T2_FS = DT_AU * N_STEPS / FS_TO_AU;
    static constexpr double T1_AU = 0.0;
    static constexpr double T2_AU = DT_AU * N_STEPS;

    // Per-screen physics-window envelope (n_sigmas of WP). Same default as
    // coronene; jellium screen windows are placeholders pending the
    // physics re-derivation tracked in docs/plans/jellium_reorg.md §14.
    static constexpr double WP_ENVELOPE_SIGMAS = 2.0;
};

}  // namespace jellium::config
