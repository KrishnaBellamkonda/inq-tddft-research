#pragma once
// ============================================================================
// run_004 — Coronene TDDFT LEED simulation
//
// Changes from run_003:
//   - d = 0.53 Å (paper value; focused beam for LEED contrast)
//   - Lz = 47.55 Å (1.5 × 31.7 Å; longer cell to separate WP from molecule)
//   - D = Lz/2 − 5d = 21.125 Å (maximised impact distance)
//   - z_obs = Lz/2 + D = 44.900 Å (WP start = LEED screen; = Lz − 5d)
//   - T1 = D/k₀ = 10.41 a.u. = 0.252 fs (WP arrival at coronene)
//   - T2 = 3×T1 = 31.23 a.u. = 0.756 fs (matches paper T2/T1 ≈ 3.3 ratio)
//   - N_steps ≈ 1561 (~3× run_003)
//   - LEED accumulator uses background-subtracted density (see run.cpp)
//
// Rationale:
//   run_001 produced a structured LEED pattern because it used d=0.53 Å (focused
//   beam) and long propagation, but had a broken 4-corner molecule geometry.
//   run_004 restores those parameters with the correct centred geometry, so that
//   the 6-fold LEED pattern should reflect single-molecule scattering.
//
// Unchanged from run_003:
//   - Lx = Ly = 18.4 Å
//   - E_cut = 40 Ha
//   - E_kin = 200 eV
//   - Dt = 4.84e-4 fs (paper value)
//   - SCF: Broyden ndim=8, α=0.1, tol=1e-4
//
// Source: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)
// All values in atomic units unless annotated.
// ============================================================================

#include <inq/inq.hpp>
#include <cmath>

namespace cfg {

using namespace inq::magnitude;

// ── Unit conversions ─────────────────────────────────────────────────────────
inline constexpr double ANG_TO_BOHR  = 1.8897259886;
inline constexpr double BOHR_TO_ANG  = 0.529177210903;
inline constexpr double HA_TO_EV     = 27.21138625;
inline constexpr double FS_TO_AU     = 1.0 / 0.024188843;

// ── Cell dimensions ──────────────────────────────────────────────────────────
// Finite (isolated molecule).  Lz = 1.5 × 31.7 Å to allow WP long propagation.
// INQ cell origin: (0,0,0) = cell CORNER; grid runs 0 → L.
// Coronene centred at (Lx/2, Ly/2, Lz/2) — see coronene_centered.xyz.
inline constexpr double LX_ANG  = 18.4;
inline constexpr double LY_ANG  = 18.4;
inline constexpr double LZ_ANG  = 47.55;       // 1.5 × 31.7 Å

inline constexpr double LX_BOHR = LX_ANG * ANG_TO_BOHR;   // 34.771 bohr
inline constexpr double LY_BOHR = LY_ANG * ANG_TO_BOHR;   // 34.771 bohr
inline constexpr double LZ_BOHR = LZ_ANG * ANG_TO_BOHR;   // 89.856 bohr

inline auto make_cell() {
    return inq::systems::cell::orthorhombic(
        LX_BOHR*1.0_b, LY_BOHR*1.0_b, LZ_BOHR*1.0_b).finite();
}

// ── Coronene geometry ─────────────────────────────────────────────────────────
// Centred at (Lx/2, Ly/2, Lz/2) Å — all atom coords positive, within [0,L].
// Verified: z = 23.775 Å = Lz/2 for all atoms (flat molecule).
inline constexpr const char* CORONENE_XYZ = "coronene_centered.xyz";

// ── Electronic structure ─────────────────────────────────────────────────────
// C24H12: 24×C(4e) + 12×H(1e) = 108 electrons → 54 occupied KS orbitals
// Extra states: 1 WP slot + 2 SCF buffer = 3 total
// Total orbital count: 57  (indices 0–53 occupied, 54–55 buffer, 56 = WP)
inline constexpr int EXTRA_STATES = 3;

// E_cut = 40 Ha → h = π/√(2×40) ≈ 0.351 bohr = 0.186 Å
// Energy minimum from 03_ecut_convergence sweep.
inline constexpr double ECUT_HA = 40.0;

// SCF settings — Broyden ndim=8 α=0.1, tol=1e-4.
// History: tol=1e-4 matches run_003 working config. Broyden reaches de<1e-4
// at iter ~80-90, before the first restart glitch at iter ~130. For a sub-fs
// TDDFT wavepacket scattering, 1e-4 Ha accuracy is sufficient.
inline constexpr double SCF_TOL          = 1.0e-4;
inline constexpr double SCF_MIXING       = 0.1;
inline constexpr int    SCF_MIXING_NDIM  = 8;
inline constexpr int    SCF_MAX_STEPS    = 300;

// ── Wavepacket parameters ─────────────────────────────────────────────────────
// Paper Eq. 1: ψ^WP = (1/(πd²))^{3/4} exp(−|r−b|²/(2d²)) exp(ik·r)
// d = 0.53 Å (paper value): focused beam; illuminates inner coronene ring well.
// Width is narrow → WP spatial extent 2σ ≈ 1.06 Å in each direction.
inline constexpr double WP_D_ANG  = 0.53;
inline constexpr double WP_D_BOHR = WP_D_ANG * ANG_TO_BOHR;   // 1.0016 bohr

// Impact distance D = Lz/2 − 5d (maximised; 5d buffer from cell boundary)
// D = 23.775 − 2.65 = 21.125 Å = 39.921 bohr
// Compare run_003: D = 6.35 Å → 3.33× larger.
inline constexpr double WP_5D_ANG        = 5.0 * WP_D_ANG;          // 2.65 Å
inline constexpr double WP_D_IMPACT_ANG  = LZ_ANG/2.0 - WP_5D_ANG;  // 21.125 Å
inline constexpr double WP_D_IMPACT_BOHR = WP_D_IMPACT_ANG * ANG_TO_BOHR;  // 39.921 bohr

// Kinetic energy 200 eV (paper value)
inline constexpr double WP_EKIN_EV = 200.0;
inline constexpr double WP_EKIN_HA = WP_EKIN_EV / HA_TO_EV;   // 7.3499 Ha

// k₀ = √(2 E_kin)   (atomic units: m=1, ℏ=1)
inline double wp_k0()   { return std::sqrt(2.0 * WP_EKIN_HA); }   // 3.834 bohr⁻¹

// Normalisation factor (1/(π d²))^{3/4}
inline double wp_norm() { return std::pow(M_PI * WP_D_BOHR * WP_D_BOHR, -0.75); }

// WP centre — at cell centre laterally, D above molecule in +z
// WP travels in −z direction (k = (0, 0, −k₀)) toward molecule at z_flake
inline double WP_BX() { return LX_BOHR / 2.0; }
inline double WP_BY() { return LY_BOHR / 2.0; }
inline double WP_BZ() { return LZ_BOHR / 2.0 + WP_D_IMPACT_BOHR; }
// WP_BZ = 44.928 + 39.921 = 84.849 bohr = Lz − 5d = 89.856 − 5.007 bohr ✓

inline constexpr double WP_OCCUPATION = 1.0;  // singly occupied incident electron

// ── Observation and snapshot planes ──────────────────────────────────────────
// Coronene plane: z = Lz/2 (molecule centred in z)
inline double Z_FLAKE_BOHR() { return LZ_BOHR / 2.0; }                 // 44.928 bohr

// Observation plane: z = Lz/2 + D (same as WP starting position; LEED screen)
inline double Z_OBS_BOHR()   { return LZ_BOHR / 2.0 + WP_D_IMPACT_BOHR; }  // 84.849 bohr

// Mid plane: halfway between molecule and observation plane
inline double Z_MID_BOHR()   { return (Z_FLAKE_BOHR() + Z_OBS_BOHR()) / 2.0; }  // 64.889 bohr

// ── TDDFT propagation parameters ─────────────────────────────────────────────
// Δt = 4.84×10⁻⁴ fs (paper value)
inline constexpr double DT_FS  = 4.84e-4;
inline constexpr double DT_AU  = DT_FS * FS_TO_AU;               // 0.020009 a.u.

// WP arrival at coronene: t₁ = D / k₀ = 39.921 / 3.834 = 10.41 a.u. = 0.252 fs
// NOTE: T1 > run_003's T2 (0.25 fs). T2 must be > T1.
inline constexpr double T1_FS  = 0.252;
inline constexpr double T1_AU  = T1_FS * FS_TO_AU;               // 10.414 a.u.

// Total propagation time t₂ = 3×t₁ ≈ 0.756 fs
// Chosen to match paper's T2/T1 ≈ 3.3 ratio, and so that WP has time to reflect,
// travel back to z_obs, and accumulate a clean LEED pattern.
inline constexpr double T2_FS  = 0.756;
inline constexpr double T2_AU  = T2_FS * FS_TO_AU;               // 31.239 a.u.

inline constexpr int N_STEPS   = static_cast<int>(T2_AU / DT_AU);  // ~1561 steps

// Observable save intervals
inline constexpr int SNAPSHOT_INTERVAL = 10;
inline constexpr int MAX_SNAPSHOTS     = (N_STEPS / SNAPSHOT_INTERVAL) + 1;

// ── LEED accumulation ─────────────────────────────────────────────────────────
// Start accumulating when WP centre is 10σ from screen on outgoing trip.
// At 10σ distance, WP density at screen ≈ exp(-50) — negligible.
// This avoids the outgoing WP dominating the pattern at early times.
// t_leed_start = 10 * d / k₀
inline double T_LEED_START_AU() { return 10.0 * WP_D_BOHR / wp_k0(); }  // ~2.61 a.u.

// Three LEED screens: original z_obs (screen 0) plus two equally-spaced above it.
// Space between z_obs and cell wall = Lz − z_obs = 5d = 5 * WP_D_BOHR.
// Screen positions: z_obs, z_obs + 5d/3, z_obs + 10d/3
inline double Z_SCREEN0_BOHR() { return Z_OBS_BOHR(); }
inline double Z_SCREEN1_BOHR() { return Z_OBS_BOHR() + (LZ_BOHR - Z_OBS_BOHR()) / 3.0; }
inline double Z_SCREEN2_BOHR() { return Z_OBS_BOHR() + 2.0*(LZ_BOHR - Z_OBS_BOHR()) / 3.0; }

} // namespace cfg
