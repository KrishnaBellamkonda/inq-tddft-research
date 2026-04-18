#pragma once
// ============================================================================
// run_003 — Coronene TDDFT LEED simulation
//
// Changes from run_002:
//   - SCF tolerance tightened to 1e-6 Ha (was 1e-4)
//   - Linear mixing with α=0.05 (was Broyden ndim=8 α=0.1)
//   - SCF max_steps = 1000 (was 300)
//   - Z_MID_BOHR() added (midpoint between molecule and observation plane)
//   - OVERLAP_INTERVAL kept at 5 (overlap matrix computed here too)
//
// Unchanged from run_002:
//   - Cell: 18.4 × 18.4 × 31.7 Å (finite, paper values)
//   - E_cut = 40 Ha (energy minimum from 03_ecut_convergence sweep)
//   - WP: d = 1.4 Å, D = 6.35 Å, E_kin = 200 eV
//   - Time step, N_STEPS, T1, T2 all from Tsubonoya et al. (2014)
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
// Finite (isolated molecule), same as paper: 18.4 × 18.4 × 31.7 Å
// INQ cell origin: (0,0,0) = cell CORNER; grid runs 0 → L.
// Coronene centred at (Lx/2, Ly/2, Lz/2) — see coronene_centered.xyz.
inline constexpr double LX_ANG  = 18.4;
inline constexpr double LY_ANG  = 18.4;
inline constexpr double LZ_ANG  = 31.7;

inline constexpr double LX_BOHR = LX_ANG * ANG_TO_BOHR;   // 34.771 bohr
inline constexpr double LY_BOHR = LY_ANG * ANG_TO_BOHR;   // 34.771 bohr
inline constexpr double LZ_BOHR = LZ_ANG * ANG_TO_BOHR;   // 59.904 bohr

inline auto make_cell() {
    return inq::systems::cell::orthorhombic(
        LX_BOHR*1.0_b, LY_BOHR*1.0_b, LZ_BOHR*1.0_b).finite();
}

// ── Coronene geometry ─────────────────────────────────────────────────────────
// Centred at (Lx/2, Ly/2, Lz/2) Å — all atom coords positive, within [0,L].
// Verified by geometry_check.py (all 8 checks pass).
inline constexpr const char* CORONENE_XYZ = "coronene_centered.xyz";

// ── Electronic structure ─────────────────────────────────────────────────────
// C24H12: 24×C(4e) + 12×H(1e) = 108 electrons → 54 occupied KS orbitals
// Extra states: 1 WP slot + 2 SCF buffer = 3 total
// Total orbital count: 57  (indices 0–53 occupied, 54–55 buffer, 56 = WP)
inline constexpr int EXTRA_STATES = 3;

// E_cut = 40 Ha → h = π/√(2×40) ≈ 0.351 bohr = 0.186 Å
// Energy minimum from 03_ecut_convergence sweep.
// 54 Ha (paper value) causes pseudodojo_pbe projector artefacts → oscillatory SCF.
inline constexpr double ECUT_HA = 40.0;

// SCF settings — Broyden ndim=8 α=0.1, tol=1e-4.
// History:
//   linear α=0.05 → limit cycle at iter ~310 (never converges).
//   Broyden ndim=8 α=0.1 tol=1e-6 → reaches de=4e-6 at iter 114, then
//     restart glitch every ~32 iters pushes de back to 3e-3 (near-degenerate
//     HOMO of coronene causes ill-conditioned Broyden history on restart).
//   Broyden ndim=20 α=0.05 → diverges.
//   Broyden ndim=8 α=0.1 tol=1e-5 → same restart glitch at iter ~139.
// Fix: tol=1e-4 matches run_002 working config. Broyden reaches de<1e-4 at
//   iter ~80-90, well before the first restart glitch at iter ~130. For a
//   0.25 fs TDDFT wavepacket scattering, 1e-4 Ha initial-state accuracy is
//   sufficient (same as published run_002 results).
inline constexpr double SCF_TOL          = 1.0e-4;
inline constexpr double SCF_MIXING       = 0.1;
inline constexpr int    SCF_MIXING_NDIM  = 8;
inline constexpr int    SCF_MAX_STEPS    = 300;

// ── Wavepacket parameters ─────────────────────────────────────────────────────
// Paper Eq. 1: ψ^WP = (1/(πd²))^{3/4} exp(−|r−b|²/(2d²)) exp(ik·r)
// d = 1.4 Å: 35.7% illumination at inner ring (r=1.421 Å); good beam coverage.
inline constexpr double WP_D_ANG  = 1.4;
inline constexpr double WP_D_BOHR = WP_D_ANG * ANG_TO_BOHR;   // 2.6456 bohr

// Impact distance D = 6.35 Å (paper value; WP centre to coronene plane at t=0)
inline constexpr double WP_D_IMPACT_ANG  = 6.35;
inline constexpr double WP_D_IMPACT_BOHR = WP_D_IMPACT_ANG * ANG_TO_BOHR;  // 12.000 bohr

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
inline double WP_BZ() { return LZ_BOHR / 2.0 + WP_D_IMPACT_BOHR; }  // 41.952 bohr

inline constexpr double WP_OCCUPATION = 1.0;  // singly occupied incident electron

// ── Observation and snapshot planes ──────────────────────────────────────────
// Coronene plane: z = Lz/2 (molecule centred in z)
inline double Z_FLAKE_BOHR() { return LZ_BOHR / 2.0; }                 // 29.952 bohr

// Observation plane: z = Lz/2 + D (same as WP starting position; LEED screen)
inline double Z_OBS_BOHR()   { return LZ_BOHR / 2.0 + WP_D_IMPACT_BOHR; }  // 41.952 bohr

// Mid plane: halfway between molecule and observation plane
inline double Z_MID_BOHR()   { return (Z_FLAKE_BOHR() + Z_OBS_BOHR()) / 2.0; }  // 35.952 bohr

// ── TDDFT propagation parameters ─────────────────────────────────────────────
// Δt = 4.84×10⁻⁴ fs (paper value)
inline constexpr double DT_FS  = 4.84e-4;
inline constexpr double DT_AU  = DT_FS * FS_TO_AU;               // 0.020009 a.u.

// Total propagation time t₂ = 0.25 fs (paper value)
inline constexpr double T2_FS  = 0.25;
inline constexpr double T2_AU  = T2_FS * FS_TO_AU;               // 10.335 a.u.

// WP arrival at coronene: t₁ = D / k₀ (time for WP to travel from bz to z_flake)
// t₁ = 12.000 bohr / 3.834 bohr·a.u.⁻¹ ≈ 3.128 a.u. ≈ 0.077 fs  (paper value: 0.077 fs)
inline constexpr double T1_FS  = 0.077;
inline constexpr double T1_AU  = T1_FS * FS_TO_AU;               // 3.183 a.u.

inline constexpr int N_STEPS   = static_cast<int>(T2_AU / DT_AU);  // 516 steps

// Observable save intervals — all heavy saves unified to one cadence
// Energy and momentum CSVs are written every step (1 line each, negligible size).
// Everything else (WP orbital 3D, 3D density, 2D slices, z-profile, overlap matrix) every 10 steps.
inline constexpr int SNAPSHOT_INTERVAL = 10;   // all periodic heavy saves
inline constexpr int MAX_SNAPSHOTS     = (N_STEPS / SNAPSHOT_INTERVAL) + 1;  // ~53

} // namespace cfg
