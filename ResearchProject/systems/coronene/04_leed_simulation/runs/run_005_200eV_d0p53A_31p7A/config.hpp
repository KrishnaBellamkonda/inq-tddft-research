#pragma once
// ============================================================================
// run_005 — Coronene TDDFT LEED simulation (paper parameters)
//
// Changes from run_004:
//   - Lz = 31.7 Å  (paper value; was 47.55 Å)
//   - D  = 6.35 Å  (paper value; was 21.125 Å)
//   - T1 = 0.077 fs (paper value: WP arrival at flake; was 0.252 fs)
//   - T2 = 0.25 fs  (paper value: WP reaches bottom boundary; was 0.756 fs)
//   - N_steps ≈ 517 (was 1561)
//   - LEED: ∫_{t1}^{t2} n_total(z_obs) dt  (paper Eq.5; no GS subtraction)
//   - One LEED screen at z_obs (removed screens 1 and 2)
//   - coronene z = 15.850 Å = Lz/2 (was 23.775 Å)
//
// Rationale for D=6.35 Å:
//   At D=21.125 Å (run_004), the WP disperses from 0.53 Å to ~5.5 Å by the time
//   it reaches the flake — much larger than the C-C bond length (1.42 Å) — blurring
//   all diffraction features. At D=6.35 Å, the WP spread at the flake is ~1.7 Å,
//   comparable to the hexagonal lattice scale → 6-fold spots are resolved.
//
// Rationale for T2=0.25 fs:
//   T2 is chosen so the simulation ends just as the transmitted WP reaches z=0
//   (the bottom boundary). During [t1, t2], only the backscattered component
//   contributes at z_obs — the transmitted WP is below z_flake and the finite-BC
//   boundary reflection has not yet returned to z_obs. This is the paper's
//   "WP exits box" condition.
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
// Paper cell: 18.4 × 18.4 × 31.7 Å³  (finite, isolated)
// INQ cell origin: (0,0,0) = cell CORNER; grid runs 0 → L.


//inline constexpr double LX_ANG  = 18.4;
//inline constexpr double LY_ANG  = 18.4;
//inline constexpr double LZ_ANG  = 31.7;    // paper value

// Changing the z value for the size of the box
// as observing that the wavepacket is dispersing before it reaches
// the lattice
inline constexpr double LX_ANG  = 18.4;
inline constexpr double LY_ANG  = 18.4;
inline constexpr double LZ_ANG  = 5;    // paper value


inline constexpr double LX_BOHR = LX_ANG * ANG_TO_BOHR;   // 34.771 bohr
inline constexpr double LY_BOHR = LY_ANG * ANG_TO_BOHR;   // 34.771 bohr
inline constexpr double LZ_BOHR = LZ_ANG * ANG_TO_BOHR;   // 59.904 bohr

inline auto make_cell() {
    return inq::systems::cell::orthorhombic(
        LX_BOHR*1.0_b, LY_BOHR*1.0_b, LZ_BOHR*1.0_b).finite();
}

// ── Coronene geometry ─────────────────────────────────────────────────────────
// Centred at (Lx/2, Ly/2, Lz/2) Å = (9.2, 9.2, 15.85) Å
inline constexpr const char* CORONENE_XYZ = "coronene_centered.xyz";

// ── Electronic structure ─────────────────────────────────────────────────────
inline constexpr int    EXTRA_STATES     = 3;
inline constexpr double ECUT_HA          = 40.0;
inline constexpr double SCF_TOL          = 1.0e-4;
inline constexpr double SCF_MIXING       = 0.1;
inline constexpr int    SCF_MIXING_NDIM  = 8;
inline constexpr int    SCF_MAX_STEPS    = 300;

// ── Wavepacket parameters ─────────────────────────────────────────────────────
// d = 0.53 Å = 1.0016 bohr  (paper value)
inline constexpr double WP_D_ANG  = 0.53;
inline constexpr double WP_D_BOHR = WP_D_ANG * ANG_TO_BOHR;   // 1.0016 bohr

// D = 6.35 Å  (paper value)
inline constexpr double WP_D_IMPACT_ANG  = 6.35;
inline constexpr double WP_D_IMPACT_BOHR = WP_D_IMPACT_ANG * ANG_TO_BOHR;  // 12.000 bohr

inline constexpr double WP_EKIN_EV = 200.0;
inline constexpr double WP_EKIN_HA = WP_EKIN_EV / HA_TO_EV;   // 7.3499 Ha

inline double wp_k0()   { return std::sqrt(2.0 * WP_EKIN_HA); }   // 3.834 bohr⁻¹
inline double wp_norm() { return std::pow(M_PI * WP_D_BOHR * WP_D_BOHR, -0.75); }

// WP centre: (Lx/2, Ly/2, z_flake + D)  — WP travels in −z toward molecule
inline double WP_BX() { return LX_BOHR / 2.0; }
inline double WP_BY() { return LY_BOHR / 2.0; }
inline double WP_BZ() { return LZ_BOHR / 2.0 + WP_D_IMPACT_BOHR; }
// WP_BZ = 29.952 + 12.000 = 41.952 bohr = 22.200 Å ✓

inline constexpr double WP_OCCUPATION = 1.0;

// ── Observation and snapshot planes ──────────────────────────────────────────
// Molecule at z = Lz/2 = 15.85 Å = 29.952 bohr
inline double Z_FLAKE_BOHR() { return LZ_BOHR / 2.0; }

// Observation plane = WP start = z_flake + D = 22.20 Å = 41.952 bohr
inline double Z_OBS_BOHR()   { return LZ_BOHR / 2.0 + WP_D_IMPACT_BOHR; }

// Mid plane
inline double Z_MID_BOHR()   { return (Z_FLAKE_BOHR() + Z_OBS_BOHR()) / 2.0; }

// ── TDDFT propagation parameters ─────────────────────────────────────────────
// Δt = 4.84×10⁻⁴ fs  (paper value)
inline constexpr double DT_FS  = 4.84e-4;
inline constexpr double DT_AU  = DT_FS * FS_TO_AU;                   // 0.020009 a.u.

// t1 = D/k0 = 12.000/3.834 = 3.129 a.u. = 0.0757 fs  (paper: 0.077 fs)
// This is when the WP centre arrives at z_flake. LEED accumulation begins here.
inline constexpr double T1_FS  = 0.077;
inline constexpr double T1_AU  = T1_FS * FS_TO_AU;                   // 3.184 a.u.

// t2 = 0.25 fs  (paper value: WP reaches z=0 boundary)
// At T2, WP centre is at z_obs - k0*T2 = 41.952 - 3.834×10.336 ≈ 2.3 bohr
// — just at the bottom boundary. Simulation ends here.
inline constexpr double T2_FS  = 0.25;
inline constexpr double T2_AU  = T2_FS * FS_TO_AU;                   // 10.336 a.u.

inline constexpr int N_STEPS   = static_cast<int>(T2_AU / DT_AU);    // ~517 steps

inline constexpr int SNAPSHOT_INTERVAL = 10;
inline constexpr int MAX_SNAPSHOTS     = (N_STEPS / SNAPSHOT_INTERVAL) + 1;

// LEED screen: z_obs (WP starting position)
inline double Z_SCREEN_BOHR()  { return Z_OBS_BOHR(); }

} // namespace cfg
