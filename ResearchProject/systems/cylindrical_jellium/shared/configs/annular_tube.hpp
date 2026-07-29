// shared/configs/annular_tube.hpp
// ----------------------------------------------------------------------------
// LOCKED geometry table for the cylindrical-jellium S(v) campaign
// (docs/campaigns/cylindrical_jellium/cylindrical_jellium_projectile.md,
//  design frozen by the user 2026-06-27).
//
// A PERIODIC annular jellium tube (positive background between concentric
// cylinders R_in < d < R_out, axis ∥ z, hollow bore), with a classical electron
// projectile gliding on-axis. Transverse box + radii are FIXED across densities
// for comparability; L_z is sized PER DENSITY to ≥ 2× the v=0.45 wake length
// λ=2πv/ω_p (ω_p=√(3/r_s³)). S is per-unit-length on a z-uniform tube, so a
// per-density L_z does not break S(r_s) comparability.
//
// N is rounded to the nearest EVEN integer; the run then sets n0 = N/V_annulus so
// ∫n₊ = N EXACTLY (exact neutrality; the G=0 cancellation requirement). This file
// is DOCUMENTATION + compile-time reference — the run binaries are pure-env
// (the Python orchestrator passes every value), mirroring campaign_autorun.
//
//   density   r_s   L_z(Bohr)   N(even)   n0=N/V_annulus (a0^-3)
//   ------    ---   ---------   -------   ----------------------
//   dilute     6       48          24      ~1.105e-3
//   mid        4       28          48      ~3.79e-3
//   dense      2       10         136      ~3.01e-2
//
// Fixed across all densities:
//   R_in = 5 Bohr, R_out = 13 Bohr (8 Bohr wall), L_xy = 40 Bohr, dx = 0.5 Bohr,
//   edge_width w = 1.0 Bohr, fully periodic cell, tube axis = z (slab_axis=2).
//
// Projectile (classical electron, EHRENFEST): charge −1 carried as fictitious
// "H", mass = m_e (1/1822.8885 amu), Gaussian erf-smoothed UPF
// electron_gaussian_wpsigma0p5.upf (VERIFIED 2026-06-28: V(0)=+2.257 Ha
// repulsive, Z=1.000, σ_pot=0.3536 = σ_WP/√2 for σ_WP=0.5). Velocities
// v = {0.15, 0.30, 0.45} a.u.; dt = 0.020 a.u.; real-time LDA, ETRS electrons +
// Ehrenfest ion.
// ----------------------------------------------------------------------------
#pragma once

#include <cmath>

namespace cylindrical_jellium::config {

// Fixed transverse geometry (Bohr).
inline constexpr double R_IN      = 5.0;
inline constexpr double R_OUT     = 13.0;
inline constexpr double L_XY      = 40.0;
inline constexpr double SPACING   = 0.5;
inline constexpr double EDGE_W    = 1.0;
inline constexpr int    TUBE_AXIS = 2;   // z

// Projectile + propagation.
inline constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;  // m_e
inline constexpr double DT_AU         = 0.020;
inline constexpr const char* PROJ_SPECIES_SYMBOL = "H";
inline constexpr const char* PROJ_PSEUDO_PATH =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/cylindrical_jellium/"
    "shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

// Per-density axial length + electron count (see table above).
struct Density { double r_s; double L_z; int N; };
inline constexpr Density RS6{6.0, 48.0, 24};
inline constexpr Density RS4{4.0, 28.0, 48};
inline constexpr Density RS2{2.0, 10.0, 136};

// V_annulus = π (R_out² − R_in²) L_z ; n0 set so ∫n₊ = N exactly.
inline constexpr double annulus_volume(double L_z) {
    return M_PI * (R_OUT * R_OUT - R_IN * R_IN) * L_z;
}
inline constexpr double n0_for(int N, double L_z) {
    return double(N) / annulus_volume(L_z);
}

}  // namespace cylindrical_jellium::config
