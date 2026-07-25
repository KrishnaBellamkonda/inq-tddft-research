/* inqkit::absorbers — mask-function shape helpers (INQ-free, host-only).
 *
 * The pure sin^2 mask profiles of De Giovannini, Larsen & Rubio,
 * "Modeling electron dynamics coupled to continuum states in finite volumes",
 * arXiv:1409.1689 (2014), Eq. 13 — split out of mask_absorber.hpp so they carry
 * NO <inq/...> dependency and can be unit-tested in the pure (no-engine) tier.
 * The GPU kernels in mask_absorber.hpp inline the identical expressions with
 * device-safe math.
 */
#pragma once

#include <cmath>

namespace inqkit {
namespace absorbers {

// Single-sided sin^2 mask M(s) (Eq. 13): M=1 for s<=z_abs0, ramps to 0 across
// [z_abs0, z_abs0+L]. z_abs0 = absorber start, L = absorber width.
inline double sin2_mask_value(double s, double z_abs0, double L) {
  if (s <= z_abs0) return 1.0;
  if (s >= z_abs0 + L) return 0.0;
  double sn = std::sin(M_PI * (s - z_abs0) / (2.0 * L));
  return 1.0 - sn * sn;
}

// Two-sided sin^2 mask, symmetric about the box centre (s=0): the single-sided
// Eq.13 shape applied to |s|. M=1 in the inner region |s|<z_in, ramps to 0 across
// |s| in [z_in, z_in+Lhalf] at BOTH ends, M=0 beyond. z_in = inner-region
// half-width, Lhalf = per-end width (= L_total/2). A packet leaving either
// boundary sees the same ramp.
inline double sin2_mask_value_twosided(double s, double z_in, double Lhalf) {
  double a = std::fabs(s);
  if (a <= z_in) return 1.0;
  if (a >= z_in + Lhalf) return 0.0;
  double sn = std::sin(M_PI * (a - z_in) / (2.0 * Lhalf));
  return 1.0 - sn * sn;
}

} // namespace absorbers
} // namespace inqkit
