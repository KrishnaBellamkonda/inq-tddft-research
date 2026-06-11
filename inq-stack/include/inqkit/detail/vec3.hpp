/*
 * Vec3 — a tiny pure 3-component double vector for inqkit.
 *
 * Purpose: give COD / current / dipole / centre coordinates a single vector
 * *unit* instead of loose x/y/z scalars (TODOs T07, T09), WITHOUT pulling
 * <inq> into otherwise-pure headers (center_of_density.hpp, observables_writer.
 * hpp). The API deliberately mirrors a subset of inq::vector3 (dot/norm and the
 * usual arithmetic) so a later switch to INQ's type — if those headers ever
 * become engine-tier — is mechanical.
 *
 * Pure: depends only on <cmath>. Tested by tests/cpp/test_vec3.cpp.
 */
#pragma once

#include <cmath>

namespace inqkit::detail {

struct Vec3 {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;

  // Index access (0=x, 1=y, 2=z) so call sites that used loose x/y/z arrays
  // (observables_writer's ctx.current[0..2]) work unchanged after the switch.
  double &operator[](int i) { return i == 0 ? x : (i == 1 ? y : z); }
  double operator[](int i) const { return i == 0 ? x : (i == 1 ? y : z); }

  double dot(Vec3 const &o) const { return x * o.x + y * o.y + z * o.z; }
  double norm2() const { return dot(*this); }
  double norm() const { return std::sqrt(norm2()); }

  Vec3 &operator+=(Vec3 const &o) { x += o.x; y += o.y; z += o.z; return *this; }
  Vec3 &operator-=(Vec3 const &o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
  Vec3 &operator*=(double s) { x *= s; y *= s; z *= s; return *this; }
};

inline Vec3 operator+(Vec3 a, Vec3 const &b) { a += b; return a; }
inline Vec3 operator-(Vec3 a, Vec3 const &b) { a -= b; return a; }
inline Vec3 operator*(double s, Vec3 v) { v *= s; return v; }
inline Vec3 operator*(Vec3 v, double s) { v *= s; return v; }

inline bool operator==(Vec3 const &a, Vec3 const &b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(Vec3 const &a, Vec3 const &b) { return !(a == b); }

} // namespace inqkit::detail
