// Pure-tier characterization + validation test for
// inqkit::observables::center_of_density (supports T07/T09/T11/T12).
//
// Locks the CURRENT behaviour (struct CenterOfDensityResult{x_bohr,y_bohr,
// z_bohr,total_weight}) so the later Vec3 struct-swap can be proven
// behaviour-preserving. Also validates analytically:
//   * coordinates are in Bohr via x = origin + (ix+0.5)*dx (cell-centre)  [T11]
//   * total_weight = integral f dV                                        [T12]
//   * the w>0 guard keeps a zero field at the origin                      [T12]
//   * axes are not transposed (centroid with three DISTINCT coords)

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inqkit/observables/center_of_density.hpp>
#include <inqkit/detail/grid_layout.hpp>

#include <vector>

using inqkit::fields::RealField3D;
using inqkit::observables::center_of_density;
using inqkit::detail::grid_layout::flatten_index;
using Catch::Approx;

namespace {
RealField3D make_field(int nx, int ny, int nz,
                       double ox, double oy, double oz,
                       double dx, double dy, double dz,
                       std::vector<double> vals) {
  RealField3D f;
  f.nx = nx; f.ny = ny; f.nz = nz;
  f.origin_x_bohr = ox; f.origin_y_bohr = oy; f.origin_z_bohr = oz;
  f.dx_bohr = dx; f.dy_bohr = dy; f.dz_bohr = dz;
  f.values = std::move(vals);
  return f;
}
} // namespace

TEST_CASE("center_of_density: uniform 1D field centroid sits at geometric centre", "[center_of_density][pure]") {
  // 3 cells in x, origin -1.5, spacing 1 -> cell centres x = -1, 0, +1; y=z=0.
  auto f = make_field(3, 1, 1, -1.5, -0.5, -0.5, 1, 1, 1, {1.0, 1.0, 1.0});
  auto r = center_of_density(f);
  CHECK(r.center_bohr.x == Approx(0.0));
  CHECK(r.center_bohr.y == Approx(0.0));
  CHECK(r.center_bohr.z == Approx(0.0));
  CHECK(r.total_weight == Approx(3.0));  // integral f dV = 3 * 1
}

TEST_CASE("center_of_density: density-weighted centroid (Bohr units)", "[center_of_density][pure]") {
  // f = {1,0,3} over x = {-1,0,+1}: <x> = (1*-1 + 3*+1)/4 = 0.5 Bohr.
  auto f = make_field(3, 1, 1, -1.5, -0.5, -0.5, 1, 1, 1, {1.0, 0.0, 3.0});
  auto r = center_of_density(f);
  CHECK(r.center_bohr.x == Approx(0.5));
  CHECK(r.total_weight == Approx(4.0));
}

TEST_CASE("center_of_density: zero field stays at origin via the w>0 guard", "[center_of_density][pure]") {
  auto f = make_field(3, 1, 1, -1.5, -0.5, -0.5, 1, 1, 1, {0.0, 0.0, 0.0});
  auto r = center_of_density(f);
  CHECK(r.total_weight == Approx(0.0));
  CHECK(r.center_bohr.x == Approx(0.0));  // guard: no divide-by-zero, default 0
  CHECK(r.center_bohr.y == Approx(0.0));
  CHECK(r.center_bohr.z == Approx(0.0));
}

TEST_CASE("center_of_density: axes are not transposed (three distinct coords)", "[center_of_density][pure]") {
  // 2x2x2 grid with per-axis offset origins so cell centres differ per axis:
  //   x centres {0,1}, y centres {1,2}, z centres {2,3}.
  // Single unit mass at (ix,iy,iz) = (1,1,1) -> centroid (1,2,3), all distinct.
  std::vector<double> vals(8, 0.0);
  vals[flatten_index(1, 1, 1, /*ny=*/2, /*nz=*/2)] = 1.0;  // flat = 7
  auto f = make_field(2, 2, 2, -0.5, 0.5, 1.5, 1, 1, 1, std::move(vals));
  auto r = center_of_density(f);
  CHECK(r.center_bohr.x == Approx(1.0));
  CHECK(r.center_bohr.y == Approx(2.0));
  CHECK(r.center_bohr.z == Approx(3.0));
  CHECK(r.total_weight == Approx(1.0));
}

TEST_CASE("center_of_density: empty field throws", "[center_of_density][pure]") {
  RealField3D f;  // default: nx=ny=nz=0, empty values
  CHECK_THROWS(center_of_density(f));
}
