// Pure-tier test for inqkit::observables::DensityDelta (T15).
//
// Question (TODO at density_delta.hpp:135): does delta at step n use step n-1 as
// the base (rolling), or is t=0 the fixed base for ALL steps?
//
// The code computes delta = current - ref_, where ref_ is captured ONCE
// (set_reference or lazy capture on the first snapshot). This test pins that
// behaviour: it picks values so a rolling base (current - previous) would yield
// a DIFFERENT L2 than the fixed-t0 base (current - reference). The returned L2
// is sigma^2 = integral |delta|^2 dV.
//
// Pure: with emit_raw_vti/emit_coarse_vti = false, snapshot() performs no file
// I/O (the writer's create_directories lives in write(), never called here).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inqkit/observables/density_delta.hpp>

#include <vector>

using inqkit::fields::RealField3D;
using inqkit::observables::DensityDelta;
using inqkit::observables::DensityDeltaConfig;
using Catch::Approx;

namespace {

// Build a RealField3D with given spacing and flat values (nx = values.size()).
RealField3D make_field(std::vector<double> vals, double dx, double dy, double dz) {
  RealField3D f;
  f.nx = static_cast<int>(vals.size());
  f.ny = 1;
  f.nz = 1;
  f.dx_bohr = dx;
  f.dy_bohr = dy;
  f.dz_bohr = dz;
  f.values = std::move(vals);
  return f;
}

// A DensityDelta that never writes files: pure compute path only.
DensityDelta make_compute_only() {
  DensityDeltaConfig cfg;
  cfg.emit_raw_vti = false;
  cfg.emit_coarse_vti = false;
  cfg.compute_l2 = true;
  return DensityDelta("/dev/null/raw_unused", "/dev/null/coarse_unused", cfg);
}

} // namespace

TEST_CASE("DensityDelta: first snapshot lazy-captures reference and returns zero", "[density_delta][pure]") {
  auto dd = make_compute_only();
  CHECK_FALSE(dd.has_reference());

  auto a = make_field({1.0}, 1.0, 1.0, 1.0);
  double l2 = dd.snapshot(a, /*time_au=*/0.0, /*step=*/0);

  CHECK(dd.has_reference());
  CHECK(l2 == Approx(0.0));  // t=t0 frame is zero by construction
}

TEST_CASE("DensityDelta: t=0 is the fixed base for ALL steps (not rolling)", "[density_delta][pure]") {
  auto dd = make_compute_only();

  // dV = 1 (unit spacing), single cell, so L2 = |current - ref|^2.
  auto A = make_field({1.0}, 1.0, 1.0, 1.0);  // reference (t0)
  auto B = make_field({3.0}, 1.0, 1.0, 1.0);
  auto C = make_field({7.0}, 1.0, 1.0, 1.0);

  CHECK(dd.snapshot(A, 0.0, 0) == Approx(0.0));   // captures ref = 1
  CHECK(dd.snapshot(B, 1.0, 1) == Approx(4.0));   // |3 - 1|^2 = 4

  // The discriminator: fixed-t0 base => |7 - 1|^2 = 36.
  //                    rolling base   => |7 - 3|^2 = 16.
  double l2_C = dd.snapshot(C, 2.0, 2);
  CHECK(l2_C == Approx(36.0));   // fixed t0 base
  CHECK(l2_C != Approx(16.0));   // NOT a rolling t->t+1 base
}

TEST_CASE("DensityDelta: explicit set_reference matches lazy capture", "[density_delta][pure]") {
  auto dd = make_compute_only();
  auto A = make_field({1.0}, 1.0, 1.0, 1.0);
  dd.set_reference(A);
  CHECK(dd.has_reference());

  auto C = make_field({7.0}, 1.0, 1.0, 1.0);
  // No lazy zero-frame consumed: first snapshot already deltas against A.
  CHECK(dd.snapshot(C, 0.0, 0) == Approx(36.0));
}

TEST_CASE("DensityDelta: L2 sums cells and scales with dV", "[density_delta][pure]") {
  auto dd = make_compute_only();

  // Two cells, dV = 0.5 * 1 * 1 = 0.5.
  auto A = make_field({1.0, 1.0}, 0.5, 1.0, 1.0);
  auto C = make_field({7.0, 4.0}, 0.5, 1.0, 1.0);

  dd.snapshot(A, 0.0, 0);                    // ref = {1,1}
  // delta = {6, 3}; sum |delta|^2 = 36 + 9 = 45; * dV(0.5) = 22.5.
  CHECK(dd.snapshot(C, 1.0, 1) == Approx(22.5));
}

TEST_CASE("DensityDelta: grid mismatch against the reference throws", "[density_delta][pure]") {
  auto dd = make_compute_only();
  auto A = make_field({1.0, 1.0}, 1.0, 1.0, 1.0);   // 2 cells
  auto B = make_field({3.0}, 1.0, 1.0, 1.0);        // 1 cell

  dd.snapshot(A, 0.0, 0);                            // ref has 2 cells
  CHECK_THROWS(dd.snapshot(B, 1.0, 1));
}
