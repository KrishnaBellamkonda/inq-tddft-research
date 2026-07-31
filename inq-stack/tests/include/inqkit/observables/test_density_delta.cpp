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

#include <cstddef>
#include <filesystem>
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

// REGRESSION (2026-07-31). snapshot() is normally called EVERY step because the
// L2 it returns is an every-step scalar, but each call also wrote a full delta
// FIELD. In the wp_highdensity_sv campaign that was 3624 x 18 MB per run instead
// of 302 -- ~66 GB each, ~490 GB across the campaign. It exhausted the 1.0 TB
// /rds quota and killed three sigma=3 runs mid-flight, their vacuum controls and
// the notebook job. cfg.emit_every gates the FIELD write only.
TEST_CASE("DensityDelta: emit_every does not gate the L2", "[density_delta][pure]") {
  DensityDeltaConfig cfg;
  cfg.emit_raw_vti = false;      // compute path only
  cfg.emit_coarse_vti = false;
  cfg.compute_l2 = true;
  cfg.emit_every = 12;           // would suppress writes, must not suppress L2
  DensityDelta dd("/dev/null/raw_unused", "/dev/null/coarse_unused", cfg);

  auto A = make_field({1.0}, 1.0, 1.0, 1.0);
  auto C = make_field({7.0}, 1.0, 1.0, 1.0);
  dd.snapshot(A, 0.0, 0);                       // reference

  // Steps 1..11 are all OFF the emit cadence; every one must still return the L2.
  for (int step = 1; step < 12; ++step) {
    CHECK(dd.snapshot(C, 0.04 * step, step) == Approx(36.0));  // (7-1)^2 * dV=1
  }
  CHECK(dd.snapshot(C, 0.48, 12) == Approx(36.0));             // and on-cadence
}

TEST_CASE("DensityDelta: emit_every gates the VTI field writes", "[density_delta][pure]") {
  namespace fs = std::filesystem;
  const auto root = fs::temp_directory_path() / "inqkit_dd_emit_every_test";
  fs::remove_all(root);
  const auto raw = root / "raw", coarse = root / "coarse";

  auto count_vti = [](fs::path const& d) {
    std::size_t n = 0;
    if (fs::is_directory(d))
      for (auto const& e : fs::directory_iterator(d))
        if (e.path().extension() == ".vti") ++n;
    return n;
  };

  DensityDeltaConfig cfg;
  cfg.emit_raw_vti = true;
  cfg.emit_coarse_vti = false;   // coarse-graining needs a larger grid; raw suffices
  cfg.compute_l2 = true;
  cfg.emit_every = 5;
  DensityDelta dd(raw.string(), coarse.string(), cfg);

  auto A = make_field({1.0}, 1.0, 1.0, 1.0);
  auto C = make_field({7.0}, 1.0, 1.0, 1.0);
  for (int step = 0; step <= 20; ++step) dd.snapshot(step == 0 ? A : C, 0.04 * step, step);

  // Steps 0,5,10,15,20 emit; the other 16 do not.
  CHECK(count_vti(raw) == 5u);

  fs::remove_all(root);
}

TEST_CASE("DensityDelta: emit_every defaults to writing every call", "[density_delta][pure]") {
  namespace fs = std::filesystem;
  const auto root = fs::temp_directory_path() / "inqkit_dd_emit_default_test";
  fs::remove_all(root);
  const auto raw = root / "raw", coarse = root / "coarse";

  // The default MUST stay 1: ~130 existing run.cpp callers rely on it.
  DensityDeltaConfig cfg;
  cfg.emit_raw_vti = true;
  cfg.emit_coarse_vti = false;
  CHECK(cfg.emit_every == 1);

  DensityDelta dd(raw.string(), coarse.string(), cfg);
  auto A = make_field({1.0}, 1.0, 1.0, 1.0);
  auto C = make_field({7.0}, 1.0, 1.0, 1.0);
  for (int step = 0; step <= 6; ++step) dd.snapshot(step == 0 ? A : C, 0.04 * step, step);

  std::size_t n = 0;
  for (auto const& e : fs::directory_iterator(raw))
    if (e.path().extension() == ".vti") ++n;
  CHECK(n == 7u);

  fs::remove_all(root);
}

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
