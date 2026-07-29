// Engine-tier single-rank characterization of inqkit::screens::PlaneScreen::
// extract: the z=0 plane of a He atom at the origin. Asserts the slice has the
// right shape, is a non-negative density, and is non-trivial (peak > 0).
// (The multi-rank E01 Allreduce bug is covered separately by the [!shouldfail]
// cross-rank test launched under mpirun.)

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/screens/plane_screen.hpp>

#include <algorithm>
#include <utility>
#include <vector>

using namespace inq;
using namespace inq::magnitude;

TEST_CASE("PlaneScreen::extract: single-rank z=0 slice of He", "[screens][plane_screen][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(ions, options::electrons{}.spacing(0.5 * 1.0_bohr));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  inqkit::screens::PlaneScreen screen(0.0, "z0");
  auto slice = screen.extract(electrons);

  const int nx = screen.nx(electrons);
  const int ny = screen.ny(electrons);
  REQUIRE(static_cast<int>(slice.size()) == ny);
  REQUIRE(static_cast<int>(slice[0].size()) == nx);

  double peak = 0.0, total = 0.0;
  for (auto const &row : slice)
    for (double v : row) {
      CHECK(v >= 0.0);          // density slice is non-negative
      peak = std::max(peak, v);
      total += v;
    }
  CHECK(peak > 0.0);            // non-trivial
  CHECK(total > 0.0);
}

// D2: the generalised axis. He at the origin is spherically symmetric, so the
// x=0, y=0 and z=0 planes carry the same density (same peak and total) — a clean
// analytic check that extract() handles every axis, not just z.
TEST_CASE("PlaneScreen::extract: x/y/z slices of a symmetric atom agree (D2)",
          "[screens][plane_screen][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(ions, options::electrons{}.spacing(0.5 * 1.0_bohr));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto stats = [](std::vector<std::vector<double>> const &s) {
    double peak = 0.0, total = 0.0;
    for (auto const &row : s) for (double v : row) { peak = std::max(peak, v); total += v; }
    return std::pair<double, double>{peak, total};
  };

  auto [px, tx] = stats(inqkit::screens::PlaneScreen(0.0, "x0", 0).extract(electrons));
  auto [py, ty] = stats(inqkit::screens::PlaneScreen(0.0, "y0", 1).extract(electrons));
  auto [pz, tz] = stats(inqkit::screens::PlaneScreen(0.0, "z0", 2).extract(electrons));

  CHECK(px > 0.0);
  CHECK(py == Catch::Approx(px).epsilon(1e-9));    // spherical symmetry (peaks exact)
  CHECK(pz == Catch::Approx(px).epsilon(1e-9));
  // Plane totals agree to ~8 sig figs; the residual is FP summation-order noise
  // across the three per-axis traversals, not a physics difference.
  CHECK(ty == Catch::Approx(tx).epsilon(1e-6));
  CHECK(tz == Catch::Approx(tx).epsilon(1e-6));
}

// D2: time-averaged screen. ⟨ρ⟩ of N identical frames (each dt) equals the frame.
TEST_CASE("TimeAveragedScreen: average of constant frames is the frame (D2)",
          "[screens][plane_screen][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(ions, options::electrons{}.spacing(0.5 * 1.0_bohr));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto frame = inqkit::screens::PlaneScreen(0.0, "z0").extract(electrons);
  inqkit::screens::TimeAveragedScreen avg;
  for (int i = 0; i < 4; ++i) avg.add(frame, 0.25);     // 4 frames, dt=0.25
  CHECK(avg.total_time() == Catch::Approx(1.0));
  auto mean = avg.average();
  REQUIRE(mean.size() == frame.size());
  for (std::size_t i = 0; i < frame.size(); ++i)
    for (std::size_t j = 0; j < frame[i].size(); ++j)
      CHECK(mean[i][j] == Catch::Approx(frame[i][j]).margin(1e-12));
}
