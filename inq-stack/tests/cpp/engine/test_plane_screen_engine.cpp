// Engine-tier single-rank characterization of inqkit::screens::PlaneScreen::
// extract: the z=0 plane of a He atom at the origin. Asserts the slice has the
// right shape, is a non-negative density, and is non-trivial (peak > 0).
// (The multi-rank E01 Allreduce bug is covered separately by the [!shouldfail]
// cross-rank test launched under mpirun.)

#include <catch2/catch_test_macros.hpp>

#include <inq/inq.hpp>
#include <inqkit/screens/plane_screen.hpp>

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
