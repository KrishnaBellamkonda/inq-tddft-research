// Engine-tier test of inqkit::screens::LeedPatternAccumulator: pattern(x,y) +=
// slice(x,y)·dt each accumulate() call. With no propagation between calls the
// slice is constant, so after two dt=0.5 accumulations the pattern equals the
// single plane_screen slice (Σdt = 1.0).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/screens/plane_screen.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

TEST_CASE("LeedPatternAccumulator: pattern = Σ slice·dt", "[screens][leed][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(ions, options::electrons{}.spacing(0.5 * 1.0_bohr));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  inqkit::screens::PlaneScreen screen(0.0, "z0");
  auto slice = screen.extract(electrons);

  inqkit::screens::LeedPatternAccumulator leed{screen};
  leed.accumulate(electrons, 0.5);
  leed.accumulate(electrons, 0.5);   // Σdt = 1.0, constant slice
  auto const &pat = leed.pattern();

  REQUIRE(pat.size() == slice.size());
  REQUIRE(pat[0].size() == slice[0].size());
  for (std::size_t iy = 0; iy < pat.size(); ++iy)
    for (std::size_t ix = 0; ix < pat[iy].size(); ++ix)
      CHECK(pat[iy][ix] == Approx(slice[iy][ix]).margin(1e-12));
}
