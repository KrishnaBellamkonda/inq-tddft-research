// Engine-tier orbital test (#5): confirms fields/orbital.hpp compiles + runs
// after the T05 change (it no longer includes fields/density.hpp, getting
// fft_shift from detail/grid_layout.hpp), and that orbital::wavefunction
// extracts the COMPLEX KS orbital with the correct shape and normalisation.
//
// A KS orbital is normalised: integral |psi|^2 dV = 1.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/fields/orbital.hpp>

#include <complex>
#include <cstddef>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

TEST_CASE("engine: orbital::wavefunction extracts a normalised complex field", "[fields][orbital][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});

  systems::electrons electrons(ions, options::electrons{}.spacing(0.5 * 1.0_bohr));
  ground_state::initial_guess(ions, electrons);

  auto psi = inqkit::fields::orbital::wavefunction(electrons, /*orbital_index=*/0);

  CHECK(psi.nx > 0);
  CHECK(psi.values.size() ==
        static_cast<std::size_t>(psi.nx) * psi.ny * psi.nz);

  // integral |psi|^2 dV = 1 for a normalised KS orbital.
  const double dV = psi.dx_bohr * psi.dy_bohr * psi.dz_bohr;
  long double norm2 = 0.0L;
  for (auto const &c : psi.values) norm2 += std::norm(c);  // |c|^2 = re^2+im^2
  CHECK(static_cast<double>(norm2 * dV) == Approx(1.0).margin(0.02));
}
