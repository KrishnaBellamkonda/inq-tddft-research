// Engine-tier smoke test: proves the INQ-linking path compiles, links and runs,
// AND that the T01/T05 fft_shift relocation still compiles inside an
// INQ-dependent inqkit header (fields/density.hpp now calls
// grid_layout::fft_shift_index).
//
// Minimal system: one He atom (closed shell, 2 electrons) in a small cubic box.
// We only need an initial guess density, not a converged SCF, to exercise
// density::total. Validates: field dims match the basis, and the density
// integrates to the electron count (a first empirical touch of T02 as well).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

TEST_CASE("engine: density::total matches basis dims and integrates to N", "[fields][density][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});

  systems::electrons electrons(ions, options::electrons{}.cutoff(30.0_Ry));
  ground_state::initial_guess(ions, electrons);

  auto rho = inqkit::fields::density::total(electrons);

  // Field shape is well-formed and consistent.
  CHECK(rho.nx > 0);
  CHECK(rho.ny > 0);
  CHECK(rho.nz > 0);
  CHECK(rho.values.size() ==
        static_cast<std::size_t>(rho.nx) * rho.ny * rho.nz);
  CHECK(rho.dx_bohr > 0.0);

  // Integral of the density = number of electrons (He -> 2).
  const double dV = rho.dx_bohr * rho.dy_bohr * rho.dz_bohr;
  long double sum = 0.0L;
  for (double v : rho.values) sum += v;
  const double integral = static_cast<double>(sum * dV);

  CHECK(integral == Approx(2.0).margin(0.05));
}
