// Engine-tier Test B (flagship Θ-coord, overall comment #12 + T06).
//
// Places a single He atom at an ASYMMETRIC off-centre position with three
// distinct signed coordinates, then recovers the centre of the density field
// via argmax and checks it lands on the nucleus. Three distinct signed coords
// make an axis transposition (x<->z) or a sign flip detectable: those produce
// errors of order L, far outside the one-cell tolerance.
//
// Run at BOTH even (L=10.0, 20^3) and odd (L=10.5, 21^3) grids at dx=0.5, the
// only difference that exercises the (size+1)/2 half-cell branch physically
// (BL-coord-1a / BL-coord-1b).
//
// The pure tests already pin the exact +0.5 cell-centre convention
// (test_center_of_density) and the fft_shift parity (test_fft_shift); Test B's
// job is the END-TO-END axis mapping through electrons.density() ->
// fields::density::total. Tolerance = one cell (dx) accordingly.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>

#include <cstddef>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {

struct Recovered { double x, y, z; double dx, dy, dz; };

// Build He at (px,py,pz) Bohr in an L-Bohr cubic box at spacing 0.5 Bohr,
// initial-guess the density, and return the argmax cell's physical coordinate.
Recovered recover_he_peak(double L, double px, double py, double pz) {
  systems::ions ions(systems::cell::cubic(L * 1.0_bohr).finite());
  ions.insert("He", {px * 1.0_bohr, py * 1.0_bohr, pz * 1.0_bohr});

  systems::electrons electrons(ions,
                               options::electrons{}.spacing(0.5 * 1.0_bohr));
  ground_state::initial_guess(ions, electrons);

  auto rho = inqkit::fields::density::total(electrons);

  // argmax over the flat field.
  std::size_t best = 0;
  double best_v = rho.values.empty() ? 0.0 : rho.values[0];
  for (std::size_t i = 1; i < rho.values.size(); ++i) {
    if (rho.values[i] > best_v) { best_v = rho.values[i]; best = i; }
  }

  // flat = ((ix*ny)+iy)*nz + iz  -> decode (ix,iy,iz).
  const std::size_t nz = static_cast<std::size_t>(rho.nz);
  const std::size_t ny = static_cast<std::size_t>(rho.ny);
  const int iz = static_cast<int>(best % nz);
  const int iy = static_cast<int>((best / nz) % ny);
  const int ix = static_cast<int>(best / (ny * nz));

  return {rho.origin_x_bohr + (ix + 0.5) * rho.dx_bohr,
          rho.origin_y_bohr + (iy + 0.5) * rho.dy_bohr,
          rho.origin_z_bohr + (iz + 0.5) * rho.dz_bohr,
          rho.dx_bohr, rho.dy_bohr, rho.dz_bohr};
}

} // namespace

TEST_CASE("engine Test B: off-centre He, EVEN grid (L=10.0)", "[fields][density][coord][engine]") {
  const double px = 1.5, py = -2.0, pz = 1.0;
  auto r = recover_he_peak(10.0, px, py, pz);

  // Recovered peak within one cell of the nucleus on every axis (catches
  // transposition/sign-flip which would be off by ~L).
  CHECK(r.x == Approx(px).margin(r.dx));
  CHECK(r.y == Approx(py).margin(r.dy));
  CHECK(r.z == Approx(pz).margin(r.dz));
}

TEST_CASE("engine Test B: off-centre He, ODD grid (L=10.5)", "[fields][density][coord][engine]") {
  const double px = 1.5, py = -2.0, pz = 1.0;
  auto r = recover_he_peak(10.5, px, py, pz);

  CHECK(r.x == Approx(px).margin(r.dx));
  CHECK(r.y == Approx(py).margin(r.dy));
  CHECK(r.z == Approx(pz).margin(r.dz));
}
