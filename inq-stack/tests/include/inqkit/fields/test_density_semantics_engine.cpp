// Engine-tier T02/T03 (WP-inclusion + bath subtraction, finding E02) and
// T29/T31 (orthogonalisation rigor, finding E03). Baseline BL-dens-1: He bath
// (N=2) + one WP injected (occ=1).
//
// KEY EMPIRICAL FINDING (E02, see docs/validation/inqkit-errors.md):
// inq electrons.density() returns a CACHED field (spin_density_) that is NOT
// refreshed by a manual WavePacket injection. So:
//   * right after injection (NO propagation): density::total EXCLUDES the WP
//     (stale post-SCF bath); total_excluding_orbital then double-subtracts.
//   * after real_time::propagate refreshes spin_density_: density::total
//     INCLUDES the WP (full), and total_excluding_orbital gives the bath.
// This test pins BOTH states so the behaviour is locked and any future change
// (e.g. forcing a refresh on injection) is caught.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <cstddef>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {
double integrate(inqkit::fields::RealField3D const &f) {
  const double dV = f.dx_bohr * f.dy_bohr * f.dz_bohr;
  long double s = 0.0L;
  for (double v : f.values) s += v;
  return static_cast<double>(s * dV);
}

systems::electrons make_he_with_wp_slots() {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(4));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());
  return electrons;
}
} // namespace

TEST_CASE("engine T02/E02: density() cache is stale after injection, refreshed by propagation", "[fields][density][wavepacket][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(4));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  // Bath-only (He) integrates to the electron count.
  CHECK(integrate(inqkit::fields::density::total(electrons)) == Approx(2.0).margin(0.05));

  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0)
                    .sigma(1.0)
                    .k0(0.0, 0.0, 0.0)
                    .inject_into_last_extra_state(electrons, 1.0);

  // ---- E02 confirmed: BEFORE propagation, density() is the stale bath. ----
  CHECK(integrate(inqkit::fields::density::total(electrons)) == Approx(2.0).margin(0.05));
  // ... so total_excluding_orbital DOUBLE-subtracts here (documents the trap):
  {
    auto bath_stale = inqkit::fields::density::total_excluding_orbital(
        electrons, report.state_index, 1.0);
    CHECK(integrate(bath_stale) == Approx(1.0).margin(0.1));  // 2 - 1, wrong-by-construction
  }

  // ---- Refresh spin_density_ via a few propagation steps. ----
  real_time::propagate(
      ions, electrons, [](auto const &) {},
      options::theory{}.lda(),
      options::real_time{}.num_steps(4).dt(0.02 * 1.0_atomictime));

  // ---- AFTER propagation: density() now INCLUDES the WP (full = 3). ----
  const double total_after = integrate(inqkit::fields::density::total(electrons));
  CHECK(total_after == Approx(3.0).margin(0.1));   // 2 bath + 1 WP

  // ---- bath subtraction is now correct: full - wp = 2. ----
  auto bath = inqkit::fields::density::total_excluding_orbital(
      electrons, report.state_index, 1.0);
  CHECK(integrate(bath) == Approx(2.0).margin(0.15));
}

TEST_CASE("engine T29/E03: orthogonalisation, single-pass limitation at strong overlap", "[wavepacket][orthogonalisation][engine]") {
  auto electrons = make_he_with_wp_slots();

  // WP centred on the He (sigma 1): deep inside the occupied 1s -> strong overlap.
  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0)
                    .sigma(1.0)
                    .k0(0.0, 0.0, 0.0)
                    .orthogonalise_against_occupied(electrons)
                    .inject_into_last_extra_state(electrons, 1.0);

  // Non-vacuous: the WP really did overlap the occupied subspace.
  CHECK(report.max_overlap > 1.0e-2);
  // Renormalised after projection.
  CHECK(report.norm_after == Approx(1.0).margin(0.03));
  // E03 FIXED: iterated (2-pass) Gram-Schmidt + residual measurement now reaches
  // tolerance even for this strong-overlap WP. (Pre-fix this was CHECK_FALSE,
  // documenting the single-pass limitation.)
  CHECK(report.passed_tolerance);
}
