// Engine-tier MULTI-RANK test documenting E01 (plane_screen missing Allreduce).
// Launched under `mpirun -np 2`. Each rank computes plane_screen.extract(); the
// returned slice should be IDENTICAL on every rank (parallelisation must not
// change the answer). Under E01 it is NOT — each rank accumulates only its local
// states (or its domain slab), so the per-rank slices differ.
//
// We measure cross-rank disagreement as the variance of the per-rank slice-sum
// signature (computed with PLUS all_reduce only). After the E01 fix (two
// all_reduce calls added to plane_screen::extract) every rank returns the
// identical complete slice, so the variance is ~0 and this is a normal green
// test. (Before the fix it was tagged [!shouldfail] to document the bug.)

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/screens/plane_screen.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <functional>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

TEST_CASE("E01 fixed: plane_screen slice agrees across ranks", "[screens][plane_screen][parallel][engine]") {
  inq::parallel::communicator world{boost::mpi3::environment::get_world_instance()};

  // NOTE: no WP injection — WavePacket injection throws under multi-rank
  // (separate limitation). He alone gives a non-trivial slice, and INQ still
  // decomposes the states/domain across the 2 ranks, which is enough to expose
  // (pre-fix) and verify (post-fix) the plane_screen all_reduce.
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(5));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  inqkit::screens::PlaneScreen screen(0.0, "z0");
  auto slice = screen.extract(electrons);

  double sig = 0.0;
  for (auto const &row : slice)
    for (double v : row) sig += v;

  // Cross-rank agreement via variance of the per-rank signature (PLUS-only).
  double sum_sig = sig;
  world.all_reduce_in_place_n(&sum_sig, 1, std::plus<>{});
  const double mean = sum_sig / world.size();
  double dev = (sig - mean) * (sig - mean);
  world.all_reduce_in_place_n(&dev, 1, std::plus<>{});

  // ~0: all ranks hold the identical full slice (post-E01-fix all_reduce).
  CHECK(dev == Approx(0.0).margin(1e-10));
}
