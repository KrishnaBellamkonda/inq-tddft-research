// Engine-tier test of the extracted WPRealSpaceStats::compute() method (stats
// testability refactor). Calls compute() DIRECTLY on an injected WP and checks
// the real-space moments: ⟨r⟩ = injected centre, N ≈ 1 (normalised), variance
// ≈ σ²/2 (|ψ|² ∝ exp(-r²/σ²)). Uses INQ's rvector_cartesian (node convention),
// so unlike center_of_density it carries no half-cell offset (E04).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <filesystem>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
namespace fs = std::filesystem;

TEST_CASE("WPRealSpaceStats::compute(): ⟨r⟩ = injected centre, N ≈ 1", "[observables][wp_real_space][compute][engine]") {
  const double CX = 1.0, CY = -1.0, CZ = 0.5, SIGMA = 1.5;

  systems::ions ions(systems::cell::cubic(16.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(2));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto report = inqkit::WavePacket{}
                    .center(CX, CY, CZ).sigma(SIGMA).k0(0.0, 0.0, 0.0)
                    .inject_into_last_extra_state(electrons, 1.0);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_wprs";
  fs::remove_all(dir);
  inqkit::observables::WPRealSpaceStats stats((dir / "wprs.csv").string(),
                                              report.state_index);

  auto m = stats.compute(electrons);

  CHECK(m.N == Approx(1.0).margin(0.05));     // normalised WP, ∫|ψ|² dV ≈ 1
  CHECK(m.x == Approx(CX).margin(0.05));
  CHECK(m.y == Approx(CY).margin(0.05));
  CHECK(m.z == Approx(CZ).margin(0.05));
  CHECK(m.sx2 == Approx(SIGMA * SIGMA / 2.0).margin(0.15));  // Var = σ²/2

  fs::remove_all(dir);
}
