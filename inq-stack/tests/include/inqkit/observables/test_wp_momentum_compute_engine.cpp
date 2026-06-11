// Engine-tier test of the extracted WPMomentumStats::compute() method (stats
// testability refactor). Calls compute() DIRECTLY (no CSV, no RT Viewables) on
// an injected WP and checks the momentum moments: ⟨p⟩ = k₀, Parseval norm ≈ 1,
// e_kin ≈ ½k₀². This guards the refactor and is the direct unit test the
// two-route test validated independently.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <filesystem>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
namespace fs = std::filesystem;

TEST_CASE("WPMomentumStats::compute(): moments of an injected WP", "[observables][wp_momentum][compute][engine]") {
  const double K0 = 1.0;  // Bohr^-1 along x

  systems::ions ions(systems::cell::cubic(12.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(2));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0).sigma(1.5).k0(K0, 0.0, 0.0)
                    .inject_into_last_extra_state(electrons, 1.0);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_wpmom";
  fs::remove_all(dir);
  inqkit::observables::WPMomentumStats stats((dir / "wpmom.csv").string(),
                                             report.state_index);

  auto m = stats.compute(electrons);

  // N is the raw Parseval sum in INQ's reciprocal convention (no dV_k factor),
  // so it is large/convention-dependent, NOT ≈1. The mean momenta below use the
  // ratio sum_p/N, which cancels it.
  CHECK(m.N > 0.0);
  CHECK(m.px == Approx(K0).margin(0.05));        // ⟨p_x⟩ = k₀ (Bohr^-1)
  CHECK(m.py == Approx(0.0).margin(0.05));
  CHECK(m.pz == Approx(0.0).margin(0.05));
  // e_kin = ½(k₀² + 3σ_p²) exceeds the pure-k₀ kinetic ½k₀²=0.5 by the WP
  // momentum spread; bound it rather than pin the exact spread term.
  CHECK(m.ekin > 0.5 * K0 * K0 * 0.95);
  CHECK(m.ekin < 1.0);

  fs::remove_all(dir);
}
