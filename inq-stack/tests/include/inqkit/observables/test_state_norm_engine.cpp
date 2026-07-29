// Engine-tier test of inqkit::observables::StateNormWriter (norm-per-state
// diagnostic). Every KS orbital is normalised, so ∫|ψ_i|² dV ≈ 1 for the He
// occupied orbital, the extra states, AND the injected WP. Asserts compute()
// returns one entry per state, each norm ≈ 1.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/state_norm_writer.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <filesystem>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
namespace fs = std::filesystem;

TEST_CASE("StateNormWriter: every orbital norm ≈ 1 (incl. injected WP)",
          "[observables][state_norm][engine]") {
  const int n_extra = 2;

  systems::ions ions(systems::cell::cubic(12.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(n_extra));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0).sigma(1.5).k0(0.5, 0.0, 0.0)
                    .inject_into_last_extra_state(electrons, 1.0);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_state_norm";
  fs::remove_all(dir);
  inqkit::observables::StateNormWriter writer((dir / "norms.csv").string());

  auto norms = writer.compute(electrons);

  // He (1 occupied orbital) + n_extra states.
  REQUIRE(static_cast<int>(norms.size()) == 1 + n_extra);

  for (auto const& s : norms) {
    CHECK(s.norm == Approx(1.0).margin(0.02));      // every orbital normalised
  }
  // the injected WP slot is also ≈ 1
  CHECK(norms.at(report.state_index).norm == Approx(1.0).margin(0.05));

  fs::remove_all(dir);
}
