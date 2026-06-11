// Engine-tier t=0 identity test for inqkit::observables::OrbitalOverlapMatrix.
// At t=0 the evolved orbitals equal the GS reference orbitals, so the n_ref×n_ref
// block of O_ij = |<ψ_i^GS|ψ_j>|² must be the identity (KS orbitals orthonormal).
// Column n_ref is the WP overlap |<ψ_i|ψ_wp>|² (valid in [0,1]).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
namespace fs = std::filesystem;

namespace {
std::vector<std::string> split(std::string const &s, char sep) {
  std::vector<std::string> o; std::stringstream ss(s); std::string x;
  while (std::getline(ss, x, sep)) o.push_back(x);
  return o;
}
} // namespace

TEST_CASE("OrbitalOverlapMatrix: t=0 identity block + WP column", "[observables][orbital_overlap][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(3));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0).sigma(1.0).k0(0.5, 0.0, 0.0)
                    .orthogonalise_against_occupied(electrons)
                    .inject_into_last_extra_state(electrons, 1.0);
  const int n_ref = report.state_index;  // states 0..n_ref-1 are GS refs

  fs::path dir = fs::temp_directory_path() / "inqkit_test_overlap";
  fs::remove_all(dir);
  {
    inqkit::observables::OrbitalOverlapMatrix obs(electrons, n_ref, dir.string());
    obs.snapshot(electrons, /*time_au=*/0.0, /*step=*/0);
  }

  std::ifstream in(dir / "overlap_000000.csv");
  REQUIRE(in.good());
  std::string line;
  std::getline(in, line);                 // skip "# step=..." comment
  REQUIRE(line[0] == '#');

  std::vector<std::vector<double>> O;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::vector<double> row;
    for (auto const &t : split(line, ',')) row.push_back(std::stod(t));
    O.push_back(row);
  }

  REQUIRE(static_cast<int>(O.size()) == n_ref);
  REQUIRE(static_cast<int>(O[0].size()) == n_ref + 1);  // n_evolved = n_ref+1

  for (int i = 0; i < n_ref; ++i) {
    CHECK(O[i][i] == Approx(1.0).margin(0.02));          // diagonal identity
    for (int j = 0; j < n_ref; ++j)
      if (i != j) CHECK(O[i][j] == Approx(0.0).margin(0.02));  // off-diagonal
    CHECK(O[i][n_ref] >= 0.0);                            // WP column valid
    CHECK(O[i][n_ref] <= 1.0 + 1e-6);
  }

  fs::remove_all(dir);
}
