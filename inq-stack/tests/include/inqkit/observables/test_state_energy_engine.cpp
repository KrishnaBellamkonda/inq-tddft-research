// Engine-tier test of inqkit::observables::StateEnergyWriter, driven through a
// REAL short real_time::propagate (its snapshot() needs data.ham(), only
// available inside the propagator). Asserts the CSV has per-state E_expect_ha
// rows that are finite and number n_states per recorded step. CSV columns:
// step,time_au,kpoint_index,state_index,weight,occupation,E_expect_ha.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/state_energy_writer.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace inq;
using namespace inq::magnitude;
namespace fs = std::filesystem;

namespace {
std::vector<std::string> split(std::string const &s, char sep) {
  std::vector<std::string> o; std::stringstream ss(s); std::string x;
  while (std::getline(ss, x, sep)) o.push_back(x);
  return o;
}
} // namespace

TEST_CASE("StateEnergyWriter: per-state energies via short propagation", "[observables][state_energy][engine]") {
  const int n_extra = 2;

  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(n_extra));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  fs::path dir = fs::temp_directory_path() / "inqkit_test_stateE";
  fs::remove_all(dir);
  fs::path csv = dir / "state_energies.csv";
  {
    inqkit::observables::StateEnergyWriter sew(csv.string(), /*emit_variance=*/false);
    real_time::propagate(
        ions, electrons, [&](auto const &data) { sew.snapshot(data); },
        options::theory{}.lda(),
        options::real_time{}.num_steps(2).dt(0.02 * 1.0_atomictime));
  }

  std::ifstream in(csv);
  std::string header; std::getline(in, header);
  auto hcols = split(header, ',');
  REQUIRE(hcols.size() >= 7);
  CHECK(hcols[6] == "E_expect_ha");

  std::map<int, int> rows_per_step;
  int finite_count = 0, total = 0;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto c = split(line, ',');
    if (c.size() < 7) continue;
    rows_per_step[std::stoi(c[0])]++;
    double e = std::stod(c[6]);
    if (std::isfinite(e)) finite_count++;
    total++;
  }

  REQUIRE(total > 0);
  CHECK(finite_count == total);            // all E_expect_ha finite
  // Each recorded step writes one row per state (He occupied + n_extra).
  for (auto const &[step, n] : rows_per_step)
    CHECK(n == 1 + n_extra);
}
