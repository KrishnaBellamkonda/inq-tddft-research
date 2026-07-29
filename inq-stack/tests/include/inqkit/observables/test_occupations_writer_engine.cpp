// Engine-tier test of inqkit::observables::OccupationsWriter via a MOCK
// Viewables (duck-types electrons()/iter()/time()). Because snapshot() is a
// template, the mock exercises the real observable code with no RT propagation.
// He bath + injected WP → CSV occupations: state 0 = 2 (He), WP slot = 1, rest 0.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
namespace fs = std::filesystem;

namespace {
struct MockView {
  inq::systems::electrons const &e;
  int it;
  double t;
  inq::systems::electrons const &electrons() const { return e; }
  int iter() const { return it; }
  double time() const { return t; }
};
std::vector<std::string> split(std::string const &s, char sep) {
  std::vector<std::string> o; std::stringstream ss(s); std::string x;
  while (std::getline(ss, x, sep)) o.push_back(x);
  return o;
}
} // namespace

TEST_CASE("OccupationsWriter: He+WP occupations via mock Viewables", "[observables][occupations][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(3));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0).sigma(1.0).k0(0.0, 0.0, 0.0)
                    .inject_into_last_extra_state(electrons, 1.0);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_occ";
  fs::remove_all(dir);
  fs::path csv = dir / "occupations.csv";
  {
    inqkit::observables::OccupationsWriter writer(csv.string());
    writer.snapshot(MockView{electrons, 0, 0.0});
  }

  std::ifstream in(csv);
  std::string header; std::getline(in, header);
  CHECK(header == "step,time_au,kpoint_index,state_index,occupation");

  std::map<int, double> occ;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto c = split(line, ',');
    REQUIRE(c.size() == 5);
    occ[std::stoi(c[3])] = std::stod(c[4]);
  }

  CHECK(occ.at(0) == Approx(2.0).margin(0.02));            // He ground orbital
  CHECK(occ.at(report.state_index) == Approx(1.0).margin(0.02)); // WP slot
  // an intermediate (empty) extra state
  if (report.state_index > 1) CHECK(occ.at(1) == Approx(0.0).margin(0.02));

  fs::remove_all(dir);
}
