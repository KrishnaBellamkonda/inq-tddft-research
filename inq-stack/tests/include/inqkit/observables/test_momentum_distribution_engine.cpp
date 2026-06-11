// Engine-tier test of inqkit::observables::MomentumDistribution via a MOCK
// Viewables. Injects a WP at k₀ along x and checks the binned WP momentum
// distribution n_wp(|k|) peaks in the bin containing |k₀|. CSV columns:
// step,time_au,k_bohr_inv,n_total,n_wp.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
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

TEST_CASE("MomentumDistribution: WP spectrum peaks at |k0|", "[observables][momentum_distribution][engine]") {
  const double K0 = 1.0;
  const double L = 12.0;

  systems::ions ions(systems::cell::cubic(L * 1.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(2));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0).sigma(1.5).k0(K0, 0.0, 0.0)
                    .inject_into_last_extra_state(electrons, 1.0);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_momdist";
  fs::remove_all(dir);
  fs::path csv = dir / "momdist.csv";
  {
    inqkit::observables::MomentumDistribution md(csv.string(), report.state_index, L);
    md.accumulate(MockView{electrons, 0, 0.0});
  }

  std::ifstream in(csv);
  std::string line;
  double peak_k = -1.0, peak_n = -1.0;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#' || line[0] == 's') continue;  // skip header/comment
    auto c = split(line, ',');
    if (c.size() != 5) continue;
    double k = std::stod(c[2]);
    double n_wp = std::stod(c[4]);
    if (n_wp > peak_n) { peak_n = n_wp; peak_k = k; }
  }

  REQUIRE(peak_n > 0.0);
  CHECK(peak_k == Approx(K0).margin(0.2));   // within ~2 bins of |k0|
}
