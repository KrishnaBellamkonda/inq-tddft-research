// Engine-tier characterization of inqkit::observables::dump_eigenvalues.
// He GS (gamma-only) → dump CSV → assert one row per state, ascending
// eigenvalues, He ground occupation (2,0,0,…), Ha→eV unit consistency, gamma kpt.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/eigenvalue_dump.hpp>

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
  std::vector<std::string> out; std::stringstream ss(s); std::string t;
  while (std::getline(ss, t, sep)) out.push_back(t);
  return out;
}
} // namespace

TEST_CASE("dump_eigenvalues: He GS eigenvalue table", "[observables][eigenvalue_dump][engine]") {
  systems::ions ions(systems::cell::cubic(10.0_bohr).finite());
  ions.insert("He", {0.0_bohr, 0.0_bohr, 0.0_bohr});
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.5 * 1.0_bohr).extra_states(3));
  ground_state::initial_guess(ions, electrons);
  ground_state::calculate(ions, electrons, options::theory{}.lda());

  fs::path dir = fs::temp_directory_path() / "inqkit_test_eig";
  fs::remove_all(dir);
  fs::path csv = dir / "eigenvalues.csv";
  inqkit::observables::dump_eigenvalues(electrons, csv.string());

  std::ifstream in(csv);
  std::string header; std::getline(in, header);
  auto hcols = split(header, ',');
  REQUIRE(hcols.size() == 9);
  CHECK(hcols[6] == "eigenvalue_ha");
  CHECK(hcols[8] == "occupation");

  std::vector<double> eig_ha, eig_ev, occ;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto c = split(line, ',');
    REQUIRE(c.size() == 9);
    // gamma point.
    CHECK(std::stod(c[1]) == Approx(0.0));
    eig_ha.push_back(std::stod(c[6]));
    eig_ev.push_back(std::stod(c[7]));
    occ.push_back(std::stod(c[8]));
  }

  REQUIRE(eig_ha.size() >= 4);                 // 1 occupied + >=3 extra
  // Eigenvalues ascending (INQ sorts by energy).
  for (std::size_t i = 1; i < eig_ha.size(); ++i)
    CHECK(eig_ha[i] >= eig_ha[i - 1] - 1e-9);
  // He closed shell: lowest state doubly occupied, rest empty.
  CHECK(occ[0] == Approx(2.0).margin(0.05));
  CHECK(occ[1] == Approx(0.0).margin(0.05));
  // Ha -> eV consistency (27.2114).
  for (std::size_t i = 0; i < eig_ha.size(); ++i)
    CHECK(eig_ev[i] == Approx(eig_ha[i] * 27.211386).epsilon(1e-4));

  fs::remove_all(dir);
}
