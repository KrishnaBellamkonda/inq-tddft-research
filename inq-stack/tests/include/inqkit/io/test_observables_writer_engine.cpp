// Engine-tier (links INQ for StepContext's inq::vector3 + systems types, but
// runs NO GPU/SCF) characterization of inqkit::io::ObservablesWriter: the CSV
// header reflects exactly the enabled columns, and append() writes the matching
// StepContext fields in order.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/step_context.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using Catch::Approx;

namespace {
std::vector<std::string> split(std::string const &s, char sep) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, sep)) out.push_back(tok);
  return out;
}
} // namespace

TEST_CASE("ObservablesWriter: header + row reflect the selection", "[io][observables_writer][engine]") {
  // NOTE: several fields default to TRUE (step, time_au, energy_total,
  // energy_kinetic, current_x/y/z). Set the exact selection we assert on.
  inqkit::io::ObservableSelection sel;
  sel.step = true;
  sel.time_au = true;
  sel.energy_total = true;
  sel.energy_kinetic = false;   // default true — disable to match expected header
  sel.current_x = sel.current_y = sel.current_z = true;
  sel.cod_x = sel.cod_y = sel.cod_z = true;

  fs::path dir = fs::temp_directory_path() / "inqkit_test_obs";
  fs::remove_all(dir); fs::create_directories(dir);
  fs::path csv = dir / "observables.csv";

  {
    inqkit::io::ObservablesWriter writer(csv.string(), sel);
    writer.write_header();

    inqkit::StepContext ctx;
    ctx.step = 7;
    ctx.time_au = 1.25;
    ctx.energy_total = -3.5;
    ctx.current = inqkit::detail::Vec3{0.1, 0.2, 0.3};   // Vec3 unit (D1)
    ctx.wp_center = inq::vector3<double>{1.0, -2.0, 3.0};
    writer.append(ctx);
    writer.finish();
  }

  std::ifstream in(csv);
  std::string header, row;
  std::getline(in, header);
  std::getline(in, row);

  auto hcols = split(header, ',');
  auto rcols = split(row, ',');
  const std::vector<std::string> expected_header = {
      "step", "time_au", "energy_total",
      "current_x", "current_y", "current_z",
      "cod_x_bohr", "cod_y_bohr", "cod_z_bohr"};
  REQUIRE(hcols == expected_header);
  REQUIRE(rcols.size() == expected_header.size());

  CHECK(std::stoi(rcols[0]) == 7);
  CHECK(std::stod(rcols[1]) == Approx(1.25));
  CHECK(std::stod(rcols[2]) == Approx(-3.5));
  CHECK(std::stod(rcols[3]) == Approx(0.1));
  CHECK(std::stod(rcols[5]) == Approx(0.3));
  CHECK(std::stod(rcols[6]) == Approx(1.0));
  CHECK(std::stod(rcols[8]) == Approx(3.0));

  fs::remove_all(dir);
}
