// Pure-tier round-trip characterization of inqkit::io::RealField3DWriter.
// Writes a known RealField3D and reads back the .raw (flat float64 binary) and
// .meta.txt sidecar, asserting the round-trip preserves values + grid metadata.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/fields/real_field_3d.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using inqkit::fields::RealField3D;
using Catch::Approx;
namespace fs = std::filesystem;

namespace {
std::map<std::string, std::string> parse_meta(fs::path const &p) {
  std::map<std::string, std::string> m;
  std::ifstream in(p);
  std::string line;
  while (std::getline(in, line)) {
    auto eq = line.find('=');
    if (eq == std::string::npos) continue;
    std::string key = line.substr(0, eq);
    std::string val = line.substr(eq + 1);
    auto trim = [](std::string &s) {
      while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.erase(s.begin());
      while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) s.pop_back();
    };
    trim(key); trim(val);
    m[key] = val;
  }
  return m;
}
} // namespace

TEST_CASE("RealField3DWriter: .raw + .meta round-trip preserves field", "[io][real_field_writer][pure]") {
  RealField3D f;
  f.nx = 2; f.ny = 3; f.nz = 4;
  f.origin_x_bohr = -1.0; f.origin_y_bohr = -1.5; f.origin_z_bohr = -2.0;
  f.dx_bohr = 0.5; f.dy_bohr = 0.25; f.dz_bohr = 0.125;
  f.values.resize(2 * 3 * 4);
  for (std::size_t i = 0; i < f.values.size(); ++i)
    f.values[i] = 0.1 * static_cast<double>(i) - 1.0;  // distinct values

  fs::path dir = fs::temp_directory_path() / "inqkit_test_realfield";
  fs::remove_all(dir);

  inqkit::io::RealField3DWriter writer(dir.string(), {.field_name = "rho"});
  writer.write(f, "frame0");

  // ---- .raw : flat float64, ((ix*ny)+iy)*nz+iz order ----
  fs::path raw = dir / "frame0.raw";
  REQUIRE(fs::exists(raw));
  REQUIRE(fs::file_size(raw) == f.values.size() * sizeof(double));
  std::vector<double> back(f.values.size());
  {
    std::ifstream in(raw, std::ios::binary);
    in.read(reinterpret_cast<char *>(back.data()),
            static_cast<std::streamsize>(back.size() * sizeof(double)));
  }
  for (std::size_t i = 0; i < f.values.size(); ++i)
    CHECK(back[i] == Approx(f.values[i]));

  // ---- .meta.txt : grid metadata ----
  fs::path meta = dir / "frame0.meta.txt";
  REQUIRE(fs::exists(meta));
  auto m = parse_meta(meta);
  CHECK(m["type"] == "real_field_3d");
  CHECK(m["dtype"] == "float64");
  CHECK(m["field_name"] == "rho");
  CHECK(m["nx"] == "2");
  CHECK(m["ny"] == "3");
  CHECK(m["nz"] == "4");
  CHECK(m["layout"] == "x_slowest_z_fastest");
  CHECK(m["value_file"] == "frame0.raw");
  // origin / spacing parse to the right numbers.
  {
    std::istringstream os(m["origin_bohr"]);
    double ox, oy, oz; os >> ox >> oy >> oz;
    CHECK(ox == Approx(-1.0)); CHECK(oy == Approx(-1.5)); CHECK(oz == Approx(-2.0));
    std::istringstream ss(m["spacing_bohr"]);
    double dx, dy, dz; ss >> dx >> dy >> dz;
    CHECK(dx == Approx(0.5)); CHECK(dy == Approx(0.25)); CHECK(dz == Approx(0.125));
  }

  fs::remove_all(dir);
}
