// Pure-tier round-trip characterization of inqkit::io::ComplexField3DWriter.
// Writes a known ComplexField3D and reads back _real.raw + _imag.raw + meta,
// asserting the real/imag split and grid metadata round-trip.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/fields/complex_field_3d.hpp>

#include <complex>
#include <filesystem>
#include <fstream>
#include <vector>

using inqkit::fields::ComplexField3D;
using Catch::Approx;
namespace fs = std::filesystem;

namespace {
std::vector<double> read_doubles(fs::path const &p, std::size_t n) {
  std::vector<double> v(n);
  std::ifstream in(p, std::ios::binary);
  in.read(reinterpret_cast<char *>(v.data()),
          static_cast<std::streamsize>(n * sizeof(double)));
  return v;
}
} // namespace

TEST_CASE("ComplexField3DWriter: _real/_imag raw round-trip", "[io][complex_field_writer][pure]") {
  ComplexField3D f;
  f.nx = 2; f.ny = 2; f.nz = 3;
  f.origin_x_bohr = 0.0; f.origin_y_bohr = 0.0; f.origin_z_bohr = 0.0;
  f.dx_bohr = f.dy_bohr = f.dz_bohr = 1.0;
  f.values.resize(2 * 2 * 3);
  for (std::size_t i = 0; i < f.values.size(); ++i)
    f.values[i] = std::complex<double>(0.5 * i, -0.25 * i);  // distinct re/im

  fs::path dir = fs::temp_directory_path() / "inqkit_test_complexfield";
  fs::remove_all(dir);

  inqkit::io::ComplexField3DWriter writer(dir.string(), {.field_name = "psi"});
  writer.write(f, "wf0");

  fs::path re = dir / "wf0_real.raw";
  fs::path im = dir / "wf0_imag.raw";
  REQUIRE(fs::exists(re));
  REQUIRE(fs::exists(im));
  REQUIRE(fs::file_size(re) == f.values.size() * sizeof(double));
  REQUIRE(fs::file_size(im) == f.values.size() * sizeof(double));

  auto reals = read_doubles(re, f.values.size());
  auto imags = read_doubles(im, f.values.size());
  for (std::size_t i = 0; i < f.values.size(); ++i) {
    CHECK(reals[i] == Approx(f.values[i].real()));
    CHECK(imags[i] == Approx(f.values[i].imag()));
  }

  REQUIRE(fs::exists(dir / "wf0.meta.txt"));
  fs::remove_all(dir);
}
