// Pure-tier T06 test: VTIImageDataWriter must convert inqkit's x-slowest/
// z-fastest buffer into VTK's x-fastest PointData order WITHOUT transposing axes.
// Each cell value encodes its own (ix,iy,iz); we parse the ascii VTI back and
// assert position j in the data stream holds encode(decode_xfastest(j)).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inqkit/io/vti_image_data_writer.hpp>
#include <inqkit/fields/real_field_3d.hpp>
#include <inqkit/detail/grid_layout.hpp>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using inqkit::fields::RealField3D;
using inqkit::detail::grid_layout::flatten_index;
using Catch::Approx;
namespace fs = std::filesystem;

namespace {
double encode(int ix, int iy, int iz) { return ix * 100.0 + iy * 10.0 + iz; }

std::string slurp(fs::path const &p) {
  std::ifstream in(p);
  std::ostringstream ss; ss << in.rdbuf();
  return ss.str();
}
} // namespace

TEST_CASE("VTIImageDataWriter: x-fastest reorder preserves axis mapping (T06)", "[io][vti][pure]") {
  const int nx = 2, ny = 3, nz = 4;  // distinct dims -> catches transposition
  RealField3D f;
  f.nx = nx; f.ny = ny; f.nz = nz;
  f.origin_x_bohr = -1.0; f.origin_y_bohr = -2.0; f.origin_z_bohr = -3.0;
  f.dx_bohr = 0.5; f.dy_bohr = 0.25; f.dz_bohr = 0.125;
  f.values.resize(nx * ny * nz);
  for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
      for (int iz = 0; iz < nz; ++iz)
        f.values[flatten_index(ix, iy, iz, ny, nz)] = encode(ix, iy, iz);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_vti";
  fs::remove_all(dir); fs::create_directories(dir);
  fs::path vti = dir / "field.vti";

  inqkit::io::VTIImageDataWriter writer({.format = inqkit::io::VTIWriteOptions::Format::ascii});
  writer.write_real(f, vti.string(), "density");

  std::string xml = slurp(vti);
  REQUIRE_FALSE(xml.empty());

  // Extract the ascii DataArray payload.
  auto da = xml.find("format=\"ascii\"");
  REQUIRE(da != std::string::npos);
  auto gt = xml.find('>', da);
  auto end = xml.find("</DataArray>", gt);
  REQUIRE(end != std::string::npos);
  std::string payload = xml.substr(gt + 1, end - gt - 1);

  std::istringstream ps(payload);
  std::vector<double> seq;
  double v;
  while (ps >> v) seq.push_back(v);
  REQUIRE(seq.size() == static_cast<std::size_t>(nx * ny * nz));

  // VTK x-fastest: j = ix + nx*(iy + ny*iz). Decode and check the value.
  for (int iz = 0; iz < nz; ++iz)
    for (int iy = 0; iy < ny; ++iy)
      for (int ix = 0; ix < nx; ++ix) {
        std::size_t j = ix + nx * (iy + ny * iz);
        CHECK(seq[j] == Approx(encode(ix, iy, iz)));
      }

  // Origin + Spacing attributes survive.
  CHECK(xml.find("Origin=\"-1") != std::string::npos);
  CHECK(xml.find("Spacing=\"0.5") != std::string::npos);

  fs::remove_all(dir);
}
