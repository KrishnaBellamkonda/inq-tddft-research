// Pure-tier smoke test for inqkit::detail::grid_layout.
//
// Purpose: prove the Catch2 + CMake + ctest PURE pipeline compiles, links and
// runs against header-only inqkit code with NO INQ dependency. Exercises the
// already-present pure helpers (flatten_index, step_suffix). The fft_shift
// permutation test (Test A) is added once fft_shift_index is moved here.

#include <catch2/catch_test_macros.hpp>

#include <inqkit/detail/grid_layout.hpp>

using inqkit::detail::grid_layout::flatten_index;
using inqkit::detail::grid_layout::step_suffix;

TEST_CASE("flatten_index: documented x-slowest z-fastest ordering", "[grid_layout][pure]") {
  // flat = ((ix * ny) + iy) * nz + iz, with ny=4, nz=5.
  const int ny = 4, nz = 5;

  // Origin.
  CHECK(flatten_index(0, 0, 0, ny, nz) == 0u);
  // z is fastest: +1 in iz -> +1 in flat.
  CHECK(flatten_index(0, 0, 1, ny, nz) == 1u);
  // y is middle: +1 in iy -> +nz in flat.
  CHECK(flatten_index(0, 1, 0, ny, nz) == 5u);
  // x is slowest: +1 in ix -> +ny*nz in flat.
  CHECK(flatten_index(1, 0, 0, ny, nz) == 20u);
  // Mixed: ((2*4)+3)*5 + 4 = 59.
  CHECK(flatten_index(2, 3, 4, ny, nz) == 59u);
}

TEST_CASE("flatten_index: every cell maps to a unique contiguous index", "[grid_layout][pure]") {
  const int nx = 3, ny = 4, nz = 5;
  std::size_t expected = 0;
  for (int ix = 0; ix < nx; ++ix)
    for (int iy = 0; iy < ny; ++iy)
      for (int iz = 0; iz < nz; ++iz)
        // Iterating in x-slowest..z-fastest order must yield 0,1,2,... densely.
        CHECK(flatten_index(ix, iy, iz, ny, nz) == expected++);
  CHECK(expected == static_cast<std::size_t>(nx) * ny * nz);
}

TEST_CASE("step_suffix: zero-padded 6-digit step tag", "[grid_layout][pure]") {
  CHECK(step_suffix(0) == "_t000000");
  CHECK(step_suffix(100) == "_t000100");
  CHECK(step_suffix(123456) == "_t123456");
}
