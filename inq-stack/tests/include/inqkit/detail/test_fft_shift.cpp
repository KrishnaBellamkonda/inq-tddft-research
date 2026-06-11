// Pure-tier Test A: fft_shift_index permutation (gates the T01/T05 move).
//
// fft_shift_index converts a contiguous ("human") output index — running
// left-to-right from -L/2 — into INQ's FFT-natural array index, where physical
// origin sits at array index 0. Oracle: the documented table in the original
// density.hpp docstring, plus parity and permutation invariants.
//
// This test targets the POST-MOVE location detail/grid_layout.hpp. It fails to
// compile against the pre-move tree (symbol lives in fields/density.hpp, which
// is INQ-dependent and cannot be included in a pure test) — that RED state is
// the gate the move must turn GREEN.

#include <catch2/catch_test_macros.hpp>

#include <inqkit/detail/grid_layout.hpp>

#include <vector>

using inqkit::detail::grid_layout::fft_shift_index;

TEST_CASE("fft_shift_index: documented even-size table (size 6)", "[fft_shift][pure]") {
  // position | output idx | FFT idx
  //   -3dx   |     0      |   3
  //   -2dx   |     1      |   4
  //   -1dx   |     2      |   5
  //    0     |     3      |   0   <- origin at FFT index 0
  //   +1dx   |     4      |   1
  //   +2dx   |     5      |   2
  CHECK(fft_shift_index(0, 6) == 3);
  CHECK(fft_shift_index(1, 6) == 4);
  CHECK(fft_shift_index(2, 6) == 5);
  CHECK(fft_shift_index(3, 6) == 0);
  CHECK(fft_shift_index(4, 6) == 1);
  CHECK(fft_shift_index(5, 6) == 2);
}

TEST_CASE("fft_shift_index: odd sizes keep physical origin at FFT index 0", "[fft_shift][pure]") {
  // For odd size the origin output index is size/2 (floor); (size+1)/2 rounding
  // must map it to FFT index 0.
  CHECK(fft_shift_index(2, 5) == 0);  // size 5: positions -2..+2, origin at idx 2
  CHECK(fft_shift_index(3, 7) == 0);  // size 7: origin at idx 3

  // Full size-5 mapping: [3,4,0,1,2].
  CHECK(fft_shift_index(0, 5) == 3);
  CHECK(fft_shift_index(1, 5) == 4);
  CHECK(fft_shift_index(3, 5) == 1);
  CHECK(fft_shift_index(4, 5) == 2);
}

TEST_CASE("fft_shift_index: origin maps to 0 for all parities", "[fft_shift][pure]") {
  for (int size = 1; size <= 16; ++size) {
    // Output index size/2 (floor) is the physical origin; it must land on 0.
    CHECK(fft_shift_index(size / 2, size) == 0);
  }
}

TEST_CASE("fft_shift_index: is a bijection of [0,size) (no collisions)", "[fft_shift][pure]") {
  for (int size = 1; size <= 16; ++size) {
    std::vector<int> seen(size, 0);
    for (int out = 0; out < size; ++out) {
      int s = fft_shift_index(out, size);
      REQUIRE(s >= 0);
      REQUIRE(s < size);
      seen[s]++;
    }
    for (int i = 0; i < size; ++i) {
      CHECK(seen[i] == 1);  // every FFT index hit exactly once
    }
  }
}
