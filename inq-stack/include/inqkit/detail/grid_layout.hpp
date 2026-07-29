/*
 * This file contains the schema for the files that are written
 * in the complex and real field writers. It defines the following -
 *
 * 1. meta file details (type, dtype, layout)
 * 2. suffixes of the files to be written
 *    a. Real field (.meta.txt, .raw)
 *    b. Complex field (.meta.txt, _real.raw, _imag.raw)
 *
 *
 * Changes: Perhaps, the filename must be changed to something like
 * schema, or something along the same lines.
 * */

#pragma once

#include <string>

namespace inqkit::detail::grid_layout {

/*
 * This file defines the shared schema contract for complex 3D raw fields.
 *
 * The writer follows this contract.
 * The Python reader should read the metadata and assume exactly this layout.
 */
struct ComplexField3DRawSchema {
  std::string type = "complex_field_3d";
  std::string dtype = "float64";

  /*
   * Flattening convention used by the writer:
   *
   *   flat = ((ix * ny) + iy) * nz + iz
   *
   * In a 3D grid, this specific flattening convention is such that
   * for a constant value of (ix, iy), all of the iz coordinates are
   * exhausted first. When done with iz values, iy -> iy+1. When all
   * the iy values are exhausted (for a constant ix, upto ny -1), then
   * ix -> ix+1. Can be thought of as dividing the system into
   * columns along z direction for a given ix and iy.
   *
   * Therefore:
   *   x = slowest varying index
   *   z = fastest varying index
   */
  std::string layout = "x_slowest_z_fastest";

  std::string real_suffix = "_real.raw";
  std::string imag_suffix = "_imag.raw";
  std::string meta_suffix = ".meta.txt";
};

struct RealField3DRawSchema {
  std::string type = "real_field_3d";
  std::string dtype = "float64";

  // Flat ordering:
  //   flat = ((ix * ny) + iy) * nz + iz
  // so x is slowest, z is fastest
  std::string layout = "x_slowest_z_fastest";

  std::string value_suffix = ".raw";
  std::string meta_suffix = ".meta.txt";
};


/* Schema initialising functions */
inline ComplexField3DRawSchema complex_field_3d_raw_schema() { return {}; }

inline RealField3DRawSchema real_field_3d_raw_schema() { return {}; }

/*
 * Shared flattening convention for all 3D fields:
 *   flat = ((ix * ny) + iy) * nz + iz
 *
 * x = slowest varying index
 * z = fastest varying index
 */
inline std::size_t flatten_index(int ix, int iy, int iz, int ny, int nz) {
  return ((static_cast<std::size_t>(ix) * ny) + iy) * nz + iz;
}

/*
 * Convert a contiguous ("human") output index into the FFT-natural array index
 * used internally by INQ.
 *
 * For a 1D grid with size = 6 cells spanning [-L/2, +L/2):
 *
 *   position | usual idx | FFT idx
 *   ---------+-----------+--------
 *    -3*dx   |     0     |    3
 *    -2*dx   |     1     |    4
 *    -1*dx   |     2     |    5
 *     0      |     3     |    0   <- origin sits at index 0 in FFT layout
 *    +1*dx   |     4     |    1
 *    +2*dx   |     5     |    2
 *
 * The usual layout sorts cells from most-negative to most-positive. The FFT
 * layout wraps around: it starts at the origin, runs through the positive half,
 * then continues with the negative half.
 *
 * To go from usual to FFT, shift by size/2 with wrap around:
 *   fft_idx = (output_idx + size/2) % size
 *
 * (size+1)/2 instead of size/2 is used to handle odd sizes correctly: integer
 * division would round down and shift the origin off by one cell, while
 * (size+1)/2 rounds up and keeps the origin at FFT index 0. For even sizes both
 * expressions are equal.
 *
 * This is the index-space equivalent of numpy.fft.ifftshift.
 * See also: inq/src/basis/grid.hpp:78-95 (to/from_symmetric_range).
 *
 * Lives here (the pure index-convention header) so both fields/density.hpp and
 * fields/orbital.hpp can share it without one including the other. Verified by
 * the pure test inq-stack/tests/cpp/test_fft_shift.cpp.
 */
inline int fft_shift_index(int output_idx, int size) {
  return (output_idx + (size + 1) / 2) % size;
}

// Returns a step-index suffix like "_t000100" for use in time-series filenames.
inline std::string step_suffix(int step) {
  char buf[16];
  // snprintf() is a function that helps format the
  // string as if it were to be printed but stores it
  // in a variable.
  std::snprintf(buf, sizeof(buf), "_t%06d", step);
  return std::string(buf);
}

} // namespace inqkit::detail::grid_layout
