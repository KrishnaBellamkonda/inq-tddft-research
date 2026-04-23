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
