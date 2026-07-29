/*
 * Defines a structure containing a complex field in three
 * dimensions, for example, the wavefunction of
 * a single orbital. Stores important information
 * such as -
 * 1. origin of the coordinates
 * 2. number of points in each dimension
 * 3. spacing between points in each dimension
 * 4. the field values stored as an array, stored as vector<complex>
 *
 * Data in 1-3 is required to recreate the coordinates of the grid
 * on which the real field is evaluated.
 *
 * The writer class and the post process (Python) interfaces must depend
 * on this abstraction.
 * */

#pragma once

#include <complex>
#include <cstddef>
#include <vector>

namespace inqkit::fields {

/*
 * The writers (wave function writer) and the post processing scripts
 * visualising these files must be based on this abstraction.
 */
struct ComplexField3D {
  int nx = 0;
  int ny = 0;
  int nz = 0;

  double origin_x_bohr = 0.0;
  double origin_y_bohr = 0.0;
  double origin_z_bohr = 0.0;

  double dx_bohr = 0.0;
  double dy_bohr = 0.0;
  double dz_bohr = 0.0;

  std::vector<std::complex<double>> values;

  std::size_t size() const { return values.size(); }

  bool empty() const { return values.empty(); }
};

} // namespace inqkit::fields
