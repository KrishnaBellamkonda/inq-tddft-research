#pragma once

#include <cstddef>
#include <vector>

namespace inqkit::fields {

/*
 * Generic owning container for any scalar field sampled on a 3D grid.
 *
 * This is intentionally independent of INQ's internal field storage.
 * Writers and post-processing interfaces should depend on this abstraction.
 */
struct RealField3D {
  int nx = 0;
  int ny = 0;
  int nz = 0;

  double origin_x_bohr = 0.0;
  double origin_y_bohr = 0.0;
  double origin_z_bohr = 0.0;

  double dx_bohr = 0.0;
  double dy_bohr = 0.0;
  double dz_bohr = 0.0;

  std::vector<double> values;

  std::size_t size() const { return values.size(); }

  bool empty() const { return values.empty(); }
};

} // namespace inqkit::fields
