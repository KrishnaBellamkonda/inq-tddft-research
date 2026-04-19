#pragma once

#include <inq/inq.hpp>

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/complex_field_3d.hpp>

#include <complex>
#include <stdexcept>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#ifndef INQKIT_GPU_SYNC
#define INQKIT_GPU_SYNC() cudaDeviceSynchronize()
#endif
#else
#ifndef INQKIT_GPU_SYNC
#define INQKIT_GPU_SYNC() ((void)0)
#endif
#endif

namespace inqkit::fields::orbital {

/*
 * Build the complex wavefunction field of one selected orbital:
 *
 *   psi(r)  for a chosen orbital index and k-point index
 *
 * This is the complex-field analogue of density::orbital(...).
 */
inline ComplexField3D wavefunction(inq::systems::electrons const &electrons,
                                   int orbital_index, int kpoint_index = 0) {
  INQKIT_GPU_SYNC();

  if (kpoint_index < 0 ||
      kpoint_index >= static_cast<int>(electrons.kpin().size())) {
    throw std::runtime_error(
        "inqkit::fields::orbital::wavefunction: kpoint_index is out of range.");
  }

  auto const &phi = electrons.kpin()[kpoint_index];

  if (orbital_index < 0 || orbital_index >= phi.spinor_set_size()) {
    throw std::runtime_error("inqkit::fields::orbital::wavefunction: "
                             "orbital_index is out of range.");
  }

  if (phi.basis().comm().size() != 1) {
    throw std::runtime_error("inqkit::fields::orbital::wavefunction: "
                             "multi-rank basis export is not implemented yet.");
  }

  if (phi.set_comm().size() != 1) {
    throw std::runtime_error("inqkit::fields::orbital::wavefunction: "
                             "multi-rank state export is not implemented yet.");
  }

  if (phi.spinor_dim() != 1) {
    throw std::runtime_error("inqkit::fields::orbital::wavefunction: "
                             "spinor_dim != 1 is not implemented yet.");
  }

  auto const global_orbital = inq::parallel::global_index(orbital_index);

  if (!phi.spinor_set_part().contains(global_orbital)) {
    throw std::runtime_error("inqkit::fields::orbital::wavefunction: requested "
                             "orbital is not present on this rank.");
  }

  auto const local_orbital =
      phi.spinor_set_part().global_to_local(global_orbital);

  auto const &basis = phi.basis();
  auto const nx = basis.sizes()[0];
  auto const ny = basis.sizes()[1];
  auto const nz = basis.sizes()[2];
  auto const spacing = basis.rspacing();

  auto hc = phi.hypercubic();

  ComplexField3D field;
  field.nx = nx;
  field.ny = ny;
  field.nz = nz;

  field.origin_x_bohr = basis.symmetric_range_begin(0) * spacing[0];
  field.origin_y_bohr = basis.symmetric_range_begin(1) * spacing[1];
  field.origin_z_bohr = basis.symmetric_range_begin(2) * spacing[2];

  field.dx_bohr = spacing[0];
  field.dy_bohr = spacing[1];
  field.dz_bohr = spacing[2];
  field.values.resize(static_cast<std::size_t>(nx) * ny * nz);

  for (int ix = 0; ix < nx; ++ix) {
    for (int iy = 0; iy < ny; ++iy) {
      for (int iz = 0; iz < nz; ++iz) {
        auto flat =
            inqkit::detail::grid_layout::flatten_index(ix, iy, iz, ny, nz);
        auto psi = hc[ix][iy][iz][local_orbital];

        field.values[flat] =
            std::complex<double>(inq::real(psi), inq::imag(psi));
      }
    }
  }

  return field;
}

} // namespace inqkit::fields::orbital
