/*
 * This file handles the storage of the complex wavefunction or the entire
 * electronic system or just an orbital (identified using k_point_index and
 * orbital_index - a sub index that identifies the correct orbital within
 * the set of wavefunctions at a given k point. The API of this module is
 * almost identical to the density.hpp (fields::density) class.
 *
 * It is important to note that the complex numbers are handled using
 * inq::complex instead of std::complex. This is to have consistency with
 * the storage of complex numbers in the INQ module.
 * */


 /* 
 
  TODO: 1.  I want to check if in making of the new observables such as momentum distribution
  and others, if the total complex number of the wavefunction has been used. This is important
  to ensure that sensible observables are being produced.  
 
  2. We currently have it that this orbtial.hpp file imports
  the density.hpp file only for the fft_shift function. 
  Perhaps, we can move fft_shift to a better place? 
 
  */

#pragma once

#include <inq/inq.hpp>

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/complex_field_3d.hpp>
// fft_shift_index now lives in detail/grid_layout.hpp (included above), so this
// file no longer depends on fields/density.hpp (resolves the cross-include).

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

  // Output index ix runs left-to-right (ix = 0 at -L/2, ix = nx-1 near +L/2).
  // INQ's hypercubic is FFT-natural, so we read hc[fft_shift_index(ix), ...].
  // Without this shift the exported wavefunction is spatially scrambled
  // relative to the metadata origin (-L/2). Same convention as
  // density.hpp::total / density.hpp::orbital.
  for (int ix = 0; ix < nx; ++ix) {
    int sx = inqkit::detail::grid_layout::fft_shift_index(ix, nx);
    for (int iy = 0; iy < ny; ++iy) {
      int sy = inqkit::detail::grid_layout::fft_shift_index(iy, ny);
      for (int iz = 0; iz < nz; ++iz) {
        int sz = inqkit::detail::grid_layout::fft_shift_index(iz, nz);
        auto flat =
            inqkit::detail::grid_layout::flatten_index(ix, iy, iz, ny, nz);
        auto psi = hc[sx][sy][sz][local_orbital];

        field.values[flat] =
            std::complex<double>(inq::real(psi), inq::imag(psi));
      }
    }
  }

  return field;
}

} // namespace inqkit::fields::orbital
