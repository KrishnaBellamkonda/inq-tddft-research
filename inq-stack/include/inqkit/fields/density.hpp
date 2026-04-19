#pragma once

#include <inq/inq.hpp>

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/real_field_3d.hpp>

#include <stdexcept>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#define INQKIT_GPU_SYNC() cudaDeviceSynchronize()
#else
#define INQKIT_GPU_SYNC() ((void)0)
#endif

namespace inqkit::fields::density {

/*
 * Build the total electronic density field:
 *   rho(r) = electrons.density()
 */
inline RealField3D total(inq::systems::electrons const &electrons) {
  INQKIT_GPU_SYNC();

  auto density = electrons.density();
  auto const &basis = density.basis();

  if (basis.comm().size() != 1) {
    throw std::runtime_error("inqkit::fields::density::total: multi-rank basis "
                             "export is not implemented yet.");
  }

  auto const nx = basis.sizes()[0];
  auto const ny = basis.sizes()[1];
  auto const nz = basis.sizes()[2];
  auto const spacing = basis.rspacing();

  auto hc = density.cubic();

  RealField3D field;
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
        field.values[flat] = hc[ix][iy][iz];
      }
    }
  }

  return field;
}

/*
 * Build the density of one selected orbital:
 *   rho_orbital(r) = |psi(r)|^2
 *
 * This is the right abstraction for a single KS orbital or a wavepacket
 * stored in one orbital slot.
 */
inline RealField3D orbital(inq::systems::electrons const &electrons,
                           int orbital_index, int kpoint_index = 0) {
  INQKIT_GPU_SYNC();

  if (kpoint_index < 0 ||
      kpoint_index >= static_cast<int>(electrons.kpin().size())) {
    throw std::runtime_error(
        "inqkit::fields::density::orbital: kpoint_index is out of range.");
  }

  auto const &phi = electrons.kpin()[kpoint_index];

  if (orbital_index < 0 || orbital_index >= phi.spinor_set_size()) {
    throw std::runtime_error(
        "inqkit::fields::density::orbital: orbital_index is out of range.");
  }

  if (phi.basis().comm().size() != 1) {
    throw std::runtime_error("inqkit::fields::density::orbital: multi-rank "
                             "basis export is not implemented yet.");
  }

  if (phi.set_comm().size() != 1) {
    throw std::runtime_error("inqkit::fields::density::orbital: multi-rank "
                             "state export is not implemented yet.");
  }

  if (phi.spinor_dim() != 1) {
    throw std::runtime_error("inqkit::fields::density::orbital: spinor_dim != "
                             "1 is not implemented yet.");
  }

  if (!phi.spinor_set_part().contains(orbital_index)) {
    throw std::runtime_error("inqkit::fields::density::orbital: requested "
                             "orbital is not present on this rank.");
  }

  auto const local_orbital = phi.spinor_set_part().global_to_local(
      inq::parallel::global_index(orbital_index));

  auto const &basis = phi.basis();
  auto const nx = basis.sizes()[0];
  auto const ny = basis.sizes()[1];
  auto const nz = basis.sizes()[2];
  auto const spacing = basis.rspacing();

  auto hc = phi.hypercubic();

  RealField3D field;
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
        auto re = inq::real(psi);
        auto im = inq::imag(psi);
        field.values[flat] = re * re + im * im;
      }
    }
  }

  return field;
}

} // namespace inqkit::fields::density
