#pragma once

#include <inq/inq.hpp>

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/real_field_3d.hpp>

#include <stdexcept>

// This is a simple if condition that checks 
// if being compiled by CUDA (using the lines __CUDAC__)
// If yes, imports the library and defines a function cudaDeviceSynchronize()
// that blocks all CPU progress until all the GPU work is done. 
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
  // Basis has the x, y and z direction basis vectors
  // and their corresponding sizes. 
  auto const &basis = density.basis();

  
  if (basis.comm().size() != 1) {
    throw std::runtime_error("inqkit::fields::density::total: multi-rank basis "
                             "export is not implemented yet.");
  }

  // Number of grid points in each dimension
  // and the spacing between each of them
  auto const nx = basis.sizes()[0];
  auto const ny = basis.sizes()[1];
  auto const nz = basis.sizes()[2];
  auto const spacing = basis.rspacing();

  // Turns the density object in to 3D vector
  // phi.cubic()[ix][iy][iz]
  // where ix, iy and iz are indices of the 3D
  // grid. 
  auto hc = density.cubic();

  // Using the RealField class to house in all of the
  // information about the field in one object
  RealField3D field;
  field.nx = nx;
  field.ny = ny;
  field.nz = nz;

  // The origin is found out using the undertanding that INQ
  // always has its origin at the (0,0,0) point (as the bottom-most
  // and the left-most point)
  field.origin_x_bohr = basis.symmetric_range_begin(0) * spacing[0];
  field.origin_y_bohr = basis.symmetric_range_begin(1) * spacing[1];
  field.origin_z_bohr = basis.symmetric_range_begin(2) * spacing[2];

  field.dx_bohr = spacing[0];
  field.dy_bohr = spacing[1];
  field.dz_bohr = spacing[2];
  field.values.resize(static_cast<std::size_t>(nx) * ny * nz);

  // Writing the density in the field values using 
  // the flattened index. 
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
 * stored in one orbital slot. Takes in two arguments to locate the right
 * orbital 
 *  1. kpoint_index: 
 *  2. orbital_index
 *
 *  INQ stores its orbitals clubebd together per kpoint. The format 
 *  is electrons.kpin()[k_point_index][orbital_index]
 */
inline RealField3D orbital(inq::systems::electrons const &electrons,
                           int orbital_index, int kpoint_index = 0) {
  INQKIT_GPU_SYNC();

  // Sanity checks 

  //  
  if (kpoint_index < 0 ||
      kpoint_index >= static_cast<int>(electrons.kpin().size())) {
    throw std::runtime_error(
        "inqkit::fields::density::orbital: kpoint_index is out of range.");
  }

  // phi contains all of the orbitals contained at a k point.
  auto const &phi = electrons.kpin()[kpoint_index];

  // To obtain the number of orbitals at a given k point, 
  // we write electrons.kpoin()[k_point_index].spinot_set_size()

  if (orbital_index < 0 || orbital_index >= phi.spinor_set_size()) {
    throw std::runtime_error(
        "inqkit::fields::density::orbital: orbital_index is out of range.");
  }

  // Below are errors that alert the user that multi rank functionality
  // for this density class is not implemented yet. Can be implemented
  // if the necessity arises.
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

  // phi.spinor_set_part() gets the distribution or the parittion object that
  // knows how the spinor orbitals are split across parallel tasks. The command
  // inq::parallel::global_index(orbital_index) gets the index of the said orbital
  // in this parallel sub-process.
  auto const local_orbital = phi.spinor_set_part().global_to_local(
      inq::parallel::global_index(orbital_index));

  auto const &basis = phi.basis();
  auto const nx = basis.sizes()[0];
  auto const ny = basis.sizes()[1];
  auto const nz = basis.sizes()[2];
  auto const spacing = basis.rspacing();

  // As opposed to the phi.cubic() defined above
  // this method turns the density into a 4D object
  // phi.hypercubic()[ix][iy][iz][ist], 
  // with ix, iy and iz having the same definition
  // while ist is the local orbital index at this given
  // point. Meaning, .cubic() gives the entire electronic
  // system's density, while hypercubic gives per orbital
  // density.
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
