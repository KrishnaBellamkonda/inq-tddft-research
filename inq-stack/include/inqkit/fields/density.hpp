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

// FFT-shift along one axis: convert a contiguous output index
// (0 -> -L/2, ..., size-1 -> +L/2 - dx) into the FFT-natural array
// index that INQ uses internally (0 -> cell centre, size/2 -> -L/2).
// See inq/src/basis/grid.hpp:78-95 (to/from_symmetric_range) and the
// run_06_centred_writer_check diagnostic. The formula
// (idx + (size+1)/2) % size handles both even and odd sizes.
inline int fft_shift_index(int output_idx, int size) {
  return (output_idx + (size + 1) / 2) % size;
}

/*
 * Build the total electronic density field:
 *   rho(r) = electrons.density()
 *
 * INQ stores the real-space grid in FFT-natural order (array index 0
 * corresponds to physical position 0, the cell centre). The metadata
 * we publish has Origin = -L/2, so when iterating output index ix from
 * 0 to nx-1 we must read INQ array element at fft_shift_index(ix, nx)
 * to get a contiguous left-to-right physical layout.
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

  // INQ's density.cubic() returns a 3D view of the GPU-resident density.
  // The original implementation iterated hc[sx][sy][sz] in a host loop,
  // which triggers ONE GPU->host element fetch per cell — at 4.7M cells
  // that was ~30 minutes per call. Bulk-copy the entire field to a host
  // multi::array first, then do the FFT-shift loop on host memory.
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

  // Single GPU->host bulk copy via boost::multi's array copy constructor.
  // host_hc is a host-allocated 3D array of the same shape.
  boost::multi::array<double, 3> host_hc{density.cubic()};

  // Output index ix runs left-to-right (ix=0 at -L/2, ix=nx-1 near +L/2).
  // INQ's hc is FFT-natural, so we read host_hc[fft_shift_index(ix)].
  for (int ix = 0; ix < nx; ++ix) {
    int sx = fft_shift_index(ix, nx);
    for (int iy = 0; iy < ny; ++iy) {
      int sy = fft_shift_index(iy, ny);
      for (int iz = 0; iz < nz; ++iz) {
        int sz = fft_shift_index(iz, nz);
        auto flat =
            inqkit::detail::grid_layout::flatten_index(ix, iy, iz, ny, nz);
        field.values[flat] = host_hc[sx][sy][sz];
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

  // INQ's phi.hypercubic() is a 4D GPU array [nx][ny][nz][nstates_local].
  // We only want one orbital, so slice first then bulk-copy that 3D slice
  // to host. Slicing a multi::array returns a non-owning view into GPU
  // memory; the copy constructor pulls it to host in one cudaMemcpy.
  // (Original implementation did per-element hc[sx][sy][sz][local_orbital]
  //  inside a host loop, costing ~30 min per 4.7M-cell call.)
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

  // NOTE: this is the ORIGINAL per-element loop (slow: each access is
  // a synchronous GPU->host fetch, ~30 min per 4.7M-cell call). Kept
  // for one-shot initialisation calls (e.g. writing density_wp_initial
  // before propagation starts). Per-step propagation callbacks should
  // AVOID calling this; instead use density::total which is bulk-copied.
  auto hc = phi.hypercubic();
  for (int ix = 0; ix < nx; ++ix) {
    int sx = fft_shift_index(ix, nx);
    for (int iy = 0; iy < ny; ++iy) {
      int sy = fft_shift_index(iy, ny);
      for (int iz = 0; iz < nz; ++iz) {
        int sz = fft_shift_index(iz, nz);
        auto flat =
            inqkit::detail::grid_layout::flatten_index(ix, iy, iz, ny, nz);
        auto psi = hc[sx][sy][sz][local_orbital];
        auto re = inq::real(psi);
        auto im = inq::imag(psi);
        field.values[flat] = re * re + im * im;
      }
    }
  }

  return field;
}

} // namespace inqkit::fields::density
