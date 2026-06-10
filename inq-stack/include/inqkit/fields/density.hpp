/* 
* This file contains the logic to find the total, orbital (specific) and system
* densities from the INQ simulation. Importantly, this file has the fft_shift_index
* function, which converts indices of a flattened array in the usual ordering 
* (for definitions of usual and fft order, check below) to the fft ordering.
* 
*
*/

#pragma once

#include <inq/inq.hpp>

#include <inqkit/detail/grid_layout.hpp>
#include <inqkit/fields/real_field_3d.hpp>

#include <stdexcept>

/*
* This is a simple if condition that checks 
* if being compiled by CUDA (using the lines __CUDAC__)
* If yes, imports the library and defines a function cudaDeviceSynchronize()
* that blocks all CPU progress until all the GPU work is done.
*/  
#ifdef __CUDACC__
#include <cuda_runtime.h>
#define INQKIT_GPU_SYNC() cudaDeviceSynchronize()
#else
#define INQKIT_GPU_SYNC() ((void)0)
#endif

namespace inqkit::fields::density {

// fft_shift_index moved to the shared pure header
// inqkit::detail::grid_layout (detail/grid_layout.hpp), so density.hpp and
// orbital.hpp can share it without one including the other. Used below via
// grid_layout::fft_shift_index. Tested by tests/cpp/test_fft_shift.cpp.

/*
 * Build the total electronic density field:
 *   rho(r) = electrons.density()
 *
 * TODO: IMPORTANT: This function utilises the density
 * function that is shipped by INQ. However, a validation
 * test must be performed here to check if it output the 
 * total or electornic system's density. I define the system as
 * target (which does not include the wavepacket). Need to make
 * a validation check to confirm this.  

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
    int sx = inqkit::detail::grid_layout::fft_shift_index(ix, nx);
    for (int iy = 0; iy < ny; ++iy) {
      int sy = inqkit::detail::grid_layout::fft_shift_index(iy, ny);
      for (int iz = 0; iz < nz; ++iz) {
        int sz = inqkit::detail::grid_layout::fft_shift_index(iz, nz);
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
  // we write electrons.kpoint()[k_point_index].spinot_set_size()

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
  // a synchronous GPU to host fetch, ~30 min per 4.7M-cell call). Kept
  // for one-shot initialisation calls (e.g. writing density_wp_initial
  // before propagation starts). Per-step propagation callbacks should
  // AVOID calling this; instead use density::total which is bulk-copied.
  auto hc = phi.hypercubic();
  for (int ix = 0; ix < nx; ++ix) {
    int sx = inqkit::detail::grid_layout::fft_shift_index(ix, nx);
    for (int iy = 0; iy < ny; ++iy) {
      int sy = inqkit::detail::grid_layout::fft_shift_index(iy, ny);
      for (int iz = 0; iz < nz; ++iz) {
        int sz = inqkit::detail::grid_layout::fft_shift_index(iz, nz);
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

/*
 * System density: the full electronic density MINUS one orbital's contribution.
 *  
 *   The system is referred to as bath in this code. 
 * 
 *   rho_system(r) = rho_total(r) - occupation * |psi_exclude(r)|^2
 *
 * Motivation: when a wave-packet (WP) projectile is injected into an extra
 * Kohn-Sham state at a fixed occupation (typically 1.0), electrons.density()
 * (and hence density::total) ALREADY INCLUDES that WP orbital. For studying the
 * *target system's* response (the jellium system wake, induced density, etc.) the
 * WP orbital must be removed so the saved "system" density integrates to the
 * bath electron count (N_electrons), not N_electrons + 1.
 *
 * Verified (run_wp_n162_L50_E20_sigma1_v2): total() integrates to 163 e
 * (162 bath + 1 WP); total_excluding_orbital(wp_idx, 1.0) integrates to 162 e.
 *
 * `occupation` defaults to 1.0 (the standard WP injection occupation, set in
 * wavepacket.hpp). Pass it explicitly rather than reading
 * electrons.occupations()[k][i] to avoid a synchronous GPU element fetch.
 *
 * Cost: one bulk-copied total() + one per-element orbital() (the orbital() call
 * is the slow per-element loop; acceptable when a WP density is already being
 * written each frame — reuse that field instead if available, see overload).
 */
inline RealField3D total_excluding_orbital(
    inq::systems::electrons const &electrons, int exclude_index,
    double occupation = 1.0, int kpoint_index = 0) {
  RealField3D bath = total(electrons);
  RealField3D orb  = orbital(electrons, exclude_index, kpoint_index);

  // Check if the number of grid points for the bath and the orbital
  // density are the same 
  if (bath.values.size() != orb.values.size() ||
      bath.nx != orb.nx || bath.ny != orb.ny || bath.nz != orb.nz) {
    throw std::runtime_error(
        "inqkit::fields::density::total_excluding_orbital: grid mismatch "
        "between total and orbital densities.");
  }
  for (std::size_t i = 0; i < bath.values.size(); ++i) {
    bath.values[i] -= occupation * orb.values[i];
  }
  return bath;
}

/*
 * Overload: subtract an ALREADY-COMPUTED orbital density from an
 * already-computed total density. Use this in per-step callbacks that already
 * write the WP density that frame, to avoid a second expensive orbital() call.
 *   rho_bath = total_field - occupation * orbital_field
 */

// TODO: Find where this is used. This function essentially uses
// orbital and total densities outside of this function
// and performs the subtraction to get the system density. 
// I feel like this function should not be used. Need to firstly document
// exactly where such as function is being used.  
inline RealField3D total_excluding_orbital(
    RealField3D const &total_field, RealField3D const &orbital_field,
    double occupation = 1.0) {
  if (total_field.values.size() != orbital_field.values.size() ||
      total_field.nx != orbital_field.nx ||
      total_field.ny != orbital_field.ny ||
      total_field.nz != orbital_field.nz) {
    throw std::runtime_error(
        "inqkit::fields::density::total_excluding_orbital: grid mismatch "
        "between total and orbital fields.");
  }
  RealField3D bath = total_field;  // copy metadata + values
  for (std::size_t i = 0; i < bath.values.size(); ++i) {
    bath.values[i] -= occupation * orbital_field.values[i];
  }
  return bath;
}

} // namespace inqkit::fields::density
