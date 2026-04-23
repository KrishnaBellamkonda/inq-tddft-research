/*
 *
 *  Limitations:
 *  1. Currently only work along the z axis. TODO: Generalise the class to all dims.
 *  2. TODO: Can implement time averaged total density. 
 *
 * */

#pragma once

// inqkit::screens::PlaneScreen
//
// Extracts the 2D electron density on a constant-z plane:
//   ρ(x, y, z_screen) = ∑_i occ_i |ψ_i(x, y, z_screen)|²
//
// Implementation ported from the proven extract_density_slice /
// save_density_slice functions in
// ResearchProject/systems/coronene/04_leed_simulation/utils.hpp. GPU_SYNC() is
// called at the start of extract() because INQ orbital data lives in CUDA
// managed memory prefetched to the device; CPU reads without a prior sync
// return stale values.

#include <cassert>
#include <cmath>
#include <fstream>
#include <inq/inq.hpp>
#include <iomanip>
#include <string>
#include <vector>

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

namespace inqkit {
namespace screens {

class PlaneScreen {
  double z_bohr_;
  std::string label_;

  // ── Nearest-z-index helper ──────────────────────────────────────────────
  static int iz_nearest(inq::systems::electrons const &electrons,
                        double z_target) {
    auto const &basis = electrons.states_basis();
    int Nz = basis.sizes()[2];
    double dz = basis.rspacing()[2];
    int iz = static_cast<int>(std::round(z_target / dz));
    return std::max(0, std::min(iz, Nz - 1));
  }

public:
  PlaneScreen() = default;

  PlaneScreen(double z_bohr, std::string label = "")
      : z_bohr_(z_bohr), label_(std::move(label)) {}

  double z_bohr() const { return z_bohr_; }
  std::string const &label() const { return label_; }

  // ── Grid dimensions for this cell ──────────────────────────────────────
  int nx(inq::systems::electrons const &electrons) const {
    return electrons.states_basis().sizes()[0];
  }
  int ny(inq::systems::electrons const &electrons) const {
    return electrons.states_basis().sizes()[1];
  }
  double dx(inq::systems::electrons const &electrons) const {
    return electrons.states_basis().rspacing()[0];
  }
  double dy(inq::systems::electrons const &electrons) const {
    return electrons.states_basis().rspacing()[1];
  }

  // ── Extract 2D density slice ────────────────────────────────────────────
  // Returns [Ny_g][Nx_g] array: ∑_i occ_i |ψ_i(x,y,z_screen)|²
  //
  // Occupation weights are read directly from electrons.occupations()[0].
  // WavePacket::inject_into_last_extra_state() sets occ[WP] = 1.0, so the
  // WP is included automatically via the normal occupation path. States with
  // occ == 0.0 (unused extra states) are skipped. This handles fractional
  // occupations (metallic smearing, open-shell) without special-casing.
  //
  // Note: electrons.density() (INQ's built-in cached field) does not include
  // the WP density even after injection. The orbital wavefunction data in
  // kpin()[0] is authoritative; reading occupations directly from
  // electrons.occupations()[0] is the correct approach here.
  //
  // TODO (MPI): no Allreduce on slice after the loop. In a multi-rank run
  // with state-dimension parallelism each process only accumulates its local
  // states. Add an Allreduce before returning if multi-rank support is needed.
  std::vector<std::vector<double>>
  extract(inq::systems::electrons const &electrons) const {
    INQKIT_GPU_SYNC(); // flush device-resident orbital data before CPU reads

    auto const &basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int iz_tgt = iz_nearest(electrons, z_bohr_);

    auto const &phi = electrons.kpin()[0];
    int Nx = phi.basis().cubic_part(0).local_size();
    int Ny = phi.basis().cubic_part(1).local_size();
    int Nz = phi.basis().cubic_part(2).local_size();

    std::vector<std::vector<double>> slice(Ny_g,
                                           std::vector<double>(Nx_g, 0.0));
    auto const &occ = electrons.occupations()[0];
    auto hc = phi.hypercubic();

    for (int ist = 0; ist < phi.set_part().local_size(); ist++) {
      double f = occ[ist];
      if (f == 0.0) continue;  // skip unused extra states
      for (int iz = 0; iz < Nz; iz++) {
        auto iz_g = phi.basis().cubic_part(2).local_to_global({iz}).value();
        if (iz_g != iz_tgt)
          continue;
        for (int ix = 0; ix < Nx; ix++) {
          auto ix_g = phi.basis().cubic_part(0).local_to_global({ix}).value();
          for (int iy = 0; iy < Ny; iy++) {
            auto iy_g = phi.basis().cubic_part(1).local_to_global({iy}).value();
            auto w = hc[ix][iy][iz][ist];
            slice[iy_g][ix_g] +=
                f * (w.real() * w.real() + w.imag() * w.imag());
          }
        }
      }
    }
    return slice;
  }

  // ── Write slice to file ─────────────────────────────────────────────────
  // Header line: "# label=LABEL z=Z_BOHR t=T_AU"
  // Data: space-separated rows (scientific notation).
  void save(std::vector<std::vector<double>> const &slice, double time_au,
            std::string const &filename) const {
    std::ofstream f(filename);
    f << "# label=" << label_ << " z=" << std::fixed << std::setprecision(6)
      << z_bohr_ << " t=" << time_au << "\n";
    for (auto const &row : slice) {
      for (std::size_t ix = 0; ix < row.size(); ix++) {
        f << std::scientific << std::setprecision(6) << row[ix];
        if (ix + 1 < row.size())
          f << " ";
      }
      f << "\n";
    }
  }
};

} // namespace screens
} // namespace inqkit
