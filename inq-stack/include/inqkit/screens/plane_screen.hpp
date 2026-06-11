/*
 * inqkit::screens::PlaneScreen — 2D electron-density slice on a plane normal to
 * a chosen axis (D2: generalised from z-only to x/y/z; default axis=2=z is
 * byte-identical to the previous behaviour). A time-averaged accumulator
 * (TimeAveragedScreen) is provided for ⟨ρ⟩ = Σ_t ρ·dt / T.
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
#include <functional>
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
  double z_bohr_;        // position (Bohr) along axis_ of the screen plane
  std::string label_;
  int axis_ = 2;         // plane NORMAL axis: 0=x, 1=y, 2=z (default z)

  // ── Nearest-index helper (along `axis`) ─────────────────────────────────
  // INQ stores fields in FFT-natural order: array index 0 is physical 0 (cell
  // centre); indices in (N/2, N-1] map to negative physical coords. The wrap
  // `((i % N) + N) % N` is the inverse of grid::to_symmetric_range (fixes the
  // old clamp that pinned every negative-coord screen to index 0).
  static int index_nearest(inq::systems::electrons const &electrons,
                           double target, int axis) {
    auto const &basis = electrons.states_basis();
    int N = basis.sizes()[axis];
    double d = basis.rspacing()[axis];
    int i = static_cast<int>(std::round(target / d));
    return ((i % N) + N) % N;               // FFT-natural wrap into [0, N)
  }

public:
  PlaneScreen() = default;

  PlaneScreen(double z_bohr, std::string label = "", int axis = 2)
      : z_bohr_(z_bohr), label_(std::move(label)), axis_(axis) {}

  double z_bohr() const { return z_bohr_; }   // position along axis_ (legacy name)
  double position() const { return z_bohr_; }
  int axis() const { return axis_; }
  std::string const &label() const { return label_; }

  // The two in-plane axes (ascending) for a given normal axis: {0,1,2}\{axis}.
  // axis 0 -> (y,z); axis 1 -> (x,z); axis 2 -> (x,y). Explicit lookup so the
  // ascending-order invariant (extract() returns [size(a1)][size(a0)]) is obvious.
  static void inplane_axes(int axis, int &a0, int &a1) {
    static constexpr int IP[3][2] = {{1, 2}, {0, 2}, {0, 1}};
    a0 = IP[axis][0];
    a1 = IP[axis][1];
  }

  // ── In-plane grid dimensions/spacings (axis-aware) ──────────────────────
  // nx()/ny()/dx()/dy() are the in-plane sizes/spacings of THIS screen's plane:
  // (a0, a1) for the normal axis_. For axis_=2 (z) this is (x, y) — unchanged
  // back-compat — and matches extract()'s [ny][nx] = [size(a1)][size(a0)] output.
  int nx(inq::systems::electrons const &electrons) const {
    int a0, a1; inplane_axes(axis_, a0, a1);
    return electrons.states_basis().sizes()[a0];
  }
  int ny(inq::systems::electrons const &electrons) const {
    int a0, a1; inplane_axes(axis_, a0, a1);
    return electrons.states_basis().sizes()[a1];
  }
  double dx(inq::systems::electrons const &electrons) const {
    int a0, a1; inplane_axes(axis_, a0, a1);
    return electrons.states_basis().rspacing()[a0];
  }
  double dy(inq::systems::electrons const &electrons) const {
    int a0, a1; inplane_axes(axis_, a0, a1);
    return electrons.states_basis().rspacing()[a1];
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
  // MPI (E01 fix): the slice is accumulated from LOCAL states (set
  // decomposition) into GLOBAL (x,y) cells (basis/domain decomposition), so in a
  // multi-rank run each rank holds only a partial slice. It is reduced across
  // BOTH the basis communicator and the state communicator before returning,
  // exactly as wp_real_space_stats does for its moment sums. Cross-rank
  // agreement is verified by tests/cpp/engine/test_plane_screen_parallel_engine.
  std::vector<std::vector<double>>
  extract(inq::systems::electrons const &electrons) const {
    INQKIT_GPU_SYNC(); // flush device-resident orbital data before CPU reads

    auto const &basis = electrons.states_basis();
    int a0, a1;                               // in-plane axes (ascending)
    inplane_axes(axis_, a0, a1);
    int S0 = basis.sizes()[a0];               // output plane: [S1][S0]
    int S1 = basis.sizes()[a1];
    int it_tgt = index_nearest(electrons, z_bohr_, axis_);

    auto const &phi = electrons.kpin()[0];
    int La = phi.basis().cubic_part(axis_).local_size();
    int Lb0 = phi.basis().cubic_part(a0).local_size();
    int Lb1 = phi.basis().cubic_part(a1).local_size();

    std::vector<std::vector<double>> slice(S1, std::vector<double>(S0, 0.0));
    auto const &occ = electrons.occupations()[0];
    auto hc = phi.hypercubic();

    // The target plane's local index is STATE-independent: find it once. If this
    // rank does not own the plane, la_match stays -1 and the slice is all-zero
    // here (the cross-rank reduce below fills it). Also precompute the in-plane
    // global indices once instead of per state.
    int la_match = -1;
    for (int la = 0; la < La; ++la)
      if (static_cast<int>(phi.basis().cubic_part(axis_).local_to_global({la}).value())
          == it_tgt) { la_match = la; break; }

    if (la_match >= 0) {
      std::vector<int> g0v(Lb0), g1v(Lb1);
      for (int l0 = 0; l0 < Lb0; ++l0)
        g0v[l0] = phi.basis().cubic_part(a0).local_to_global({l0}).value();
      for (int l1 = 0; l1 < Lb1; ++l1)
        g1v[l1] = phi.basis().cubic_part(a1).local_to_global({l1}).value();

      for (int ist = 0; ist < phi.set_part().local_size(); ist++) {
        double f = occ[ist];
        if (f == 0.0) continue;  // skip unused extra states
        for (int l0 = 0; l0 < Lb0; l0++) {
          for (int l1 = 0; l1 < Lb1; l1++) {
            int loc[3];
            loc[axis_] = la_match; loc[a0] = l0; loc[a1] = l1;  // (ix,iy,iz) local
            auto w = hc[loc[0]][loc[1]][loc[2]][ist];
            slice[g1v[l1]][g0v[l0]] += f * (w.real() * w.real() + w.imag() * w.imag());
          }
        }
      }
    }

    // E01 fix: reduce partial [S1][S0] slices across the basis (domain) and
    // state communicators so every rank returns the complete slice.
    if (phi.basis().comm().size() > 1 || phi.set_comm().size() > 1) {
      std::vector<double> buf(static_cast<std::size_t>(S1) * S0);
      for (int i1 = 0; i1 < S1; ++i1)
        for (int i0 = 0; i0 < S0; ++i0)
          buf[static_cast<std::size_t>(i1) * S0 + i0] = slice[i1][i0];

      if (phi.basis().comm().size() > 1)
        phi.basis().comm().all_reduce_in_place_n(buf.data(), buf.size(),
                                                 std::plus<>{});
      if (phi.set_comm().size() > 1)
        phi.set_comm().all_reduce_in_place_n(buf.data(), buf.size(),
                                             std::plus<>{});

      for (int i1 = 0; i1 < S1; ++i1)
        for (int i0 = 0; i0 < S0; ++i0)
          slice[i1][i0] = buf[static_cast<std::size_t>(i1) * S0 + i0];
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

// Time-averaged density screen: ⟨ρ⟩ = Σ_t ρ(t)·dt / Σ_t dt (D2). Accumulate a
// PlaneScreen::extract() slice per step with its dt; average() returns the
// dt-weighted mean. For constant frames this returns the frame exactly.
class TimeAveragedScreen {
  std::vector<std::vector<double>> accum_;
  double total_dt_ = 0.0;

public:
  void add(std::vector<std::vector<double>> const &slice, double dt) {
    if (accum_.empty())
      accum_.assign(slice.size(),
                    std::vector<double>(slice.empty() ? 0 : slice[0].size(), 0.0));
    for (std::size_t i = 0; i < slice.size(); ++i)
      for (std::size_t j = 0; j < slice[i].size(); ++j)
        accum_[i][j] += slice[i][j] * dt;
    total_dt_ += dt;
  }

  std::vector<std::vector<double>> average() const {
    auto out = accum_;
    if (total_dt_ > 0.0)
      for (auto &row : out)
        for (auto &v : row)
          v /= total_dt_;
    return out;
  }

  double total_time() const { return total_dt_; }
};

} // namespace screens
} // namespace inqkit
