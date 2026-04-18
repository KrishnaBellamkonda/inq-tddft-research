#pragma once
// ============================================================================
// run_003/utils.hpp
//
// Utility functions for the coronene LEED simulation.
//
// Critical fix over run_002:
//   cudaDeviceSynchronize() is called at the top of every function that reads
//   from phi.hypercubic(). In run_002, the TDDFT propagation kernel ran
//   asynchronously on GPU. The CPU callback read stale pre-propagation values
//   from the WP orbital (values were 10^11× too small and showed wrong sign).
//   The sync ensures all GPU writes are complete before any CPU read.
//
// New functions vs run_002:
//   save_orbital_3d        — write single KS orbital as human-readable text
//   save_density_3d        — write total electron density (3D) as text
//   compute_overlap_matrix — S_ij(t) = <φ_i^GS | φ_j(t)> for all i,j
//   save_overlap_matrix    — append overlap matrix block to text file
//   save_grid_coords       — write x/y/z grid coordinates once after SCF
//
// Text formats (all human-readable):
//
//   Orbital file (complex):
//     # ist=I t=T_AU step=N Nx=NX Ny=NY Nz=NZ
//     # Format: one complex value per line, "real imag", C-order (ix slowest)
//     -1.234567e-03  5.678901e-04
//     ...
//
//   Density file (real):
//     # t=T_AU step=N Nx=NX Ny=NY Nz=NZ
//     # Format: one float per line, C-order (ix slowest)
//     1.234567e-03
//     ...
//
//   Overlap matrix block (appended to one file):
//     # step=N t=T_AU n_states=57
//     # S_ij = <phi_i_GS | phi_j(t)>, i=row j=col, pairs "re im"
//     -9.99e-01  0.00e+00   4.56e-07 -7.89e-08 ...
//     ...
//
// Source: Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)
// ============================================================================

#include "config.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cassert>
#include <cerrno>
#include <cstring>
#include <sys/stat.h>

// GPU sync — required before any hypercubic() access from CPU
#ifdef __CUDACC__
#include <cuda_runtime.h>
#define GPU_SYNC() cudaDeviceSynchronize()
#else
#define GPU_SYNC() ((void)0)
#endif

namespace leed_utils {

// ── Helper: zero-pad integer to width W ──────────────────────────────────────
inline std::string pad(int n, int w) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(w) << n;
    return ss.str();
}

// ── Helper: mkdir -p (recursive, like `mkdir -p`) ────────────────────────────
// Creates all missing intermediate directories in path.
// Uses POSIX mkdir; silently ignores EEXIST at each level.
inline void mkdir_p(std::string const & path) {
    std::string tmp = path;
    // Strip trailing slashes
    while(!tmp.empty() && tmp.back() == '/') tmp.pop_back();
    for(std::size_t i = 1; i <= tmp.size(); ++i) {
        if(i == tmp.size() || tmp[i] == '/') {
            std::string sub = tmp.substr(0, i);
            int r = mkdir(sub.c_str(), 0755);
            if(r != 0 && errno != EEXIST) {
                // Not a fatal error for simulation — print warning and continue
                std::fprintf(stderr, "mkdir_p: cannot create '%s': %s\n",
                             sub.c_str(), std::strerror(errno));
            }
        }
    }
}

// ── Nearest grid index for a given z value ────────────────────────────────────
// Maps z_bohr → nearest iz in the global grid [0, Nz-1].
// For the finite cell (corner at origin), grid point iz is at z = iz × dz.
inline int iz_nearest(inq::systems::electrons const & electrons, double z_bohr) {
    auto const & basis = electrons.states_basis();
    int    Nz_g = basis.sizes()[2];
    double dz   = basis.rspacing()[2];
    int    iz   = static_cast<int>(std::round(z_bohr / dz));
    return std::max(0, std::min(iz, Nz_g - 1));
}

// ── WP injection ─────────────────────────────────────────────────────────────
// ψ^WP(r) = (1/(πd²))^{3/4} exp(−|r−b|²/(2d²)) exp(ik·r)
// Written into the LAST KS orbital slot (ist_wp = set_size - 1).
//
// Uses gpu::run (GPU kernel) to write orbital values — required because INQ
// allocates orbital data in CUDA managed memory prefetched to device. CPU loops
// writing to device-resident UVM pages fail silently (values do not persist).
// All INQ orbital writes (randomize, kick, etc.) use gpu::run for this reason.
// After the kernel, GPU_SYNC() ensures the write is complete before any CPU read.
inline void inject_wp(inq::systems::electrons & electrons,
                      double bx, double by, double bz,
                      double kx, double ky, double kz)
{
    using namespace inq;
    using complex = inq::complex;

    auto & phi   = electrons.kpin()[0];
    auto & basis = phi.basis();

    int ist_wp         = phi.set_part().local_size() - 1;
    const double d     = cfg::WP_D_BOHR;
    const double norm_ = cfg::wp_norm();

    // Capture grid spacings and partition offsets as plain scalars — avoids
    // capturing point_operator (which contains parallel::partition with non-trivial
    // GPU-capture semantics that caused rvector_cartesian to return wrong coords).
    // For a single-rank orthogonal cell, r = (ix+x0)*dx, (iy+y0)*dy, (iz+z0)*dz.
    double dx_sp = basis.rspacing()[0];
    double dy_sp = basis.rspacing()[1];
    double dz_sp = basis.rspacing()[2];
    int x0 = basis.cubic_part(0).start();
    int y0 = basis.cubic_part(1).start();
    int z0 = basis.cubic_part(2).start();

    auto phicub_ = begin(phi.hypercubic());

    gpu::run(basis.local_sizes()[2],   // Nz
             basis.local_sizes()[1],   // Ny
             basis.local_sizes()[0],   // Nx
        [=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
            double rx  = (ix + x0) * dx_sp;
            double ry  = (iy + y0) * dy_sp;
            double rz  = (iz + z0) * dz_sp;
            double dx_ = rx - bx;
            double dy_ = ry - by;
            double dz_ = rz - bz;
            double r2  = dx_*dx_ + dy_*dy_ + dz_*dz_;
            double amp = norm_ * exp(-r2 / (2.0 * d * d));
            double phase_ = kx*rx + ky*ry + kz*rz;
            phicub_[ix][iy][iz][ist_wp] = complex(amp * cos(phase_), amp * sin(phase_));
        });
    GPU_SYNC();  // ensure GPU write complete before any CPU read
}

// ── WP validation ─────────────────────────────────────────────────────────────
// Returns {norm, ke_nominal}. Norm should be ≈ 1.0 (±3% at 40 Ha).
// Uses a GPU reduction kernel so the norm is computed on-device — avoids UVM
// page-migration issues when trying to CPU-read GPU-resident orbital data.
inline std::pair<double,double> validate_wp(inq::systems::electrons const & electrons)
{
    using namespace inq;

    auto const & phi  = electrons.kpin()[0];
    double dV = phi.basis().volume_element();
    int ist_wp = phi.set_part().local_size() - 1;
    int n_pts  = phi.basis().local_size();

    // Use same pattern as operations/overlap_diagonal_impl: begin(phi.matrix())
    // captured in GPU_LAMBDA; mat1[ip][ist] accesses element (ip, ist).
    auto mat_    = begin(phi.matrix());
    int  n_st_   = phi.set_part().local_size();

    // Compute ∑_ip |ψ_wp(ip)|² * dV on GPU via reduction
    auto result = gpu::run(1, gpu::reduce(n_pts), 0.0,
        [dV, mat_, ist_wp_=ist_wp] GPU_LAMBDA (auto /*ist*/, auto ip) {
            auto v = mat_[ip][ist_wp_];
            return dV * (v.real()*v.real() + v.imag()*v.imag());
        });
    GPU_SYNC();
    double norm_val = result[0];  // UVM read: single element after GPU reduce
    return {norm_val, cfg::WP_EKIN_HA};
}

// ── Total density 2D slice ────────────────────────────────────────────────────
// ∑_i f_i |ψ_i(x,y,z_target)|² — all states, at nearest z grid plane.
// Returns [Ny_g][Nx_g] array.
inline std::vector<std::vector<double>>
extract_density_slice(inq::systems::electrons const & electrons, double z_target_bohr)
{
    GPU_SYNC();

    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int iz_target = iz_nearest(electrons, z_target_bohr);

    auto const & phi = electrons.kpin()[0];
    int Nx = phi.basis().cubic_part(0).local_size();
    int Ny = phi.basis().cubic_part(1).local_size();
    int Nz = phi.basis().cubic_part(2).local_size();

    std::vector<std::vector<double>> slice(Ny_g, std::vector<double>(Nx_g, 0.0));
    auto const & occ = electrons.occupations()[0];
    auto hc = phi.hypercubic();

    for(int ist = 0; ist < phi.set_part().local_size(); ist++){
        double f = occ[ist];
        if(f == 0.0) continue;
        for(int iz = 0; iz < Nz; iz++){
            auto iz_g = phi.basis().cubic_part(2).local_to_global({iz}).value();
            if(iz_g != iz_target) continue;
            for(int ix = 0; ix < Nx; ix++){
                auto ix_g = phi.basis().cubic_part(0).local_to_global({ix}).value();
                for(int iy = 0; iy < Ny; iy++){
                    auto iy_g = phi.basis().cubic_part(1).local_to_global({iy}).value();
                    auto w = hc[ix][iy][iz][ist];
                    slice[iy_g][ix_g] += f * (w.real()*w.real() + w.imag()*w.imag());
                }
            }
        }
    }
    return slice;
}

// ── Save 2D density slice ─────────────────────────────────────────────────────
// Header: "# t=T z=Z\n". Data: scientific, space-separated rows.
inline void save_density_slice(std::vector<std::vector<double>> const & slice,
                                double time_au, double z_bohr,
                                std::string const & filename)
{
    std::ofstream f(filename);
    f << "# t=" << std::fixed << std::setprecision(6) << time_au
      << " z=" << z_bohr << "\n";
    for(auto const & row : slice){
        for(size_t ix = 0; ix < row.size(); ix++){
            f << std::scientific << std::setprecision(6) << row[ix];
            if(ix + 1 < row.size()) f << " ";
        }
        f << "\n";
    }
}

// ── 1D z-profile at cell centre ───────────────────────────────────────────────
// Returns density[iz_g] at (ix_c=Nx/2, iy_c=Ny/2).
inline std::vector<double>
extract_z_profile(inq::systems::electrons const & electrons)
{
    GPU_SYNC();

    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int Nz_g = basis.sizes()[2];
    int ix_c = Nx_g / 2;
    int iy_c = Ny_g / 2;

    auto const & phi = electrons.kpin()[0];
    int Nx = phi.basis().cubic_part(0).local_size();
    int Ny = phi.basis().cubic_part(1).local_size();
    int Nz = phi.basis().cubic_part(2).local_size();
    auto const & occ = electrons.occupations()[0];
    auto hc = phi.hypercubic();

    std::vector<double> profile(Nz_g, 0.0);
    for(int ist = 0; ist < phi.set_part().local_size(); ist++){
        double f = occ[ist];
        if(f == 0.0) continue;
        for(int ix = 0; ix < Nx; ix++){
            auto ix_g = phi.basis().cubic_part(0).local_to_global({ix}).value();
            if(ix_g != ix_c) continue;
            for(int iy = 0; iy < Ny; iy++){
                auto iy_g = phi.basis().cubic_part(1).local_to_global({iy}).value();
                if(iy_g != iy_c) continue;
                for(int iz = 0; iz < Nz; iz++){
                    auto iz_g = phi.basis().cubic_part(2).local_to_global({iz}).value();
                    auto w = hc[ix][iy][iz][ist];
                    profile[iz_g] += f * (w.real()*w.real() + w.imag()*w.imag());
                }
            }
        }
    }
    return profile;
}

// ── Save a single KS orbital as human-readable text (3D) ─────────────────────
// Format:
//   # ist=I t=T_AU step=N Nx=NX Ny=NY Nz=NZ
//   # Format: one complex value per line, "real imag", C-order (ix slowest)
//   re im
//   ...
// Grid coordinates reconstructed via grid_metadata.txt.
inline void save_orbital_3d(inq::systems::electrons const & electrons,
                             int ist_global,
                             double time_au,
                             int step,
                             std::string const & filename)
{
    GPU_SYNC();

    auto const & phi  = electrons.kpin()[0];
    auto const & basis = phi.basis();

    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int Nz_g = basis.sizes()[2];
    int Nx   = basis.cubic_part(0).local_size();
    int Ny   = basis.cubic_part(1).local_size();
    int Nz   = basis.cubic_part(2).local_size();

    // Map global state index to local (single-rank assumed; ist_global == ist_local)
    int ist_local = ist_global;
    if(ist_local >= phi.set_part().local_size()) return;  // state not on this rank

    // Gather into a global array [Nx_g][Ny_g][Nz_g]
    // Single-rank path: local == global
    using cpx = std::complex<double>;
    std::vector<cpx> global(Nx_g * Ny_g * Nz_g, cpx(0.0, 0.0));

    auto hc = phi.hypercubic();
    for(int ix = 0; ix < Nx; ix++){
        auto ix_g = basis.cubic_part(0).local_to_global({ix}).value();
        for(int iy = 0; iy < Ny; iy++){
            auto iy_g = basis.cubic_part(1).local_to_global({iy}).value();
            for(int iz = 0; iz < Nz; iz++){
                auto iz_g = basis.cubic_part(2).local_to_global({iz}).value();
                auto v = hc[ix][iy][iz][ist_local];
                global[ix_g * Ny_g * Nz_g + iy_g * Nz_g + iz_g] = cpx(v.real(), v.imag());
            }
        }
    }

    std::ofstream f(filename);
    f << "# ist=" << ist_global
      << " t=" << std::fixed << std::setprecision(8) << time_au
      << " step=" << step
      << " Nx=" << Nx_g << " Ny=" << Ny_g << " Nz=" << Nz_g << "\n";
    f << "# Format: one complex value per line, real imag, C-order (ix slowest)\n";
    f << std::scientific << std::setprecision(8);
    for(int idx = 0; idx < Nx_g * Ny_g * Nz_g; idx++){
        f << global[idx].real() << "  " << global[idx].imag() << "\n";
    }
}

// ── Save total electron density as 3D text ────────────────────────────────────
// Format:
//   # t=T_AU step=N Nx=NX Ny=NY Nz=NZ
//   # Format: one float per line, C-order (ix slowest)
//   value
//   ...
inline void save_density_3d(inq::systems::electrons const & electrons,
                             double time_au,
                             int step,
                             std::string const & filename)
{
    GPU_SYNC();

    auto const & phi  = electrons.kpin()[0];
    auto const & basis = phi.basis();

    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int Nz_g = basis.sizes()[2];
    int Nx   = basis.cubic_part(0).local_size();
    int Ny   = basis.cubic_part(1).local_size();
    int Nz   = basis.cubic_part(2).local_size();

    auto const & occ = electrons.occupations()[0];

    std::vector<double> density(Nx_g * Ny_g * Nz_g, 0.0);
    auto hc = phi.hypercubic();

    for(int ist = 0; ist < phi.set_part().local_size(); ist++){
        double f = occ[ist];
        if(f == 0.0) continue;
        for(int ix = 0; ix < Nx; ix++){
            auto ix_g = basis.cubic_part(0).local_to_global({ix}).value();
            for(int iy = 0; iy < Ny; iy++){
                auto iy_g = basis.cubic_part(1).local_to_global({iy}).value();
                for(int iz = 0; iz < Nz; iz++){
                    auto iz_g = basis.cubic_part(2).local_to_global({iz}).value();
                    auto v = hc[ix][iy][iz][ist];
                    density[ix_g * Ny_g * Nz_g + iy_g * Nz_g + iz_g]
                        += f * (v.real()*v.real() + v.imag()*v.imag());
                }
            }
        }
    }

    std::ofstream f(filename);
    f << "# t=" << std::fixed << std::setprecision(8) << time_au
      << " step=" << step
      << " Nx=" << Nx_g << " Ny=" << Ny_g << " Nz=" << Nz_g << "\n";
    f << "# Format: one float per line, C-order (ix slowest)\n";
    f << std::scientific << std::setprecision(8);
    for(double v : density)
        f << v << "\n";
}

// ── Compute full overlap matrix S_ij(t) = <φ_i^GS | φ_j(t)> ─────────────────
// i indexes GS orbitals (0..n_states-1), j indexes time-evolved orbitals.
// Returns n×n complex matrix.
// GPU sync is applied once at entry.
// Computational cost: n²×N_grid ≈ 57²×1.68M ≈ 5.5G FLOPs (CPU).
// Called every OVERLAP_INTERVAL=5 steps → ~100×, ~9 min total overhead.
inline std::vector<std::vector<std::complex<double>>>
compute_overlap_matrix(inq::systems::electrons const & gs_electrons,
                       inq::systems::electrons const & electrons)
{
    GPU_SYNC();

    auto const & phi_gs = gs_electrons.kpin()[0];
    auto const & phi_t  = electrons.kpin()[0];
    auto const & basis  = phi_t.basis();

    int n_states = phi_t.set_part().local_size();
    double dV    = basis.volume_element();

    int Nx = basis.cubic_part(0).local_size();
    int Ny = basis.cubic_part(1).local_size();
    int Nz = basis.cubic_part(2).local_size();

    auto hc_gs = phi_gs.hypercubic();
    auto hc_t  = phi_t.hypercubic();

    // S[i][j] = <φ_i^GS | φ_j(t)>
    using cpx = std::complex<double>;
    std::vector<std::vector<cpx>> S(n_states, std::vector<cpx>(n_states, cpx(0.0, 0.0)));

    for(int i = 0; i < n_states; i++){
        for(int j = 0; j < n_states; j++){
            cpx acc(0.0, 0.0);
            for(int ix = 0; ix < Nx; ix++){
                for(int iy = 0; iy < Ny; iy++){
                    for(int iz = 0; iz < Nz; iz++){
                        auto vi = hc_gs[ix][iy][iz][i];
                        auto vj = hc_t [ix][iy][iz][j];
                        // conj(vi) * vj
                        acc += cpx( vi.real()*vj.real() + vi.imag()*vj.imag(),
                                    vi.real()*vj.imag() - vi.imag()*vj.real() ) * dV;
                    }
                }
            }
            S[i][j] = acc;
        }
    }
    return S;
}

// ── Save overlap matrix block to text file (appended) ─────────────────────────
inline void save_overlap_matrix(
    std::vector<std::vector<std::complex<double>>> const & S,
    double time_au, int step,
    std::ofstream & f)
{
    int n = static_cast<int>(S.size());
    f << "# step=" << step
      << " t=" << std::fixed << std::setprecision(8) << time_au
      << " n_states=" << n << "\n";
    f << "# S_ij = <phi_i_GS | phi_j(t)>  i=row j=col\n";
    f << "# Format: n rows of n pairs \"re im\" separated by spaces\n";
    f << std::scientific << std::setprecision(6);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            f << S[i][j].real() << "  " << S[i][j].imag();
            if(j + 1 < n) f << "   ";
        }
        f << "\n";
    }
    f << "\n";  // blank line between blocks
}

// ── Save grid coordinates (called once after SCF) ────────────────────────────
// Writes results/grid/grid_x.txt, grid_y.txt, grid_z.txt, grid_metadata.txt.
// Each coordinate file: one bohr value per line.
inline void save_grid_coords(inq::systems::electrons const & electrons,
                             std::string const & outdir)
{
    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int Nz_g = basis.sizes()[2];
    double dx = basis.rspacing()[0];
    double dy = basis.rspacing()[1];
    double dz = basis.rspacing()[2];

    // x-coords
    {
        std::ofstream f(outdir + "/grid_x.txt");
        f << "# x-coordinates of grid points in bohr (ix=0..Nx-1, x = ix*dx)\n";
        f << std::scientific << std::setprecision(8);
        for(int ix = 0; ix < Nx_g; ix++) f << ix * dx << "\n";
    }
    // y-coords
    {
        std::ofstream f(outdir + "/grid_y.txt");
        f << "# y-coordinates of grid points in bohr\n";
        f << std::scientific << std::setprecision(8);
        for(int iy = 0; iy < Ny_g; iy++) f << iy * dy << "\n";
    }
    // z-coords
    {
        std::ofstream f(outdir + "/grid_z.txt");
        f << "# z-coordinates of grid points in bohr\n";
        f << std::scientific << std::setprecision(8);
        for(int iz = 0; iz < Nz_g; iz++) f << iz * dz << "\n";
    }
    // metadata
    {
        std::ofstream f(outdir + "/grid_metadata.txt");
        f << "# Grid metadata for run_003\n";
        f << "# Columns: Nx Ny Nz dx_bohr dy_bohr dz_bohr Lx_bohr Ly_bohr Lz_bohr\n";
        f << Nx_g << " " << Ny_g << " " << Nz_g << " "
          << std::scientific << std::setprecision(8)
          << dx << " " << dy << " " << dz << " "
          << Nx_g * dx << " " << Ny_g * dy << " " << Nz_g * dz << "\n";
    }
}

} // namespace leed_utils
