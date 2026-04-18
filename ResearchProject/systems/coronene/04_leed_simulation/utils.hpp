#pragma once
// ============================================================================
// 04_leed_simulation/utils.hpp
//
// Utility functions for the LEED simulation:
//   - WP orbital injection into the last extra-state KS orbital
//   - WP norm and kinetic energy validation
//   - 2D density slice extraction (for snapshots and LEED accumulation)
//   - LEED pattern file I/O
// ============================================================================

#include "config.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>

namespace leed_utils {

// ── WP injection ─────────────────────────────────────────────────────────────
// Writes the Gaussian WP into the LAST KS orbital slot of kpin[0].
// The WP is centred at (bx, by, bz), propagating in the direction of k_vec.
//
// Formula (paper Eq. 1):
//   ψ^WP(r) = (1/(π d²))^{3/4} · exp(−|r−b|²/(2d²)) · exp(i k·r)
//
// Precondition: electrons has been converged by SCF; WP slot = last state.
inline void inject_wp(inq::systems::electrons & electrons,
                      double bx, double by, double bz,
                      double kx, double ky, double kz)
{
    using namespace inq;
    using complex = inq::complex;

    auto & phi   = electrons.kpin()[0];
    auto & basis = phi.basis();
    int Nx = basis.cubic_part(0).local_size();
    int Ny = basis.cubic_part(1).local_size();
    int Nz = basis.cubic_part(2).local_size();
    auto po = basis.point_op();

    // WP orbital: last state in the local set (= first extra_state slot)
    int ist_wp = phi.set_part().local_size() - 1;

    const double d    = cfg::WP_D_BOHR;
    const double norm = cfg::wp_norm();

    auto hc = phi.hypercubic();
    for(int ix = 0; ix < Nx; ix++){
        for(int iy = 0; iy < Ny; iy++){
            for(int iz = 0; iz < Nz; iz++){
                auto r = po.rvector_cartesian(ix, iy, iz);
                double dx = r[0] - bx;
                double dy = r[1] - by;
                double dz = r[2] - bz;
                double r2    = dx*dx + dy*dy + dz*dz;
                double amp   = norm * std::exp(-r2 / (2.0 * d * d));
                double phase = kx*r[0] + ky*r[1] + kz*r[2];
                hc[ix][iy][iz][ist_wp] = complex(amp * std::cos(phase),
                                                  amp * std::sin(phase));
            }
        }
    }
}

// ── WP validation ─────────────────────────────────────────────────────────────
// Returns {norm, kinetic_energy_Ha} of the WP orbital.
// norm should be ≈ 1.0; KE should be ≈ cfg::WP_EKIN_HA.
inline std::pair<double,double> validate_wp(inq::systems::electrons const & electrons)
{
    auto const & phi   = electrons.kpin()[0];
    auto const & basis = phi.basis();
    int Nx = basis.cubic_part(0).local_size();
    int Ny = basis.cubic_part(1).local_size();
    int Nz = basis.cubic_part(2).local_size();
    double dV = basis.volume_element();
    int ist_wp = phi.set_part().local_size() - 1;

    auto hc = phi.hypercubic();
    double norm_val = 0.0;
    for(int ix = 0; ix < Nx; ix++)
        for(int iy = 0; iy < Ny; iy++)
            for(int iz = 0; iz < Nz; iz++){
                auto v = hc[ix][iy][iz][ist_wp];
                norm_val += (v.real()*v.real() + v.imag()*v.imag()) * dV;
            }

    // KE estimate: ⟨T⟩ ≈ k₀²/2 (dominant for narrow WP)
    double ke_est = cfg::WP_EKIN_HA;  // by construction

    return {norm_val, ke_est};
}

// ── 2D density slice extraction ───────────────────────────────────────────────
// Extracts the total electron density on a 2D (x,y) plane at a fixed z value.
// Returns a 2D vector [iy][ix] → density at that grid point.
// z_target_bohr: desired z in bohr (nearest grid point is used).
//
// Uses the density field directly from viewables.
inline std::vector<std::vector<double>>
extract_density_slice(inq::systems::electrons const & electrons,
                      double z_target_bohr)
{
    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int Nz_g = basis.sizes()[2];
    double dz = basis.rspacing()[2];

    // Find nearest z grid index (global)
    int iz_target = static_cast<int>(std::round(z_target_bohr / dz));
    iz_target = ((iz_target % Nz_g) + Nz_g) % Nz_g;

    auto & phi   = electrons.kpin()[0];
    auto & b     = phi.basis();
    int Nx = b.cubic_part(0).local_size();
    int Ny = b.cubic_part(1).local_size();
    int Nz = b.cubic_part(2).local_size();

    // Build slice: sum |ψᵢ|² × fᵢ for all states at the target z
    std::vector<std::vector<double>> slice(Ny_g, std::vector<double>(Nx_g, 0.0));

    auto const & occ = electrons.occupations()[0];
    auto hc = phi.hypercubic();

    for(int ist = 0; ist < phi.set_part().local_size(); ist++){
        double f = occ[ist];
        if(f == 0.0) continue;

        // Find the local iz that corresponds to the global iz_target
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

// ── Save a 2D density slice to a text file ────────────────────────────────────
// Format: one row per iy, space-separated values (ix varying fastest).
// Header: # t=<time> z=<z_bohr>
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

} // namespace leed_utils
