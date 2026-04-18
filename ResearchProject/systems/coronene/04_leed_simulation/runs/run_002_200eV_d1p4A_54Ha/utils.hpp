#pragma once
// ============================================================================
// run_002/utils.hpp
//
// Utility functions for the LEED simulation:
//   leed_utils::inject_wp          — write Gaussian WP into last KS orbital
//   leed_utils::validate_wp        — check WP norm
//   leed_utils::extract_density_slice — total density at fixed z (all orbitals)
//   leed_utils::extract_wp_slice   — WP orbital density at fixed z (ist_wp only)
//   leed_utils::extract_z_profile  — 1D density profile along z at (ix_c, iy_c)
//   leed_utils::save_density_slice — write 2D slice to text file
//   leed_utils::iz_nearest         — map z_bohr to nearest grid index
// ============================================================================

#include "config.hpp"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>

namespace leed_utils {

// ── Nearest grid index for a given z value ────────────────────────────────────
inline int iz_nearest(inq::systems::electrons const & electrons, double z_bohr) {
    auto const & basis = electrons.states_basis();
    int Nz_g = basis.sizes()[2];
    double dz = basis.rspacing()[2];
    int iz = static_cast<int>(std::round(z_bohr / dz));
    return ((iz % Nz_g) + Nz_g) % Nz_g;
}

// ── WP injection ─────────────────────────────────────────────────────────────
// Writes ψ^WP(r) = (1/(πd²))^{3/4} exp(−|r−b|²/(2d²)) exp(ik·r)
// into the LAST KS orbital slot of kpin[0].
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

    int ist_wp = phi.set_part().local_size() - 1;
    const double d    = cfg::WP_D_BOHR;
    const double norm = cfg::wp_norm();

    auto hc = phi.hypercubic();
    for(int ix = 0; ix < Nx; ix++){
        for(int iy = 0; iy < Ny; iy++){
            for(int iz = 0; iz < Nz; iz++){
                auto r   = po.rvector_cartesian(ix, iy, iz);
                double dx = r[0] - bx;
                double dy = r[1] - by;
                double dz = r[2] - bz;
                double r2    = dx*dx + dy*dy + dz*dz;
                double amp   = norm * std::exp(-r2 / (2.0 * d * d));
                double phase = kx*r[0] + ky*r[1] + kz*r[2];
                hc[ix][iy][iz][ist_wp] = complex(amp*std::cos(phase),
                                                  amp*std::sin(phase));
            }
        }
    }
}

// ── WP validation ─────────────────────────────────────────────────────────────
// Returns {norm, ke_estimate} of the WP orbital.
inline std::pair<double,double> validate_wp(inq::systems::electrons const & electrons)
{
    auto const & phi   = electrons.kpin()[0];
    auto const & basis = phi.basis();
    int Nx = basis.cubic_part(0).local_size();
    int Ny = basis.cubic_part(1).local_size();
    int Nz = basis.cubic_part(2).local_size();
    double dV   = basis.volume_element();
    int ist_wp  = phi.set_part().local_size() - 1;

    auto hc = phi.hypercubic();
    double norm_val = 0.0;
    for(int ix = 0; ix < Nx; ix++)
        for(int iy = 0; iy < Ny; iy++)
            for(int iz = 0; iz < Nz; iz++){
                auto v = hc[ix][iy][iz][ist_wp];
                norm_val += (v.real()*v.real() + v.imag()*v.imag()) * dV;
            }
    return {norm_val, cfg::WP_EKIN_HA};
}

// ── Total density 2D slice ────────────────────────────────────────────────────
// Sum |ψᵢ|² × fᵢ for ALL states at the nearest grid plane to z_target.
inline std::vector<std::vector<double>>
extract_density_slice(inq::systems::electrons const & electrons, double z_target_bohr)
{
    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int iz_target = iz_nearest(electrons, z_target_bohr);

    auto & phi = electrons.kpin()[0];
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

// ── WP-only orbital 2D slice ─────────────────────────────────────────────────
// Returns |ψ_WP|² at the nearest grid plane to z_target (ist_wp only, f=1).
inline std::vector<std::vector<double>>
extract_wp_slice(inq::systems::electrons const & electrons, double z_target_bohr)
{
    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int iz_target = iz_nearest(electrons, z_target_bohr);

    auto & phi = electrons.kpin()[0];
    int Nx = phi.basis().cubic_part(0).local_size();
    int Ny = phi.basis().cubic_part(1).local_size();
    int Nz = phi.basis().cubic_part(2).local_size();
    int ist_wp = phi.set_part().local_size() - 1;

    std::vector<std::vector<double>> slice(Ny_g, std::vector<double>(Nx_g, 0.0));
    auto hc = phi.hypercubic();

    for(int iz = 0; iz < Nz; iz++){
        auto iz_g = phi.basis().cubic_part(2).local_to_global({iz}).value();
        if(iz_g != iz_target) continue;
        for(int ix = 0; ix < Nx; ix++){
            auto ix_g = phi.basis().cubic_part(0).local_to_global({ix}).value();
            for(int iy = 0; iy < Ny; iy++){
                auto iy_g = phi.basis().cubic_part(1).local_to_global({iy}).value();
                auto w = hc[ix][iy][iz][ist_wp];
                slice[iy_g][ix_g] = w.real()*w.real() + w.imag()*w.imag();
            }
        }
    }
    return slice;
}

// ── 1D z-profile at cell centre (ix_c, iy_c) ─────────────────────────────────
// Returns density[iz_global] at the grid point (ix_c, iy_c) for each z.
// Used to track density along the WP trajectory.
inline std::vector<double>
extract_z_profile(inq::systems::electrons const & electrons)
{
    auto const & basis = electrons.states_basis();
    int Nx_g = basis.sizes()[0];
    int Ny_g = basis.sizes()[1];
    int Nz_g = basis.sizes()[2];
    int ix_c  = Nx_g / 2;  // cell centre in x
    int iy_c  = Ny_g / 2;  // cell centre in y

    auto & phi = electrons.kpin()[0];
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

// ── Save 2D density slice ─────────────────────────────────────────────────────
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
