/*
 * Writes ground-state Kohn-Sham eigenvalues, occupations, and k-point
 * coordinates / weights to a long-format CSV file. Intended to be called
 * once after electrons.load() or SCF convergence at the start of a TDDFT
 * run, so the time-zero band structure is always recoverable.
 *
 * CSV schema (one row per (k-point, state) pair)
 * -----------------------------------------------
 *   kpoint_index   Global k-point index.
 *   kx, ky, kz     K-point coordinates in units of 2π/a, matching the
 *                  convention in inq's brillouin.hpp.
 *   weight         Brillouin zone integration weight for this k-point.
 *   state_index    Band index (zero-based).
 *   eigenvalue_ha  Kohn-Sham eigenvalue in Hartree.
 *   eigenvalue_ev  Kohn-Sham eigenvalue in eV  (× 27.21138625).
 *   occupation     State occupation number.
 *
 * Usage
 * -----
 *   inqkit::observables::dump_eigenvalues(electrons, "/output/eigenvalues.csv");
 *
 * Multi-k-point behaviour
 * -----------------------
 * Each MPI rank writes its local kpin slice (kpin_part().start() to
 * kpin_part().start() + kpin_size()). For single-rank runs the local
 * slice covers the full Brillouin zone, producing the complete eigenvalue
 * table in one file. Parallel multi-rank output is not yet merged
 * automatically.
 *
 * Note: single-rank only for production use, consistent with the existing
 * inqkit writers.
 */
#pragma once

#include <inq/inq.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace inqkit::observables {

inline constexpr double HA_TO_EV = 27.21138625;


/* 
* All the states per k point have their eigen-energy dumped to  
* a CSV file. 
*/

inline void dump_eigenvalues(inq::systems::electrons const& electrons,
                             std::string const& csv_path) {
    namespace fs = std::filesystem;
    if (auto parent = fs::path(csv_path).parent_path(); !parent.empty())
        fs::create_directories(parent);

    // Path to offload the csv into 
    std::ofstream f(csv_path);
    if (!f) throw std::runtime_error("dump_eigenvalues: cannot open '" + csv_path + "'");

    f << "kpoint_index,kx,ky,kz,weight,state_index,"
         "eigenvalue_ha,eigenvalue_ev,occupation\n";
    f << std::setprecision(15);

        
    // TODO: Understand the data structure type of these different
    // properties of electrons
    auto const& bz       = electrons.brillouin_zone(); // BZ k point mesh
    auto const& weights  = electrons.kpin_weights(); // Integration weights of each point
    auto const& eigs_arr = electrons.eigenvalues(); // Eigenvalues per k point
    auto const& occs_arr = electrons.occupations(); // occupations (of each state) per k point

    // 
    const long kp_offset = electrons.kpin_part().start();
    const long nk_local  = electrons.kpin_size();

    for (long ilot = 0; ilot < nk_local; ++ilot) {
        const long ik_global = kp_offset + ilot;
        auto kp = bz.kpoint(static_cast<int>(ik_global));
        // Match brillouin.hpp print convention: divide by 2*pi for the
        // tabulated kx/ky/kz so values are in "fractions of reciprocal
        // lattice vectors" (cubic cell only — for non-cubic cells the
        // user should consult cell-relative kpoint() output).
        const double kx = kp[0] / (2.0 * M_PI);
        const double ky = kp[1] / (2.0 * M_PI);
        const double kz = kp[2] / (2.0 * M_PI);
        const double w  = static_cast<double>(weights[ilot]);

        auto const& eigs = eigs_arr[ilot];
        auto const& occs = occs_arr[ilot];
        const long nst   = static_cast<long>(eigs.size());

        for (long ist = 0; ist < nst; ++ist) {
            const double e_ha = static_cast<double>(eigs[ist]);
            const double occ  = (ist < static_cast<long>(occs.size()))
                                  ? static_cast<double>(occs[ist]) : 0.0;
            f << ik_global << ',' << kx << ',' << ky << ',' << kz << ','
              << w << ',' << ist << ','
              << e_ha << ',' << e_ha * HA_TO_EV << ','
              << occ << '\n';
        }
    }
    f.flush();
}

} // namespace inqkit::observables
