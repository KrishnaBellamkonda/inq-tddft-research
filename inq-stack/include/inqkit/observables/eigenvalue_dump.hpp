// ============================================================================
// inqkit/observables/eigenvalue_dump.hpp
//
// Writes ground-state Kohn-Sham eigenvalues + occupations + k-point
// coordinates / weights for every k-point in the Brillouin zone, in
// long-format CSV. Intended to be called once after `electrons.load(...)`
// (or after SCF) at the start of a TDDFT run, so the time-zero band
// structure is always recoverable.
//
// CSV schema (one row per (k, state)):
//
//   kpoint_index, kx, ky, kz, weight,
//   state_index, eigenvalue_ha, eigenvalue_ev, occupation
//
// Where (kx, ky, kz) is in units of 2*pi / a_lattice, matching the
// convention used by inq's brillouin printer (brillouin.hpp:139).
//
// Single-rank assumption: same as state_energy_writer.hpp — each MPI rank
// writes its local kpin slice; for the Li 54-atom 2x2x2 single-GPU runs
// there is only one rank, so kpin_part().start() == 0 and we get the
// global eigenvalue table.
//
// Adapted from coronene::eigenvalues::dump
// (ResearchProject/systems/coronene/shared/cpp/eigenvalues_writer.hpp),
// extended to multi-kpoint.
// ============================================================================
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

inline void dump_eigenvalues(inq::systems::electrons const& electrons,
                             std::string const& csv_path) {
    namespace fs = std::filesystem;
    if (auto parent = fs::path(csv_path).parent_path(); !parent.empty())
        fs::create_directories(parent);

    std::ofstream f(csv_path);
    if (!f) throw std::runtime_error("dump_eigenvalues: cannot open '" + csv_path + "'");

    f << "kpoint_index,kx,ky,kz,weight,state_index,"
         "eigenvalue_ha,eigenvalue_ev,occupation\n";
    f << std::setprecision(15);

    auto const& bz       = electrons.brillouin_zone();
    auto const& weights  = electrons.kpin_weights();
    auto const& eigs_arr = electrons.eigenvalues();
    auto const& occs_arr = electrons.occupations();

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
