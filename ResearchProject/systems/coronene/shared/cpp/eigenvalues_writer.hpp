// ============================================================================
// shared/cpp/eigenvalues_writer.hpp
//
// Writes the gamma-point KS orbital eigenvalues + occupations to a CSV
// pair (eigenvalues.csv, occupations.csv). Called from:
//
//   * each save_gs/<sig>/run.cpp, immediately after the SCF + electrons.save(),
//     writing into both <CHECKPOINT_DIR>/ and the run's results/ subtree.
//   * run_template.hpp's run_propagation, which invokes the
//     copy_eigenvalues_from_checkpoint helper to mirror the canonical files
//     from the checkpoint into results/raw/observables/eigenvalues/.
//
// Single-kpoint (Gamma-only) only — matches the rest of the framework.
// ============================================================================
#pragma once

#include <inq/inq.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <string>

namespace coronene::eigenvalues {

inline constexpr double HA_TO_EV = 27.21138625;

// Write eigenvalues.csv + occupations.csv into ``out_dir``. Schema:
//   eigenvalues.csv: state_index,eigenvalue_ha,eigenvalue_ev
//   occupations.csv: state_index,occupation
inline void dump(inq::systems::electrons const &electrons,
                 std::string const &out_dir) {
    std::filesystem::create_directories(out_dir);

    auto const &eigs = electrons.eigenvalues()[0];
    auto const &occs = electrons.occupations()[0];
    int n = static_cast<int>(eigs.size());

    {
        std::ofstream f(out_dir + "/eigenvalues.csv");
        f << "state_index,eigenvalue_ha,eigenvalue_ev\n";
        f << std::setprecision(16);
        for (int i = 0; i < n; ++i) {
            f << i << "," << eigs[i] << "," << eigs[i] * HA_TO_EV << "\n";
        }
    }
    {
        std::ofstream f(out_dir + "/occupations.csv");
        f << "state_index,occupation\n";
        f << std::setprecision(16);
        int n_occ = static_cast<int>(occs.size());
        for (int i = 0; i < n_occ; ++i) {
            f << i << "," << occs[i] << "\n";
        }
    }
}

// Copy eigenvalues.csv + occupations.csv from <checkpoint_dir>/ into
// <results_eig_dir>/ if they exist there. Silent no-op if the
// checkpoint doesn't carry them (e.g. legacy checkpoints from before
// Phase 4); the retrofit script
// scripts/extract_eigenvalues_from_log.py is the fallback.
inline void copy_from_checkpoint(std::string const &checkpoint_dir,
                                 std::string const &results_eig_dir) {
    namespace fs = std::filesystem;
    fs::create_directories(results_eig_dir);
    for (auto const &name : {"eigenvalues.csv", "occupations.csv"}) {
        fs::path src = fs::path(checkpoint_dir) / name;
        if (!fs::exists(src)) continue;
        fs::path dst = fs::path(results_eig_dir) / name;
        fs::copy_file(src, dst, fs::copy_options::overwrite_existing);
    }
}

}  // namespace coronene::eigenvalues
