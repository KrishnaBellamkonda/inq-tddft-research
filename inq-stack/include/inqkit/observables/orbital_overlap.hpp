#pragma once

// inqkit::observables::OrbitalOverlapMatrix
//
// Tracks the time-evolution of KS orbital character by computing:
//
//   O_ij(t) = |<ψ_i^GS | ψ_j(t)>|²
//
// where i indexes the reference (ground-state) orbital set (0..n_ref-1)
// and j indexes evolved orbitals (0..n_ref, inclusive — column n_ref is the WP).
//
// Call snapshot() at each RT step to record the n_ref × (n_ref+1) matrix.
// Results are written to output_dir/overlap_XXXXXX.csv with an index file.
//
// Validation: at t=0 (step 0), columns 0..n_ref-1 should be identity-like
// (O_jj≈1, rest≈0) and column n_ref (WP) should be ≈0 everywhere
// (WP was orthogonalised against occupied states before injection).
//
// Performance note: snapshot() extracts (n_ref+1) complex wavefunctions from
// GPU memory per call.  For a 38-electron system (n_ref≈19) and 80³ grid,
// this is ~21 GPU syncs and ~380 complex dot products of length 512K per step.
// On a modern CPU this is ~50–500ms per step; disable by commenting out the
// snapshot() call if run time is critical.

#include <inqkit/fields/complex_field_3d.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inq/inq.hpp>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace inqkit::observables {

class OrbitalOverlapMatrix {
    std::vector<fields::ComplexField3D> ref_wfns_;
    int    n_ref_;       // number of GS orbitals = wp_idx
    int    n_evolved_;   // n_ref_ + 1 (last column = WP orbital)
    double dv_;          // grid-cell volume dx*dy*dz (bohr³)
    std::string output_dir_;
    std::ofstream index_file_;

public:
    // Construct from electrons state *after* WP injection.
    // n_ref = wp_idx (the state index of the WP orbital).
    // Reference wavefunctions 0..n_ref-1 are extracted at construction time;
    // they are the GS orbitals and are unchanged by the injection.
    OrbitalOverlapMatrix(inq::systems::electrons const& electrons,
                         int n_ref,
                         std::string const& output_dir)
        : n_ref_(n_ref), n_evolved_(n_ref + 1), output_dir_(output_dir)
    {
        std::filesystem::create_directories(output_dir_);

        ref_wfns_.reserve(n_ref_);
        for (int i = 0; i < n_ref_; ++i)
            ref_wfns_.push_back(fields::orbital::wavefunction(electrons, i));

        if (!ref_wfns_.empty()) {
            auto const& f = ref_wfns_[0];
            dv_ = f.dx_bohr * f.dy_bohr * f.dz_bohr;
        }

        index_file_.open(output_dir_ + "/index.csv");
        if (!index_file_)
            throw std::runtime_error("OrbitalOverlapMatrix: cannot open index file in " + output_dir_);
        index_file_ << "step,time_au,file\n";
    }

    ~OrbitalOverlapMatrix() {
        if (index_file_.is_open()) index_file_.close();
    }

    // Compute O_ij(t) and write to output_dir/overlap_XXXXXX.csv.
    // Rows = GS reference orbital index i (0..n_ref-1).
    // Columns = evolved orbital index j (0..n_ref); last column is the WP.
    void snapshot(inq::systems::electrons const& electrons,
                  double time_au, int step)
    {
        // Extract evolved wavefunctions 0..n_evolved_-1
        std::vector<fields::ComplexField3D> evolved;
        evolved.reserve(n_evolved_);
        for (int j = 0; j < n_evolved_; ++j)
            evolved.push_back(fields::orbital::wavefunction(electrons, j));

        std::size_t n_pts = ref_wfns_.empty() ? 0 : ref_wfns_[0].values.size();

        // Compute overlap matrix
        std::vector<std::vector<double>> O(n_ref_, std::vector<double>(n_evolved_, 0.0));
        for (int i = 0; i < n_ref_; ++i) {
            for (int j = 0; j < n_evolved_; ++j) {
                std::complex<double> inner(0.0, 0.0);
                auto const& ri = ref_wfns_[i].values;
                auto const& ej = evolved[j].values;
                for (std::size_t r = 0; r < n_pts; ++r)
                    inner += std::conj(ri[r]) * ej[r];
                inner *= dv_;
                O[i][j] = std::norm(inner);  // |inner|^2
            }
        }

        // Build filename
        std::ostringstream ss;
        ss << output_dir_ << "/overlap_"
           << std::setfill('0') << std::setw(6) << step << ".csv";
        std::string fname = ss.str();

        std::ofstream f(fname);
        if (!f)
            throw std::runtime_error("OrbitalOverlapMatrix: cannot open " + fname);

        f << "# step=" << step
          << " time_au=" << std::fixed << std::setprecision(6) << time_au
          << " n_ref=" << n_ref_ << " n_evolved=" << n_evolved_ << "\n";
        for (int i = 0; i < n_ref_; ++i) {
            for (int j = 0; j < n_evolved_; ++j) {
                f << std::fixed << std::setprecision(8) << O[i][j];
                if (j + 1 < n_evolved_) f << ",";
            }
            f << "\n";
        }

        // Update index
        std::string basename = std::filesystem::path(fname).filename().string();
        index_file_ << step << ","
                    << std::fixed << std::setprecision(6) << time_au << ","
                    << basename << "\n";
        index_file_.flush();
    }
};

} // namespace inqkit::observables
