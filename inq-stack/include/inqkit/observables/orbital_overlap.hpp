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

    // WP-only variant: compute O_i,wp(t) = |<psi_i^GS | psi_wp(t)>|^2 for
    // every reference orbital i in 0..n_ref-1. Cost is n_ref complex dot
    // products per call instead of n_ref*(n_ref+1) for the full matrix —
    // ~62x cheaper for the coronene paper-replica setup.
    //
    // Output format:
    //   overlap_XXXXXX.csv contains a header line plus a single row of
    //   n_ref comma-separated |overlap|^2 values, in column order i=0..n_ref-1.
    //   index.csv entries match snapshot()'s schema (step,time_au,file).
    //
    // Use this when you only need to track decoherence / decay of the WP
    // into the occupied subspace, not the full KS-orbital cross-talk matrix.
    void snapshot_wp_only(inq::systems::electrons const& electrons,
                          double time_au, int step)
    {
        // Pull only the WP wavefunction (column j = n_ref) from the device.
        auto wp = fields::orbital::wavefunction(electrons, n_ref_);

        std::size_t n_pts = wp.values.size();
        std::vector<double> O(n_ref_, 0.0);
        for (int i = 0; i < n_ref_; ++i) {
            std::complex<double> inner(0.0, 0.0);
            auto const& ri = ref_wfns_[i].values;
            auto const& ej = wp.values;
            for (std::size_t r = 0; r < n_pts; ++r)
                inner += std::conj(ri[r]) * ej[r];
            inner *= dv_;
            O[i] = std::norm(inner);
        }

        std::ostringstream ss;
        ss << output_dir_ << "/overlap_"
           << std::setfill('0') << std::setw(6) << step << ".csv";
        std::string fname = ss.str();

        std::ofstream f(fname);
        if (!f)
            throw std::runtime_error("OrbitalOverlapMatrix: cannot open " + fname);

        f << "# step=" << step
          << " time_au=" << std::fixed << std::setprecision(6) << time_au
          << " n_ref=" << n_ref_ << " mode=wp_only\n";
        for (int i = 0; i < n_ref_; ++i) {
            f << std::fixed << std::setprecision(8) << O[i];
            if (i + 1 < n_ref_) f << ",";
        }
        f << "\n";

        std::string basename = std::filesystem::path(fname).filename().string();
        index_file_ << step << ","
                    << std::fixed << std::setprecision(6) << time_au << ","
                    << basename << "\n";
        index_file_.flush();
    }

    // ----- Proxy variant -------------------------------------------------
    // Compute O_{i, j} for fixed i in 0..n_ref-1 and j in proxy_indices
    // (typically 2 orbitals per shell of degeneracy). This is the scalable
    // form for the GS-projected occupation observable:
    //
    //   n_i^GS(t) ≈ sum_shells (g_s / N_proxies(s)) sum_{j in P_s}
    //                f_j(0) * |<psi_i^GS | psi_j(t)>|^2
    //
    // Cost: |proxy_indices| evolved-wavefunction extractions and
    //       n_ref * |proxy_indices| dot products per call.
    //
    // Output format:
    //   overlap_XXXXXX.csv:
    //     # step=N time_au=T n_ref=R n_proxies=P proxy_indices=p1,p2,...
    //     row 0: O_{0,p1}, O_{0,p2}, ..., O_{0,pP}
    //     row 1: O_{1,p1}, O_{1,p2}, ..., O_{1,pP}
    //     ...
    //   index.csv: step,time_au,file (shared with snapshot/snapshot_wp_only)
    //
    // The caller is expected to also write a sibling shells.csv mapping
    // proxy_indices to (shell_id, degeneracy) so that the postprocess can
    // do the shell-averaged projection — this method does not write it.
    void snapshot_proxies(inq::systems::electrons const& electrons,
                          std::vector<int> const& proxy_indices,
                          double time_au, int step)
    {
        const int P = static_cast<int>(proxy_indices.size());

        // Extract evolved wavefunctions for proxies only.
        std::vector<fields::ComplexField3D> evolved;
        evolved.reserve(P);
        for (int p = 0; p < P; ++p)
            evolved.push_back(
                fields::orbital::wavefunction(electrons, proxy_indices[p]));

        std::size_t n_pts = ref_wfns_.empty() ? 0 : ref_wfns_[0].values.size();

        std::vector<std::vector<double>> O(n_ref_, std::vector<double>(P, 0.0));
        for (int i = 0; i < n_ref_; ++i) {
            for (int p = 0; p < P; ++p) {
                std::complex<double> inner(0.0, 0.0);
                auto const& ri = ref_wfns_[i].values;
                auto const& ej = evolved[p].values;
                for (std::size_t r = 0; r < n_pts; ++r)
                    inner += std::conj(ri[r]) * ej[r];
                inner *= dv_;
                O[i][p] = std::norm(inner);
            }
        }

        std::ostringstream ss;
        ss << output_dir_ << "/overlap_"
           << std::setfill('0') << std::setw(6) << step << ".csv";
        std::string fname = ss.str();

        std::ofstream f(fname);
        if (!f)
            throw std::runtime_error("OrbitalOverlapMatrix: cannot open " + fname);

        f << "# step=" << step
          << " time_au=" << std::fixed << std::setprecision(6) << time_au
          << " n_ref=" << n_ref_ << " n_proxies=" << P
          << " proxy_indices=";
        for (int p = 0; p < P; ++p) {
            f << proxy_indices[p];
            if (p + 1 < P) f << ",";
        }
        f << "\n";
        for (int i = 0; i < n_ref_; ++i) {
            for (int p = 0; p < P; ++p) {
                f << std::fixed << std::setprecision(8) << O[i][p];
                if (p + 1 < P) f << ",";
            }
            f << "\n";
        }

        std::string basename = std::filesystem::path(fname).filename().string();
        index_file_ << step << ","
                    << std::fixed << std::setprecision(6) << time_au << ","
                    << basename << "\n";
        index_file_.flush();
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
