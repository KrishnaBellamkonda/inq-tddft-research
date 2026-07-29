// ============================================================================
// inqkit/observables/momentum_distribution.hpp
//
// On-the-fly one-body momentum distribution for a real-time TDDFT run.
//
//     n(|k|, t) = sum_i f_i |phi-tilde_i(k, t)|^2          (total)
//     n_wp(|k|, t) =     |phi-tilde_wp(k, t)|^2 * f_wp     (wave-packet)
//
// where phi-tilde_i is the forward FFT of the i-th KS orbital, f_i is its
// occupation, and the sum runs over every orbital in every k-point. For
// gamma-only jellium, this collapses to a single k-point. Cross-terms
// vanish under unitary propagation because the KS orbitals stay
// orthonormal in real space, so this is the diagonal of the one-body
// reduced density matrix in momentum space.
//
// We bin |k| (Bohr^-1) into n_bins uniform-width bins on [0, k_max]. For
// jellium L = 60 Bohr the natural Delta-k = 2*pi/L ~= 0.105 Bohr^-1, and
// k_max = pi/spacing ~= pi/0.50 = 6.28 Bohr^-1, so 64 bins gives
// dk ~= 0.10 Bohr^-1, comparable to Delta-k.
//
// Per-state normalisation: each orbital is renormalised to 1 in the
// histogram (sum of bins = f_i for a normalised orbital), so the total
// integrates to sum_i f_i = total electron count, irrespective of FFT
// pre-factor conventions.
//
// CSV layout (long format, one row per (step, |k|-bin)):
//     step,time_au,k_bohr_inv,n_total,n_wp
// Header line includes box length L_BOHR and dk so the postprocess
// can label axes.
// ============================================================================
#pragma once

#include <inq/inq.hpp>
#include <operations/transform.hpp>
#include <basis/fourier_space.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace inqkit::observables {

struct MomentumDistributionConfig {
    int    n_bins      = 64;
    double k_max_bohr_inv = 0.0;   // 0 -> auto (pi / spacing)
    int    write_every = 1;        // call accumulate every Nth iteration
};

class MomentumDistribution {
public:
    MomentumDistribution(std::string csv_path,
                         int wp_state_index,
                         double l_bohr,
                         MomentumDistributionConfig cfg = {})
        : csv_path_(std::move(csv_path)),
          wp_idx_(wp_state_index),
          l_bohr_(l_bohr),
          cfg_(cfg) {
        namespace fs = std::filesystem;
        if (auto parent = fs::path(csv_path_).parent_path(); !parent.empty())
            fs::create_directories(parent);
        file_.open(csv_path_);
        if (!file_)
            throw std::runtime_error(
                "MomentumDistribution: cannot open '" + csv_path_ + "'");
        // Header
        file_ << "# l_bohr=" << l_bohr_
              << "  n_bins=" << cfg_.n_bins
              << "  wp_idx=" << wp_idx_ << '\n';
        file_ << "step,time_au,k_bohr_inv,n_total,n_wp\n";
    }

    template <typename Viewables>
    void maybe_accumulate(Viewables const& data) {
        if (cfg_.write_every <= 0) return;
        if (data.iter() % cfg_.write_every != 0) return;
        accumulate(data);
    }

    template <typename Viewables>
    void accumulate(Viewables const& data) {
        using namespace inq;
        auto const& electrons = data.electrons();

        // Determine k_max once.
        if (k_max_ <= 0.0) {
            const int nk = electrons.kpin_size();
            if (nk == 0)
                throw std::runtime_error(
                    "MomentumDistribution: no local k-points");
            // Spacing comes from the real-space basis of the first kpin set.
            auto const& rs_basis = electrons.kpin()[0].basis();
            const double dx = rs_basis.rspacing()[0];
            k_max_ = (cfg_.k_max_bohr_inv > 0.0)
                ? cfg_.k_max_bohr_inv
                : M_PI / dx;
            dk_ = k_max_ / cfg_.n_bins;
        }

        std::vector<double> total(cfg_.n_bins, 0.0);
        std::vector<double> wp_only(cfg_.n_bins, 0.0);

        const int nk_local = electrons.kpin_size();
        for (int ik = 0; ik < nk_local; ++ik) {
            auto const& phi = electrons.kpin()[ik];
            auto fphi = inq::operations::transform::to_fourier(phi);

            // Bulk-copy the full hypercubic [ix][iy][iz][ist] view into a
            // contiguous host array — one transfer per snapshot rather than
            // 14M page faults. boost::multi's unary `+` builds a host copy.
            auto host_hc = +fphi.hypercubic();

            auto const& fbasis = fphi.basis();
            auto point_op = fbasis.point_op();
            const auto sizes = fbasis.local_sizes();
            const int nx = sizes[0];
            const int ny = sizes[1];
            const int nz = sizes[2];

            const long state_start = phi.set_part().start();
            const int nst = static_cast<int>(phi.set_part().local_size());

            auto const& occs = electrons.occupations()[ik];
            const double weight = electrons.kpin_weights()[ik];

            // Per-state Parseval norm (constant for all orbitals if they
            // are normalised the same way in real space; computed empirically
            // once on the first snapshot for safety against future changes).
            std::vector<long double> per_state_sum(nst, 0.0L);
            for (int ist = 0; ist < nst; ++ist) {
                long double s = 0.0L;
                for (int ix = 0; ix < nx; ++ix)
                    for (int iy = 0; iy < ny; ++iy)
                        for (int iz = 0; iz < nz; ++iz) {
                            auto z = host_hc[ix][iy][iz][ist];
                            const double v = inq::real(z) * inq::real(z)
                                           + inq::imag(z) * inq::imag(z);
                            s += v;
                        }
                per_state_sum[ist] = s > 0.0L ? s : 1.0L;
            }

            // Bin-accumulate.
            for (int ist = 0; ist < nst; ++ist) {
                const double f_i = (static_cast<std::size_t>(ist) <
                                    occs.size())
                    ? static_cast<double>(occs[ist])
                    : 0.0;
                if (f_i == 0.0) continue;
                const double norm = static_cast<double>(per_state_sum[ist]);
                const long ist_global = state_start + ist;
                const bool is_wp = (ist_global == wp_idx_);
                for (int ix = 0; ix < nx; ++ix)
                    for (int iy = 0; iy < ny; ++iy)
                        for (int iz = 0; iz < nz; ++iz) {
                            const double k_mag = std::sqrt(
                                point_op.g2(ix, iy, iz));
                            const int bin = std::min(
                                cfg_.n_bins - 1,
                                static_cast<int>(std::floor(k_mag / dk_)));
                            if (bin < 0) continue;
                            auto z = host_hc[ix][iy][iz][ist];
                            const double v = inq::real(z) * inq::real(z)
                                           + inq::imag(z) * inq::imag(z);
                            const double contrib = weight * f_i * v / norm;
                            total[bin] += contrib;
                            if (is_wp) wp_only[bin] += contrib;
                        }
            }
        }

        // Write rows.
        const int step = data.iter();
        const double time = data.time();
        file_ << std::setprecision(12);
        for (int b = 0; b < cfg_.n_bins; ++b) {
            const double k_centre = (b + 0.5) * dk_;
            file_ << step << ',' << time << ',' << k_centre << ','
                  << total[b] << ',' << wp_only[b] << '\n';
        }
        file_.flush();
    }

    ~MomentumDistribution() {
        if (file_.is_open()) file_.close();
    }

private:
    std::string csv_path_;
    int wp_idx_;
    double l_bohr_;
    MomentumDistributionConfig cfg_;
    std::ofstream file_;
    double k_max_ = 0.0;
    double dk_    = 0.0;
};

} // namespace inqkit::observables
