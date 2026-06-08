// ============================================================================
// inqkit/observables/wp_momentum_stats.hpp
//
// Per-step momentum statistics of a single wave-packet orbital in a real-time
// TDDFT run. Computes (in atomic units):
//
//     <p_d>(t)   = integral k_d |phi-tilde_wp(k,t)|^2 dV_k / N
//     <p_d^2>(t) = integral k_d^2 |phi-tilde_wp(k,t)|^2 dV_k / N
//     sigma_p_d^2(t) = <p_d^2> - <p_d>^2
//     E_kin(t)   = 0.5 (<p_x^2> + <p_y^2> + <p_z^2>)            (Hartree)
//
// where N = integral |phi-tilde_wp(k,t)|^2 dV_k is the Parseval normalisation
// of the WP orbital in Fourier space (kept as a diagnostic, written to CSV as
// `norm_check`). Dividing every moment by N makes the result independent of
// INQ's FFT prefactor convention.
//
// Design (mirrors inq/src/observables/dipole.hpp):
//   - On-device 3D GPU reduction over the Fourier-space grid points.
//   - Host-side MPI all_reduce_in_place_n on the small fixed-size summary,
//     across both the basis communicator (FFT-grid decomposition) and the
//     set/state communicator (state decomposition — only the rank holding the
//     WP contributes a non-zero local sum). This is the host-after-reduction
//     pattern the campaign requires; the legacy MomentumDistribution observable
//     does not all-reduce and is single-rank-safe only.
//
// CSV layout (one row per recorded step):
//   step,time_au,px_mean,py_mean,pz_mean,
//                px2_mean,py2_mean,pz2_mean,
//                sigma_px2,sigma_py2,sigma_pz2,
//                e_kin_ha,norm_check
//
// Known-case validation: Tutorial/wp-momentum-stats-test/. A Gaussian WP
// injected with sigma_r = 5 Bohr, k0 = (0, 0, 2.711) Bohr^-1 must give
// sigma_p_d = 1/(2 sigma_r) = 0.1 Bohr^-1 in every cartesian direction,
// <p_z> = 2.711 Bohr^-1, and E_kin = 0.5 (k0^2 + 3 sigma_p^2) ~= 3.69 Ha
// (~= 100.4 eV).
//
// Single k-point only (matches inqkit::WavePacket).
// ============================================================================
#pragma once

#include <inq/inq.hpp>
#include <operations/transform.hpp>
#include <basis/fourier_space.hpp>
#include <math/vector3.hpp>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

namespace inqkit::observables {

struct WPMomentumStatsConfig {
    int write_every = 1;   // accumulate every Nth iteration; <=0 disables
};

class WPMomentumStats {
public:
    WPMomentumStats(std::string csv_path,
                    int wp_state_index,
                    WPMomentumStatsConfig cfg = {})
        : csv_path_(std::move(csv_path)),
          wp_idx_(wp_state_index),
          cfg_(cfg) {
        namespace fs = std::filesystem;
        if (auto parent = fs::path(csv_path_).parent_path(); !parent.empty())
            fs::create_directories(parent);
        file_.open(csv_path_);
        if (!file_)
            throw std::runtime_error(
                "WPMomentumStats: cannot open '" + csv_path_ + "'");
        file_ << "# wp_state_index=" << wp_idx_
              << "  write_every=" << cfg_.write_every << '\n';
        file_ << "step,time_au,"
                 "px_mean,py_mean,pz_mean,"
                 "px2_mean,py2_mean,pz2_mean,"
                 "sigma_px2,sigma_py2,sigma_pz2,"
                 "e_kin_ha,norm_check\n";
    }

    ~WPMomentumStats() {
        if (file_.is_open()) file_.close();
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

        if (electrons.kpin_size() != 1)
            throw std::runtime_error(
                "WPMomentumStats: only single-kpoint (gamma-only) "
                "runs are supported.");

        auto fphi = inq::operations::transform::to_fourier(electrons.kpin()[0]);
        auto const& fbasis = fphi.basis();

        // Is the requested WP state local to this rank's set partition?
        const long st_start = fphi.set_part().start();
        const long st_size  = fphi.set_part().local_size();
        const bool wp_local = (wp_idx_ >= st_start &&
                               wp_idx_ <  st_start + st_size);
        const int  ist_l    = wp_local
            ? static_cast<int>(wp_idx_ - st_start)
            : 0;

        // 7 scalar sums: [N, <px>*N, <py>*N, <pz>*N, <px^2>*N, <py^2>*N, <pz^2>*N]
        double sum_n  = 0.0;
        inq::vector3<double> sum_p {0.0, 0.0, 0.0};
        inq::vector3<double> sum_p2{0.0, 0.0, 0.0};

        if (wp_local) {
            auto const sizes = fbasis.local_sizes();
            auto fhc       = begin(fphi.hypercubic());
            auto point_op  = fbasis.point_op();

            // --- N = sum |phi-tilde(k)|^2 -------------------------------------
            sum_n = gpu::run(
                gpu::reduce(sizes[2]),
                gpu::reduce(sizes[1]),
                gpu::reduce(sizes[0]),
                0.0,
                [fhc, ist_l] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto v = fhc[ix][iy][iz][ist_l];
                    return inq::real(v) * inq::real(v)
                         + inq::imag(v) * inq::imag(v);
                });

            // --- <p_d> * N = sum k_d |phi-tilde(k)|^2 -------------------------
            sum_p = gpu::run(
                gpu::reduce(sizes[2]),
                gpu::reduce(sizes[1]),
                gpu::reduce(sizes[0]),
                inq::vector3<double>{0.0, 0.0, 0.0},
                [fhc, ist_l, point_op] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto v = fhc[ix][iy][iz][ist_l];
                    double w = inq::real(v) * inq::real(v)
                             + inq::imag(v) * inq::imag(v);
                    auto k = point_op.gvector_cartesian(ix, iy, iz);
                    return inq::vector3<double>{k[0]*w, k[1]*w, k[2]*w};
                });

            // --- <p_d^2> * N = sum k_d^2 |phi-tilde(k)|^2 ---------------------
            sum_p2 = gpu::run(
                gpu::reduce(sizes[2]),
                gpu::reduce(sizes[1]),
                gpu::reduce(sizes[0]),
                inq::vector3<double>{0.0, 0.0, 0.0},
                [fhc, ist_l, point_op] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto v = fhc[ix][iy][iz][ist_l];
                    double w = inq::real(v) * inq::real(v)
                             + inq::imag(v) * inq::imag(v);
                    auto k = point_op.gvector_cartesian(ix, iy, iz);
                    return inq::vector3<double>{
                        k[0]*k[0]*w, k[1]*k[1]*w, k[2]*k[2]*w};
                });
        }

        // Host-side MPI reductions: basis-decomposition and set-decomposition.
        // The WP only contributes on the rank holding its state; all_reduce
        // sums the per-rank zeros with the one non-zero contribution.
        double host_buf[7] = {sum_n,
                              sum_p[0],  sum_p[1],  sum_p[2],
                              sum_p2[0], sum_p2[1], sum_p2[2]};

        if (fbasis.comm().size() > 1)
            fbasis.comm().all_reduce_in_place_n(host_buf, 7, std::plus<>{});
        if (fphi.set_comm().size() > 1)
            fphi.set_comm().all_reduce_in_place_n(host_buf, 7, std::plus<>{});

        const double N = host_buf[0];
        if (!(N > 0.0))
            throw std::runtime_error(
                "WPMomentumStats: non-positive Parseval norm for WP orbital "
                "(state " + std::to_string(wp_idx_) + ").");

        const double px  = host_buf[1] / N;
        const double py  = host_buf[2] / N;
        const double pz  = host_buf[3] / N;
        const double px2 = host_buf[4] / N;
        const double py2 = host_buf[5] / N;
        const double pz2 = host_buf[6] / N;
        const double sx2 = px2 - px*px;
        const double sy2 = py2 - py*py;
        const double sz2 = pz2 - pz*pz;
        const double ekin = 0.5 * (px2 + py2 + pz2);

        const int    step = data.iter();
        const double t_au = data.time();
        file_ << std::setprecision(12);
        file_ << step << ',' << t_au << ','
              << px  << ',' << py  << ',' << pz  << ','
              << px2 << ',' << py2 << ',' << pz2 << ','
              << sx2 << ',' << sy2 << ',' << sz2 << ','
              << ekin << ',' << N << '\n';
        file_.flush();
    }

private:
    std::string csv_path_;
    int wp_idx_;
    WPMomentumStatsConfig cfg_;
    std::ofstream file_;
};

} // namespace inqkit::observables
