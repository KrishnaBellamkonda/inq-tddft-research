// ============================================================================
// inqkit/observables/wp_real_space_stats.hpp
//
// Per-step real-space statistics of a single wave-packet orbital in a
// real-time TDDFT run. Computes (in Bohr):
//
//     <x_d>(t)   = integral x_d |psi_wp(r,t)|^2 dV / N
//     <x_d^2>(t) = integral x_d^2 |psi_wp(r,t)|^2 dV / N
//     sigma_r_d^2(t) = <x_d^2> - <x_d>^2
//
// where N = integral |psi_wp(r,t)|^2 dV is the real-space norm (written to
// CSV as `norm_check`; should remain ~1 throughout a unitary propagation).
//
// Design (mirrors inq/src/observables/dipole.hpp):
//   - On-device 3D GPU reduction over the real-space grid points.
//   - Host-side MPI all_reduce_in_place_n on the small fixed-size summary
//     across basis().comm() (FFT-grid decomposition) and set_comm() (state
//     decomposition).
//
// CSV layout (one row per recorded step):
//   step,time_au,x_mean,y_mean,z_mean,
//                x2_mean,y2_mean,z2_mean,
//                sigma_x2,sigma_y2,sigma_z2,
//                norm_check
//
// Known-case validation: Tutorial/wp-real-space-stats-test/. A Gaussian WP
// injected with inqkit::WavePacket sigma = 5 Bohr and centre = (0, 0, 0)
// must give <x_d> ~= 0 and sigma_r_d = 5 / sqrt(2) ~= 3.5355 Bohr (the
// injector writes psi ~ exp(-r^2 / (2 sigma^2)), so the density is a
// Gaussian with std-dev sigma / sqrt(2) — see
// wp_momentum_stats.hpp header for the convention discussion).
//
// Single k-point only (matches inqkit::WavePacket).
// ============================================================================
#pragma once

#include <inq/inq.hpp>
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

struct WPRealSpaceStatsConfig {
    int write_every = 1;   // accumulate every Nth iteration; <=0 disables
};

class WPRealSpaceStats {
public:
    WPRealSpaceStats(std::string csv_path,
                     int wp_state_index,
                     WPRealSpaceStatsConfig cfg = {})
        : csv_path_(std::move(csv_path)),
          wp_idx_(wp_state_index),
          cfg_(cfg) {
        namespace fs = std::filesystem;
        if (auto parent = fs::path(csv_path_).parent_path(); !parent.empty())
            fs::create_directories(parent);
        file_.open(csv_path_);
        if (!file_)
            throw std::runtime_error(
                "WPRealSpaceStats: cannot open '" + csv_path_ + "'");
        file_ << "# wp_state_index=" << wp_idx_
              << "  write_every=" << cfg_.write_every << '\n';
        file_ << "step,time_au,"
                 "x_mean,y_mean,z_mean,"
                 "x2_mean,y2_mean,z2_mean,"
                 "sigma_x2,sigma_y2,sigma_z2,"
                 "norm_check\n";
    }

    ~WPRealSpaceStats() {
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
                "WPRealSpaceStats: only single-kpoint (gamma-only) "
                "runs are supported.");

        auto const& phi   = electrons.kpin()[0];
        auto const& basis = phi.basis();
        const double dV   = basis.volume_element();

        // Is the requested WP state local to this rank's set partition?
        const long st_start = phi.set_part().start();
        const long st_size  = phi.set_part().local_size();
        const bool wp_local = (wp_idx_ >= st_start &&
                               wp_idx_ <  st_start + st_size);
        const int  ist_l    = wp_local
            ? static_cast<int>(wp_idx_ - st_start)
            : 0;

        // 7 sums in real-space units: [N, <x>*N, <y>*N, <z>*N,
        //                              <x^2>*N, <y^2>*N, <z^2>*N]
        // All carry the volume element dV; positions are Cartesian Bohr.
        double sum_n  = 0.0;
        inq::vector3<double> sum_r {0.0, 0.0, 0.0};
        inq::vector3<double> sum_r2{0.0, 0.0, 0.0};

        if (wp_local) {
            auto const sizes = basis.local_sizes();
            auto phic      = begin(phi.hypercubic());
            auto point_op  = basis.point_op();

            // --- N = sum |psi(r)|^2 dV ---------------------------------------
            sum_n = gpu::run(
                gpu::reduce(sizes[2]),
                gpu::reduce(sizes[1]),
                gpu::reduce(sizes[0]),
                0.0,
                [phic, ist_l, dV] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto v = phic[ix][iy][iz][ist_l];
                    return dV * (inq::real(v) * inq::real(v)
                               + inq::imag(v) * inq::imag(v));
                });

            // --- <x_d> * N = sum x_d |psi(r)|^2 dV ---------------------------
            sum_r = gpu::run(
                gpu::reduce(sizes[2]),
                gpu::reduce(sizes[1]),
                gpu::reduce(sizes[0]),
                inq::vector3<double>{0.0, 0.0, 0.0},
                [phic, ist_l, dV, point_op] GPU_LAMBDA (
                    auto iz, auto iy, auto ix) {
                    auto v = phic[ix][iy][iz][ist_l];
                    double w = dV * (inq::real(v) * inq::real(v)
                                   + inq::imag(v) * inq::imag(v));
                    auto r = point_op.rvector_cartesian(ix, iy, iz);
                    return inq::vector3<double>{r[0]*w, r[1]*w, r[2]*w};
                });

            // --- <x_d^2> * N = sum x_d^2 |psi(r)|^2 dV -----------------------
            sum_r2 = gpu::run(
                gpu::reduce(sizes[2]),
                gpu::reduce(sizes[1]),
                gpu::reduce(sizes[0]),
                inq::vector3<double>{0.0, 0.0, 0.0},
                [phic, ist_l, dV, point_op] GPU_LAMBDA (
                    auto iz, auto iy, auto ix) {
                    auto v = phic[ix][iy][iz][ist_l];
                    double w = dV * (inq::real(v) * inq::real(v)
                                   + inq::imag(v) * inq::imag(v));
                    auto r = point_op.rvector_cartesian(ix, iy, iz);
                    return inq::vector3<double>{
                        r[0]*r[0]*w, r[1]*r[1]*w, r[2]*r[2]*w};
                });
        }

        double host_buf[7] = {sum_n,
                              sum_r[0],  sum_r[1],  sum_r[2],
                              sum_r2[0], sum_r2[1], sum_r2[2]};

        if (basis.comm().size() > 1)
            basis.comm().all_reduce_in_place_n(host_buf, 7, std::plus<>{});
        if (phi.set_comm().size() > 1)
            phi.set_comm().all_reduce_in_place_n(host_buf, 7, std::plus<>{});

        const double N = host_buf[0];
        if (!(N > 0.0))
            throw std::runtime_error(
                "WPRealSpaceStats: non-positive norm for WP orbital "
                "(state " + std::to_string(wp_idx_) + ").");

        const double x  = host_buf[1] / N;
        const double y  = host_buf[2] / N;
        const double z  = host_buf[3] / N;
        const double x2 = host_buf[4] / N;
        const double y2 = host_buf[5] / N;
        const double z2 = host_buf[6] / N;
        const double sx2 = x2 - x*x;
        const double sy2 = y2 - y*y;
        const double sz2 = z2 - z*z;

        const int    step = data.iter();
        const double t_au = data.time();
        file_ << std::setprecision(12);
        file_ << step << ',' << t_au << ','
              << x  << ',' << y  << ',' << z  << ','
              << x2 << ',' << y2 << ',' << z2 << ','
              << sx2 << ',' << sy2 << ',' << sz2 << ','
              << N << '\n';
        file_.flush();
    }

private:
    std::string csv_path_;
    int wp_idx_;
    WPRealSpaceStatsConfig cfg_;
    std::ofstream file_;
};

} // namespace inqkit::observables
