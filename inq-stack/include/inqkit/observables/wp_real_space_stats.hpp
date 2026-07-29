/*
 * This file tracks the real-space statistical moments of a single wave-packet
 * orbital at each time step of a TDDFT run. Given the orbital ψ_wp(r, t), it
 * computes the density-weighted integrals (in Bohr):
 *
 *   N(t)         = ∫ |ψ_wp(r, t)|² dV                 (norm; should stay ≈ 1)
 *
 *   <x_d>(t)     = (1/N) ∫ x_d |ψ_wp(r, t)|² dV       (centroid, d = x, y, z)
 *
 *   <x_d²>(t)    = (1/N) ∫ x_d² |ψ_wp(r, t)|² dV      (second moment)
 *
 *   σ²_d(t)      = <x_d²> − <x_d>²                    (spatial variance)
 *
 * where x_d denotes the d-th Cartesian coordinate of a real-space grid point.
 * The centroid tracks the classical trajectory of the packet; the variance
 * σ²_d measures its spread along each axis and grows with dispersion.
 *
 * CSV layout (one row per recorded step):
 * ----------------------------------------
 *   step, time_au,
 *   x_mean, y_mean, z_mean,
 *   x2_mean, y2_mean, z2_mean,
 *   sigma_x2, sigma_y2, sigma_z2,
 *   norm_check
 *
 * norm_check is N(t) and should remain close to 1.0 throughout a unitary
 * propagation. Significant drift indicates a numerical problem upstream.
 *
 * Parallelism
 * -----------
 * The three reductions (N, <x_d>·N, <x_d²>·N) are performed on-device as
 * 3D GPU reductions over the local real-space partition. The resulting 7
 * partial sums are then reduced across MPI ranks via two all_reduce_in_place
 * calls: one over basis().comm() (the FFT-grid decomposition) and one over
 * set_comm() (the state decomposition). Only gamma-point (single k-point)
 * runs are supported, matching the inqkit::WavePacket injector.
 *
 * Known-case validation
 * ---------------------
 * A Gaussian wave packet injected with inqkit::WavePacket at centre (0,0,0)
 * and spread σ (Bohr) writes ψ ∝ exp(−r²/(2σ²)), so the density is a
 * Gaussian with standard deviation σ/√2. The expected observables at t = 0
 * are therefore:
 *
 *   <x_d>   = 0
 *   σ²_d    = σ² / 2
 *
 * For σ = 5 Bohr this gives σ_d = 5/√2 ≈ 3.5355 Bohr. See
 * Tutorial/wp-real-space-stats-test/ for the reference run.
 *
 * Usage
 * -----
 *   WPRealSpaceStatsConfig cfg;
 *   cfg.write_every = 5;           // record every 5th propagation step
 *
 *   WPRealSpaceStats obs("output/wp_rs_stats.csv", wp_state_index, cfg);
 *
 *   // inside the real-time callback:
 *   obs.maybe_accumulate(data);    // no-op on skipped steps
 */
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

// TODO: Is the write_every from the real time session passed on here?
// TODO: Again, we should avoid w, and other such variable names, only ix, iy, iz,
// in some cases, x, y, z (and such like neames wx, wy, wz). However, is it worth
// it to fix this convention throughout?
struct WPRealSpaceStatsConfig {
    int write_every = 1;   // accumulate every Nth iteration; <=0 disables
};

// WP real-space moments — exactly the values accumulate() writes to one CSV row.
struct WPRealSpaceMoments {
    double x = 0, y = 0, z = 0;          // mean position ⟨r⟩ (Bohr, node convention)
    double x2 = 0, y2 = 0, z2 = 0;       // second moments
    double sx2 = 0, sy2 = 0, sz2 = 0;    // variances
    double N = 0;                        // norm ∫|ψ|² dV (≈ 1 for a normalised WP)
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

    // TODO: What does this syntax really do?
    ~WPRealSpaceStats() {
        if (file_.is_open()) file_.close();
    }

    template <typename Viewables>
    void maybe_accumulate(Viewables const& data) {
        if (cfg_.write_every <= 0) return;
        if (data.iter() % cfg_.write_every != 0) return;
        accumulate(data);
    }

    // Compute the WP real-space moments from the current electrons state. Split
    // from accumulate() so it is unit-testable directly (no CSV, no RT
    // Viewables) and reusable. accumulate() = compute() + one CSV row.
    WPRealSpaceMoments compute(inq::systems::electrons const& electrons) const {
        using namespace inq;

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
        // TODO: Need to test that the paralellisation is working as expected. Come up
        // with a test case to test this. 
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
            
            /* Each GPU run runs a highly parallelised loop over every single
            grid space point. */
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
        
        // TODO: MPI reduction here. Is this GPU communication? Need to check? 
        // How are GPU and MPI compatible? Is GPU being paralellised here?
        if (basis.comm().size() > 1)
            basis.comm().all_reduce_in_place_n(host_buf, 7, std::plus<>{});
        if (phi.set_comm().size() > 1)
            phi.set_comm().all_reduce_in_place_n(host_buf, 7, std::plus<>{});

        const double N = host_buf[0];
        if (!(N > 0.0))
            throw std::runtime_error(
                "WPRealSpaceStats: non-positive norm for WP orbital "
                "(state " + std::to_string(wp_idx_) + ").");

        WPRealSpaceMoments m;
        m.N  = N;
        m.x  = host_buf[1] / N;
        m.y  = host_buf[2] / N;
        m.z  = host_buf[3] / N;
        m.x2 = host_buf[4] / N;
        m.y2 = host_buf[5] / N;
        m.z2 = host_buf[6] / N;
        m.sx2 = m.x2 - m.x*m.x;
        m.sy2 = m.y2 - m.y*m.y;
        m.sz2 = m.z2 - m.z*m.z;
        return m;
    }

    // accumulate() = compute() + one CSV row (format unchanged).
    template <typename Viewables>
    void accumulate(Viewables const& data) {
        auto const m = compute(data.electrons());
        const int    step = data.iter();
        const double t_au = data.time();
        file_ << std::setprecision(12);
        file_ << step << ',' << t_au << ','
              << m.x  << ',' << m.y  << ',' << m.z  << ','
              << m.x2 << ',' << m.y2 << ',' << m.z2 << ','
              << m.sx2 << ',' << m.sy2 << ',' << m.sz2 << ','
              << m.N << '\n';
        file_.flush();
    }

private:
    std::string csv_path_;
    int wp_idx_;
    WPRealSpaceStatsConfig cfg_;
    std::ofstream file_;
};

} // namespace inqkit::observables
