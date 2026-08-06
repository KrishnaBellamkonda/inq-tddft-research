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
 * PERIODIC-AWARE MOMENTS (circular statistics)
 * --------------------------------------------
 * The moments above use the raw Cartesian coordinate x_d, which is DISCONTINUOUS
 * across a periodic boundary. Once a packet straddles a cell face, ⟨x_d⟩ and
 * σ²_d are meaningless — the estimator averages the two halves of a split packet
 * and returns a point somewhere in the middle of the cell, which looks plausible
 * and is completely wrong.
 *
 * The periodic-safe replacement is the circular (phase) estimator. Working in
 * fractional coordinates f_d = (b_d · r) / 2π (b_d = reciprocal lattice vector,
 * so f_d advances by exactly 1 across the cell), define the resultant
 *
 *   C_d = (1/N) ∫ cos(2π f_d) |ψ|² dV ,  S_d = (1/N) ∫ sin(2π f_d) |ψ|² dV
 *
 * and then
 *
 *   f̄_d   = atan2(S_d, C_d) / 2π              (mean fractional position)
 *   ⟨r⟩_circ = Σ_d f̄_d a_d                     (mean position, Cartesian Bohr)
 *   R_d    = sqrt(C_d² + S_d²)                (resultant length, 0 … 1)
 *   σ_d,circ = (L_d/2π) · sqrt(−2 ln R_d)     (circular spread, Bohr)
 *
 * ⟨r⟩_circ is exact in a periodic cell and agrees with the naive ⟨r⟩ to
 * machine precision for a packet well away from the faces, so the two are a
 * mutual cross-check. R_d → 1 for a tightly localised packet and R_d → 0 for
 * one spread uniformly over the cell; σ_d,circ inverts the wrapped-Gaussian
 * relation R = exp(−σ_θ²/2) with σ_θ = 2π σ / L, so it recovers the ordinary
 * standard deviation when σ ≪ L and saturates gracefully when the packet fills
 * the cell (where the naive σ_d would keep growing without bound).
 *
 * L_d here is the length of the d-th lattice vector. For a non-orthogonal cell
 * σ_d,circ is measured along that lattice direction, not along a Cartesian axis.
 *
 * Reference: R. Resta, "Quantum-Mechanical Position Operator in Extended
 * Systems", Phys. Rev. Lett. 80, 1800 (1998) — the same phase construction,
 * introduced there for the many-body polarisation. The single-orbital form used
 * here is the standard directional-statistics circular mean.
 *
 * CSV layout (one row per recorded step):
 * ----------------------------------------
 *   step, time_au,
 *   x_mean, y_mean, z_mean,
 *   x2_mean, y2_mean, z2_mean,
 *   sigma_x2, sigma_y2, sigma_z2,
 *   norm_check,
 *   x_mean_circ, y_mean_circ, z_mean_circ,
 *   R_x, R_y, R_z,
 *   sigma_x_circ, sigma_y_circ, sigma_z_circ
 *
 * The nine circular columns are APPENDED after norm_check so that the legacy
 * schema (through norm_check) is byte-compatible with every pre-2026-07-30 run
 * and its post-processing.
 *
 * f̄_d is returned in (−0.5, +0.5], matching the node convention of
 * rvector_cartesian, so ⟨r⟩_circ lands in the same cell-centred window as the
 * naive ⟨r⟩. A trajectory that crosses a face therefore shows a single jump of
 * one cell length in ⟨r⟩_circ, which post-processing UNWRAPS (np.unwrap on
 * 2π f̄) to recover a continuous path — as opposed to the naive centroid, whose
 * boundary behaviour is not a jump but a smooth slide to the wrong answer and
 * cannot be repaired after the fact.
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

    // Periodic-aware (circular) moments — see the header comment. These are the
    // ONLY position estimates valid once the packet straddles a cell face.
    double xc = 0, yc = 0, zc = 0;       // ⟨r⟩_circ (Bohr, node convention)
    double Rx = 0, Ry = 0, Rz = 0;       // resultant length, 0 (spread) … 1 (sharp)
    double sxc = 0, syc = 0, szc = 0;    // circular spread along each lattice dir (Bohr)
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
                 "norm_check,"
                 "x_mean_circ,y_mean_circ,z_mean_circ,"
                 "R_x,R_y,R_z,"
                 "sigma_x_circ,sigma_y_circ,sigma_z_circ\n";
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
        // Circular (periodic-aware) accumulators: Σ cos(2π f_d)|ψ|² dV and
        // Σ sin(2π f_d)|ψ|² dV, with f_d the fractional coordinate.
        inq::vector3<double> sum_cos{0.0, 0.0, 0.0};
        inq::vector3<double> sum_sin{0.0, 0.0, 0.0};

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

            // --- circular moments: sum cos(2*pi*f_d) |psi|^2 dV --------------
            // f_d = fractional coordinate = (b_d . r) / 2pi, obtained from the
            // cell itself so this is correct for ANY lattice, not just
            // orthorhombic. cos/sin are evaluated on the phase 2*pi*f_d.
            sum_cos = gpu::run(
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
                    auto f = point_op.cell().to_contravariant(r);
                    return inq::vector3<double>{
                        w * cos(2.0*M_PI*f[0]),
                        w * cos(2.0*M_PI*f[1]),
                        w * cos(2.0*M_PI*f[2])};
                });

            // --- circular moments: sum sin(2*pi*f_d) |psi|^2 dV --------------
            sum_sin = gpu::run(
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
                    auto f = point_op.cell().to_contravariant(r);
                    return inq::vector3<double>{
                        w * sin(2.0*M_PI*f[0]),
                        w * sin(2.0*M_PI*f[1]),
                        w * sin(2.0*M_PI*f[2])};
                });
        }

        double host_buf[13] = {sum_n,
                               sum_r[0],   sum_r[1],   sum_r[2],
                               sum_r2[0],  sum_r2[1],  sum_r2[2],
                               sum_cos[0], sum_cos[1], sum_cos[2],
                               sum_sin[0], sum_sin[1], sum_sin[2]};

        // TODO: MPI reduction here. Is this GPU communication? Need to check?
        // How are GPU and MPI compatible? Is GPU being paralellised here?
        if (basis.comm().size() > 1)
            basis.comm().all_reduce_in_place_n(host_buf, 13, std::plus<>{});
        if (phi.set_comm().size() > 1)
            phi.set_comm().all_reduce_in_place_n(host_buf, 13, std::plus<>{});

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

        // ----- circular (periodic-aware) position and spread ----------------
        // Mean fractional coordinate per lattice direction, then back to
        // Cartesian through the lattice vectors.
        const inq::vector3<double> C{host_buf[7] / N, host_buf[8] / N,
                                     host_buf[9] / N};
        const inq::vector3<double> S{host_buf[10] / N, host_buf[11] / N,
                                     host_buf[12] / N};

        auto const& cell = basis.cell();
        inq::vector3<double, inq::contravariant> fbar{0.0, 0.0, 0.0};
        double Rd[3];
        for (int d = 0; d < 3; ++d) {
            fbar[d] = std::atan2(S[d], C[d]) / (2.0 * M_PI);   // in (-1/2, 1/2]
            Rd[d]   = std::sqrt(C[d]*C[d] + S[d]*S[d]);        // in [0, 1]
        }
        auto rc = cell.to_cartesian(fbar);
        m.xc = rc[0];  m.yc = rc[1];  m.zc = rc[2];
        m.Rx = Rd[0];  m.Ry = Rd[1];  m.Rz = Rd[2];

        // Circular spread: invert the wrapped-Gaussian relation R = exp(-s^2/2)
        // with s = 2*pi*sigma/L, i.e. sigma = (L/2pi) sqrt(-2 ln R). R -> 0 means
        // the packet fills the cell and sigma is unbounded; clamp to the uniform
        // -distribution value so the column stays finite and monotone.
        constexpr double R_FLOOR = 1.0e-12;
        for (int d = 0; d < 3; ++d) {
            const double L_d = cell.lattice(d).length();
            const double R   = Rd[d] > R_FLOOR ? Rd[d] : R_FLOOR;
            const double s   = (L_d / (2.0*M_PI)) * std::sqrt(-2.0*std::log(R));
            if (d == 0) m.sxc = s;
            else if (d == 1) m.syc = s;
            else m.szc = s;
        }
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
              << m.N << ','
              << m.xc  << ',' << m.yc  << ',' << m.zc  << ','
              << m.Rx  << ',' << m.Ry  << ',' << m.Rz  << ','
              << m.sxc << ',' << m.syc << ',' << m.szc << '\n';
        file_.flush();
    }

private:
    std::string csv_path_;
    int wp_idx_;
    WPRealSpaceStatsConfig cfg_;
    std::ofstream file_;
};

} // namespace inqkit::observables
