/*
 * This file tracks the momentum-space statistical moments of a single
 * wave-packet orbital at each time step of a TDDFT run. The orbital is
 * Fourier-transformed on-the-fly at each recorded step, and the following
 * density-weighted integrals are computed over the resulting k-space grid
 * (all quantities in atomic units):
 *
 *   N(t)         = ∫ |ψ̃_wp(k, t)|² dV_k              (Parseval norm; see note)
 *
 *   <p_d>(t)     = (1/N) ∫ k_d |ψ̃_wp(k, t)|² dV_k    (mean momentum, d = x,y,z)
 *
 *   <p_d²>(t)    = (1/N) ∫ k_d² |ψ̃_wp(k, t)|² dV_k   (second moment)
 *
 *   σ²_pd(t)     = <p_d²> − <p_d>²                    (momentum variance)
 *
 *   E_kin(t)     = ½ (<p_x²> + <p_y²> + <p_z²>)       (kinetic energy, Hartree)
 *
 * where k_d is the d-th Cartesian component of the reciprocal-space grid
 * point. Dividing every moment by N renders the result independent of INQ's
 * FFT prefactor convention; N itself is written to CSV as norm_check.
 *
 * NOTE on the magnitude of N (clarified 2026-07-30): N is NOT ~ 1. Unlike the
 * real-space companion it carries no explicit dV — INQ's FFT prefactor sets its
 * scale, so it is large and grid-dependent (~4.9e7 on a 120x120x200 grid for a
 * normalised packet). That is harmless: every moment is divided by it. What
 * matters is that N stays CONSTANT in time; a drift signals real norm loss from
 * the orbital. Do NOT gate on N ~ 1 — gate on the real-space norm from
 * WPRealSpaceStats, which does carry dV and is ~ 1.
 *
 * Note: unlike the real-space companion (WPRealSpaceStats), no explicit
 * volume element dV appears inside the GPU reductions — the FFT prefactor
 * absorbed into ψ̃ already encodes it. N therefore plays the same
 * normalisation role that the real-space norm plays in WPRealSpaceStats.
 *
 * CSV layout (one row per recorded step):
 * ----------------------------------------
 *   step, time_au,
 *   px_mean, py_mean, pz_mean,
 *   px2_mean, py2_mean, pz2_mean,
 *   sigma_px2, sigma_py2, sigma_pz2,
 *   e_kin_ha, norm_check
 *
 * Parallelism
 * -----------
 * The real-space orbital is transformed to Fourier space with
 * operations::transform::to_fourier() at each recorded step. The three
 * k-space reductions (N, <p_d>·N, <p_d²>·N) are then performed on-device
 * as 3D GPU reductions over the local Fourier-space partition. The resulting
 * 7 partial sums are reduced across MPI ranks via two all_reduce_in_place
 * calls: one over fbasis.comm() (the FFT-grid decomposition) and one over
 * set_comm() (the state decomposition). Ranks that do not hold the WP
 * orbital contribute zero to all sums. Only gamma-point (single k-point)
 * runs are supported, matching the inqkit::WavePacket injector.
 *
 * Known-case validation
 * ---------------------
 * For a Gaussian wave packet injected with inqkit::WavePacket at real-space
 * spread σ_r (Bohr) and initial momentum k₀ — i.e. ψ ∝ exp(−r²/2σ_r²) e^{ik₀·r} —
 * the Fourier transform is also Gaussian,
 *
 *   ψ̃(k) ∝ exp(−σ_r² (k − k₀)² / 2)   ⇒   |ψ̃(k)|² ∝ exp(−σ_r² (k − k₀)²)
 *
 * so the momentum DENSITY has variance σ²_pd = 1/(2 σ_r²). The expected
 * observables at t = 0 are therefore:
 *
 *   <p_d>   = k₀_d
 *   σ²_pd   = 1 / (2 σ_r²)
 *   E_kin   = ½ (|k₀|² + 3 σ_p²) = ½ (|k₀|² + 3/(2σ_r²))
 *
 * CORRECTED 2026-07-30. This block previously claimed σ_p = 1/(2σ_r) and
 * σ²_pd = 1/(4σ_r²). That is WRONG by a factor of 2 in the variance, and it is
 * wrong in a way that is easy to check: the real-space density |ψ|² ∝
 * exp(−r²/σ_r²) has standard deviation σ_d = σ_r/√2 (see the companion
 * WPRealSpaceStats), so the old value would give
 *
 *   σ_d · σ_p = (σ_r/√2)·(1/2σ_r) = 0.354 < ½,
 *
 * i.e. it VIOLATES the Heisenberg bound. The corrected value gives
 * σ_d · σ_p = (σ_r/√2)·(1/(√2 σ_r)) = ½ exactly — a minimum-uncertainty packet,
 * which is what a Gaussian must be. The CODE was always right; only this
 * docstring was wrong. It was copied into a run's t=0 gate on 2026-07-30 and
 * aborted a production job (see docs/handovers/bulk-jellium-ks-stopping.md).
 *
 * Worked example, σ_r = 2 Bohr, k₀ = (0, 0, 2.7111) Bohr⁻¹:
 *   σ²_pd  = 1/(2·4) = 0.125          (σ_p = 0.3536 Bohr⁻¹)
 *   <p_z>  = 2.7111 Bohr⁻¹
 *   E_kin  = ½(7.3499 + 0.375) = 3.8624 Ha = 105.10 eV
 *   E_kin − |k₀|²/2 = 3/(4σ_r²) = 0.1875 Ha = 5.102 eV   (localisation energy)
 * All five verified against a live run (job 32401321, 2026-07-30).
 *
 * For σ_r = 5 Bohr and k₀ = (0, 0, 2.711) Bohr⁻¹: σ_p = 0.1414 Bohr⁻¹,
 * <p_z> = 2.711 Bohr⁻¹, E_kin ≈ 3.705 Ha (≈ 100.8 eV).
 * See Tutorial/wp-momentum-stats-test/ for the reference run.
 *
 * Usage
 * -----
 *   WPMomentumStatsConfig cfg;
 *   cfg.write_every = 5;           // record every 5th propagation step
 *
 *   WPMomentumStats obs("output/wp_mom_stats.csv", wp_state_index, cfg);
 *
 *   // inside the real-time callback:
 *   obs.maybe_accumulate(data);    // no-op on skipped steps
 */
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

// WP momentum moments — exactly the values accumulate() writes to one CSV row.
struct WPMomentumMoments {
    double px = 0, py = 0, pz = 0;       // mean momentum components (Bohr^-1)
    double px2 = 0, py2 = 0, pz2 = 0;    // second moments
    double sx2 = 0, sy2 = 0, sz2 = 0;    // variances
    double ekin = 0;                     // kinetic energy (Ha)
    double N = 0;                        // Parseval norm (norm_check)
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

    // Compute the WP momentum moments from the current electrons state. Split
    // from accumulate() so it is unit-testable directly (no CSV, no RT
    // Viewables) and reusable. accumulate() = compute() + one CSV row.
    WPMomentumMoments compute(inq::systems::electrons const& electrons) const {
        using namespace inq;

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

        WPMomentumMoments m;
        m.N   = N;
        m.px  = host_buf[1] / N;
        m.py  = host_buf[2] / N;
        m.pz  = host_buf[3] / N;
        m.px2 = host_buf[4] / N;
        m.py2 = host_buf[5] / N;
        m.pz2 = host_buf[6] / N;
        m.sx2 = m.px2 - m.px*m.px;
        m.sy2 = m.py2 - m.py*m.py;
        m.sz2 = m.pz2 - m.pz*m.pz;
        m.ekin = 0.5 * (m.px2 + m.py2 + m.pz2);
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
              << m.px  << ',' << m.py  << ',' << m.pz  << ','
              << m.px2 << ',' << m.py2 << ',' << m.pz2 << ','
              << m.sx2 << ',' << m.sy2 << ',' << m.sz2 << ','
              << m.ekin << ',' << m.N << '\n';
        file_.flush();
    }

private:
    std::string csv_path_;
    int wp_idx_;
    WPMomentumStatsConfig cfg_;
    std::ofstream file_;
};

} // namespace inqkit::observables
