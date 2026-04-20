#pragma once
// ============================================================================
// inqkit::WavePacket
//
// Builder-pattern class for constructing and injecting a Gaussian wavepacket
// into the last extra-state slot of an INQ electrons object.
//
// Usage:
//   auto report = inqkit::WavePacket{}
//       .center(cx_bohr, cy_bohr, cz_bohr)
//       .sigma(sigma_bohr)
//       .k0(kx, ky, kz)
//       .orthogonalise_against_occupied(electrons)
//       .inject_into_last_extra_state(electrons, 1.0);
//
// GPU injection pattern adapted from:
//   ResearchProject/systems/coronene/04_leed_simulation/runs/run_005/utils.hpp
//   Tsubonoya, Hu, Watanabe, PRB 90, 035416 (2014)
//
// Single-rank only (multi-rank basis/set partitioning is not supported).
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/wavepacket/injection_report.hpp>

#include <cmath>
#include <stdexcept>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#ifndef INQKIT_GPU_SYNC
#define INQKIT_GPU_SYNC() cudaDeviceSynchronize()
#endif
#else
#ifndef INQKIT_GPU_SYNC
#define INQKIT_GPU_SYNC() ((void)0)
#endif
#endif

namespace inqkit {

class WavePacket {
    double cx_ = 0, cy_ = 0, cz_ = 0;
    double sigma_ = 1.0;
    double kx_ = 0, ky_ = 0, kz_ = 0;
    bool   do_ortho_ = false;
    double ortho_tol_ = 1e-6;

public:
    WavePacket& center(double x_bohr, double y_bohr, double z_bohr) {
        cx_ = x_bohr; cy_ = y_bohr; cz_ = z_bohr;
        return *this;
    }

    WavePacket& sigma(double sigma_bohr) {
        sigma_ = sigma_bohr;
        return *this;
    }

    WavePacket& k0(double kx, double ky, double kz) {
        kx_ = kx; ky_ = ky; kz_ = kz;
        return *this;
    }

    // Mark the WP for orthogonalisation against all occupied states during injection.
    // Actual projection is performed inside inject_into_last_extra_state using
    // GPU kernels (modified Gram-Schmidt). Passing electrons here is required by the
    // API for consistency with future pre-computation of overlaps; currently unused.
    WavePacket& orthogonalise_against_occupied(
        inq::systems::electrons const& /*electrons*/,
        double tolerance = 1e-6)
    {
        do_ortho_  = true;
        ortho_tol_ = tolerance;
        return *this;
    }

    InjectionReport inject_into_last_extra_state(
        inq::systems::electrons& electrons,
        double occupation = 1.0) const;
};

// ────────────────────────────────────────────────────────────────────────────

inline InjectionReport WavePacket::inject_into_last_extra_state(
    inq::systems::electrons& electrons,
    double occupation) const
{
    using complex = inq::complex;

    if (electrons.kpin().size() != 1) {
        throw std::runtime_error(
            "inqkit::WavePacket: only single-kpoint (gamma-only) runs are supported.");
    }

    auto& phi   = electrons.kpin()[0];
    auto& basis = phi.basis();

    if (phi.basis().comm().size() != 1 || phi.set_comm().size() != 1) {
        throw std::runtime_error(
            "inqkit::WavePacket: multi-rank basis/set partitioning is not supported.");
    }

    int  ist_wp = phi.set_part().local_size() - 1;
    int  n_pts  = phi.basis().local_size();
    double dV   = basis.volume_element();

    InjectionReport report;
    report.kpoint_index = 0;
    report.state_index  = ist_wp;

    // ── 1. Norm of existing last slot (before injection) ───────────────────────
    INQKIT_GPU_SYNC();
    {
        auto mat_ = begin(phi.matrix());
        auto res = gpu::run(1, gpu::reduce(n_pts), 0.0,
            [dV, mat_, ist_wp_=ist_wp] GPU_LAMBDA (auto, auto ip) {
                auto v = mat_[ip][ist_wp_];
                return dV * (inq::real(v)*inq::real(v) + inq::imag(v)*inq::imag(v));
            });
        INQKIT_GPU_SYNC();
        report.norm_before = std::sqrt(res[0]);
    }

    // ── 2. Inject raw Gaussian wavepacket ──────────────────────────────────────
    // psi_wp(r) = (pi sigma^2)^{-3/4} exp(-|r-b|^2 / (2 sigma^2)) exp(i k.r)
    double norm_fac = std::pow(M_PI * sigma_ * sigma_, -0.75);
    double sig  = sigma_;
    double bx   = cx_,  by  = cy_,  bz  = cz_;
    double kxv  = kx_,  kyv = ky_,  kzv = kz_;
    double dx_sp = basis.rspacing()[0];
    double dy_sp = basis.rspacing()[1];
    double dz_sp = basis.rspacing()[2];
    int    x0    = basis.cubic_part(0).start();
    int    y0    = basis.cubic_part(1).start();
    int    z0    = basis.cubic_part(2).start();
    int    ist_w = ist_wp;

    auto phicub_ = begin(phi.hypercubic());

    gpu::run(basis.local_sizes()[2],
             basis.local_sizes()[1],
             basis.local_sizes()[0],
        [=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
            double rx  = (ix + x0) * dx_sp;
            double ry  = (iy + y0) * dy_sp;
            double rz  = (iz + z0) * dz_sp;
            double dx_ = rx - bx,  dy_ = ry - by,  dz_ = rz - bz;
            double r2  = dx_*dx_ + dy_*dy_ + dz_*dz_;
            double amp = norm_fac * exp(-r2 / (2.0 * sig * sig));
            double ph  = kxv*rx + kyv*ry + kzv*rz;
            phicub_[ix][iy][iz][ist_w] = complex(amp * cos(ph), amp * sin(ph));
        });
    INQKIT_GPU_SYNC();

    // ── 3. Orthogonalise against occupied states (modified Gram-Schmidt on GPU) ─
    if (do_ortho_) {
        double max_ov = 0.0;
        auto   mat_   = begin(phi.matrix());

        for (int i = 0; i < ist_wp; ++i) {
            // Real part of <psi_i | psi_wp>
            auto res_re = gpu::run(1, gpu::reduce(n_pts), 0.0,
                [dV, mat_, i_=i, ist_=ist_wp] GPU_LAMBDA (auto, auto ip) {
                    auto vi = mat_[ip][i_];
                    auto vw = mat_[ip][ist_];
                    return dV * (inq::real(vi)*inq::real(vw) + inq::imag(vi)*inq::imag(vw));
                });
            INQKIT_GPU_SYNC();
            double ov_re = res_re[0];

            // Imaginary part of <psi_i | psi_wp>
            auto res_im = gpu::run(1, gpu::reduce(n_pts), 0.0,
                [dV, mat_, i_=i, ist_=ist_wp] GPU_LAMBDA (auto, auto ip) {
                    auto vi = mat_[ip][i_];
                    auto vw = mat_[ip][ist_];
                    return dV * (inq::real(vi)*inq::imag(vw) - inq::imag(vi)*inq::real(vw));
                });
            INQKIT_GPU_SYNC();
            double ov_im = res_im[0];

            max_ov = std::max(max_ov, std::sqrt(ov_re*ov_re + ov_im*ov_im));

            // Subtract projection: psi_wp -= (ov_re + i*ov_im) * psi_i
            double re_ = ov_re, im_ = ov_im;
            int    i_  = i;
            gpu::run(basis.local_sizes()[2],
                     basis.local_sizes()[1],
                     basis.local_sizes()[0],
                [=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto vi  = phicub_[ix][iy][iz][i_];
                    auto vw  = phicub_[ix][iy][iz][ist_w];
                    // (ov_re + i*ov_im) * vi
                    double sub_re = re_ * inq::real(vi) - im_ * inq::imag(vi);
                    double sub_im = re_ * inq::imag(vi) + im_ * inq::real(vi);
                    phicub_[ix][iy][iz][ist_w] = complex(
                        inq::real(vw) - sub_re,
                        inq::imag(vw) - sub_im);
                });
            INQKIT_GPU_SYNC();
        }

        report.max_overlap    = max_ov;
        report.orthogonalised = true;

        // Renormalise after projection
        {
            auto res = gpu::run(1, gpu::reduce(n_pts), 0.0,
                [dV, mat_, ist_=ist_wp] GPU_LAMBDA (auto, auto ip) {
                    auto v = mat_[ip][ist_];
                    return dV * (inq::real(v)*inq::real(v) + inq::imag(v)*inq::imag(v));
                });
            INQKIT_GPU_SYNC();
            double scale = 1.0 / std::sqrt(res[0]);

            gpu::run(basis.local_sizes()[2],
                     basis.local_sizes()[1],
                     basis.local_sizes()[0],
                [=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
                    auto v = phicub_[ix][iy][iz][ist_w];
                    phicub_[ix][iy][iz][ist_w] = complex(
                        inq::real(v) * scale, inq::imag(v) * scale);
                });
            INQKIT_GPU_SYNC();
        }

        report.passed_tolerance = (max_ov < ortho_tol_ * 10.0);
    } else {
        report.passed_tolerance = true;
    }

    // ── 4. Norm after injection (and orthogonalisation if applied) ─────────────
    {
        auto mat_ = begin(phi.matrix());
        auto res = gpu::run(1, gpu::reduce(n_pts), 0.0,
            [dV, mat_, ist_=ist_wp] GPU_LAMBDA (auto, auto ip) {
                auto v = mat_[ip][ist_];
                return dV * (inq::real(v)*inq::real(v) + inq::imag(v)*inq::imag(v));
            });
        INQKIT_GPU_SYNC();
        report.norm_after = std::sqrt(res[0]);
    }

    // ── 5. Set occupation ──────────────────────────────────────────────────────
    electrons.occupations()[0][ist_wp] = occupation;

    return report;
}

} // namespace inqkit
