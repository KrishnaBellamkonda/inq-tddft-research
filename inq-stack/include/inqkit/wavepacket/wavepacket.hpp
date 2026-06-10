/*
 * This file provides a builder-pattern class for constructing a Gaussian wave
 * packet and injecting it into the last extra-state slot of an INQ electrons
 * object. The injected orbital takes the form:
 *
 *   ψ_wp(r) = (π σ²)^{-3/4} exp(−|r − b|² / (2σ²)) exp(i k₀·r)
 *
 * where b is the packet centre (Bohr), σ is the real-space spread (Bohr), and
 * k₀ is the initial momentum (Bohr⁻¹). The prefactor (π σ²)^{-3/4} ensures
 * ∫ |ψ_wp|² dV = 1. The resulting density is a Gaussian with standard
 * deviation σ/√2 along each axis — see WPRealSpaceStats for the
 * corresponding validation.
 *
 * Injection sequence
 * ------------------
 * inject_into_last_extra_state() performs the following steps in order:
 *
 *   1. Norm before    Records ‖ψ_last‖ of the slot before overwriting, so
 *                     the caller can verify the slot was previously empty.
 *
 *   2. Raw injection  Writes ψ_wp(r) into the last extra-state slot on the
 *                     GPU, using basis.point_op().rvector_cartesian() for
 *                     the physical coordinates. This is essential for
 *                     orthorhombic cells: the symmetric-range convention
 *                     maps grid index ig to (ig ≤ N/2 ? ig : ig − N)·dr,
 *                     so naively computing ig·dr would silently corrupt the
 *                     exp(i k·r) phase across half the cell.
 *
 *   3. Orthogonalisation (optional)
 *                     If orthogonalise_against_occupied() was called,
 *                     modified Gram-Schmidt is applied against all occupied
 *                     KS orbitals. The outer loop over occupied states is
 *                     serial on the CPU (each projection depends on the
 *                     previous subtraction); the inner overlap integrals and
 *                     subtraction steps run as GPU kernels. The packet is
 *                     renormalised after all projections. The maximum overlap
 *                     before subtraction is recorded in the report.
 *
 *   4. Norm after     Records ‖ψ_wp‖ of the injected (and optionally
 *                     orthogonalised) state. Should be ≈ 1.0.
 *
 *   5. Occupation     Sets electrons.occupations()[0][ist_wp] to the
 *                     requested value (default 1.0).
 *
 * Builder usage
 * -------------
 *   auto report = inqkit::WavePacket{}
 *       .center(cx_bohr, cy_bohr, cz_bohr)
 *       .sigma(sigma_bohr)
 *       .k0(kx, ky, kz)
 *       .orthogonalise_against_occupied(electrons)   // optional
 *       .inject_into_last_extra_state(electrons, 1.0);
 *
 * The returned InjectionReport carries norm_before, norm_after,
 * max_overlap (if orthogonalised), and passed_tolerance, and is
 * suitable for logging or assertion in the calling run script.
 *
 * Current limitations
 * -------------------
 *   - Single k-point (gamma-only) runs only.
 *   - Single MPI rank only; multi-rank basis/set partitioning will throw.
 *   - Spherical Gaussians only (one σ shared across all three axes).
 *   - Occupation values other than 1.0 are untested.
 */


// TODO: Make a momentum space version of Gram-Schmidt and see if there is any difference
// in a given simulations observables. If there is, then consider it, and make a decision. 

// TODO: Need to check the parallelisation code, and ensure that all is working as
// expected. Need to write tests to prove the understanding of GPU::run gained, and 
// reduce functions. The same must be done with MPI. 

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


// TODO: Can have double orthogonalisation using GS algorithm. 

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
  // All of the input coordinates and values in bohrs
  double cx_ = 0, cy_ = 0, cz_ = 0;
  double sigma_ = 1.0;
  /* TODO: Do not know what units the k values are defined in. Best guess
  is 1/bohr to perform a sanity check using the simulations (free propagation o
  of a wavepacket if it works, or jellium wave packet propagation) to determine
  if the k coordinate is working as expected.
  */
  double kx_ = 0, ky_ = 0, kz_ = 0;
  bool do_ortho_ = false;
  double ortho_tol_ = 1e-6;

public:
  WavePacket &center(double x_bohr, double y_bohr, double z_bohr) {
    cx_ = x_bohr;
    cy_ = y_bohr;
    cz_ = z_bohr;
    return *this;
  }

  WavePacket &sigma(double sigma_bohr) {
    sigma_ = sigma_bohr;
    return *this;
  }

  WavePacket &k0(double kx, double ky, double kz) {
    kx_ = kx;
    ky_ = ky;
    kz_ = kz;
    return *this;
  }

  // Mark the WP for orthogonalisation against all occupied states during
  // injection. Actual projection is performed inside
  // inject_into_last_extra_state using GPU kernels (modified Gram-Schmidt).
  // Passing electrons here is required by the API for consistency with future
  // pre-computation of overlaps; currently unused.

  // TODO: Test the orthogonalisation rigorously. 
  WavePacket &
  orthogonalise_against_occupied(inq::systems::electrons const & /*electrons*/,
                                 double tolerance = 1e-6) {
    do_ortho_ = true;
    ortho_tol_ = tolerance;
    return *this;
  }

  InjectionReport
  inject_into_last_extra_state(inq::systems::electrons &electrons,
                               double occupation = 1.0) const;
};

// ────────────────────────────────────────────────────────────────────────────
// This line injects the gaussian wave packet into the last extra state that is defined
// in the system
inline InjectionReport
WavePacket::inject_into_last_extra_state(inq::systems::electrons &electrons,
                                         double occupation) const {
  using complex = inq::complex;

  // TODO: Update the code such that it would work with any number of kpoint
  // configurations
  if (electrons.kpin().size() != 1) {
    throw std::runtime_error("inqkit::WavePacket: only single-kpoint "
                             "(gamma-only) runs are supported.");
  }

  auto &phi = electrons.kpin()[0];
  auto &basis = phi.basis();

  if (phi.basis().comm().size() != 1 || phi.set_comm().size() != 1) {
    throw std::runtime_error("inqkit::WavePacket: multi-rank basis/set "
                             "partitioning is not supported.");
  }

  // Defining the last electron state (extra state) as the
  // index of the wave packet.
  int ist_wp = phi.set_part().local_size() - 1;
  int n_pts = phi.basis().local_size();
  double dV = basis.volume_element();

  InjectionReport report;
  report.kpoint_index = 0;
  report.state_index = ist_wp;

  // ── 1. Norm of existing last slot (before injection) ───────────────────────
 
  // ip runs over all of the grid points. ist_wp selects the wavefunction
  // value at the given coorindate. Then, the density at that point
  // dV* (phi)^2 is calculated and added up over the entire grid. 
  // The summation should equal one for a self normalised function
  INQKIT_GPU_SYNC();
  {
    auto mat_ = begin(phi.matrix());
    auto res = gpu::run(1, gpu::reduce(n_pts), 0.0,
                        [dV, mat_, ist_wp_ = ist_wp] GPU_LAMBDA(auto, auto ip) {
                          auto v = mat_[ip][ist_wp_];
                          return dV * (inq::real(v) * inq::real(v) +
                                       inq::imag(v) * inq::imag(v));
                        });
    INQKIT_GPU_SYNC();
    report.norm_before = std::sqrt(res[0]);
  }

  // ── 2. Inject raw Gaussian wavepacket ──────────────────────────────────────
  // psi_wp(r) = (pi sigma^2)^{-3/4} exp(-|r-b|^2 / (2 sigma^2)) exp(i k.r)
  //
  // Coordinate convention: INQ stores centred orthorhombic cells with the
  // origin at array index 0 and uses to_symmetric_range so that the physical
  // coordinate at global index ig is (ig <= N/2 ? ig : ig - N) * dr — i.e.
  // r ∈ [-L/2, +L/2]. We therefore obtain the Cartesian position from
  // basis.point_op().rvector_cartesian(...) instead of recomputing it as
  // ig*dr (which would be wrong for ig > N/2 and silently corrupt the
  // e^{i k·r} phase across half the cell).
  double norm_fac = std::pow(M_PI * sigma_ * sigma_, -0.75);
  double sig = sigma_;
  double bx = cx_, by = cy_, bz = cz_; // initial position of the wp
  double kxv = kx_, kyv = ky_, kzv = kz_;
  int ist_w = ist_wp; // index of the wave-packet

  // Gets access to the first element in the 4D vector
  // Used to get the pointer that can be manipulated in GPU accelerated
  // functions.
  auto phicub_ = begin(phi.hypercubic());
  auto point_op_ = basis.point_op();

  gpu::run(basis.local_sizes()[2], basis.local_sizes()[1],
           basis.local_sizes()[0], [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
             auto rvec = point_op_.rvector_cartesian(ix, iy, iz);
             double rx = rvec[0];
             double ry = rvec[1];
             double rz = rvec[2];
             double dx_ = rx - bx, dy_ = ry - by, dz_ = rz - bz;
             double r2 = dx_ * dx_ + dy_ * dy_ + dz_ * dz_;
             double amp = norm_fac * exp(-r2 / (2.0 * sig * sig));
             double ph = kxv * rx + kyv * ry + kzv * rz;
             phicub_[ix][iy][iz][ist_w] = complex(amp * cos(ph), amp * sin(ph));
           });
  INQKIT_GPU_SYNC();

  // ── 3. Orthogonalise against occupied states (modified Gram-Schmidt on GPU)
  // ─
  if (do_ortho_) {
    double max_ov_initial = 0.0;   // pre-ortho overlap (first pass) — reported
    double max_ov_residual = 0.0;  // residual overlap (final pass) — gates tol
    auto mat_ = begin(phi.matrix());

    /* Iterated modified Gram-Schmidt (E03 fix). Each pass runs over all KS
     * states below the WP slot, measures <psi_i|psi_wp> and subtracts the
     * projection. KS orbitals are mutually orthonormal, so a single pass
     * orthogonalises in EXACT arithmetic; a SECOND pass removes the
     * finite-precision residual a single pass leaves. Crucially, measuring the
     * overlap on the FINAL pass lets passed_tolerance reflect the TRUE
     * post-orthogonalisation residual rather than the (large) pre-ortho overlap
     * the first pass sees. The overlap reduction is GPU-accelerated; the loop
     * over states is serial (each subtraction depends on the previous). */
    const int n_ortho_passes = 2;
    for (int pass = 0; pass < n_ortho_passes; ++pass) {
    double pass_max = 0.0;
    for (int i = 0; i < ist_wp; ++i) {
      // Real part of <psi_i | psi_wp>
      auto res_re =
          gpu::run(1, gpu::reduce(n_pts), 0.0,
                   [dV, mat_, i_ = i, ist_ = ist_wp] GPU_LAMBDA(auto, auto ip) {
                     auto vi = mat_[ip][i_];
                     auto vw = mat_[ip][ist_];
                     return dV * (inq::real(vi) * inq::real(vw) +
                                  inq::imag(vi) * inq::imag(vw));
                   });
      INQKIT_GPU_SYNC();
      double ov_re = res_re[0];

      // Imaginary part of <psi_i | psi_wp>
      auto res_im =
          gpu::run(1, gpu::reduce(n_pts), 0.0,
                   [dV, mat_, i_ = i, ist_ = ist_wp] GPU_LAMBDA(auto, auto ip) {
                     auto vi = mat_[ip][i_];
                     auto vw = mat_[ip][ist_];
                     return dV * (inq::real(vi) * inq::imag(vw) -
                                  inq::imag(vi) * inq::real(vw));
                   });
      INQKIT_GPU_SYNC();
      double ov_im = res_im[0];

      pass_max = std::max(pass_max, std::sqrt(ov_re * ov_re + ov_im * ov_im));

      // Subtract projection: psi_wp -= (ov_re + i*ov_im) * psi_i
      double re_ = ov_re, im_ = ov_im;
      int i_ = i;
      gpu::run(basis.local_sizes()[2], basis.local_sizes()[1],
               basis.local_sizes()[0],
               [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
                 auto vi = phicub_[ix][iy][iz][i_];
                 auto vw = phicub_[ix][iy][iz][ist_w];
                 // (ov_re + i*ov_im) * vi
                 double sub_re = re_ * inq::real(vi) - im_ * inq::imag(vi);
                 double sub_im = re_ * inq::imag(vi) + im_ * inq::real(vi);
                 phicub_[ix][iy][iz][ist_w] =
                     complex(inq::real(vw) - sub_re, inq::imag(vw) - sub_im);
               });
      INQKIT_GPU_SYNC();
    }
    if (pass == 0) max_ov_initial = pass_max;
    max_ov_residual = pass_max;
    }  // end iterated-GS pass loop

    report.max_overlap = max_ov_initial;
    report.orthogonalised = true;

    // Renormalise after projection
    {
      auto res = gpu::run(1, gpu::reduce(n_pts), 0.0,
                          [dV, mat_, ist_ = ist_wp] GPU_LAMBDA(auto, auto ip) {
                            auto v = mat_[ip][ist_];
                            return dV * (inq::real(v) * inq::real(v) +
                                         inq::imag(v) * inq::imag(v));
                          });
      INQKIT_GPU_SYNC();
      double scale = 1.0 / std::sqrt(res[0]);

      gpu::run(basis.local_sizes()[2], basis.local_sizes()[1],
               basis.local_sizes()[0],
               [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
                 auto v = phicub_[ix][iy][iz][ist_w];
                 phicub_[ix][iy][iz][ist_w] =
                     complex(inq::real(v) * scale, inq::imag(v) * scale);
               });
      INQKIT_GPU_SYNC();
    }

    report.passed_tolerance = (max_ov_residual < ortho_tol_ * 10.0);
  } else {
    report.passed_tolerance = true;
  }

  // ── 4. Norm after injection (and orthogonalisation if applied) ─────────────
  {
    auto mat_ = begin(phi.matrix());
    auto res = gpu::run(1, gpu::reduce(n_pts), 0.0,
                        [dV, mat_, ist_ = ist_wp] GPU_LAMBDA(auto, auto ip) {
                          auto v = mat_[ip][ist_];
                          return dV * (inq::real(v) * inq::real(v) +
                                       inq::imag(v) * inq::imag(v));
                        });
    INQKIT_GPU_SYNC();
    report.norm_after = std::sqrt(res[0]);
  }

  // ── 5. Set occupation ──────────────────────────────────────────────────────
  electrons.occupations()[0][ist_wp] = occupation;

  return report;
}

} // namespace inqkit
