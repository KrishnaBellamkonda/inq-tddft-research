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
  bool minimum_image_ = false;
  double ortho_tol_ = 1e-6;
  // Longitudinal (z) focusing: launch a wider, converging packet whose waist
  // (density std sigma_/sqrt2) forms a focal distance ahead (e.g. the slab face).
  bool   do_focus_z_ = false;
  double focus_dist_ = 0.0;   // Bohr, launch -> focal point along +z
  double focus_mass_ = 1.0;   // effective mass of the WP (sets v = kz/m)

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

  // Build the packet from the MINIMUM-IMAGE displacement, so it WRAPS around the
  // cell faces instead of being CLIPPED by them. Defaults to false: every
  // previously published run keeps its exact behaviour.
  //
  // WHY THIS EXISTS (2026-08-01, channeling twin). It is often said that a
  // wavepacket needs no special boundary treatment because a KS orbital lives on
  // a plain 3-D FFT basis and wraps exactly. That is true of the PROPAGATION and
  // false of the INJECTION: the kernel below builds the Gaussian from a plain
  // Cartesian displacement, so a packet launched within a couple of sigma of a
  // face is TRUNCATED, and normalisation then hides it in the norm.
  //
  // The damage is not subtle. A sigma_WP = 4 packet launched 2 Bohr (0.71 density
  // sigma) from the -z face failed six of its own t=0 analytic gates, and every
  // failure was in z alone while x and y were perfect:
  //     <p_z>        1.882 vs 1.917   (-1.8 %)
  //     var(p_z)     0.473 vs 0.0313  (+1413 %)  <- the sharp real-space edge
  //     T1 - T2      0.268 vs 0.0469  (+471 %)
  //     centroid z  -26.97 vs -28     (truncating the left tail pulls the mean right)
  //     sigma_z      2.133 vs 2.828   (-24.6 %)  (a narrower, truncated packet)
  // i.e. the run would have measured a packet that was not the one it claimed,
  // with a momentum spread fifteen times too large -- fatal for any observable
  // built on var(p).
  //
  // In a classical/wavepacket TWIN this also breaks the pair: the classical half
  // uses gaussian_density_minimum_image, so a clipped wavepacket differs from its
  // twin precisely at the boundary the study introduces on purpose.
  //
  // The PHASE is built from the same minimum-image displacement, not from the raw
  // coordinate. Wrapping the amplitude while leaving the phase as exp(i k.r) puts
  // a jump of exp(i k.L) across the seam whenever k.L is not a multiple of 2 pi
  // (here k0*L_z = 115.0 rad = 18.3 x 2 pi -- a 1.9 rad discontinuity). Using
  // exp(i k.d) makes the local momentum k everywhere and differs from the old
  // form only by the global constant exp(i k.b).
  WavePacket &minimum_image(bool on = true) {
    minimum_image_ = on;
    return *this;
  }

  // Launch a LONGITUDINALLY FOCUSING packet: instead of a minimum-width Gaussian
  // at t=0, inject a wider, converging packet along z whose waist (density std
  // sigma_/sqrt2) forms after travelling `focal_distance_bohr` (launch -> slab
  // face). `effective_mass` is the WP mass (v = kz/m). Transverse (x,y) stay at
  // sigma_. Width + chirp derived from time-reversed free propagation of a
  // min-width Gaussian (validated 1D known-case). Requires k0(...,kz>0).
  WavePacket &focus_z(double focal_distance_bohr, double effective_mass) {
    do_focus_z_ = true;
    focus_dist_ = focal_distance_bohr;
    focus_mass_ = effective_mass;
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

  // Longitudinal focusing setup. Default (no focus): sigz = sig, chirp = 0 — the
  // ordinary spherical min-width packet. With focus_z() set, launch a wider,
  // converging z-packet: q^2 = a0^4 + (tau/m)^2 with a0 = sigma_, tau = D/v,
  // v = kz/m; launch z-width sigz = |q|/a0, chirp = -(tau/m)/(2 q^2). Derived from
  // time-reversed free propagation (1D-validated); the anisotropic norm_fac is
  // approximate and corrected by the post-orthogonalisation renormalise.
  double sigz = sig;
  double chirp = 0.0;
  if (do_focus_z_) {
    double a0 = sigma_;
    double vgz = kzv / focus_mass_;
    double tau = focus_dist_ / vgz;
    double tm = tau / focus_mass_;
    double q2 = a0 * a0 * a0 * a0 + tm * tm;
    sigz = std::sqrt(q2) / a0;
    chirp = -tm / (2.0 * q2);
  }

  // Gets access to the first element in the 4D vector
  // Used to get the pointer that can be manipulated in GPU accelerated
  // functions.
  auto phicub_ = begin(phi.hypercubic());
  auto point_op_ = basis.point_op();
  const bool min_img = minimum_image_;

  gpu::run(basis.local_sizes()[2], basis.local_sizes()[1],
           basis.local_sizes()[0], [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
             auto rvec = point_op_.rvector_cartesian(ix, iy, iz);
             double rx = rvec[0];
             double ry = rvec[1];
             double rz = rvec[2];
             double dx_ = rx - bx, dy_ = ry - by, dz_ = rz - bz;
             if(min_img) {
               // Fold the separation into [-L/2, L/2) per lattice direction (the
               // same window as systems::cell::position_in_cell), so the packet
               // wraps around every face instead of being clipped by it.
               inq::vector3<double> dsep{dx_, dy_, dz_};
               auto fr = point_op_.cell().to_contravariant(dsep);
               for(int idir = 0; idir < 3; idir++) {
                 fr[idir] -= floor(fr[idir]);
                 if(fr[idir] >= 0.5) fr[idir] -= 1.0;
               }
               auto dc = point_op_.cell().to_cartesian(fr);
               dx_ = dc[0]; dy_ = dc[1]; dz_ = dc[2];
             }
             // anisotropic when focusing (sigz != sig); spherical otherwise
             double amp = norm_fac *
                 exp(-(dx_ * dx_ + dy_ * dy_) / (2.0 * sig * sig) -
                     dz_ * dz_ / (2.0 * sigz * sigz));
             // converging quadratic phase on z (chirp = 0 when not focusing).
             // In minimum-image mode the phase MUST use the same wrapped
             // displacement, or the wrapped lobe carries a exp(i k.L) jump.
             double ph = min_img ? (kxv * dx_ + kyv * dy_ + kzv * dz_ + chirp * dz_ * dz_)
                                 : (kxv * rx  + kyv * ry  + kzv * rz  + chirp * dz_ * dz_);
             phicub_[ix][iy][iz][ist_w] = complex(amp * cos(ph), amp * sin(ph));
           });
  INQKIT_GPU_SYNC();

  // ── 2b. Norm of the RAW Gaussian, before any orthogonalisation ────────────
  // The analytic norm_fac normalises the CONTINUUM Gaussian; on a finite grid
  // the discrete norm is only ~1. removed_weight must therefore be a RATIO
  // against this measured value, not against a hard 1.0, or the discretisation
  // error would masquerade as orthogonalisation loss.
  {
    auto mat_ = begin(phi.matrix());
    auto res = gpu::run(1, gpu::reduce(n_pts), 0.0,
                        [dV, mat_, ist_ = ist_wp] GPU_LAMBDA(auto, auto ip) {
                          auto v = mat_[ip][ist_];
                          return dV * (inq::real(v) * inq::real(v) +
                                       inq::imag(v) * inq::imag(v));
                        });
    INQKIT_GPU_SYNC();
    report.norm_pre_ortho = std::sqrt(res[0]);
  }

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
    double pass_sum_sq = 0.0;   // sum_i |<psi_i|psi_wp>|^2, FIRST pass only
    for (int pass = 0; pass < n_ortho_passes; ++pass) {
    double pass_max = 0.0;
    double this_pass_sum_sq = 0.0;
    if (pass == 0) report.overlap_by_state.assign(ist_wp, 0.0);
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

      const double ov_sq = ov_re * ov_re + ov_im * ov_im;
      pass_max = std::max(pass_max, std::sqrt(ov_sq));
      this_pass_sum_sq += ov_sq;
      if (pass == 0) report.overlap_by_state[i] = std::sqrt(ov_sq);

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
    if (pass == 0) { max_ov_initial = pass_max; pass_sum_sq = this_pass_sum_sq; }
    max_ov_residual = pass_max;
    }  // end iterated-GS pass loop

    report.max_overlap = max_ov_initial;
    report.sum_overlap_sq = pass_sum_sq;
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
      // res[0] = ||psi||^2 AFTER the projection but BEFORE rescaling — the one
      // moment at which the orthogonalisation loss is still visible. Capture it
      // here; one line later it is scaled away for good.
      report.norm_pre_renorm = std::sqrt(res[0]);
      if (report.norm_pre_ortho > 0.0) {
        const double ratio = report.norm_pre_renorm / report.norm_pre_ortho;
        report.removed_weight = 1.0 - ratio * ratio;
      }
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
    // No orthogonalisation requested: nothing was removed, and the packet is
    // left exactly as constructed (not even renormalised).
    report.norm_pre_renorm = report.norm_pre_ortho;
    report.removed_weight  = 0.0;
    report.sum_overlap_sq  = 0.0;
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
