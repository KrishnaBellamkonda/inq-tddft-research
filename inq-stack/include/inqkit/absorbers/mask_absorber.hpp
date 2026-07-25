/* inqkit::absorbers — mask function absorber (MFA)
 *
 * Implements the sin² mask absorber of De Giovannini, Larsen & Rubio,
 * "Modeling electron dynamics coupled to continuum states in finite volumes",
 * arXiv:1409.1689 (2014), Eq. 12–13:
 *
 *    psi(t+dt) = M(z) * U(t+dt,t) * psi(t)            (Eq. 12)
 *
 *   M(z) = 1                          for z <= z_abs0          (Eq. 13, x<0)
 *        = 1 - sin^2( (z-z_abs0) pi / (2L) )   for z_abs0 < z < z_abs0 + L
 *        = 0                          for z >= z_abs0 + L
 *
 * The absorber acts along ONE Cartesian axis (default z, index 2) on a single
 * orbital (the wavepacket). It is applied ENTIRELY in this wrapper, inside the
 * per-step `real_time::propagate` callback, by multiplying the orbital in place
 * on the GPU. INQ's `viewables` observer is const; the callback instead mutates
 * the captured non-const `electrons` object that propagate() holds by reference,
 * so the masked orbital is what the next ETRS step propagates. This realises
 * Eq. 12 exactly and leaves inq/ and inq-study byte-identical. (Mechanism
 * verified: fidelity PASS, feedthrough PASS — see
 * ResearchProject/systems/vacuum/tests/mask_mechanism_check/.)
 *
 * `inner_region_norm()` computes the reflection error epsilon (paper Eq. 7) as
 * the surviving norm of the orbital in the inner region z < z_abs0 at the final
 * time tau:  epsilon = ∫_{z<z_abs0} |psi_wp|^2 dV.
 *
 * GPU loops mirror the established inqkit patterns:
 *   - mask multiply  → WavePacket::inject_into_last_extra_state (wavepacket.hpp)
 *   - region reduce  → WPRealSpaceStats::compute (wp_real_space_stats.hpp)
 *
 * Limitations (inherited from WavePacket / WPRealSpaceStats): single-kpoint
 * (gamma-only), single MPI rank.
 */
#pragma once

#include <inq/inq.hpp>

#include <inqkit/absorbers/mask_shape.hpp>  // sin2_mask_value[_twosided] (INQ-free)

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
namespace absorbers {

// sin2_mask_value / sin2_mask_value_twosided are defined in mask_shape.hpp
// (INQ-free, pure-tier tested). The GPU kernels below inline the same maths.

class MaskAbsorber {
public:
  // axis: 0=x, 1=y, 2=z (propagation axis). z_abs0 = absorber start (Bohr),
  // L = absorber width (Bohr). wp_idx = global state index of the orbital to
  // mask (default -1 → the last extra state, where WavePacket injects).
  MaskAbsorber(int axis, double z_abs0, double L, long wp_idx = -1)
      : axis_(axis), z_abs0_(z_abs0), L_(L), wp_idx_(wp_idx) {
    if (axis_ < 0 || axis_ > 2)
      throw std::runtime_error("MaskAbsorber: axis must be 0, 1 or 2.");
    if (L_ <= 0.0)
      throw std::runtime_error("MaskAbsorber: absorber width L must be > 0.");
  }

  // Multiply the masked orbital by M(coordinate) in place on the GPU. Call this
  // once per step, AFTER the ETRS step (inside the propagate callback), passing
  // the same non-const `electrons` propagate() holds by reference.
  void apply(inq::systems::electrons &electrons) const {
    using namespace inq;
    auto &phi = electrons.kpin()[0];
    auto &basis = phi.basis();
    const int ist = resolve_local_index(phi);
    const int ax = axis_;
    const double z0 = z_abs0_;
    const double LL = L_;

    auto phic = begin(phi.hypercubic());
    auto point_op = basis.point_op();
    auto const sizes = basis.local_sizes();

    gpu::run(sizes[2], sizes[1], sizes[0],
             [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
               auto r = point_op.rvector_cartesian(ix, iy, iz);
               double s = r[ax];
               double M;
               if (s <= z0) {
                 M = 1.0;
               } else if (s >= z0 + LL) {
                 M = 0.0;
               } else {
                 double sn = sin(M_PI * (s - z0) / (2.0 * LL));
                 M = 1.0 - sn * sn;
               }
               phic[ix][iy][iz][ist] = phic[ix][iy][iz][ist] * M;
             });
    INQKIT_GPU_SYNC();
  }

private:
  int axis_;
  double z_abs0_;
  double L_;
  long wp_idx_;

  template <typename Phi> int resolve_local_index(Phi const &phi) const {
    if (phi.basis().comm().size() != 1 || phi.set_comm().size() != 1)
      throw std::runtime_error(
          "MaskAbsorber: multi-rank basis/set partitioning is not supported.");
    const int last = static_cast<int>(phi.set_part().local_size()) - 1;
    if (wp_idx_ < 0) return last;
    return static_cast<int>(wp_idx_);
  }
};

// Two-sided mask absorber: a symmetric sin^2 mask on BOTH boundaries of `axis`,
// applied once per step after the ETRS step (same in-callback mutation as
// MaskAbsorber). z_in = inner-region half-width (Bohr); Lhalf = per-end width
// (= L_total/2). The absorber spans |z| in [z_in, z_in+Lhalf] at each end. Used
// by the two-sided CAP-vs-mask study (docs/plans/twosided-cap-vs-mask.md); the
// CAP counterpart is `absorbing(+end) + absorbing(-end)` via perturbations::sum.
class TwoSidedMaskAbsorber {
public:
  TwoSidedMaskAbsorber(int axis, double z_in, double Lhalf, long wp_idx = -1)
      : axis_(axis), z_in_(z_in), Lhalf_(Lhalf), wp_idx_(wp_idx) {
    if (axis_ < 0 || axis_ > 2)
      throw std::runtime_error("TwoSidedMaskAbsorber: axis must be 0, 1 or 2.");
    if (Lhalf_ <= 0.0)
      throw std::runtime_error("TwoSidedMaskAbsorber: per-end width Lhalf must be > 0.");
    if (z_in_ <= 0.0)
      throw std::runtime_error("TwoSidedMaskAbsorber: inner half-width z_in must be > 0.");
  }

  void apply(inq::systems::electrons &electrons) const {
    using namespace inq;
    auto &phi = electrons.kpin()[0];
    auto &basis = phi.basis();
    const int ist = resolve_local_index(phi);
    const int ax = axis_;
    const double zin = z_in_;
    const double Lh = Lhalf_;

    auto phic = begin(phi.hypercubic());
    auto point_op = basis.point_op();
    auto const sizes = basis.local_sizes();

    gpu::run(sizes[2], sizes[1], sizes[0],
             [=] GPU_LAMBDA(auto iz, auto iy, auto ix) {
               auto r = point_op.rvector_cartesian(ix, iy, iz);
               double a = fabs(r[ax]);
               double M;
               if (a <= zin) {
                 M = 1.0;
               } else if (a >= zin + Lh) {
                 M = 0.0;
               } else {
                 double sn = sin(M_PI * (a - zin) / (2.0 * Lh));
                 M = 1.0 - sn * sn;
               }
               phic[ix][iy][iz][ist] = phic[ix][iy][iz][ist] * M;
             });
    INQKIT_GPU_SYNC();
  }

private:
  int axis_;
  double z_in_;
  double Lhalf_;
  long wp_idx_;

  template <typename Phi> int resolve_local_index(Phi const &phi) const {
    if (phi.basis().comm().size() != 1 || phi.set_comm().size() != 1)
      throw std::runtime_error(
          "TwoSidedMaskAbsorber: multi-rank basis/set partitioning is not supported.");
    const int last = static_cast<int>(phi.set_part().local_size()) - 1;
    if (wp_idx_ < 0) return last;
    return static_cast<int>(wp_idx_);
  }
};

// Symmetric inner-region norm for the two-sided geometry: the orbital's surviving
// norm in |coordinate| < z_in (the inner region between the two absorbers). This
// is the two-sided reflection error ε (un-absorbed fraction); divide by N0.
inline double inner_region_norm_twosided(inq::systems::electrons const &electrons,
                                         int axis, double z_in, long wp_idx = -1) {
  using namespace inq;
  if (electrons.kpin_size() != 1)
    throw std::runtime_error(
        "inner_region_norm_twosided: only single-kpoint (gamma-only) runs supported.");

  auto const &phi = electrons.kpin()[0];
  auto const &basis = phi.basis();
  const double dV = basis.volume_element();

  const long st_start = phi.set_part().start();
  const long st_size = phi.set_part().local_size();
  const long idx_global = (wp_idx < 0) ? (st_start + st_size - 1) : wp_idx;
  const bool wp_local =
      (idx_global >= st_start && idx_global < st_start + st_size);
  if (!wp_local) return 0.0;
  const int ist = static_cast<int>(idx_global - st_start);
  const int ax = axis;
  const double zin = z_in;

  auto const sizes = basis.local_sizes();
  auto phic = begin(phi.hypercubic());
  auto point_op = basis.point_op();

  double eps = gpu::run(
      gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]), 0.0,
      [phic, ist, dV, ax, zin, point_op] GPU_LAMBDA(auto iz, auto iy, auto ix) {
        auto r = point_op.rvector_cartesian(ix, iy, iz);
        if (fabs(r[ax]) >= zin) return 0.0;
        auto v = phic[ix][iy][iz][ist];
        return dV * (inq::real(v) * inq::real(v) + inq::imag(v) * inq::imag(v));
      });
  return eps;
}

// Reflection error epsilon (paper Eq. 7): the orbital's surviving norm in the
// inner region (coordinate along `axis` < z_abs0). For a normalised WP this is
// the reflected fraction at t = tau: ~0 for a perfect absorber, ~1 for a hard
// wall. const-correct (read-only): safe to call from a viewables observer.
inline double inner_region_norm(inq::systems::electrons const &electrons,
                                int axis, double z_abs0, long wp_idx = -1) {
  using namespace inq;
  if (electrons.kpin_size() != 1)
    throw std::runtime_error(
        "inner_region_norm: only single-kpoint (gamma-only) runs supported.");

  auto const &phi = electrons.kpin()[0];
  auto const &basis = phi.basis();
  const double dV = basis.volume_element();

  const long st_start = phi.set_part().start();
  const long st_size = phi.set_part().local_size();
  const long idx_global = (wp_idx < 0) ? (st_start + st_size - 1) : wp_idx;
  const bool wp_local =
      (idx_global >= st_start && idx_global < st_start + st_size);
  if (!wp_local) return 0.0;
  const int ist = static_cast<int>(idx_global - st_start);
  const int ax = axis;
  const double z0 = z_abs0;

  auto const sizes = basis.local_sizes();
  auto phic = begin(phi.hypercubic());
  auto point_op = basis.point_op();

  double eps = gpu::run(
      gpu::reduce(sizes[2]), gpu::reduce(sizes[1]), gpu::reduce(sizes[0]), 0.0,
      [phic, ist, dV, ax, z0, point_op] GPU_LAMBDA(auto iz, auto iy, auto ix) {
        auto r = point_op.rvector_cartesian(ix, iy, iz);
        if (r[ax] >= z0) return 0.0;
        auto v = phic[ix][iy][iz][ist];
        return dV * (inq::real(v) * inq::real(v) + inq::imag(v) * inq::imag(v));
      });
  return eps;
}

} // namespace absorbers
} // namespace inqkit
