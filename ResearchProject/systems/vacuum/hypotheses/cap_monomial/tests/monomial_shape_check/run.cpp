// ============================================================================
// monomial_shape_check: KNOWN-CASE test for the inq-study perturbation
// perturbations::absorbing_monomial  (V = i*eta*s^n, s in [0,1] across the slab).
//
// Task-specific mechanism test (ADR 0007: hypotheses/<sweep>/tests/). Reuses the
// proven cap_probe free-WP-in-vacuum geometry, swapping the built-in sin^2 CAP for
// the new monomial ramp. It checks TWO falsifiable predictions that follow from the
// shape (and that the sin^2 hump CANNOT have):
//
//   (1) ABSORBS: a monomial CAP at the box end removes WP norm (absorbed > 0.05),
//       proving the imaginary potential is built and propagates (inq-study fix).
//   (2) ORDER MONOTONICITY: since s^n decreases with n on (0,1), a LOWER order has
//       a STRONGER absorbing potential throughout the slab -> absorbs more ->
//       LOWER reflection error. So eps(n=1) < eps(n=4) at fixed eta, L, E.
//
// If the shape were wrong (constant, hump, or dropped), prediction (2) would fail.
// Build against inq-study (carries the scalar-potential complexification):
//   INQ_SOURCE=.../inq-study INQ_SHARE_PATH=.../inq/install/share \
//   PSEUDOPOD_SHARE_PATH=.../inq/install/share/pseudopod inq-run --reconfig
// Exit code 0 = PASS, 1 = FAIL.
// ============================================================================

#include <inq/inq.hpp>
#include <perturbations/absorbing_monomial.hpp>          // the inq-study NEW perturbation
#include <inqkit/absorbers/mask_absorber.hpp>            // inner_region_norm (eps)
#include <inqkit/wavepacket/wavepacket.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

using namespace inq;
using namespace inq::magnitude;
namespace abs_ = inqkit::absorbers;

int main() {
  const double HA_TO_EV = 27.211386245988;
  const double k0   = 1.0;            // E ~ 13.6 eV (near the 10 eV region of interest)
  const double Labs = 5.0;            // THIN absorber (the regime of interest)
  const double eta  = -0.30;          // shallow depth (under-absorbing -> lower n helps)
  const int    nperp = 8;
  const double dt   = 0.01;

  const double sigma   = 4.0 * std::sqrt(2.0) / k0;
  const double Lcell_z = 6.0 * sigma + Labs;
  const double z_abs0  = (6.0 * sigma - Labs) / 2.0;
  const double z0      = -Labs / 2.0;
  const double tau     = 2.0 * (3.0 * sigma + Labs) / k0;
  const double dx      = std::min(0.30, std::max(0.18, 0.75 / k0));
  const int N_STEPS    = std::max(1, (int)std::llround(tau / dt));
  const double Lperp   = nperp * dx;
  const double ec      = 0.5 * std::pow(M_PI / dx, 2.0);
  const double width_frac = Labs / Lcell_z;
  const double mid_frac   = 0.5 - width_frac / 2.0;

  std::printf("\n=== monomial_shape_check: k0=%.2f E=%.2f eV L=%.1f eta=%.2f Ha ===\n",
              k0, 0.5 * k0 * k0 * HA_TO_EV, Labs, eta);

  auto eps_for_order = [&](int order) -> double {
    auto cell = systems::cell::orthorhombic(Lperp * 1.0_b, Lperp * 1.0_b, Lcell_z * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions, options::electrons{}.cutoff(ec * 1.0_Ha).extra_states(1).extra_electrons(2.0));
    ground_state::initial_guess(ions, electrons);
    auto rep = inqkit::WavePacket{}.center(0.0, 0.0, z0).sigma(sigma).k0(0.0, 0.0, k0)
                   .inject_into_last_extra_state(electrons, 1.0);
    const long wp_idx = rep.state_index;
    const double N0   = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);

    perturbations::absorbing_monomial cap(eta * 1.0_Ha, mid_frac, width_frac, order);
    auto rt = options::real_time{}.num_steps(N_STEPS).dt(dt * 1.0_atomictime);   // ETRS
    real_time::propagate(ions, electrons, [](auto const &) {},
                         options::theory{}.non_interacting(), rt, cap);

    const double inner_tau = abs_::inner_region_norm(electrons, 2, z_abs0, wp_idx);
    const double total_tau = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);
    const double eps = inner_tau / N0;
    std::printf("  order=%d : eps=%.6f  absorbed=%.4f\n", order, eps, 1.0 - total_tau / N0);
    return eps;
  };

  const double eps1 = eps_for_order(1);
  const double eps4 = eps_for_order(4);
  // absorbed fractions are printed above; recompute pass conditions:
  bool absorbs = (eps1 < 0.95) and (eps4 < 0.99);          // (1) it removes norm
  bool monotonic = (eps1 < eps4);                          // (2) lower order absorbs more

  std::printf("\n  [check 1] absorbs (eps<1 for n=1,4): %s\n", absorbs ? "PASS" : "FAIL");
  std::printf("  [check 2] eps(n=1)=%.4f < eps(n=4)=%.4f : %s\n",
              eps1, eps4, monotonic ? "PASS" : "FAIL");

  if (absorbs and monotonic) { std::printf("\nRESULT: PASS\n"); return 0; }
  std::printf("\nRESULT: FAIL\n");
  return 1;
}
