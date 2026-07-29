// ============================================================================
// gate1_mask_absorber: gate-1 validation of the PRODUCTION mask absorber.
//
// Exercises inqkit::absorbers::MaskAbsorber + inner_region_norm (the code the
// sweep uses), not an inline copy. Four checks, each with a pre-accepted
// analytic expected value (NOT retrofitted to code output):
//
//   T1 mask shape (host)   sin2_mask_value: M(z0)=1, M(z0+L)=0, M(z0+L/2)=0.5,
//                          monotone decreasing on (z0, z0+L)         tol 1e-12
//   T2 epsilon known-case  symmetric Gaussian split by inner_region_norm:
//                          z_abs0 far right → eps≈1; far left → eps≈0;
//                          at the centre → eps≈0.5                   tol 1e-3
//   T3 fidelity            MaskAbsorber with M≡1 (z_abs0 beyond box) applied
//                          every step ⇒ bit-identical to no-absorber baseline
//   T4 feedthrough         a real sin² absorber ⇒ surviving norm drops a lot
//
// Exit 0 iff all pass. This is the oracle the overnight pipeline gates on.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/absorbers/mask_absorber.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <cmath>
#include <cstdio>
#include <filesystem>

using namespace inq;
using namespace inq::magnitude;
namespace fs = std::filesystem;
namespace abs_ = inqkit::absorbers;

static int failures = 0;
static void check(bool ok, const char *name, double got, double want,
                  double tol) {
  std::printf("  [%-12s] got=%+.8f want=%+.8f tol=%.0e  %s\n", name, got, want,
              tol, ok ? "PASS" : "FAIL");
  if (!ok) failures++;
}

// One free-WP propagation; optional MaskAbsorber applied every step. Returns
// final (full-box norm, z-centroid).
struct Res { double norm, z; };
static Res propagate_wp(systems::ions &ions, systems::electrons &electrons,
                        double sigma, double k0z, double z0, int n_steps,
                        double dt, const abs_::MaskAbsorber *mask) {
  auto rep = inqkit::WavePacket{}
                 .center(0.0, 0.0, z0)
                 .sigma(sigma)
                 .k0(0.0, 0.0, k0z)
                 .inject_into_last_extra_state(electrons, 1.0);

  fs::path dir = fs::temp_directory_path() / "gate1_mask";
  fs::create_directories(dir);
  // WPRealSpaceStats needs the EXPLICIT global WP index — it does NOT honour a
  // -1 sentinel (unlike MaskAbsorber/inner_region_norm). Passing -1 here zeroes
  // the norm and throws. Use the injection report's state_index.
  inqkit::observables::WPRealSpaceStats wp_rs((dir / "rs.csv").string(),
                                              rep.state_index);

  inqkit::observables::WPRealSpaceMoments mT{};
  real_time::propagate(
      ions, electrons,
      [&](auto const &data) {
        if (mask) mask->apply(electrons);
        mT = wp_rs.compute(data.electrons());
      },
      options::theory{}.non_interacting(),
      options::real_time{}.num_steps(n_steps).dt(dt * 1.0_atomictime));
  return {mT.N, mT.z};
}

int main() {
  std::printf("\n=== gate1_mask_absorber ===\n");

  // ---- T1: mask shape on the host ----------------------------------------
  std::printf("T1 mask shape (host)\n");
  const double z0 = 0.0, Lm = 10.0;
  check(std::abs(abs_::sin2_mask_value(z0 - 1.0, z0, Lm) - 1.0) < 1e-12,
        "M(z<z0)", abs_::sin2_mask_value(z0 - 1.0, z0, Lm), 1.0, 1e-12);
  check(std::abs(abs_::sin2_mask_value(z0, z0, Lm) - 1.0) < 1e-12, "M(z0)",
        abs_::sin2_mask_value(z0, z0, Lm), 1.0, 1e-12);
  check(std::abs(abs_::sin2_mask_value(z0 + Lm / 2, z0, Lm) - 0.5) < 1e-12,
        "M(z0+L/2)", abs_::sin2_mask_value(z0 + Lm / 2, z0, Lm), 0.5, 1e-12);
  check(std::abs(abs_::sin2_mask_value(z0 + Lm, z0, Lm) - 0.0) < 1e-12,
        "M(z0+L)", abs_::sin2_mask_value(z0 + Lm, z0, Lm), 0.0, 1e-12);
  check(std::abs(abs_::sin2_mask_value(z0 + 2 * Lm, z0, Lm) - 0.0) < 1e-12,
        "M(z>z0+L)", abs_::sin2_mask_value(z0 + 2 * Lm, z0, Lm), 0.0, 1e-12);
  // monotone decreasing on (z0, z0+L)
  bool mono = true;
  double prev = 1.0;
  for (int i = 1; i <= 20; i++) {
    double v = abs_::sin2_mask_value(z0 + i * Lm / 20.0, z0, Lm);
    if (v > prev + 1e-14) mono = false;
    prev = v;
  }
  check(mono, "monotone", mono ? 1.0 : 0.0, 1.0, 0.0);

  // ---- engine setup (free WP, cubic periodic) ----------------------------
  const double L = 40.0, SPACING = 0.5, SIGMA = 1.0;
  auto cell =
      systems::cell::orthorhombic(L * 1.0_b, L * 1.0_b, L * 1.0_b).periodic();
  auto ions = systems::ions(cell);
  const double ec = 0.5 * std::pow(M_PI / SPACING, 2.0);
  auto electrons = systems::electrons(
      ions, options::electrons{}.cutoff(ec * 1.0_Ha).extra_states(1).extra_electrons(2.0));
  ground_state::initial_guess(ions, electrons);

  // ---- T2: epsilon known-case (static WP, no propagation) ----------------
  // Cut at a MID-GRID position zc = dx/2 with the WP centred there, so NO grid
  // point lies on the cut and the symmetric Gaussian splits EXACTLY in half.
  // (Cutting at z=0 — a grid point on the Gaussian peak — would exclude that
  // whole plane under the strict z<z_abs0 rule, giving (1-f0)/2, not 0.5.)
  std::printf("T2 epsilon known-case (symmetric Gaussian split, mid-grid cut)\n");
  const double ZC = SPACING / 2.0; // 0.25, midway between grid points 0 and dx
  inqkit::WavePacket{}.center(0.0, 0.0, ZC).sigma(SIGMA).k0(0.0, 0.0, 0.0)
      .inject_into_last_extra_state(electrons, 1.0);
  double eps_all = abs_::inner_region_norm(electrons, 2, +1e6); // all z < z_abs0
  double eps_none = abs_::inner_region_norm(electrons, 2, -1e6); // none
  double eps_half = abs_::inner_region_norm(electrons, 2, ZC);   // split at centre
  check(std::abs(eps_all - 1.0) < 1e-3, "eps_all~1", eps_all, 1.0, 1e-3);
  check(std::abs(eps_none - 0.0) < 1e-3, "eps_none~0", eps_none, 0.0, 1e-3);
  check(std::abs(eps_half - 0.5) < 1e-3, "eps_half~0.5", eps_half, 0.5, 1e-3);

  // ---- T3/T4: fidelity + feedthrough via MaskAbsorber --------------------
  std::printf("T3/T4 fidelity + feedthrough (MaskAbsorber)\n");
  const double K0Z = 2.0, Z0 = -10.0, DT = 0.05;
  const int N = 200;
  // Each propagate_wp wrapped so a thrown exception is a COUNTED failure (oracle
  // integrity) rather than an uncounted abort.
  auto safe_prop = [&](const abs_::MaskAbsorber *m, const char *what) -> Res {
    try {
      return propagate_wp(ions, electrons, SIGMA, K0Z, Z0, N, DT, m);
    } catch (std::exception const &e) {
      std::printf("  [%-12s] threw: %s  FAIL\n", what, e.what());
      failures++;
      return {-1.0, 0.0};
    }
  };
  auto base = safe_prop(nullptr, "base");
  // M≡1: absorber placed entirely to the right of the box ⇒ M=1 everywhere
  abs_::MaskAbsorber noop(2, +1e6, 1.0);
  auto idn = safe_prop(&noop, "M==1");
  // baseline must itself be a sane free WP (norm ≈ 1) for the diffs to mean anything
  check(std::abs(base.norm - 1.0) < 0.05, "base_norm~1", base.norm, 1.0, 0.05);
  check(std::abs(idn.norm - base.norm) < 1e-9 && std::abs(idn.z - base.z) < 1e-9,
        "fidelity", std::abs(idn.norm - base.norm), 0.0, 1e-9);
  // real absorber over z in [0,15]: surviving-norm drop must be substantial AND
  // bounded above by the baseline (two-sided — a mis-targeted mask that nukes
  // the wrong state, or total annihilation, must not pass).
  abs_::MaskAbsorber absb(2, 0.0, 15.0);
  auto ab = safe_prop(&absb, "sin2-abs");
  double drop = base.norm - ab.norm;
  check(drop > 0.5 && drop < base.norm && ab.norm > 1e-4, "feedthrough", drop,
        1.0, 0.5);

  std::printf("\n=== gate1 verdict: %s (%d failure%s) ===\n",
              failures == 0 ? "PASS" : "FAIL", failures,
              failures == 1 ? "" : "s");
  return failures == 0 ? 0 : 1;
}
