// Engine-tier INTEGRATION test (IV-M11): a free Gaussian wave packet propagated
// under NON-INTERACTING theory (kinetic-only H) reproduces the exact analytic
// free-particle evolution. This is the first inqkit integration test — it
// chains WavePacket injection → real_time::propagate → WP{RealSpace,Momentum}
// Stats::compute() and checks all of them against closed-form values.
//
// Canonical free-WP recipe (mirrors run_free_wp_*): empty `ions` (no external
// potential), a "ghost" occupied orbital via extra_electrons(2.0) so INQ has
// num_electrons>0, the WP injected into the single extra state, propagated with
// options::theory{}.non_interacting().
//
// Analytic free Gaussian (atomic units, m=ℏ=1). The inqkit injector writes
// |ψ|² ∝ exp(-r²/σ²), so the density variance is σ²/2 (verified separately by
// test_wp_real_space_compute_engine). A free minimum-uncertainty Gaussian then
// obeys, per Cartesian axis:
//     ⟨r⟩(t)      = r₀ + (k₀/m) t                 (ballistic centroid)
//     Var(t)      = σ²/2 + t²/(2 σ²)              (dispersive spreading)
//     ⟨p⟩(t)      = k₀                            (momentum conserved)
//     ∫|ψ|² dV    = 1                             (norm conserved)
//     E_kin(t)    = E_kin(0)                      (free particle)
// Reference: standard wave-packet spreading, e.g. Sakurai & Napolitano,
// "Modern Quantum Mechanics", free Gaussian wave-packet section.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>

#include <cmath>
#include <filesystem>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
namespace fs = std::filesystem;

TEST_CASE("free WP: non-interacting Gaussian matches analytic spreading + drift",
          "[wavepacket][free][integration][engine]") {
  // --- geometry / WP parameters (small, fast; WP stays clear of the box edge)
  const double L       = 20.0;     // Bohr, cubic periodic
  const double SPACING = 0.5;      // Bohr
  const double SIGMA   = 2.0;      // injector sigma → density Var = σ²/2 = 2.0
  const double K0Z     = 0.8;      // Bohr⁻¹ along z (modest, no boundary wrap)
  const int    N_STEPS = 60;
  const double DT      = 0.05;     // a.u.  → T = 3.0 a.u.
  const double T       = N_STEPS * DT;

  // analytic predictions
  const double V0      = SIGMA * SIGMA / 2.0;            // 2.0
  const double VAR_T   = V0 + (T * T) / (2.0 * SIGMA * SIGMA);   // 3.125
  const double Z_T     = K0Z * T;                        // 2.4 (m = 1)

  auto cell = systems::cell::orthorhombic(L * 1.0_b, L * 1.0_b, L * 1.0_b).periodic();
  auto ions = systems::ions(cell);

  // ghost occupied orbital (state 0) so INQ has electrons; WP → the extra state.
  const double ec_ha = 0.5 * std::pow(M_PI / SPACING, 2.0);
  auto electrons = systems::electrons(
      ions, options::electrons{}.cutoff(ec_ha * 1.0_Ha)
                                .extra_states(1)
                                .extra_electrons(2.0));
  ground_state::initial_guess(ions, electrons);

  auto report = inqkit::WavePacket{}
                    .center(0.0, 0.0, 0.0).sigma(SIGMA).k0(0.0, 0.0, K0Z)
                    .inject_into_last_extra_state(electrons, 1.0);
  const int wp_idx = report.state_index;
  REQUIRE(report.norm_after == Approx(1.0).margin(0.05));

  fs::path dir = fs::temp_directory_path() / "inqkit_test_free_wp";
  fs::remove_all(dir);
  fs::create_directories(dir);
  inqkit::observables::WPRealSpaceStats wp_rs((dir / "rs.csv").string(), wp_idx);
  inqkit::observables::WPMomentumStats  wp_mom((dir / "mom.csv").string(), wp_idx);

  // --- t = 0 -------------------------------------------------------------
  auto m0 = wp_rs.compute(electrons);
  auto p0 = wp_mom.compute(electrons);

  CHECK(m0.N  == Approx(1.0).margin(0.05));
  CHECK(m0.x  == Approx(0.0).margin(0.10));
  CHECK(m0.y  == Approx(0.0).margin(0.10));
  CHECK(m0.z  == Approx(0.0).margin(0.10));
  CHECK(m0.sx2 == Approx(V0).margin(0.35));    // density Var = σ²/2
  CHECK(m0.sz2 == Approx(V0).margin(0.35));
  CHECK(p0.pz == Approx(K0Z).margin(0.04));    // injected momentum
  CHECK(p0.px == Approx(0.0).margin(0.04));
  CHECK(p0.py == Approx(0.0).margin(0.04));

  // --- propagate freely, capturing the final-step moments ----------------
  inqkit::observables::WPRealSpaceMoments mT{};
  inqkit::observables::WPMomentumMoments  pT{};
  real_time::propagate(
      ions, electrons,
      [&](auto const &data) {
        mT = wp_rs.compute(data.electrons());     // overwritten each step → last
        pT = wp_mom.compute(data.electrons());
      },
      options::theory{}.non_interacting(),
      options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime));

  // --- analytic free-particle checks at t = T ----------------------------
  CHECK(mT.N  == Approx(1.0).margin(0.05));               // norm conserved
  CHECK(mT.x  == Approx(0.0).margin(0.15));               // no transverse drift
  CHECK(mT.y  == Approx(0.0).margin(0.15));
  CHECK(mT.z  == Approx(Z_T).margin(0.20));               // ballistic centroid
  CHECK(mT.sz2 == Approx(VAR_T).margin(0.6));             // dispersive spreading
  CHECK(mT.sx2 == Approx(VAR_T).margin(0.6));             // isotropic spreading

  CHECK(pT.pz == Approx(K0Z).margin(0.04));               // momentum conserved
  CHECK(pT.px == Approx(0.0).margin(0.04));
  CHECK(pT.py == Approx(0.0).margin(0.04));
  CHECK(pT.ekin == Approx(p0.ekin).epsilon(0.05));        // kinetic E conserved

  // --- qualitative sanity: the packet really moved and really spread -----
  CHECK(mT.z  > m0.z + 1.0);                              // moved downrange
  CHECK(mT.sz2 > m0.sz2 + 0.3);                           // measurably wider

  fs::remove_all(dir);
}
