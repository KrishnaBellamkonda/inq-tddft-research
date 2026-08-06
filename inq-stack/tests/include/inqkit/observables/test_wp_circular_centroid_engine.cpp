// Known-case test for the PERIODIC-AWARE (circular) position moments added to
// WPRealSpaceStats on 2026-07-30 for the bulk-jellium KS-stopping runs
// (docs/plans/bulk-jellium-ks-stopping.md, work item W1).
//
// WHY A HAND-WRITTEN FIELD AND NOT inqkit::WavePacket:
// the injector fills psi from the RAW Cartesian displacement (r - b), not the
// minimum image (see wavepacket.hpp GPU_LAMBDA), so injecting near a face gives a
// TRUNCATED Gaussian, not a wrapped one. To test the estimator against a genuine
// periodic packet we build the minimum-image Gaussian ourselves. This also keeps
// the test non-circular: the expected values below are analytic, derived from the
// wrapped-Gaussian relations, not read back from the implementation.
//
// Analytic expectations for psi ~ exp(-|r-b|^2 / 2 sigma^2), so |psi|^2 is a
// Gaussian of standard deviation sigma_d = sigma/sqrt(2):
//
//   <r>_circ = b                                   (exact, any b in the cell)
//   R_d      = exp(-(2 pi sigma_d / L)^2 / 2)      (wrapped-Gaussian resultant)
//   sigma_circ = (L/2pi) sqrt(-2 ln R) = sigma_d   (the inverse of the above)
//
// The decisive case is C2: a packet centred 0.5 Bohr from the +z face. The
// circular estimator must recover it; the naive integral ∫z|psi|^2 dV must NOT
// (it averages the two halves of the split packet and lands mid-cell). That
// failure is the entire reason the circular columns exist.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>

#include <cmath>
#include <filesystem>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
namespace fs = std::filesystem;

namespace {

constexpr double L_BOX = 16.0;   // Bohr, cubic periodic
constexpr double SIGMA = 1.5;    // psi-width; density std = SIGMA/sqrt(2)

// Fill state `ist` with a MINIMUM-IMAGE Gaussian centred at (bx,by,bz).
// Normalisation is irrelevant: every moment in WPRealSpaceStats::compute() is
// divided by the norm N, so an unnormalised field gives identical moments.
void fill_periodic_gaussian(systems::electrons& electrons, int ist,
                            double bx, double by, double bz, double sig) {
  auto& phi        = electrons.kpin()[0];
  auto  basis      = phi.basis();
  auto  phicub     = begin(phi.hypercubic());
  auto  point_op   = basis.point_op();
  auto  sizes      = basis.local_sizes();

  gpu::run(sizes[2], sizes[1], sizes[0],
           [=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
             auto r = point_op.rvector_cartesian(ix, iy, iz);
             // Minimum image: wrap the fractional displacement into [-1/2, 1/2).
             auto d = inq::vector3<double, inq::cartesian>{r[0]-bx, r[1]-by, r[2]-bz};
             auto f = point_op.cell().to_contravariant(d);
             for (int k = 0; k < 3; ++k) f[k] -= floor(f[k] + 0.5);
             auto dmin = point_op.cell().to_cartesian(f);
             const double r2 = dmin[0]*dmin[0] + dmin[1]*dmin[1] + dmin[2]*dmin[2];
             phicub[ix][iy][iz][ist] = complex(exp(-r2 / (2.0*sig*sig)), 0.0);
           });
}

// Fill state `ist` with a constant — a maximally delocalised packet.
void fill_uniform(systems::electrons& electrons, int ist) {
  auto& phi      = electrons.kpin()[0];
  auto  phicub   = begin(phi.hypercubic());
  auto  sizes    = phi.basis().local_sizes();
  gpu::run(sizes[2], sizes[1], sizes[0],
           [=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
             phicub[ix][iy][iz][ist] = complex(1.0, 0.0);
           });
}

systems::electrons make_electrons() {
  systems::ions ions(systems::cell::cubic(L_BOX * 1.0_bohr).periodic());
  systems::electrons electrons(
      ions, options::electrons{}.spacing(0.4 * 1.0_bohr)
                                .extra_electrons(2)
                                .extra_states(2));
  ground_state::initial_guess(ions, electrons);
  return electrons;
}

}  // namespace

TEST_CASE("WPRealSpaceStats circular centroid: interior packet — circular == naive",
          "[observables][wp_real_space][circular][engine]") {
  auto electrons = make_electrons();
  REQUIRE(electrons.kpin()[0].set_comm().size() == 1);   // single-rank test
  const int ist = electrons.states().num_states() - 1;

  const double BX = 1.0, BY = -1.0, BZ = 0.5;   // well away from every face
  fill_periodic_gaussian(electrons, ist, BX, BY, BZ, SIGMA);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_wpcirc_interior";
  fs::remove_all(dir);
  inqkit::observables::WPRealSpaceStats stats((dir / "s.csv").string(), ist);
  auto m = stats.compute(electrons);

  const double sigma_d = SIGMA / std::sqrt(2.0);

  // Both estimators are correct here, and must agree closely with each other.
  CHECK(m.x  == Approx(BX).margin(0.02));
  CHECK(m.y  == Approx(BY).margin(0.02));
  CHECK(m.z  == Approx(BZ).margin(0.02));
  CHECK(m.xc == Approx(BX).margin(0.02));
  CHECK(m.yc == Approx(BY).margin(0.02));
  CHECK(m.zc == Approx(BZ).margin(0.02));
  CHECK(m.zc == Approx(m.z).margin(0.01));

  // Analytic wrapped-Gaussian resultant R = exp(-(2 pi sigma_d / L)^2 / 2).
  const double theta = 2.0 * M_PI * sigma_d / L_BOX;
  const double R_exp = std::exp(-0.5 * theta * theta);
  CHECK(m.Rx == Approx(R_exp).epsilon(0.02));
  CHECK(m.Ry == Approx(R_exp).epsilon(0.02));
  CHECK(m.Rz == Approx(R_exp).epsilon(0.02));

  // The circular spread inverts that relation, so it must return sigma_d, i.e.
  // the same width the naive variance gives for an interior packet.
  CHECK(m.sxc == Approx(sigma_d).epsilon(0.03));
  CHECK(m.szc == Approx(sigma_d).epsilon(0.03));
  CHECK(m.szc == Approx(std::sqrt(m.sz2)).epsilon(0.05));

  fs::remove_all(dir);
}

TEST_CASE("WPRealSpaceStats circular centroid: packet straddling the +z face",
          "[observables][wp_real_space][circular][engine]") {
  auto electrons = make_electrons();
  REQUIRE(electrons.kpin()[0].set_comm().size() == 1);
  const int ist = electrons.states().num_states() - 1;

  // 0.5 Bohr inside the +z face (face at +L/2 = +8). With sigma_d = 1.06 this
  // puts ~32% of the density through the boundary and into the -z end of the box.
  const double BZ = 0.5 * L_BOX - 0.5;   // = 7.5
  fill_periodic_gaussian(electrons, ist, 0.0, 0.0, BZ, SIGMA);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_wpcirc_straddle";
  fs::remove_all(dir);
  inqkit::observables::WPRealSpaceStats stats((dir / "s.csv").string(), ist);
  auto m = stats.compute(electrons);

  // THE POINT OF THE WHOLE FEATURE: the circular estimator is still exact.
  CHECK(m.zc == Approx(BZ).margin(0.05));

  // ...while the naive one is badly wrong. Asserted as a hard lower bound on the
  // error so that a future "optimisation" that silently reverts to the naive
  // integral cannot pass this test.
  CHECK(std::abs(m.z - BZ) > 2.0);

  // Transverse axes are untouched by the straddle and must still be exact.
  CHECK(m.xc == Approx(0.0).margin(0.02));
  CHECK(m.yc == Approx(0.0).margin(0.02));

  // The width is a property of the packet, not of where it sits: sigma_z_circ
  // must match the interior case even though the packet is split in two.
  CHECK(m.szc == Approx(SIGMA / std::sqrt(2.0)).epsilon(0.03));

  fs::remove_all(dir);
}

TEST_CASE("WPRealSpaceStats circular centroid: fully delocalised packet",
          "[observables][wp_real_space][circular][engine]") {
  auto electrons = make_electrons();
  REQUIRE(electrons.kpin()[0].set_comm().size() == 1);
  const int ist = electrons.states().num_states() - 1;

  fill_uniform(electrons, ist);

  fs::path dir = fs::temp_directory_path() / "inqkit_test_wpcirc_uniform";
  fs::remove_all(dir);
  inqkit::observables::WPRealSpaceStats stats((dir / "s.csv").string(), ist);
  auto m = stats.compute(electrons);

  // A uniform |psi|^2 has zero resultant: no preferred phase on any axis.
  CHECK(m.Rx == Approx(0.0).margin(1e-8));
  CHECK(m.Ry == Approx(0.0).margin(1e-8));
  CHECK(m.Rz == Approx(0.0).margin(1e-8));

  // sigma_circ must stay finite (the R floor) and report a width at least of
  // order the cell, rather than NaN or inf.
  CHECK(std::isfinite(m.sxc));
  CHECK(std::isfinite(m.szc));
  CHECK(m.szc > L_BOX / 4.0);

  fs::remove_all(dir);
}
