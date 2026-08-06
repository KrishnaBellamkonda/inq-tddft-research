// Known-case test for inqkit::jellium::gaussian_density_minimum_image, added
// 2026-07-31 for the wrap-around slab KS-stopping twin
// (docs/plans/slab-ks-orbital-stopping-wrap.md §4).
//
// WHAT IS BEING PINNED DOWN
// A classical Gaussian projectile built from a PLAIN Cartesian displacement is
// CLIPPED by the cell face: the half of the blob that falls outside the grid is
// simply absent, so its integral drops below 1 and the force it exerts is wrong
// for as long as it straddles. A KS orbital has no such problem — the
// wavefunction basis is a 3-D FFT, periodic in all directions, so a wavepacket
// wraps exactly. In a wrap-around classical-vs-wavepacket twin that difference
// lands precisely on the boundary the study introduces on purpose.
//
// The expectations below are analytic, not read back from the implementation:
//
//   interior blob   : both kernels give integral 1 and agree pointwise
//   blob at +z face : minimum image keeps integral 1;
//                     plain keeps only Phi((z_edge - b)/sigma), where z_edge is
//                     the LAST GRID POINT's half-cell edge, L/2 - dx/2 = 7.8 --
//                     NOT L/2. INQ's nodes start at exactly -L/2 (the
//                     vti-coordinate-mapping convention), so the domain's
//                     effective upper edge is half a cell short of the face.
//                     (= 0.46 for b = 7.9, sigma = 1.0, L = 16, dx = 0.4;
//                      assuming L/2 predicts 0.54 and is wrong by 17 %.)

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>

#include <cmath>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {

constexpr double L_BOX   = 16.0;   // Bohr, cubic periodic
constexpr double SPACING = 0.4;    // Bohr
constexpr double SIGMA   = 1.0;    // charge std

systems::electrons make_electrons() {
	systems::ions ions(systems::cell::cubic(L_BOX * 1.0_bohr).periodic());
	systems::electrons electrons(
		ions, options::electrons{}.spacing(SPACING * 1.0_bohr)
		                          .extra_electrons(2)
		                          .extra_states(2));
	ground_state::initial_guess(ions, electrons);
	return electrons;
}

// Standard normal CDF, for the analytic clipped fraction.
double Phi(double x) { return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0))); }

} // namespace

TEST_CASE("gaussian_density_minimum_image: interior blob matches the plain kernel",
          "[jellium][gaussian][minimum_image][engine]") {
	auto electrons = make_electrons();
	auto basis = electrons.density().basis();

	const vector3<double> centre{0.3, -0.7, 0.5};   // well away from every face

	auto n_plain = inqkit::jellium::gaussian_density(basis, centre, SIGMA);
	auto n_mimg  = inqkit::jellium::gaussian_density_minimum_image(basis, centre, SIGMA);

	// Both integrate to 1 (the Gaussian is normalised and fully contained).
	CHECK(operations::integral(n_plain) == Approx(1.0).epsilon(1e-6));
	CHECK(operations::integral(n_mimg)  == Approx(1.0).epsilon(1e-6));

	// And they are the SAME field: the minimum image is a no-op this far from
	// any face, so the difference must be at round-off.
	auto diff = operations::integral_product(n_plain, n_mimg);
	auto self = operations::integral_product(n_mimg,  n_mimg);
	CHECK(diff == Approx(self).epsilon(1e-10));
}

TEST_CASE("gaussian_density_minimum_image: blob on the +z face keeps its charge",
          "[jellium][gaussian][minimum_image][engine]") {
	auto electrons = make_electrons();
	auto basis = electrons.density().basis();

	const double bz = 0.5 * L_BOX - 0.1;            // 7.9: 0.1 Bohr inside the face
	const vector3<double> centre{0.0, 0.0, bz};

	auto n_plain = inqkit::jellium::gaussian_density(basis, centre, SIGMA);
	auto n_mimg  = inqkit::jellium::gaussian_density_minimum_image(basis, centre, SIGMA);

	const double q_plain = operations::integral(n_plain);
	const double q_mimg  = operations::integral(n_mimg);

	// The wrapped charge is conserved — this is the whole point of the function.
	CHECK(q_mimg == Approx(1.0).epsilon(1e-4));

	// The plain kernel loses everything past the last grid point. The effective
	// upper edge of the domain is NOT L/2: INQ's nodes start at exactly -L/2
	// (the same convention the vti-coordinate-mapping rule is about), so the
	// LAST node sits at L/2 - dx and the half-cell it represents ends at
	//     z_edge = L/2 - dx/2 = 8.0 - 0.2 = 7.8 Bohr.
	// Analytic survivor: Phi((z_edge - bz)/sigma) = Phi(-0.1) = 0.4602.
	// (Using L/2 here predicts 0.5398 and is wrong by 17 % — measured 0.4599.)
	const double z_edge = 0.5 * L_BOX - 0.5 * SPACING;
	const double expected_clipped = Phi((z_edge - bz) / SIGMA);
	CHECK(q_plain == Approx(expected_clipped).epsilon(0.02));
	CHECK(q_plain < 0.6);                            // decisively short of 1

	// State the failure the study would otherwise have suffered: the plain kernel
	// loses 54 % of the projectile charge as it crosses the face, which would read
	// as a spurious collapse of the drag at every wrap.
	CHECK((q_mimg - q_plain) > 0.4);
}
