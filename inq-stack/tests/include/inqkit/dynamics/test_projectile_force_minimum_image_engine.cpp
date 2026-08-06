// Known-case test for the minimum-image option of
// inqkit::dynamics::projectile_force_z / drag_energy / projectile_force_axis,
// added 2026-08-01 for the annular-tube channeling twin
// (docs/plans/cylindrical-channeling-ks-stopping.md §4).
//
// THE BUG THIS PINS. `drag_energy` builds the projectile blob with
// `gaussian_density`, which uses a PLAIN Cartesian displacement. A projectile
// within a few sigma_pot of a cell face therefore loses the part of its charge
// that falls outside the grid, and — worse for a FORCE — loses it ASYMMETRICALLY
// between the +delta and -delta finite-difference evaluations. The spurious term
// does not cancel, so the run sees a fake force exactly at launch.
//
// The channeling twin launches a sigma_pot = 2.83 Gaussian 2 Bohr from the -z
// face of a periodic tube, i.e. squarely inside that regime, and its wavepacket
// counterpart wraps EXACTLY (a KS orbital lives on a 3-D FFT basis). Driving the
// classical twin with a clipped force while its perturbation uses the minimum
// image would make the two projectiles differ at precisely the boundary the
// study relies on being identical.
//
// EXACT EXPECTATION. Take phi_drag(z) = cos(2 pi z / L), a smooth periodic field.
// For a unit-normalised periodic Gaussian of std sigma the convolution is exactly
// its Fourier coefficient at k = 2 pi / L:
//
//   E_R(Z) = integral n_proj(r - Z zhat) cos(k z) dr = exp(-k^2 sigma^2 / 2) cos(k Z)
//   F_z(Z) = -dE_R/dZ                                = exp(-k^2 sigma^2 / 2) k sin(k Z)
//
// The grid sum equals the continuum integral to ~exp(-(2 pi/dx)^2 sigma^2/2) here
// (utterly negligible), and the symmetric finite difference contributes only a
// sinc(k delta) = 0.99994 factor, so this is a genuine closed-form target rather
// than a regression baseline.
//
// TWO POSITIONS, and the contrast between them is the whole test:
//   Z = 3        blob deep inside the cell -> plain and minimum-image AGREE.
//   Z = L/2 - 1  blob straddling the +z face -> only the minimum image is right.
//
// WHY THE BOX IS 24 BOHR AND NOT 16 (this test failed once, 2026-08-01).
// "Deep inside the cell" is a claim about DISTANCE IN SIGMA, not about looking
// far from the edge on a plot. The first version used L = 16, which puts Z = 3
// only 3.57 sigma from the +z face: the plain Gaussian loses 1.8e-4 of its
// charge there, and because the minimum-image version WRAPS that tail round to
// z ~ -8 where cos(2 pi z/L) = -1 instead of dropping it, the two kernels
// disagreed by 5.3e-4 in the energy and 1.5e-3 in the force — a thousand times
// the 1e-6 the test demanded. The library was right; the test's premise was not.
//
// L = 24 puts Z = 3 at 6.43 sigma, where the clipped tail is ~1e-10 and the two
// kernels agree to 9.7e-11 (energy) / 1.8e-9 (force). Derived by quadrature
// before this edit, not tuned until it passed. The straddling case is unaffected
// by the change and still exposes the clipped kernel by 328 %.
//
// The residual 2.9e-5 between the minimum-image force and the closed form is NOT
// slack: it is the finite-difference sinc error, (k*delta)^2/6 = 2.86e-5 for
// k = 2 pi/24 and delta = 0.05. That it lands exactly there is itself a check.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/dynamics/projectile_force.hpp>

#include <cmath>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {

constexpr double L_BOX   = 24.0;   // Bohr, cubic periodic (60^3 grid)
constexpr double SPACING = 0.4;    // Bohr
constexpr double SIGMA   = 1.4;    // sigma_pot of the classical charge blob

systems::electrons make_electrons() {
	systems::ions ions(systems::cell::cubic(L_BOX * 1.0_bohr).periodic());
	systems::electrons electrons(
		ions, options::electrons{}.spacing(SPACING * 1.0_bohr)
		                          .extra_electrons(2)
		                          .extra_states(2));
	ground_state::initial_guess(ions, electrons);
	return electrons;
}

// phi_drag(r) = cos(2 pi z / L) on the density basis.
basis::field<basis::real_space, double> cosine_field(systems::electrons const & electrons) {
	auto basis_ = electrons.density().basis();
	basis::field<basis::real_space, double> f(basis_);
	auto point_op = basis_.point_op();
	auto cub      = begin(f.cubic());
	const double k = 2.0 * M_PI / L_BOX;
	gpu::run(basis_.local_sizes()[2], basis_.local_sizes()[1], basis_.local_sizes()[0],
		[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
			auto r = point_op.rvector_cartesian(ix, iy, iz);
			cub[ix][iy][iz] = cos(k * r[2]);
		});
	return f;
}

const double K_WAVE  = 2.0 * M_PI / L_BOX;
const double DAMPING = std::exp(-0.5 * K_WAVE * K_WAVE * SIGMA * SIGMA);

double exact_energy(double Z) { return DAMPING * std::cos(K_WAVE * Z); }
double exact_force (double Z) { return DAMPING * K_WAVE * std::sin(K_WAVE * Z); }

} // namespace

TEST_CASE("projectile force: deep inside the cell, clipped and minimum-image agree",
          "[dynamics][projectile_force][engine]") {
	auto electrons = make_electrons();
	auto phi = cosine_field(electrons);

	// 5 Bohr from either face is 3.6 sigma — the plain Gaussian is not clipped, so
	// the two kernels must agree and both must hit the closed form.
	const vector3<double> centre{0.0, 0.0, 3.0};

	const double e_plain = inqkit::dynamics::drag_energy(phi, centre, SIGMA, false);
	const double e_mini  = inqkit::dynamics::drag_energy(phi, centre, SIGMA, true);
	CHECK(e_plain == Approx(e_mini).epsilon(1e-6));
	CHECK(e_mini  == Approx(exact_energy(3.0)).epsilon(1e-6));

	const double f_plain = inqkit::dynamics::projectile_force_z(phi, centre, SIGMA, 0.05, false);
	const double f_mini  = inqkit::dynamics::projectile_force_z(phi, centre, SIGMA, 0.05, true);
	CHECK(f_plain == Approx(f_mini).epsilon(1e-6));
	CHECK(f_mini  == Approx(exact_force(3.0)).epsilon(2e-4));
}

TEST_CASE("projectile force: straddling the face, only the minimum image is right",
          "[dynamics][projectile_force][engine]") {
	auto electrons = make_electrons();
	auto phi = cosine_field(electrons);

	// 1 Bohr inside the +z face = 0.71 sigma: the plain Gaussian keeps only
	// Phi(1/1.4) = 76 % of its charge.
	const double Z = 0.5 * L_BOX - 1.0;   // +7.0
	const vector3<double> centre{0.0, 0.0, Z};

	const double f_mini  = inqkit::dynamics::projectile_force_z(phi, centre, SIGMA, 0.05, true);
	const double f_plain = inqkit::dynamics::projectile_force_z(phi, centre, SIGMA, 0.05, false);

	CHECK(f_mini == Approx(exact_force(Z)).epsilon(2e-3));

	// And the clipped one is NOT — that failure is the reason the flag exists, so
	// it is asserted rather than merely noted. A future change that quietly makes
	// gaussian_density periodic would trip this and should be looked at.
	CHECK(std::abs(f_plain - exact_force(Z)) > 0.05 * std::abs(exact_force(Z)));
}

TEST_CASE("projectile force: the axis-general form reproduces the z form",
          "[dynamics][projectile_force][engine]") {
	auto electrons = make_electrons();
	auto phi = cosine_field(electrons);
	const vector3<double> centre{0.0, 0.0, 2.0};

	// axis = 2 must be byte-identical to the z convenience wrapper ...
	CHECK(inqkit::dynamics::projectile_force_axis(phi, centre, SIGMA, 2, 0.05, true)
	      == Approx(inqkit::dynamics::projectile_force_z(phi, centre, SIGMA, 0.05, true)));

	// ... and phi_drag has no x-dependence, so the transverse force is zero.
	// This is the symmetry the on-axis channeling projectile relies on: with an
	// axially symmetric tube the transverse force vanishes and the projectile
	// stays on the axis without being constrained there.
	CHECK(inqkit::dynamics::projectile_force_axis(phi, centre, SIGMA, 0, 0.05, true)
	      == Approx(0.0).margin(1e-10));
}
