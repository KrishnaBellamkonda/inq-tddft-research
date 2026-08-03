// Engine-tier T0 test for the localised jellium background.
//
// Validates the n₊ BUILDER and the perturbation's sign, against analytic
// targets fixed BEFORE running (no retrofitting):
//   T0.1  slab   ∫n₊ = n₀·(2a)·Lx·Ly
//   T0.2  sphere ∫n₊ = n₀·(4/3)πR³
//   T0.3  the perturbation produces an ATTRACTIVE well: ∫ v_bg·n₊ < 0
//
// Tolerances are set by SHARP-EDGE GRID DISCRETISATION, not tuned to output:
// a Heaviside boundary is resolved to ±½ grid cell per face, so the relative
// error on the integral is ~ spacing/extent. With spacing 0.25 and slab
// thickness 6 that is ~4%; the sphere (curved surface) is looser. These are
// the honest physical tolerances; the v_bg shape itself is validated
// physically at T1 (electrons bind inside the region — VC-3).

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/analytics.hpp>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;
using namespace inqkit::jellium;

namespace {
auto make_basis_electrons(double L) {
	systems::ions ions(systems::cell::cubic(L * 1.0_bohr).periodic());
	systems::electrons electrons(ions,
		options::electrons{}.spacing(0.25_bohr).extra_electrons(2));
	return electrons;
}
}

TEST_CASE("T0.1 localised background slab integrates to N", "[jellium][localised][engine]") {
	const double L = 12.0;
	auto electrons = make_basis_electrons(L);
	auto basis = electrons.density_basis();

	const double n0 = 0.01, a = 3.0;
	localised_background_params p;
	p.shape = background_shape::slab;
	p.n0 = n0; p.half_width = a; p.slab_axis = 2;

	auto nplus = make_localised_background(basis, p);
	const double integ  = operations::integral(nplus);
	const double expect = n0 * (2.0 * a) * L * L;   // = 8.64

	CHECK(integ == Approx(expect).epsilon(0.05));
}

TEST_CASE("T0.2 localised background sphere integrates to N", "[jellium][localised][engine]") {
	const double L = 14.0;
	auto electrons = make_basis_electrons(L);
	auto basis = electrons.density_basis();

	const double n0 = 0.01, R = 4.0;
	localised_background_params p;
	p.shape = background_shape::sphere;
	p.n0 = n0; p.half_width = R;

	auto nplus = make_localised_background(basis, p);
	const double integ  = operations::integral(nplus);
	const double expect = n0 * (4.0 / 3.0) * M_PI * R * R * R; // = 2.681

	CHECK(integ == Approx(expect).epsilon(0.08));
}

TEST_CASE("T0.3 background perturbation well is attractive", "[jellium][localised][engine]") {
	const double L = 12.0;
	auto electrons = make_basis_electrons(L);
	auto basis = electrons.density_basis();

	localised_background_params p;
	p.shape = background_shape::slab;
	p.n0 = 0.01; p.half_width = 3.0; p.slab_axis = 2;

	auto nplus = make_localised_background(basis, p);

	basis::field<basis::real_space, double> v(basis);
	v.fill(0.0);

	localised_background_perturbation pert(p);
	pert.potential(0.0, v);   // v += v_bg = −poisson(n₊)

	// Electron–background interaction energy ∫ v_bg·n₊ must be NEGATIVE.
	const double e_ebg = operations::integral_product(v, nplus);
	CHECK(e_ebg < 0.0);
}

// --- annulus (hollow periodic tube) -----------------------------------------
// The annular shape underlies the cylindrical-jellium campaign. The erfc edge
// profile and axial uniformity are proven ANALYTICALLY by the formula-validation
// agent (CONFIRM, 2026-06-28); these engine tests validate the grid BUILDER via
// decomposition-safe integrals only (no point access). Tube axis = z (slab_axis=2).
namespace {
// ∫n₊ for an annulus(R_in,R_out) with sharp edges (w=0), tube ∥ z, on an L³ cell.
double annulus_integral(double L, double n0, double R_in, double R_out, double w) {
	auto electrons = make_basis_electrons(L);
	auto basis = electrons.density_basis();
	localised_background_params p;
	p.shape = background_shape::annulus;
	p.n0 = n0; p.half_width = R_out; p.inner_radius = R_in;
	p.slab_axis = 2; p.edge_width = w;
	return operations::integral(make_localised_background(basis, p));
}
}

TEST_CASE("T0.4 annulus integrates to n0·π(Rout²−Rin²)·Lz", "[jellium][localised][annulus][engine]") {
	const double L = 20.0, n0 = 0.01, R_in = 3.0, R_out = 7.0;
	const double integ  = annulus_integral(L, n0, R_in, R_out, /*w=*/0.0);
	const double expect = n0 * M_PI * (R_out*R_out - R_in*R_in) * L;  // = 25.13
	// Two curved (cylindrical) radial surfaces → sphere-like discretisation error.
	CHECK(integ == Approx(expect).epsilon(0.08));
}

TEST_CASE("T0.5 annulus = outer cylinder − inner cylinder (bore carved correctly)",
          "[jellium][localised][annulus][engine]") {
	// A solid cylinder of radius R is annulus(R_in=0, w=0) — the inner complement
	// becomes 1 everywhere (background_mask(d,0,0)=0 for physical d≥0). Then the
	// continuum identity ∫annulus(Rin,Rout) + ∫cyl(0,Rin) = ∫cyl(0,Rout) must hold;
	// the shared R_in surface is discretised identically in both terms and cancels.
	//
	// NOTE this identity is used here only as a SHARP-EDGE (w=0) algebraic check.
	// Production filled tubes must use background_shape::cylinder — at w > 0 the
	// R_in = 0 annulus is a trap (n₀/2 on the axis); see T0.7.
	const double L = 20.0, n0 = 0.01, R_in = 3.0, R_out = 7.0;
	const double ann   = annulus_integral(L, n0, R_in, R_out, 0.0);
	const double inner = annulus_integral(L, n0, 0.0,  R_in,  0.0);
	const double outer = annulus_integral(L, n0, 0.0,  R_out, 0.0);
	CHECK((ann + inner) == Approx(outer).epsilon(0.05));
}

TEST_CASE("T0.6 erfc-smoothed (production w=1) annulus conserves charge",
          "[jellium][localised][annulus][engine]") {
	// Production uses edge_width w≈1 on a thick wall (R_out−R_in=4 ≫ w). The two
	// erfc edges shift charge only by a few % (partial inner/outer cancellation),
	// so the integral stays near the sharp-edge target. Caller rescales n0 for
	// EXACT neutrality in production; here we check the builder is well-behaved.
	const double L = 20.0, n0 = 0.01, R_in = 3.0, R_out = 7.0;
	const double integ  = annulus_integral(L, n0, R_in, R_out, /*w=*/1.0);
	const double expect = n0 * M_PI * (R_out*R_out - R_in*R_in) * L;
	CHECK(integ == Approx(expect).epsilon(0.10));
}

// --- filled tube: background_shape::cylinder ---------------------------------
// Added while designing the R_in → 0 proximity ladder (2026-08-02), which needs a
// FILLED r_s=3 tube as its final rung.
//
// The filled tube is its own shape rather than annulus(R_in = 0). Reason: the erfc
// step is centred ON its nominal edge, so background_mask(0,0,w) = ½ for every
// w > 0 — an annulus with a degenerate inner edge yields n₊ = n₀/2 EXACTLY ON THE
// TUBE AXIS, relaxing to n₀ only by d ≈ 2w. That is precisely where a channeling
// projectile flies, so the error would be silent and maximal at the same point.
//
// T0.5 above cannot catch that: it probes R_in = 0 at w = 0, where
// background_mask(d,0,0) = 0 for all physical d ≥ 0 and the composition is
// accidentally right. Only the SOFTENED branch — i.e. only production, w = 0.5 —
// is affected, so the shape check below is done at w > 0.
TEST_CASE("T0.7 filled cylinder carries full n0 ON the axis (w>0)",
          "[jellium][localised][cylinder]") {
	const double R = 7.0, w = 0.5;

	// the whole point: n₀, not n₀/2, at d = 0
	CHECK(cylinder_mask(0.0, R, w) == Approx(1.0).margin(1e-12));
	// ...and the value that would have appeared there via a degenerate inner edge
	CHECK(background_mask(0.0, 0.0, w) == Approx(0.5).margin(1e-12));

	// uniform across the whole interior
	for(double d : {0.0, 0.25, 0.5, 1.0, 2.0, 5.0})
		CHECK(cylinder_mask(d, R, w) == Approx(1.0).epsilon(1e-6));

	// one radial boundary, erfc-softened: ½ on it, → 0 outside
	CHECK(cylinder_mask(R, R, w) == Approx(0.5).margin(1e-12));
	CHECK(cylinder_mask(R + 4.0*w, R, w) < 1e-6);
}

TEST_CASE("T0.8 hollow tube still carves its bore (w>0)",
          "[jellium][localised][annulus]") {
	// The new shape must not disturb the hollow production geometry (R_in=10,R_out=14).
	const double R_in = 10.0, R_out = 14.0, w = 0.5;
	CHECK(annulus_mask(0.0,   R_out, R_in, w) < 1e-12);                     // axis EMPTY
	CHECK(annulus_mask(R_in,  R_out, R_in, w) == Approx(0.5).margin(1e-9)); // on bore edge
	CHECK(annulus_mask(12.0,  R_out, R_in, w) == Approx(1.0).epsilon(1e-6));// mid-wall
	CHECK(annulus_mask(R_out, R_out, R_in, w) == Approx(0.5).margin(1e-9)); // on outer edge
}

TEST_CASE("T0.9 cylinder integrates to n0·πR²·Lz (builder branch)",
          "[jellium][localised][cylinder][engine]") {
	// Exercises the shape through make_localised_background, not just the mask:
	// proves the enum reaches the right branch and reads half_width as the RADIUS.
	const double L = 20.0, n0 = 0.01, R = 7.0;
	auto electrons = make_basis_electrons(L);
	localised_background_params p;
	p.shape = background_shape::cylinder;
	p.n0 = n0; p.half_width = R; p.slab_axis = 2; p.edge_width = 0.0;
	const double integ  = operations::integral(make_localised_background(electrons.density_basis(), p));
	const double expect = n0 * M_PI * R * R * L;                 // = 30.79
	// one curved radial surface → same discretisation class as T0.4
	CHECK(integ == Approx(expect).epsilon(0.08));

	// and it must NOT depend on inner_radius, which the shape ignores
	p.inner_radius = 3.0;
	CHECK(operations::integral(make_localised_background(electrons.density_basis(), p))
	      == Approx(integ).epsilon(1e-12));
}
