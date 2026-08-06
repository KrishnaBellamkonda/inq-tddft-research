// Known-case test for inqkit::fields::from_inq_field, added 2026-08-01.
//
// from_inq_field exports an ARBITRARY INQ real-space field as a RealField3D so it
// can be written as a VTI. The one thing it must get right is the ORDERING
// CONVENTION: INQ stores real-space fields FFT-naturally (index 0 = the origin,
// the negative half wrapped to the top) while RealField3DWriter stamps
// Origin = -L/2 and expects PHYSICAL order (index 0 = the left edge).
//
// Getting that wrong swaps the centre and the edges and produces a picture that
// looks completely plausible — the slab/cluster appears at the box faces with
// vacuum in the middle. It is not caught by eye, which is why
// .claude/rules/vti-coordinate-mapping.md exists and why this is pinned here.
//
// THE DECISIVE TEST is the last one: from_inq_field(electrons.density()) must be
// ELEMENT-WISE IDENTICAL to fields::density::total(electrons). The two are
// independent implementations of the same convention (density::total predates
// this function and is exercised by every production run ever done), so equality
// pins the new code against the proven one rather than against my own re-reading
// of the convention.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/inq_field.hpp>

#include <cmath>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {

constexpr double L_BOX   = 16.0;   // Bohr, cubic periodic
constexpr double SPACING = 0.4;    // Bohr -> 40 points per side

systems::electrons make_electrons() {
	systems::ions ions(systems::cell::cubic(L_BOX * 1.0_bohr).periodic());
	systems::electrons electrons(
		ions, options::electrons{}.spacing(SPACING * 1.0_bohr)
		                          .extra_electrons(2)
		                          .extra_states(2));
	ground_state::initial_guess(ions, electrons);
	return electrons;
}

// A Gaussian bump centred at `c`, built directly on the density basis.
basis::field<basis::real_space, double>
bump(systems::electrons const & electrons, inq::vector3<double> c, double sig) {
	auto b = electrons.density().basis();
	basis::field<basis::real_space, double> f(b);
	auto point_op = b.point_op();
	auto cub      = begin(f.cubic());
	const double cx = c[0], cy = c[1], cz = c[2], inv2s2 = 1.0/(2.0*sig*sig);
	gpu::run(b.local_sizes()[2], b.local_sizes()[1], b.local_sizes()[0],
		[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
			auto r = point_op.rvector_cartesian(ix, iy, iz);
			const double dx = r[0]-cx, dy = r[1]-cy, dz = r[2]-cz;
			cub[ix][iy][iz] = exp(-(dx*dx + dy*dy + dz*dz)*inv2s2);
		});
	return f;
}

// Index of the largest value, as (ix, iy, iz).
std::array<int,3> argmax(inqkit::fields::RealField3D const & f) {
	std::size_t best = 0;
	for(std::size_t i = 1; i < f.values.size(); ++i)
		if(f.values[i] > f.values[best]) best = i;
	const int iz = int(best % f.nz);
	const int iy = int((best / f.nz) % f.ny);
	const int ix = int(best / (std::size_t(f.ny) * f.nz));
	return {ix, iy, iz};
}

} // namespace

TEST_CASE("from_inq_field: geometry metadata matches the basis",
          "[fields][inq_field][engine]") {
	auto electrons = make_electrons();
	auto f = inqkit::fields::from_inq_field(bump(electrons, {0.0,0.0,0.0}, 1.5));

	CHECK(f.nx == 40); CHECK(f.ny == 40); CHECK(f.nz == 40);
	CHECK(f.dx_bohr == Approx(SPACING));
	CHECK(f.dz_bohr == Approx(SPACING));
	// Origin at the LEFT EDGE, not at the centre — this is what makes index 0
	// mean -L/2 for the VTI writer.
	CHECK(f.origin_x_bohr == Approx(-L_BOX/2.0).margin(SPACING));
	CHECK(f.origin_z_bohr == Approx(-L_BOX/2.0).margin(SPACING));
	CHECK(f.values.size() == std::size_t(40)*40*40);
}

TEST_CASE("from_inq_field: a feature at the box CENTRE lands at the MIDDLE index",
          "[fields][inq_field][engine]") {
	auto electrons = make_electrons();
	auto f = inqkit::fields::from_inq_field(bump(electrons, {0.0,0.0,0.0}, 1.0));

	// Without the FFT shift the peak would sit at index (0,0,0) — the classic
	// "the slab is at the edges and there is vacuum in the middle" picture.
	auto m = argmax(f);
	CHECK(m[0] == 20); CHECK(m[1] == 20); CHECK(m[2] == 20);
	CHECK(m[2] != 0);
}

TEST_CASE("from_inq_field: an OFF-CENTRE feature lands at its physical index",
          "[fields][inq_field][engine]") {
	auto electrons = make_electrons();
	// z = +4 Bohr -> index (4 - (-8))/0.4 = 30
	auto f = inqkit::fields::from_inq_field(bump(electrons, {0.0,0.0,4.0}, 1.0));

	auto m = argmax(f);
	CHECK(m[0] == 20); CHECK(m[1] == 20);
	CHECK(m[2] == 30);
	// And the reconstructed coordinate of that index is the physical position.
	CHECK(f.origin_z_bohr + m[2]*f.dz_bohr == Approx(4.0).margin(SPACING));
}

TEST_CASE("from_inq_field: identical to density::total on the same field",
          "[fields][inq_field][engine]") {
	auto electrons = make_electrons();

	// density::total(electrons) and from_inq_field(electrons.density()) are two
	// independent implementations of the same convention. density::total is
	// exercised by every production run in the project, so equality here pins the
	// new function against proven code rather than against a re-reading of the
	// convention. Bit-identical is the right expectation: same source array, same
	// shift, same flattening — only the entry point differs.
	auto a = inqkit::fields::density::total(electrons);
	auto b = inqkit::fields::from_inq_field(electrons.density());

	REQUIRE(a.nx == b.nx); REQUIRE(a.ny == b.ny); REQUIRE(a.nz == b.nz);
	REQUIRE(a.values.size() == b.values.size());
	CHECK(a.origin_x_bohr == Approx(b.origin_x_bohr));
	CHECK(a.origin_y_bohr == Approx(b.origin_y_bohr));
	CHECK(a.origin_z_bohr == Approx(b.origin_z_bohr));
	CHECK(a.dx_bohr == Approx(b.dx_bohr));

	double max_abs_diff = 0.0;
	for(std::size_t i = 0; i < a.values.size(); ++i)
		max_abs_diff = std::max(max_abs_diff, std::abs(a.values[i] - b.values[i]));
	CHECK(max_abs_diff == 0.0);
}
