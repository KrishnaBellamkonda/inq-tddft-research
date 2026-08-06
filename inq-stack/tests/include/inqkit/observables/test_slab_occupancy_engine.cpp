// Known-case test for inqkit::observables::slab_occupancy, added 2026-07-31 for
// the wrap-around slab KS-stopping study
// (docs/plans/slab-ks-orbital-stopping-wrap.md §3, Window B).
//
// f(t) = (norm of the orbital inside the slab band) / (total orbital norm) is
// what converts a centroid path into an IN-SLAB path, so that -dT/ds_slab is a
// stopping power per Bohr of medium traversed even after the wavepacket has
// spread wider than the slab. If f is wrong, every S in the study is wrong by
// the same factor, so it is pinned against closed-form values here.
//
// Expectations are analytic:
//
//   uniform orbital                 f = 2h / L                    (filling factor)
//   narrow Gaussian at band centre  f = erf(h / (sigma_d sqrt2))  (~1)
//   full-cell band                  f = 1                          (exactly)
//   Gaussian across the face, band  f = Phi((h-mu)/s) + Phi((h+mu)/s) - 1
//   across the face too                with mu the MINIMUM-IMAGE offset
//
// The last case is the decisive one. A packet centred at z = -7.5 and a band
// centred at z = +7.5 are 15 Bohr apart in raw Cartesian terms but only 1 Bohr
// apart in the periodic cell. Without the minimum image f collapses to ~0; with
// it, f = 0.998. That is exactly the situation every wrap of the production run
// puts the packet in.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/slab_occupancy.hpp>

#include <cmath>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {

constexpr double L_BOX   = 16.0;   // Bohr, cubic periodic
constexpr double SPACING = 0.4;    // Bohr
constexpr double SIGMA   = 1.5;    // psi-width; density std = SIGMA/sqrt(2)

systems::electrons make_electrons() {
	systems::ions ions(systems::cell::cubic(L_BOX * 1.0_bohr).periodic());
	systems::electrons electrons(
		ions, options::electrons{}.spacing(SPACING * 1.0_bohr)
		                          .extra_electrons(2)
		                          .extra_states(2));
	ground_state::initial_guess(ions, electrons);
	return electrons;
}

// Minimum-image Gaussian into state `ist`. Normalisation is irrelevant:
// slab_occupancy divides by the orbital norm.
void fill_periodic_gaussian(systems::electrons & electrons, int ist,
                            double bx, double by, double bz, double sig) {
	auto & phi     = electrons.kpin()[0];
	auto   basis   = phi.basis();
	auto   phicub  = begin(phi.hypercubic());
	auto   point_op = basis.point_op();
	auto   sizes   = basis.local_sizes();

	gpu::run(sizes[2], sizes[1], sizes[0],
		[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
			auto r = point_op.rvector_cartesian(ix, iy, iz);
			auto d = inq::vector3<double, inq::cartesian>{r[0]-bx, r[1]-by, r[2]-bz};
			auto f = point_op.cell().to_contravariant(d);
			for(int k = 0; k < 3; ++k) f[k] -= floor(f[k] + 0.5);
			auto dmin = point_op.cell().to_cartesian(f);
			const double r2 = dmin[0]*dmin[0] + dmin[1]*dmin[1] + dmin[2]*dmin[2];
			phicub[ix][iy][iz][ist] = complex(exp(-r2 / (2.0*sig*sig)), 0.0);
		});
}

void fill_uniform(systems::electrons & electrons, int ist) {
	auto & phi    = electrons.kpin()[0];
	auto   phicub = begin(phi.hypercubic());
	auto   sizes  = phi.basis().local_sizes();
	gpu::run(sizes[2], sizes[1], sizes[0],
		[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
			phicub[ix][iy][iz][ist] = complex(1.0, 0.0);
		});
}

double Phi(double x) { return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0))); }

} // namespace

TEST_CASE("slab_occupancy: a uniform orbital returns the geometric filling factor",
          "[observables][slab_occupancy][engine]") {
	auto electrons = make_electrons();
	REQUIRE(electrons.kpin()[0].set_comm().size() == 1);
	const int ist = electrons.states().num_states() - 1;
	fill_uniform(electrons, ist);

	// This is the delocalised limit the production runs end up in: the packet
	// fills the box, so the drag it feels is diluted by exactly 2h/L.
	for(double h : {2.0, 4.0, 6.0}) {
		auto occ = inqkit::observables::slab_occupancy(
			electrons, ist, {.axis = 2, .center = 0.0, .half_width = h});
		// Margin = one grid layer / L: the discrete sum can include or exclude a
		// single plane of points at the band edge.
		CHECK(occ.fraction == Approx(2.0 * h / L_BOX).margin(SPACING / L_BOX));
	}
}

TEST_CASE("slab_occupancy: a band covering the whole cell returns exactly 1",
          "[observables][slab_occupancy][engine]") {
	auto electrons = make_electrons();
	const int ist = electrons.states().num_states() - 1;
	fill_periodic_gaussian(electrons, ist, 0.0, 0.0, 2.0, SIGMA);

	auto occ = inqkit::observables::slab_occupancy(
		electrons, ist, {.axis = 2, .center = 0.0, .half_width = 0.5 * L_BOX});
	CHECK(occ.fraction == Approx(1.0).epsilon(1e-12));
	CHECK(occ.norm_in  == Approx(occ.norm_total).epsilon(1e-12));
}

TEST_CASE("slab_occupancy: a narrow packet at the band centre is fully inside",
          "[observables][slab_occupancy][engine]") {
	auto electrons = make_electrons();
	const int ist = electrons.states().num_states() - 1;
	fill_periodic_gaussian(electrons, ist, 0.0, 0.0, 0.0, SIGMA);

	const double sigma_d = SIGMA / std::sqrt(2.0);   // density std = 1.0607
	const double h = 4.0;
	auto occ = inqkit::observables::slab_occupancy(
		electrons, ist, {.axis = 2, .center = 0.0, .half_width = h});

	// 1-D marginal of the density is N(0, sigma_d): f = erf(h / (sigma_d sqrt2)).
	const double expected = std::erf(h / (sigma_d * std::sqrt(2.0)));
	CHECK(occ.fraction == Approx(expected).margin(2e-3));
	CHECK(occ.fraction > 0.99);
}

TEST_CASE("slab_occupancy: packet and band on opposite faces — minimum image decides",
          "[observables][slab_occupancy][engine]") {
	auto electrons = make_electrons();
	const int ist = electrons.states().num_states() - 1;

	const double b_packet = -0.5 * L_BOX + 0.5;   // -7.5
	const double c_band   =  0.5 * L_BOX - 0.5;   // +7.5
	fill_periodic_gaussian(electrons, ist, 0.0, 0.0, b_packet, SIGMA);

	const double h = 4.0;
	auto occ = inqkit::observables::slab_occupancy(
		electrons, ist, {.axis = 2, .center = c_band, .half_width = h});

	// Raw separation is 15 Bohr; minimum-image separation is 1 Bohr.
	const double mu = 1.0, s = SIGMA / std::sqrt(2.0);
	const double expected = Phi((h - mu)/s) + Phi((h + mu)/s) - 1.0;   // 0.9977

	CHECK(occ.fraction == Approx(expected).margin(3e-3));
	CHECK(occ.fraction > 0.99);   // a non-periodic implementation would give ~0
}
