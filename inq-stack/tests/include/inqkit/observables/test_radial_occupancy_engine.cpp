// Known-case test for inqkit::observables::radial_occupancy, added 2026-08-01 for
// the annular-tube channeling KS-stopping study
// (docs/plans/cylindrical-channeling-ks-stopping.md §3).
//
// f_bore(t) — the fraction of the wavepacket still inside the hollow bore — is
// the observable that decides whether the study's central premise holds. The aim
// is that a CHANNELING packet reproduces the classical stopping power because it
// barely touches the wall; if f_bore is wrong the "channeling" claim is
// unfalsifiable and the fit window is picked blind. So it is pinned here against
// closed-form values.
//
// EXPECTATIONS ARE ANALYTIC. For psi = exp(-r^2/(2 sigma^2)) the density
// |psi|^2 is a Gaussian of per-axis std sigma_d = sigma/sqrt2, so the transverse
// radius r_perp = sqrt(dx^2 + dy^2) of a packet centred ON the axis is RAYLEIGH
// distributed with scale sigma_d:
//
//   P(r_perp <  R)   = 1 - exp(-R^2 / (2 sigma_d^2))
//   <r_perp>         = sigma_d sqrt(pi/2)
//   <r_perp^2>       = 2 sigma_d^2
//   sigma_r          = sigma_d sqrt(2 - pi/2)
//
// THE DECISIVE CASE is the last one. The tube axis is placed 2 Bohr from the +x
// face, so the packet's own tail wraps around to the far side of the cell. With
// the minimum image those wrapped grid points sit at r_perp ~ 2.5 and are counted
// inside the bore; without it they sit at r_perp ~ 13.5 and are counted outside,
// and f_bore collapses from 0.98 to ~0.6. That is exactly the situation the
// production run is in for its whole duration, because the projectile is
// launched 2 Bohr from the -z face of a periodic tube.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include <inq/inq.hpp>
#include <inqkit/observables/radial_occupancy.hpp>

#include <cmath>

using namespace inq;
using namespace inq::magnitude;
using Catch::Approx;

namespace {

constexpr double L_BOX   = 16.0;   // Bohr, cubic periodic
constexpr double SPACING = 0.4;    // Bohr
constexpr double SIGMA   = 1.5;    // psi-width; density std = SIGMA/sqrt2 = 1.06066

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
// radial_occupancy divides every quantity by the orbital norm.
void fill_periodic_gaussian(systems::electrons & electrons, int ist,
                            double bx, double by, double bz, double sig) {
	auto & phi      = electrons.kpin()[0];
	auto   basis    = phi.basis();
	auto   phicub   = begin(phi.hypercubic());
	auto   point_op = basis.point_op();
	auto   sizes    = basis.local_sizes();

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

constexpr double SIGMA_D = SIGMA / 1.4142135623730951;   // 1.0606601717798212

double rayleigh_cdf(double R) { return 1.0 - std::exp(-R*R / (2.0*SIGMA_D*SIGMA_D)); }

} // namespace

TEST_CASE("radial_occupancy: the three shells partition the orbital exactly",
          "[observables][radial_occupancy][engine]") {
	auto electrons = make_electrons();
	REQUIRE(electrons.kpin()[0].set_comm().size() == 1);
	const int ist = electrons.states().num_states() - 1;
	fill_periodic_gaussian(electrons, ist, 0.0, 0.0, 1.0, SIGMA);

	auto occ = inqkit::observables::radial_occupancy(
		electrons, ist, {.axis = 2, .center = {0.0, 0.0, 0.0},
		                 .r_inner = 3.0, .r_outer = 5.0});

	// Partition identity: no grid point may be double-counted or dropped. This is
	// exact arithmetic, not physics, so it holds to machine precision.
	CHECK(occ.norm_bore + occ.norm_wall + occ.norm_outside
	      == Approx(occ.norm_total).epsilon(1e-12));
	CHECK(occ.f_bore + occ.f_wall + occ.f_outside == Approx(1.0).epsilon(1e-12));
	CHECK(occ.norm_total > 0.0);
}

TEST_CASE("radial_occupancy: an on-axis packet follows the Rayleigh law",
          "[observables][radial_occupancy][engine]") {
	auto electrons = make_electrons();
	const int ist = electrons.states().num_states() - 1;
	fill_periodic_gaussian(electrons, ist, 0.0, 0.0, 0.0, SIGMA);

	auto occ = inqkit::observables::radial_occupancy(
		electrons, ist, {.axis = 2, .center = {0.0, 0.0, 0.0},
		                 .r_inner = 3.0, .r_outer = 5.0});

	// f_bore = 1 - exp(-4) = 0.98168 ; f_bore+f_wall = 1 - exp(-11.11) = 0.99999
	CHECK(occ.f_bore == Approx(rayleigh_cdf(3.0)).margin(3e-3));
	CHECK(occ.f_bore + occ.f_wall == Approx(rayleigh_cdf(5.0)).margin(3e-3));

	// Moments of the Rayleigh distribution. These are smooth integrals (no sharp
	// indicator), so they are tighter than the shell fractions.
	CHECK(occ.r_mean  == Approx(SIGMA_D * std::sqrt(M_PI/2.0)).epsilon(5e-3));
	CHECK(occ.r2_mean == Approx(2.0 * SIGMA_D * SIGMA_D).epsilon(5e-3));
	CHECK(occ.sigma_r == Approx(SIGMA_D * std::sqrt(2.0 - M_PI/2.0)).epsilon(1e-2));
}

TEST_CASE("radial_occupancy: r_inner = 0 leaves an empty bore",
          "[observables][radial_occupancy][engine]") {
	auto electrons = make_electrons();
	const int ist = electrons.states().num_states() - 1;
	fill_periodic_gaussian(electrons, ist, 0.0, 0.0, 0.0, SIGMA);

	auto occ = inqkit::observables::radial_occupancy(
		electrons, ist, {.axis = 2, .center = {0.0, 0.0, 0.0},
		                 .r_inner = 0.0, .r_outer = 5.0});
	CHECK(occ.norm_bore == Approx(0.0).margin(1e-14));
	CHECK(occ.f_bore    == Approx(0.0).margin(1e-14));
	CHECK(occ.f_wall    == Approx(rayleigh_cdf(5.0)).margin(3e-3));
}

TEST_CASE("radial_occupancy: a packet whose tail wraps the transverse face — "
          "minimum image decides", "[observables][radial_occupancy][engine]") {
	auto electrons = make_electrons();
	const int ist = electrons.states().num_states() - 1;

	// Tube axis 2 Bohr inside the +x face (x = 6 of a [-8, 8) box), packet on it.
	// Its tail therefore lives at x > 8, i.e. wrapped round to x < -6.
	const double axis_x = 0.5 * L_BOX - 2.0;   // +6.0
	fill_periodic_gaussian(electrons, ist, axis_x, 0.0, 0.0, SIGMA);

	auto occ = inqkit::observables::radial_occupancy(
		electrons, ist, {.axis = 2, .center = {axis_x, 0.0, 0.0},
		                 .r_inner = 3.0, .r_outer = 5.0});

	// Packet and axis coincide, so the Rayleigh law must hold EXACTLY as in the
	// on-axis case — but only if the wrapped tail is folded back. A non-periodic
	// implementation puts that tail at r_perp ~ 13 and reports f_bore ~ 0.6.
	CHECK(occ.f_bore == Approx(rayleigh_cdf(3.0)).margin(3e-3));
	CHECK(occ.f_bore > 0.97);
	CHECK(occ.r_mean == Approx(SIGMA_D * std::sqrt(M_PI/2.0)).epsilon(5e-3));
}

TEST_CASE("radial_occupancy: an off-axis packet sits further out",
          "[observables][radial_occupancy][engine]") {
	auto electrons = make_electrons();
	const int ist = electrons.states().num_states() - 1;

	// Displaced 4 Bohr off the axis: <r_perp> must be close to that displacement
	// (Rice distribution, -> mu for mu >> sigma_d), and most of the packet leaves
	// a 3-Bohr bore. This is the signature of a packet that has stopped channeling.
	fill_periodic_gaussian(electrons, ist, 4.0, 0.0, 0.0, SIGMA);

	auto occ = inqkit::observables::radial_occupancy(
		electrons, ist, {.axis = 2, .center = {0.0, 0.0, 0.0},
		                 .r_inner = 3.0, .r_outer = 5.0});

	// <r_perp^2> = mu^2 + 2 sigma_d^2 is EXACT for any offset mu — the tight check.
	// The Rice mean sits ~3.5 % above mu at mu/sigma_d = 3.77, so r_mean is only
	// bounded loosely; it is here to catch a sign/axis blunder, not to be precise.
	CHECK(occ.r2_mean == Approx(16.0 + 2.0*SIGMA_D*SIGMA_D).epsilon(5e-3));
	CHECK(occ.r_mean  == Approx(4.0).epsilon(0.08));
	CHECK(occ.f_bore  < 0.30);      // mostly out of the bore
	CHECK(occ.f_wall  > 0.50);      // and into the wall
	CHECK(occ.f_bore + occ.f_wall > 0.60);
}
