// Pure host unit test for inqkit::dynamics::Projectile (velocity-Verlet).
// No INQ / GPU — the integrator is pure. Uses Catch2WithMain.
//
// Validates the integrator against closed-form motion:
//   1. zero force            -> constant velocity, linear drift
//   2. constant force (rest) -> V(t)=a t, R(t)=1/2 a t^2 (velocity-Verlet exact)
//   3. harmonic oscillator   -> total energy bounded (symplectic, no drift)
//   4. F/m acceleration      -> a = F/mass respected (mass scaling)

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <inqkit/dynamics/projectile.hpp>

#include <cmath>
#include <vector>

using inqkit::dynamics::Projectile;
using Vec3 = inqkit::detail::Vec3;
using Catch::Matchers::WithinAbs;

TEST_CASE("Projectile: zero force keeps velocity constant, drifts linearly", "[projectile]") {
	Projectile p(1.0, -1.0, Vec3{0, 0, 0}, Vec3{0, 0, 0.3});
	const double dt = 0.05;
	double z_expected = 0.0;
	for(int n = 0; n < 20; ++n) {
		CHECK_THAT(p.R().z, WithinAbs(z_expected, 1e-12));   // R_n before advance
		p.advance(Vec3{0, 0, 0}, dt);
		z_expected += 0.3 * dt;                              // R_{n+1} = R_n + V dt
	}
	CHECK_THAT(p.V().z, WithinAbs(0.3, 1e-12));
	CHECK_THAT(p.ke(), WithinAbs(0.5 * 0.09, 1e-12));
}

TEST_CASE("Projectile: constant force from rest is velocity-Verlet exact", "[projectile]") {
	const double m = 2.0, Fz = 0.4, a = Fz / m, dt = 0.02;
	Projectile p(m, -1.0, Vec3{0, 0, 0}, Vec3{0, 0, 0});
	for(int n = 0; n < 50; ++n) {
		const double t = n * dt;
		CHECK_THAT(p.R().z, WithinAbs(0.5 * a * t * t, 1e-9));   // R(n dt)=1/2 a t^2
		p.advance(Vec3{0, 0, Fz}, dt);
		CHECK_THAT(p.V().z, WithinAbs(a * t, 1e-9));             // V(n dt)=a t
	}
}

TEST_CASE("Projectile: acceleration scales as F/mass", "[projectile]") {
	Projectile light(1.0, -1.0, Vec3{0, 0, 0}, Vec3{0, 0, 0});
	Projectile heavy(4.0, -1.0, Vec3{0, 0, 0}, Vec3{0, 0, 0});
	const double dt = 0.01;
	light.advance(Vec3{0, 0, 1.0}, dt);
	heavy.advance(Vec3{0, 0, 1.0}, dt);
	// after one seed step, R = 1/2 (F/m) dt^2 -> heavy moves 4x less
	CHECK_THAT(light.R().z / heavy.R().z, WithinAbs(4.0, 1e-9));
}

TEST_CASE("Projectile: harmonic oscillator conserves energy (symplectic)", "[projectile]") {
	// F = -k x, m=1 -> omega=sqrt(k). Energy E = 1/2 v^2 + 1/2 k x^2 bounded.
	const double k = 1.0, m = 1.0, dt = 0.01;
	Projectile p(m, -1.0, Vec3{0, 0, 1.0}, Vec3{0, 0, 0});   // start at x=1, rest
	const double E0 = 0.5 * k * 1.0;                         // 0.5
	double emin = E0, emax = E0;
	for(int n = 0; n < 20000; ++n) {                        // ~30 periods
		const double Rz_n = p.R().z;                        // R_n (before advance)
		p.advance(Vec3{0, 0, -k * Rz_n}, dt);               // -> V_n (after)
		double E = p.ke() + 0.5 * k * Rz_n * Rz_n;          // both at time n dt
		emin = std::min(emin, E);
		emax = std::max(emax, E);
	}
	// velocity-Verlet: energy error is O(dt^2), bounded (no secular drift)
	CHECK((emax - emin) < 1e-3);
	CHECK_THAT(0.5 * (emax + emin), WithinAbs(E0, 2e-3));
}

// ---------------------------------------------------------------------------
// Periodic wrapping (docs/plans/slab-ks-orbital-stopping-wrap.md §4). set_R and
// wrap_into_cell exist so a classical projectile can re-enter the cell the way a
// KS orbital already does on the FFT grid.
// ---------------------------------------------------------------------------

TEST_CASE("Projectile: set_R moves the position and leaves the velocity alone", "[projectile][wrap]") {
	Projectile p(1.0, -1.0, Vec3{0, 0, 10.0}, Vec3{0, 0, 2.0});
	p.set_R(Vec3{1.0, -2.0, -30.0});
	CHECK_THAT(p.R().x, WithinAbs( 1.0, 1e-15));
	CHECK_THAT(p.R().y, WithinAbs(-2.0, 1e-15));
	CHECK_THAT(p.R().z, WithinAbs(-30.0, 1e-15));
	CHECK_THAT(p.V().z, WithinAbs(2.0, 1e-15));      // untouched
	CHECK_THAT(p.ke(),  WithinAbs(2.0, 1e-15));      // 1/2 * 1 * 2^2
}

TEST_CASE("Projectile: wrap_into_cell uses the [-L/2, +L/2) window", "[projectile][wrap]") {
	const Vec3 L{35.0, 35.0, 85.0};

	SECTION("interior point is untouched and reports no wrap") {
		Projectile p(1.0, -1.0, Vec3{0, 0, -24.0}, Vec3{0, 0, 2.0});
		CHECK(p.wrap_into_cell(L) == false);
		CHECK_THAT(p.R().z, WithinAbs(-24.0, 1e-15));
	}
	SECTION("past the +z face re-enters at -z, one lattice vector down") {
		Projectile p(1.0, -1.0, Vec3{0, 0, 43.0}, Vec3{0, 0, 2.0});
		CHECK(p.wrap_into_cell(L) == true);
		CHECK_THAT(p.R().z, WithinAbs(43.0 - 85.0, 1e-13));   // = -42.0
	}
	SECTION("the window is half-open: +L/2 wraps, -L/2 does not") {
		Projectile hi(1.0, -1.0, Vec3{0, 0,  42.5}, Vec3{0, 0, 0});
		Projectile lo(1.0, -1.0, Vec3{0, 0, -42.5}, Vec3{0, 0, 0});
		CHECK(hi.wrap_into_cell(L) == true);
		CHECK_THAT(hi.R().z, WithinAbs(-42.5, 1e-13));
		CHECK(lo.wrap_into_cell(L) == false);
		CHECK_THAT(lo.R().z, WithinAbs(-42.5, 1e-15));
	}
	SECTION("a zero length means that axis never wraps") {
		Projectile p(1.0, -1.0, Vec3{100.0, 0, 43.0}, Vec3{0, 0, 2.0});
		CHECK(p.wrap_into_cell(Vec3{0.0, 0.0, 85.0}) == true);
		CHECK_THAT(p.R().x, WithinAbs(100.0, 1e-15));         // x left alone
		CHECK_THAT(p.R().z, WithinAbs(-42.0, 1e-13));
	}
}

TEST_CASE("Projectile: wrapping every step leaves the unwrapped path intact", "[projectile][wrap]") {
	// The physical observable is the DISTANCE travelled, which must be unaffected
	// by relabelling the position across a face. Drift 400 Bohr at constant
	// velocity through a 85-Bohr cell (4+ wraps) and check the accumulated path
	// against v*t, while the wrapped coordinate always stays inside the cell.
	const double Lz = 85.0, v = 2.0, dt = 0.04;
	const int nsteps = 5000;
	// The loop reconstructs the path the way POST-PROCESSING will have to: from
	// the recorded WRAPPED coordinate alone, re-adding one lattice vector
	// whenever consecutive samples jump backwards by more than half a cell.
	Projectile p(1.0, -1.0, Vec3{0, 0, -24.0}, Vec3{0, 0, v});
	double path = 0.0, z_prev = p.R().z;
	int wraps = 0;
	for(int n = 0; n < nsteps; ++n) {
		p.advance(Vec3{0, 0, 0}, dt);
		if(p.wrap_into_cell(Vec3{0.0, 0.0, Lz})) ++wraps;
		const double z_now = p.R().z;
		double dz = z_now - z_prev;
		if(dz < -0.5 * Lz) dz += Lz;          // unwrap a backwards jump
		path += dz;
		z_prev = z_now;
		CHECK(p.R().z >= -0.5 * Lz);
		CHECK(p.R().z <   0.5 * Lz);
	}
	CHECK(wraps == 4);                                          // 400/85 = 4.7
	CHECK_THAT(path, WithinAbs(v * nsteps * dt, 1e-9));         // 400 Bohr
	CHECK_THAT(p.V().z, WithinAbs(v, 1e-12));                   // wrapping is not a force
}
