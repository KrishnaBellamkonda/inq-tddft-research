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
