// Pure host unit test for the analytic Gaussian-projectile SELF-ENERGY constant used
// by inqkit::jellium::compute_coulomb_direct (e_pp = 1/(2 sigma sqrt(pi))).
//
// No INQ / GPU. Validates the hardcoded closed form against an INDEPENDENT numerical
// self-Coulomb integral of a normalized 3D Gaussian charge:
//   n(r) = (2 pi s^2)^{-3/2} exp(-r^2 / 2 s^2),   phi(r) = erf(r/(sqrt2 s))/r
//   U_self = 1/2 * Integral n(r) phi(r) d3r
//          = 1/2 * Integral_0^inf 4 pi r^2 n(r) phi(r) dr
// Closed form (Ewald Gaussian self-energy): U_self = 1/(2 s sqrt(pi)).
// This is the value compute_coulomb_direct assigns to e_pp (constant, no clip/kink).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>

using Catch::Matchers::WithinRel;

// analytic constant (mirrors interaction_energies.hpp: 1/(2 sigma sqrt(pi)))
static double analytic_self_energy(double sigma) {
	return 1.0 / (2.0 * sigma * std::sqrt(M_PI));
}

// independent numerical self-Coulomb of the normalized Gaussian (radial quadrature)
static double numeric_self_energy(double sigma, int n = 400000, double rmax_sig = 12.0) {
	const double s2   = sigma * sigma;
	const double norm = std::pow(2.0 * M_PI * s2, -1.5);
	const double inv  = 1.0 / (std::sqrt(2.0) * sigma);
	const double rmax = rmax_sig * sigma;
	const double h    = rmax / n;
	double acc = 0.0;                                   // trapezoid on [h/2 .. rmax]
	for(int i = 0; i < n; ++i) {
		const double r  = (i + 0.5) * h;               // midpoint (avoids r=0 singularity)
		const double nr = norm * std::exp(-r * r / (2.0 * s2));
		const double phi = std::erf(r * inv) / r;      // potential of the Gaussian
		acc += 4.0 * M_PI * r * r * nr * phi * h;      // 4 pi r^2 n phi dr
	}
	return 0.5 * acc;
}

TEST_CASE("Gaussian projectile self-energy matches 1/(2 sigma sqrt(pi))", "[interaction_energies][direct]") {
	for(double sigma : {0.353553390593, 0.5, 1.0}) {   // sigma_pot(WP 0.5), and two more
		const double an = analytic_self_energy(sigma);
		const double nu = numeric_self_energy(sigma);
		CHECK_THAT(nu, WithinRel(an, 1e-4));
	}
}

TEST_CASE("Self-energy value for the sigma_pot=0.35355 projectile is ~0.798 Ha", "[interaction_energies][direct]") {
	// sigma_pot = sigma_WP/sqrt2 = 0.5/sqrt2 = 0.353553...; guards the reported ~0.80 Ha
	CHECK_THAT(analytic_self_energy(0.5 / std::sqrt(2.0)), WithinRel(0.7979, 1e-3));
}
