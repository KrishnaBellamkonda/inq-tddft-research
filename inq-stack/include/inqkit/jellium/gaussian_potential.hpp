/*
 * inqkit::jellium::gaussian_potential — the DIRECT free-space Hartree potential of a
 * unit Gaussian charge (density std `sigma`) centred at R, evaluated on the grid:
 *
 *     V(r) = erf(|r-R| / (sqrt2 * sigma)) / |r-R|        (limit sqrt(2/pi)/sigma at r->R)
 *
 * This is the potential a SINGLE isolated projectile presents — added as an external
 * field WITHOUT building a charge density, WITHOUT a Poisson solve, and WITHOUT the
 * periodic G=0 neutralizing background. Contrast poisson(gaussian_density(...)), which
 * solves the periodic Poisson of the charge and so carries the neutralizing background +
 * transverse periodic images (whose offset lurches as the clipped charge crosses the box
 * boundary — the exit transient). The gradient (hence the force) of erf/r and of
 * poisson(gaussian) agree; they differ only by the (charge-dependent) background offset.
 *
 * Uses the SAME direct Cartesian coordinate as gaussian_density (point_op.rvector_cartesian,
 * no minimum image) — for a projectile centred in the transverse plane this is
 * periodic-continuous within the cell.
 *
 * No inq/ or inq-study/ edit — wrapper-only.
 */
#ifndef INQKIT__JELLIUM__GAUSSIAN_POTENTIAL
#define INQKIT__JELLIUM__GAUSSIAN_POTENTIAL

#include <inq/inq.hpp>
#include <cmath>

namespace inqkit {
namespace jellium {

template <class Basis>
inq::basis::field<Basis, double> gaussian_potential(Basis const & basis,
                                                    inq::vector3<double> center,
                                                    double sigma) {
	inq::basis::field<Basis, double> V(basis);
	auto point_op = basis.point_op();
	auto cub = begin(V.cubic());
	const double inv = 1.0 / (std::sqrt(2.0) * sigma);      // 1/(sqrt2 sigma)
	const double vzero = std::sqrt(2.0 / M_PI) / sigma;     // V(r->0) = sqrt(2/pi)/sigma
	const double cx = center[0], cy = center[1], cz = center[2];
	gpu::run(basis.local_sizes()[2], basis.local_sizes()[1], basis.local_sizes()[0],
		[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
			auto r = point_op.rvector_cartesian(ix, iy, iz);
			double dx = r[0] - cx, dy = r[1] - cy, dz = r[2] - cz;
			double d = sqrt(dx*dx + dy*dy + dz*dz);
			cub[ix][iy][iz] = (d < 1e-8) ? vzero : erf(d * inv) / d;
		});
	return V;
}

} // namespace jellium
} // namespace inqkit

#endif
