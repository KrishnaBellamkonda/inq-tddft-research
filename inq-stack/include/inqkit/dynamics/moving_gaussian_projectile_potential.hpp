/*
 * inqkit::dynamics::moving_gaussian_projectile_potential — a projectile perturbation that
 * adds the DIRECT free-space potential of the Gaussian projectile (erf/r), tracking a live
 * Projectile. Drop-in alternative to moving_gaussian_projectile_perturbation.
 *
 * Difference: this adds inqkit::jellium::gaussian_potential(basis, R, sigma) DIRECTLY —
 * NO charge density built, NO poisson::solve, NO periodic neutralizing background. So there
 * is no charge in the cell to clip, and the exit transient (the ±hundreds-of-eV G=0-offset
 * lurch as the clipped charge crosses the box face) does not occur. The gradient (force) is
 * the same as the Poisson version; only the (charge-dependent) background offset in the
 * energy bookkeeping is removed.
 *
 * −1 projectile ⇒ +V added (repulsive to bath electrons), matching the sign convention of
 * moving_gaussian_projectile_perturbation.
 *
 * No inq/ or inq-study/ edit — wrapper-only.
 */
#ifndef INQKIT__DYNAMICS__MOVING_GAUSSIAN_PROJECTILE_POTENTIAL
#define INQKIT__DYNAMICS__MOVING_GAUSSIAN_PROJECTILE_POTENTIAL

#include <inq/inq.hpp>
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/jellium/gaussian_potential.hpp>

namespace inqkit {
namespace dynamics {

class moving_gaussian_projectile_potential : public inq::perturbations::none {

public:

	moving_gaussian_projectile_potential(Projectile const & proj, double sigma_pot)
		: proj_(&proj), sigma_(sigma_pot) {}

	auto has_potential() const { return true; }

	template <typename PotentialType>
	void potential(const double /*time*/, PotentialType & potential) const {
		auto R = proj_->R();
		auto V = inqkit::jellium::gaussian_potential(potential.basis(),
			inq::vector3<double>{R.x, R.y, R.z}, sigma_);
		auto v_cub = begin(potential.cubic());
		auto V_cub = begin(V.cubic());
		gpu::run(potential.basis().local_sizes()[2],
		         potential.basis().local_sizes()[1],
		         potential.basis().local_sizes()[0],
			[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
				v_cub[ix][iy][iz] += V_cub[ix][iy][iz];   // −1 projectile ⇒ +V
			});
	}

	double sigma() const { return sigma_; }

private:
	Projectile const * proj_;
	double sigma_;
};

} // namespace dynamics
} // namespace inqkit

#endif
