/*
 * inqkit::dynamics::moving_gaussian_projectile_perturbation — a Gaussian-charge
 * projectile whose center TRACKS a live inqkit::dynamics::Projectile.
 *
 * Identical physics to jellium::gaussian_projectile_perturbation (v_proj =
 * +poisson(n_proj), −1 projectile ⇒ repulsive well), but the center is read from
 * a Projectile* each time INQ asks for the potential, and φ_proj is re-solved when
 * the center moves. The Projectile is advanced in the RT step callback (Ehrenfest);
 * this perturbation is the read side of that shared state — the mask-absorber
 * pattern (reference_inq_propagator_mask_absorber). Holding a POINTER (not a copy
 * of the state) means the coupling survives INQ copying the perturbation by value.
 *
 * Within one RT step the center is constant, so φ is cached and only re-solved
 * across steps (one extra Poisson solve per step).
 *
 * No inq/ or inq-study/ edit — wrapper-only.
 */
#ifndef INQKIT__DYNAMICS__MOVING_GAUSSIAN_PROJECTILE_PERTURBATION
#define INQKIT__DYNAMICS__MOVING_GAUSSIAN_PROJECTILE_PERTURBATION

#include <inq/inq.hpp>
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density

#include <optional>

namespace inqkit {
namespace dynamics {

class moving_gaussian_projectile_perturbation : public inq::perturbations::none {

public:

	// proj: live projectile state (must outlive this perturbation).
	// sigma_pot: charge std = σ_WP/√2.
	// minimum_image: build the blob from the MINIMUM-IMAGE displacement, so it
	//   wraps smoothly around the cell faces instead of being clipped by them.
	//   Defaults to false — every previously published run keeps its exact
	//   behaviour. Set it for wrap-around runs, where the classical projectile
	//   re-enters the cell and must do so the way a KS orbital already does on
	//   the FFT grid (docs/plans/slab-ks-orbital-stopping-wrap.md §4).
	moving_gaussian_projectile_perturbation(Projectile const & proj, double sigma_pot,
	                                        bool minimum_image = false)
		: proj_(&proj), sigma_(sigma_pot), minimum_image_(minimum_image) {}

	auto has_potential() const { return true; }

	template <typename PotentialType>
	void potential(const double /*time*/, PotentialType & potential) const {

		auto R = proj_->R();
		inq::vector3<double> center{R.x, R.y, R.z};

		if(not phi_.has_value() or center != cached_center_) {
			auto nproj = minimum_image_
				? inqkit::jellium::gaussian_density_minimum_image(potential.basis(), center, sigma_)
				: inqkit::jellium::gaussian_density(potential.basis(), center, sigma_);
			phi_.emplace(inq::solvers::poisson::solve(nproj));   // φ_proj = poisson(n_proj)
			cached_center_ = center;
		}

		auto phi_cub = begin(phi_->cubic());
		auto vk_cub  = begin(potential.cubic());

		gpu::run(potential.basis().local_sizes()[2],
		         potential.basis().local_sizes()[1],
		         potential.basis().local_sizes()[0],
			[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
				vk_cub[ix][iy][iz] += phi_cub[ix][iy][iz];   // −1 projectile ⇒ +φ
			});
	}

	double sigma() const { return sigma_; }
	bool minimum_image() const { return minimum_image_; }

private:
	Projectile const * proj_;
	double sigma_;
	bool   minimum_image_ = false;
	mutable std::optional<inq::basis::field<inq::basis::real_space, double>> phi_;
	mutable inq::vector3<double> cached_center_{0.0, 0.0, 0.0};
};

} // namespace dynamics
} // namespace inqkit

#endif
