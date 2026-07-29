/*
 * inqkit::jellium::gaussian_projectile_perturbation — represent a CLASSICAL
 * projectile as a STATIONARY Gaussian charge distribution added to the Kohn–Sham
 * potential via its Poisson potential, instead of a UPF pseudopotential ghost.
 *
 * WHY (campaign localised-jellium-dynamics-analysis): the UPF ghost's local
 * potential erf(r/0.5)/r ≈ 1/r is long-range; with Z_valence=0 INQ places the
 * whole tail as a "short-range" local potential truncated at r_cut, so the
 * as-implemented projectile↔background term aliases / diverges with r_cut (see
 * reference_ghost_upf_tail_aliasing). Representing the projectile as a genuine
 * Gaussian CHARGE and taking its potential from the SAME Poisson solver as the
 * background removes the pseudopotential entirely: no r_cut, no aliasing, and the
 * projectile↔electron term (captured by INQ in E_external) and projectile↔background
 * term are computed in one consistent convention.
 *
 * Sign: the projectile is a −1 electron (like the wavepacket). An electron near a
 * −1 charge is REPELLED, so the KS potential gains v_proj = +poisson(n_proj), the
 * opposite sign to the +background well (background: v_bg = −poisson(n₊)). Mirrors
 * localised_background_perturbation exactly (pointwise `vk += φ`, complex-safe).
 *
 * Composable with the background via inq::perturbations::sum(background, projectile).
 * The projectile is an EXTERNAL potential, not an INQ ion and not an electron — so
 * the cell charge count is unchanged (still neutral for a slab GS), and INQ adds
 * ∫n_e·v_proj to E_external automatically. Its self-energy and its interaction with
 * the background are NOT in INQ's total (same as any classical projectile) — the
 * background term is the diagnostic U_proj_bg = −∫n_proj·φ₊.
 *
 * No inq/ or inq-study/ edit — wrapper-only, on top of the existing perturbation.
 */
#ifndef INQKIT__JELLIUM__GAUSSIAN_PROJECTILE_PERTURBATION
#define INQKIT__JELLIUM__GAUSSIAN_PROJECTILE_PERTURBATION

#include <inq/inq.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density

#include <optional>

namespace inqkit {
namespace jellium {

class gaussian_projectile_perturbation : public inq::perturbations::none {

public:

	// center: projectile position (Bohr); sigma_pot: charge std = σ_WP/√2.
	gaussian_projectile_perturbation(inq::vector3<double> center, double sigma_pot)
		: center_(center), sigma_(sigma_pot) {}

	auto has_potential() const { return true; }

	// Add v_proj = +poisson(n_proj) to the (real or complex) KS potential field
	// (−1 projectile ⇒ repulsive well). φ_proj computed once and cached.
	template <typename PotentialType>
	void potential(const double /*time*/, PotentialType & potential) const {

		if(not phi_.has_value()) {
			auto nproj = gaussian_density(potential.basis(), center_, sigma_);
			phi_.emplace(inq::solvers::poisson::solve(nproj));   // φ_proj = poisson(n_proj)
		}

		auto phi_cub = begin(phi_->cubic());
		auto vk_cub  = begin(potential.cubic());

		gpu::run(potential.basis().local_sizes()[2],
		         potential.basis().local_sizes()[1],
		         potential.basis().local_sizes()[0],
			[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
				vk_cub[ix][iy][iz] += phi_cub[ix][iy][iz];   // −1 projectile ⇒ +φ (repulsive)
			});
	}

	// Projectile charge density n_proj (∫ = 1) on a given basis (for diagnostics).
	template <class Basis>
	auto charge_density(Basis const & basis) const {
		return gaussian_density(basis, center_, sigma_);
	}

	inq::vector3<double> center() const { return center_; }
	double sigma() const { return sigma_; }

private:
	inq::vector3<double> center_;
	double sigma_;
	mutable std::optional<inq::basis::field<inq::basis::real_space, double>> phi_;
};

} // namespace jellium
} // namespace inqkit

#endif
