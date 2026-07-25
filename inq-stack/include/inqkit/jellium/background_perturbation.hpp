/*
 * inqkit::jellium::localised_background_perturbation — inject a confined jellium
 * background into the Kohn–Sham potential as a STATIC external potential.
 *
 * Mechanism (see docs/notes/localised-jellium-theory.md, ADR-0008)
 * ----------------------------------------------------------------
 * INQ exposes a Perturbation hook used by BOTH ground_state::calculate and
 * real_time::propagate; self_consistency::update_hamiltonian calls
 * pert_.potential(t, vscalar) while assembling the KS potential. We use it to add
 * the electrostatic well of a localised positive background:
 *
 *     n₊(r)            built by make_localised_background()  (positive charge)
 *     φ(r) = poisson(n₊)                                     (its potential, +hump)
 *     v_bg(r) = −φ(r)  added to the KS potential             (electron attraction)
 *
 * Total electrostatics the electrons feel = poisson(n_elec) − φ(n₊) =
 * poisson(n_elec − n₊), which is charge-neutral and exact when ∫n₊ = N (the
 * caller's responsibility; the dropped G=0 cancels — theory Part 2.4).
 *
 * Because the same object is handed to the SCF and the propagator, the well
 * confines the electrons in the ground state AND persists while a projectile
 * flies. The potential is time-INDEPENDENT (static well); φ is computed once and
 * cached.
 *
 * Real-into-complex: on inq-study the KS scalar potential is complexified (so the
 * sin² CAP's imaginary part propagates). We add a REAL well via an explicit
 * pointwise loop (`vk -= φ`), which is `complex += double` on inq-study and
 * `double += double` on stock inq — no operations::increment type clash.
 *
 * No `inq/` or `inq-study/` edit: this is a wrapper-only perturbation, composable
 * with the CAP via perturbations::sum(background, CAP_−z, CAP_+z).
 *
 * Status: UNVERIFIED — pending compile + T1 engine test (SCF binds electrons
 * inside the region; e–background energy < 0).
 */

#ifndef INQKIT__JELLIUM__BACKGROUND_PERTURBATION
#define INQKIT__JELLIUM__BACKGROUND_PERTURBATION

#include <inq/inq.hpp>
#include <inqkit/jellium/localised_background.hpp>

#include <optional>

namespace inqkit {
namespace jellium {

class localised_background_perturbation : public inq::perturbations::none {

public:

	explicit localised_background_perturbation(localised_background_params params)
		: params_(params) {}

	auto has_potential() const { return true; }

	// EXPLAIN: What is phi? Where does it come from?
	// What does phi_.emplace do?
	// Where does is say that we are adding to the potential? 

	// Add v_bg = −poisson(n₊) to the (real or complex) KS potential field.
	template <typename PotentialType>
	void potential(const double /*time*/, PotentialType & potential) const {

		if(not phi_.has_value()) {
			auto nplus = make_localised_background(potential.basis(), params_);
			phi_.emplace(inq::solvers::poisson::solve(nplus)); // φ = poisson(n₊), +hump
		}

		auto phi_cub = begin(phi_->cubic());
		auto vk_cub  = begin(potential.cubic());

		gpu::run(potential.basis().local_sizes()[2],
		         potential.basis().local_sizes()[1],
		         potential.basis().local_sizes()[0],
			[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
				vk_cub[ix][iy][iz] -= phi_cub[ix][iy][iz]; // electron well = −φ
			});
	}

	// Background charge density n₊ on a given basis (for diagnostics / energy
	// bookkeeping: ∫v_bg·n is captured by INQ automatically, but E_self of the
	// background must be added separately — see analytics.hpp).
	template <class Basis>
	auto background_density(Basis const & basis) const {
		return make_localised_background(basis, params_);
	}

	localised_background_params const & params() const { return params_; }

private:
	localised_background_params params_;
	mutable std::optional<inq::basis::field<inq::basis::real_space, double>> phi_;
};

} // namespace jellium
} // namespace inqkit

#endif
