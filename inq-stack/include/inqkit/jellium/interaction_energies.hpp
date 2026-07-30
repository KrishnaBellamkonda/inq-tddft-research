/*
 * inqkit::jellium::interaction_energies — decompose the electrostatic energy into
 * the physical pairwise terms of the three charge groups P (projectile), S (slab
 * electrons), B (positive background), using INQ's OWN boundary-matched Poisson so
 * the terms sum EXACTLY back to the broadcast INQ scalars (the closure check).
 *
 * Given n_slab, n_P (both INQ fields) and φ₊ = poisson(n₊):
 *   E_SS = ½∫n_slab·φ_slab            slab-slab           (= INQ E_hartree, classical)
 *   E_PP = ½∫n_P·φ_P                  projectile self-Hartree (the quantum residual)
 *   E_PS =  ∫n_slab·φ_P               projectile-slab
 *   E_SB = -∫n_slab·φ₊  (= ∫n_slab·v_bg)   slab-background
 *   E_PB = -∫n_P·φ₊     (= ∫n_P·v_bg)      projectile-background (= U_proj_bg)
 * plus norms (∫n_slab, ∫n_P) as sanity.  E_BB = ½∫n₊·φ₊ is a constant (background
 * self, G=0-gauge) computed once via background_self_energy().
 *
 * Closure (see docs/plans/twin-run-rung2-dynamic-spec.md):
 *   classical: E_hartree = E_SS ; E_external = E_SB + E_PS
 *   WP:        E_hartree = E_SS + E_PS + E_PP ; E_external = E_SB + E_PB
 *
 * Absolute E_SB/E_PB/E_BB carry the charged-cell G=0 gauge; only closure sums and
 * WP-classical differences are gauge-clean (reference_charged_cell_hartree_convention).
 *
 * No inq/ or inq-study/ edit — wrapper-only.
 */
#ifndef INQKIT__JELLIUM__INTERACTION_ENERGIES
#define INQKIT__JELLIUM__INTERACTION_ENERGIES

#include <inq/inq.hpp>
#include <cmath>

namespace inqkit {
namespace jellium {

struct coulomb_terms {
	double e_ss = 0.0;   // slab-slab (½∫n_slab·φ_slab)
	double e_pp = 0.0;   // projectile self-Hartree (½∫n_P·φ_P)
	double e_ps = 0.0;   // projectile-slab (∫n_slab·φ_P)
	double e_sb = 0.0;   // slab-background (-∫n_slab·φ₊)
	double e_pb = 0.0;   // projectile-background (-∫n_P·φ₊)
	double norm_slab = 0.0;
	double norm_p = 0.0;
};

// n_slab, n_P: INQ real-space fields; phiplus = poisson(n₊). All on one basis.
template <class Field>
coulomb_terms compute_coulomb(Field const & n_slab, Field const & n_P, Field const & phiplus) {
	auto phi_slab = inq::solvers::poisson::solve(n_slab);
	auto phi_P    = inq::solvers::poisson::solve(n_P);
	coulomb_terms t;
	t.e_ss = 0.5 * inq::operations::integral_product(n_slab, phi_slab);
	t.e_pp = 0.5 * inq::operations::integral_product(n_P, phi_P);
	t.e_ps =       inq::operations::integral_product(n_slab, phi_P);
	t.e_sb = -     inq::operations::integral_product(n_slab, phiplus);
	t.e_pb = -     inq::operations::integral_product(n_P, phiplus);
	t.norm_slab =  inq::operations::integral(n_slab);
	t.norm_p    =  inq::operations::integral(n_P);
	return t;
}

// DIRECT-potential pairwise decomposition (companion to the direct erf/r projectile
// perturbation, moving_gaussian_projectile_potential). The projectile terms are formed
// from the DIRECT free-space potential v_proj = gaussian_potential(basis,center,σ_pot)
// = erf(|r-R|/(√2σ))/|r-R| — NO charge on the grid, NO Poisson, NO periodic
// neutralizing background. This removes BOTH charge-based artifacts of compute_coulomb:
//   (i)  the exit KINK (the charge no longer clips at the z-wall), and
//   (ii) the linear e_ps DRIFT (the localised potential is not replicated into an
//        x,y-periodic charged sheet), so e_ps is POSITIVE and →0 as the projectile
//        recedes — the physical trend (verified against ∫n·erf/r on saved frames).
//
//   E_PP = 1/(2·σ_pot·√π)            projectile self-Hartree (analytic, CONSTANT for a
//                                    rigid Gaussian — no clip, no kink)
//   E_PS = +∫ n_slab·v_proj          projectile-slab   (repulsion: −proj, −electrons > 0)
//   E_SB = −∫ n_slab·φ₊              slab-background   (Poisson, unchanged)
//   E_PB = −∫ n₊·v_proj              projectile-background (reciprocity; attraction < 0)
//   E_SS = ½∫ n_slab·φ_slab          slab-slab (= INQ E_hartree, classical)
// Closure (classical, direct rep): E_hartree = E_SS ; E_external = E_SB + E_PS
// (the applied external field is v_bg + v_proj, so ∫n·v_ext = e_sb + e_ps EXACTLY).
// norm_p is set to 1.0 (the direct potential carries no grid charge to integrate).
template <class Field>
coulomb_terms compute_coulomb_direct(Field const & n_slab, Field const & v_proj,
                                     Field const & nplus, Field const & phiplus,
                                     double sigma_pot) {
	auto phi_slab = inq::solvers::poisson::solve(n_slab);
	coulomb_terms t;
	t.e_ss = 0.5 * inq::operations::integral_product(n_slab, phi_slab);
	t.e_pp = 1.0 / (2.0 * sigma_pot * std::sqrt(M_PI));            // analytic self-energy (constant)
	t.e_ps =        inq::operations::integral_product(n_slab, v_proj);   // + , →0
	t.e_sb = -      inq::operations::integral_product(n_slab, phiplus);
	t.e_pb = -      inq::operations::integral_product(nplus,  v_proj);   // attraction, →0
	t.norm_slab =   inq::operations::integral(n_slab);
	t.norm_p    =   1.0;   // direct potential: no charge density on the grid
	return t;
}

// Background self-Coulomb E_BB = ½∫n₊·φ₊ (constant; a gauge estimate). Compute once.
template <class Field>
double background_self_energy(Field const & nplus, Field const & phiplus) {
	return 0.5 * inq::operations::integral_product(nplus, phiplus);
}

// Build the WP orbital's contribution to the density (occ·|ψ_idx|²) as an INQ
// field, so it can be Poisson-solved (single-rank, gamma point). Mirrors INQ's
// observables::density::calculate_add but for ONE orbital.
inline inq::basis::field<inq::basis::real_space, double>
orbital_density_field(inq::systems::electrons const & electrons, int idx, int ik = 0) {
	auto const & phi = electrons.kpin()[ik];
	inq::basis::field<inq::basis::real_space, double> n(phi.basis());
	n.fill(0.0);
	gpu::run(phi.basis().part().local_size(),
		[idx, occ = electrons.occupations()[ik].cbegin(),
		 ph = phi.matrix().cbegin(), nn = n.linear().begin()] GPU_LAMBDA (auto ip) {
			nn[ip] = occ[idx] * norm(ph[ip][idx]);
		});
	return n;
}

// WP-run decomposition WITHOUT forming n_slab, via Poisson linearity
// (φ_total = φ_slab + φ_wp) and reciprocity, so the terms close EXACTLY:
//   E_PP = ½∫n_wp·φ_wp
//   E_PS = ∫n_wp·φ_total − 2·E_PP
//   E_SS = ½∫n_total·φ_total − ∫n_wp·φ_total + E_PP     (== E_hartree − E_PS − E_PP)
//   E_PB = −∫n_wp·φ₊ ;  E_SB = −∫n_total·φ₊ − E_PB
// e_hartree_check = ½∫n_total·φ_total must equal INQ E_hartree; e_external_check
// = −∫n_total·φ₊ must equal INQ E_external (closure gates).
struct coulomb_terms_wp {
	double e_ss = 0.0, e_pp = 0.0, e_ps = 0.0, e_sb = 0.0, e_pb = 0.0;
	double e_hartree_check = 0.0, e_external_check = 0.0;
	double norm_wp = 0.0, norm_total = 0.0;
};

template <class Field>
coulomb_terms_wp compute_coulomb_wp(Field const & n_total, Field const & n_wp, Field const & phiplus) {
	auto phi_total = inq::solvers::poisson::solve(n_total);
	auto phi_wp    = inq::solvers::poisson::solve(n_wp);
	coulomb_terms_wp t;
	t.e_pp   = 0.5 * inq::operations::integral_product(n_wp, phi_wp);
	const double cross = inq::operations::integral_product(n_wp, phi_total);   // ∫n_wp·φ_total
	t.e_ps   = cross - 2.0 * t.e_pp;
	t.e_hartree_check = 0.5 * inq::operations::integral_product(n_total, phi_total);
	t.e_ss   = t.e_hartree_check - cross + t.e_pp;
	t.e_pb   = -inq::operations::integral_product(n_wp, phiplus);
	t.e_external_check = -inq::operations::integral_product(n_total, phiplus);
	t.e_sb   = t.e_external_check - t.e_pb;
	t.norm_wp    = inq::operations::integral(n_wp);
	t.norm_total = inq::operations::integral(n_total);
	return t;
}

} // namespace jellium
} // namespace inqkit

#endif
