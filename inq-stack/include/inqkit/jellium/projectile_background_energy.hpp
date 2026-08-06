/*
 * inqkit::jellium::projectile_background_energy — the Coulomb interaction energy
 * between a classical projectile and the localised positive jellium background.
 *
 * WHY (campaign localised-jellium-dynamics-analysis, Phase 1): INQ's fixed 8-term
 * total energy has NO slot for the ghost↔background interaction (the ghost is an
 * INQ ion; the background is a `localised_background_perturbation`, not an ion —
 * so no Ewald `ion` term couples them). This is the missing `U_proj_bg` ledger
 * column. It is DIAGNOSTIC ONLY — never added to energy_total or the SCF/dynamics.
 *
 * TWO reciprocal estimates are returned (user decision — track both):
 *   ideal = ∫ n_proj · v_bg
 *       n_proj = the projectile's TRUE Gaussian charge (std σ_pot = σ_WP/√2, ∫=1),
 *       v_bg   = the background potential the electrons feel = −poisson(n₊)
 *                (the SAME field `localised_background_perturbation` adds to the KS
 *                 potential; background_perturbation.hpp:65-77). Uses the projectile's
 *                real charge ⇒ INDEPENDENT of the pseudopotential radial cutoff r_cut
 *                (Phase 2: the r_cut effect then lives entirely in E_external).
 *   impl  = −∫ n₊ · v_ion
 *       v_ion = the projectile's local pseudopotential AS IMPLEMENTED
 *               = poisson(atomic_pot.ionic_density) + atomic_pot.local_potential
 *               (mirrors self_consistency.hpp:102). Truncating r_cut changes v_ion,
 *               so `impl` DOES vary with r_cut. Equals `ideal` by reciprocity for the
 *               full (untruncated) Gaussian ghost — the code-gate check.
 *
 * Both integrals are gauge-tied to the p2 open-z G=0 Poisson convention (the same
 * charged-cell caveat as E_hartree/E_external individually — reference_charged_cell
 * _hartree_convention); the ledger uses them only in the WP−CL DIFFERENCE where the
 * convention constant cancels.
 *
 * No inq/ or inq-study/ edit — wrapper-only, on top of the existing perturbation.
 */
#ifndef INQKIT__JELLIUM__PROJECTILE_BACKGROUND_ENERGY
#define INQKIT__JELLIUM__PROJECTILE_BACKGROUND_ENERGY

#include <inq/inq.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include <cmath>

namespace inqkit {
namespace jellium {

struct projectile_background_result {
	double ideal;   // ∫ n_proj · v_bg   (true Gaussian charge; r_cut-invariant)
	double impl;    // −∫ n₊ · v_ion      (as-implemented pseudopotential; r_cut-dependent)
	double n_proj_norm;  // ∫ n_proj  (sanity: should be ≈ 1)
};

// Build a normalised 3D Gaussian charge density (∫ = 1) of std `sigma` centred at
// `center`, on the given real-space basis. σ_pot ≪ box and the projectile sits on
// the z-axis interior, so lateral periodic images are negligible (single Gaussian).
template <class Basis>
inq::basis::field<Basis, double>
gaussian_density(Basis const & basis, inq::vector3<double> center, double sigma) {
	inq::basis::field<Basis, double> n(basis);
	auto point_op = basis.point_op();
	auto cub = begin(n.cubic());
	const double norm    = 1.0 / std::pow(2.0 * M_PI * sigma * sigma, 1.5);
	const double inv2s2  = 1.0 / (2.0 * sigma * sigma);
	const double cx = center[0], cy = center[1], cz = center[2];
	gpu::run(basis.local_sizes()[2], basis.local_sizes()[1], basis.local_sizes()[0],
		[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
			auto r = point_op.rvector_cartesian(ix, iy, iz);
			double dx = r[0] - cx, dy = r[1] - cy, dz = r[2] - cz;
			cub[ix][iy][iz] = norm * exp(-(dx*dx + dy*dy + dz*dz) * inv2s2);
		});
	return n;
}

// Same normalised Gaussian, but built from the MINIMUM-IMAGE displacement: the
// separation r − center is taken to fractional coordinates, wrapped into
// [−0.5, +0.5) per lattice direction (the same window as
// systems::cell::position_in_cell), and converted back. The blob therefore wraps
// smoothly around every cell face instead of being clipped by it.
//
// WHY THIS EXISTS SEPARATELY (docs/plans/slab-ks-orbital-stopping-wrap.md §4).
// `gaussian_density` above uses a plain Cartesian displacement, so a projectile
// sitting on a box face loses the half of its charge that falls outside — the
// integral drops below 1 and the force is wrong for as long as it straddles. A
// KS ORBITAL has no such problem: the wavefunction basis is a plain 3-D FFT,
// periodic in all three directions, so a wavepacket wraps exactly. In a
// wrap-around classical-vs-wavepacket twin the two projectiles would then differ
// precisely at the boundary the study deliberately introduces. This function
// gives the classical charge the same wrap the wavepacket already has.
//
// It is a SEPARATE function, not a flag on the original, so every previously
// published run keeps its exact binary behaviour.
//
// The wrap is lattice-general (it goes through the cell's contravariant basis),
// so it is correct for non-orthogonal cells too. For a blob whose σ is small
// against every lattice length and which sits well inside the cell, it agrees
// with `gaussian_density` to machine precision — that equivalence is the test.
//
// NOTE ON ELECTROSTATICS. Wrapping the CHARGE does not change the Poisson
// boundary condition. Under periodicity(2) a blob straddling the z face is
// solved as two lumps at opposite ends of a z-open box — which is exactly what
// the Poisson solver also does with the straddling wavepacket density, so the
// twins stay matched.
template <class Basis>
inq::basis::field<Basis, double>
gaussian_density_minimum_image(Basis const & basis, inq::vector3<double> center, double sigma) {
	inq::basis::field<Basis, double> n(basis);
	auto point_op = basis.point_op();
	auto cub = begin(n.cubic());
	const double norm    = 1.0 / std::pow(2.0 * M_PI * sigma * sigma, 1.5);
	const double inv2s2  = 1.0 / (2.0 * sigma * sigma);
	const double cx = center[0], cy = center[1], cz = center[2];
	gpu::run(basis.local_sizes()[2], basis.local_sizes()[1], basis.local_sizes()[0],
		[=] GPU_LAMBDA (auto iz, auto iy, auto ix) {
			auto r = point_op.rvector_cartesian(ix, iy, iz);
			inq::vector3<double> d{r[0] - cx, r[1] - cy, r[2] - cz};
			auto f = point_op.cell().to_contravariant(d);
			for(int idir = 0; idir < 3; idir++) {
				f[idir] -= floor(f[idir]);            // -> [0, 1)
				if(f[idir] >= 0.5) f[idir] -= 1.0;    // -> [-0.5, 0.5)
			}
			auto dc = point_op.cell().to_cartesian(f);
			cub[ix][iy][iz] = norm * exp(-(dc[0]*dc[0] + dc[1]*dc[1] + dc[2]*dc[2]) * inv2s2);
		});
	return n;
}

// E_proj_bg (both estimates) for a projectile at `proj_center` with charge-std
// `sigma_pot`, against the background of `pert`. `electrons` supplies the density
// basis + atomic_pot; `ions` holds the projectile (for v_ion).
template <class Perturbation, class Electrons, class Ions>
projectile_background_result
projectile_background_energy(Perturbation const & pert, Electrons & electrons,
                             Ions const & ions, inq::vector3<double> proj_center,
                             double sigma_pot) {
	auto basis = electrons.density().basis();
	auto & comm = electrons.states_comm();

	// background n₊ and its potential; v_bg (electron well) = −φ₊
	auto nplus   = pert.background_density(basis);
	auto phiplus = inq::solvers::poisson::solve(nplus);          // φ₊ = poisson(n₊)

	// projectile true Gaussian charge (∫ = 1)
	auto nproj = gaussian_density(basis, proj_center, sigma_pot);

	// projectile local potential as implemented (mirrors self_consistency.hpp:102)
	auto & atomic_pot = electrons.atomic_pot();
	auto v_ion = inq::operations::add(
		inq::solvers::poisson::solve(atomic_pot.ionic_density(comm, basis, ions)),
		atomic_pot.local_potential(comm, basis, ions));

	projectile_background_result out;
	out.n_proj_norm = inq::operations::integral(nproj);
	// ideal = ∫ n_proj · v_bg = −∫ n_proj · φ₊
	out.ideal = -inq::operations::integral_product(nproj, phiplus);
	// impl  = −∫ n₊ · v_ion
	out.impl  = -inq::operations::integral_product(nplus, v_ion);
	return out;
}

} // namespace jellium
} // namespace inqkit

#endif
