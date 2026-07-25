/*
 * inqkit::dynamics::projectile_force_z — the Hellmann-Feynman force on a rigid
 * Gaussian classical projectile, along the beam axis z.
 *
 * The only R-dependent part of the total energy is the projectile's interaction
 * with everything ELSE (electrons + background):
 *     E_R(R) = ∫ n_e·v_proj(·-R)  +  U_proj_bg(R)
 *            = ∫ n_proj(·-R)·φ_e   −  ∫ n_proj(·-R)·φ_+        [reciprocity]
 *            = ∫ n_proj(·-R)·φ_drag,     φ_drag = φ_e − φ_+ = poisson(n_e − n_+)
 * so F_z = −dE_R/dR_z. The projectile self-Hartree is R-independent for a rigid
 * symmetric Gaussian, so it contributes NO self-force (excluded automatically).
 *
 * Computed by a symmetric finite difference of the (cheap) reciprocal integral
 * ∫ n_proj(·-R±δẑ)·φ_drag — one Poisson solve for φ_drag per step (done by the
 * caller), then two integral_products. Robust and free of analytic gradients.
 *
 * NOTE on sign: the absolute Poisson G=0 gauge does not affect the GRADIENT, so
 * F_z is gauge-clean. The overall sign is validated by ENERGY CONSERVATION of the
 * run (E_elec + E_proj_KE + U_proj_bg flat) — the correctness gate in the spec.
 *
 * No inq/ or inq-study/ edit — wrapper-only.
 */
#ifndef INQKIT__DYNAMICS__PROJECTILE_FORCE
#define INQKIT__DYNAMICS__PROJECTILE_FORCE

#include <inq/inq.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density
#include <inqkit/jellium/gaussian_potential.hpp>             // gaussian_potential (direct erf/r)

namespace inqkit {
namespace dynamics {

// E_R at a shifted center: ∫ n_proj(·-center)·φ_drag.
template <class Field>
double drag_energy(Field const & phi_drag, inq::vector3<double> center, double sigma_pot) {
	auto nproj = inqkit::jellium::gaussian_density(phi_drag.basis(), center, sigma_pot);
	return inq::operations::integral_product(nproj, phi_drag);
}

// F_z = −dE_R/dR_z via a symmetric finite difference (step delta, Bohr).
// phi_drag = poisson(n_e − n_+) supplied by the caller (recomputed each step).
template <class Field>
double projectile_force_z(Field const & phi_drag, inq::vector3<double> center,
                          double sigma_pot, double delta = 0.05) {
	auto cp = center, cm = center;
	cp[2] += delta;
	cm[2] -= delta;
	const double e_plus  = drag_energy(phi_drag, cp, sigma_pot);
	const double e_minus = drag_energy(phi_drag, cm, sigma_pot);
	return -(e_plus - e_minus) / (2.0 * delta);
}

// ---------------------------------------------------------------------------
// INQ-NATIVE analytic Hellmann-Feynman force on a Gaussian-charge projectile.
//
// Byte-for-byte the integrand INQ uses for the LOCAL ionic force in
// inq/src/observables/forces_stress.hpp:182-187 — the DENSITY-GRADIENT form
//     F = − ∫ V_loc(r−R) ∇n(r) dr        (analytic ∇n via operations::gradient)
// reduced in covariant components then mapped to Cartesian with the same
// volume_element * cell.to_cartesian(...) as the engine. For a Gaussian charge the
// projectile's local potential is V_proj = poisson(n_proj) (all "long-range", no
// short-range split), so V_long+V_short → V_proj. This is NOT a finite difference;
// it is the identical HF formula INQ applies to a (local) pseudopotential ion, so a
// perturbation projectile and a ghost-UPF ion of the same V_proj feel the SAME force.
//
// density : the electron density field n(r) (basis::field<real_space,double>).
// cell    : the simulation cell (for covariant→Cartesian, = ions.cell()).
template <class DensityField, class CellType>
inq::vector3<double> projectile_force_analytic(DensityField const & density,
                                               CellType const & cell,
                                               inq::vector3<double> center,
                                               double sigma_pot) {
	using namespace inq;
	auto nproj    = inqkit::jellium::gaussian_density(density.basis(), center, sigma_pot);
	auto vproj    = solvers::poisson::solve(nproj);            // V_proj = poisson(n_proj)
	auto gdensity = operations::gradient(density);             // ∇n (covariant field)

	auto force_cov = -gpu::run(gpu::reduce(density.basis().local_size()),
		zero<vector3<double, covariant>>(),
		[vp = vproj.linear().cbegin(), gn = gdensity.linear().cbegin()] GPU_LAMBDA (auto ip) {
			return vp[ip] * gn[ip];
		});
	return density.basis().volume_element() * cell.to_cartesian(force_cov);
}

// Convenience: z-component only.
template <class DensityField, class CellType>
double projectile_force_analytic_z(DensityField const & density, CellType const & cell,
                                    inq::vector3<double> center, double sigma_pot) {
	return projectile_force_analytic(density, cell, center, sigma_pot)[2];
}

// DIRECT-potential variant: same HF integrand F = −∫ V_proj·∇n, but V_proj is the direct
// free-space erf/r potential (gaussian_potential) instead of poisson(gaussian_density) —
// consistent with the moving_gaussian_projectile_potential perturbation (no charge/Poisson/
// background). The force VALUE is essentially identical to the Poisson version (the gradient
// is insensitive to the background offset); this keeps the potential and force built the
// SAME way for the direct-perturbation runs.
template <class DensityField, class CellType>
inq::vector3<double> projectile_force_direct(DensityField const & density, CellType const & cell,
                                             inq::vector3<double> center, double sigma_pot) {
	using namespace inq;
	auto vproj    = inqkit::jellium::gaussian_potential(density.basis(), center, sigma_pot);
	auto gdensity = operations::gradient(density);
	auto force_cov = -gpu::run(gpu::reduce(density.basis().local_size()),
		zero<vector3<double, covariant>>(),
		[vp = vproj.linear().cbegin(), gn = gdensity.linear().cbegin()] GPU_LAMBDA (auto ip) {
			return vp[ip] * gn[ip];
		});
	return density.basis().volume_element() * cell.to_cartesian(force_cov);
}

template <class DensityField, class CellType>
double projectile_force_direct_z(DensityField const & density, CellType const & cell,
                                 inq::vector3<double> center, double sigma_pot) {
	return projectile_force_direct(density, cell, center, sigma_pot)[2];
}

} // namespace dynamics
} // namespace inqkit

#endif
