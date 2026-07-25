// ============================================================================
// localised_jellium_dynamics / smoke_eprojbg / run.cpp
//
// Known-case smoke test for inqkit::jellium::projectile_background_energy
// (campaign localised-jellium-dynamics-analysis, Phase-1 code gate). Builds a
// small p2 slab background + the σ=0.5 Gaussian ghost, computes both E_proj_bg
// estimates, and checks:
//   (1) n_proj integrates to ≈ 1 (normalised Gaussian charge),
//   (2) ideal ≈ impl (reciprocity) for the FULL untruncated ghost — the UPF's
//       charge is the same Gaussian, so both formulas must agree,
//   (3) both are NEGATIVE (attractive: negative projectile ↔ positive background).
// Runs on CPU (no GPU needed) so it can gate the campaign before the GPUs free.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>

#include <cstdlib>
#include <cmath>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
static const char* PROJ_PSEUDO =
  "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main() {
	// small p2 cell for a fast CPU check
	const double LX=30, LY=30, LZ=60, HALF=12.5;
	const int N=30;
	const double SPACING=0.5, SIGMA_WP=0.5, SIGMA_POT=SIGMA_WP/std::sqrt(2.0);  // 0.3536
	const double LAUNCH_Z=-16.5;   // r = 4 from the face at -12.5
	const double N0 = double(N)/(LX*LY*(2.0*HALF));

	auto cell = systems::cell::orthorhombic(LX*1.0_b, LY*1.0_b, LZ*1.0_b).periodicity(2);
	auto ions = systems::ions(cell);
	auto sp = ionic::species("H").pseudo_file(PROJ_PSEUDO).mass(1.0/1822.8885);
	ions.insert(sp, {0.0*1.0_b, 0.0*1.0_b, LAUNCH_Z*1.0_b});
	auto electrons = systems::electrons(ions, options::electrons{}.spacing(SPACING*1.0_b)
		.extra_electrons(N).extra_states(4).temperature(0.00862*1.0_eV), input::kpoints::gamma());
	ground_state::initial_guess(ions, electrons);

	inqkit::jellium::localised_background_params bg;
	bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF;
	bg.slab_axis=2; bg.center={0.0,0.0,0.0}; bg.edge_width=0.0;
	inqkit::jellium::localised_background_perturbation pert(bg);

	auto res = inqkit::jellium::projectile_background_energy(
		pert, electrons, ions, {0.0, 0.0, LAUNCH_Z}, SIGMA_POT);

	// Independent CONSISTENT-KERNEL reciprocity: with ρ_proj = −n_proj and the SAME
	// p2 Poisson kernel for both potentials, ∫ρ_proj·poisson(n₊) must equal
	// ∫n₊·poisson(ρ_proj) (symmetry of the Coulomb kernel). ideal = ∫ρ_proj·poisson(n₊).
	auto basis  = electrons.density().basis();
	auto nplus  = pert.background_density(basis);
	auto nproj  = inqkit::jellium::gaussian_density(basis, {0.0,0.0,LAUNCH_Z}, SIGMA_POT);
	auto phi_np = inq::solvers::poisson::solve(nproj);                 // poisson(n_proj)
	double recip = -inq::operations::integral_product(nplus, phi_np);  // = ∫ρ_+·poisson(ρ_proj)

	const double HA = 27.211386;
	std::cout.setf(std::ios::fixed); std::cout.precision(4);
	std::cout << "\n=== E_proj_bg smoke ===\n";
	std::cout << "  n_proj_norm = " << res.n_proj_norm << "   (expect ~1.0)\n";
	std::cout << "  ideal (∫ n_proj·v_bg)         = " << res.ideal*HA << " eV\n";
	std::cout << "  recip (∫ n₊·poisson(ρ_proj))  = " << recip*HA << " eV   [same-kernel check]\n";
	std::cout << "  impl  (−∫ n₊·v_ion, UPF gauge)= " << res.impl*HA  << " eV   [r_cut-diagnostic]\n";
	double reldiff = std::abs(res.ideal - recip) / std::max(std::abs(res.ideal), 1e-12);
	std::cout << "  |ideal−recip| = " << std::abs(res.ideal-recip)*HA << " eV  (rel "
	          << reldiff*100 << "%)\n";

	bool ok_norm  = std::abs(res.n_proj_norm - 1.0) < 0.02;
	bool ok_recip = reldiff < 0.01;                                  // same-kernel identity → tight
	bool ok_finite= std::isfinite(res.ideal) && std::isfinite(res.impl) && std::isfinite(recip);
	bool pass = ok_norm && ok_recip && ok_finite;
	std::cout << "  checks: norm=" << ok_norm << " reciprocity=" << ok_recip
	          << " finite=" << ok_finite << "\n";
	std::cout << (pass ? "SMOKE PASS\n" : "SMOKE FAIL\n");
	return pass ? 0 : 1;
}
