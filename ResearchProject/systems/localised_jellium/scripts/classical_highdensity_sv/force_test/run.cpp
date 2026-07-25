// ============================================================================
// classical_highdensity_sv / force_test / run.cpp
//
// Phase 2 Test A (independent force validation). Validates the Hellmann-Feynman
// force operator inqkit::dynamics::projectile_force_z (finite-difference of the
// on-grid Poisson interaction integral) against the CLOSED-FORM two-Gaussian
// Coulomb force — the one independent, analytic anchor for the whole classical
// projectile machinery.
//
// Setup: a FIXED Gaussian SOURCE charge (width sigma_s, unit norm) at z=0, and the
// projectile Gaussian (width sigma_pot) swept along z. phi_drag = poisson(n_source);
// for each projectile center z_c we record
//   E_num(z_c) = integral( n_proj(.-z_c) * phi_drag )        [drag_energy]
//   F_num(z_c) = -d E_num/dz_c  via projectile_force_z         [the operator tested]
// The notebook compares to the analytic two-Gaussian E(d), F(d).
//
// FINITE box (periodicity 0) => free-space Poisson => the closed form applies.
// Env: FT_L(85) FT_SPACING(0.5) FT_SIGMA_POT(0.35355) FT_SIGMA_S(0.5)
//      FT_DELTA(0.05) FT_ZMIN(1.0) FT_ZMAX(14.0) FT_ZSTEP(0.5)
// No inq/ or inq-study/ edit — wrapper-only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density
#include <inqkit/dynamics/projectile_force.hpp>              // drag_energy, projectile_force_z

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
static double env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}

int main(){
	const double L=env_d("FT_L",85), SPACING=env_d("FT_SPACING",0.5);
	const double SIGMA_POT=env_d("FT_SIGMA_POT",0.35355), SIGMA_S=env_d("FT_SIGMA_S",0.5);
	const double DELTA=env_d("FT_DELTA",0.05);
	const double ZMIN=env_d("FT_ZMIN",1.0), ZMAX=env_d("FT_ZMAX",14.0), ZSTEP=env_d("FT_ZSTEP",0.5);

	std::cout<<std::setprecision(10)<<"\n=== force_test  FINITE(periodicity 0)  L="<<L
	         <<" dx="<<SPACING<<"  sigma_pot="<<SIGMA_POT<<" sigma_s="<<SIGMA_S
	         <<" delta="<<DELTA<<" ===\n";

	// Finite (open) box => free-space Poisson, so the two-Gaussian closed form applies.
	auto cell=systems::cell::orthorhombic(L*1.0_b,L*1.0_b,L*1.0_b).finite();
	auto ions=systems::ions(cell);
	auto electrons=systems::electrons(ions,
		options::electrons{}.spacing(SPACING*1.0_b).extra_electrons(2),
		input::kpoints::gamma());
	ground_state::initial_guess(ions, electrons);
	auto basis=electrons.density().basis();

	// Fixed source Gaussian charge at origin; phi_drag = poisson(n_source).
	auto n_source=inqkit::jellium::gaussian_density(basis, {0.0,0.0,0.0}, SIGMA_S);
	const double src_norm=operations::integral(n_source);
	auto phi_drag=solvers::poisson::solve(n_source);
	std::cout<<"  source norm (should be 1) = "<<src_norm<<"\n";

	std::filesystem::create_directories("results");
	std::ofstream csv("results/force_test.csv");
	csv<<std::setprecision(12)<<"z_c,d,E_num,F_num\n";
	for(double zc=ZMIN; zc<=ZMAX+1e-9; zc+=ZSTEP){
		inq::vector3<double> center{0.0,0.0,zc};
		const double E_num=inqkit::dynamics::drag_energy(phi_drag, center, SIGMA_POT);
		const double F_num=inqkit::dynamics::projectile_force_z(phi_drag, center, SIGMA_POT, DELTA);
		csv<<zc<<","<<zc<<","<<E_num<<","<<F_num<<"\n";
		std::cout<<"  z_c="<<std::setw(6)<<zc<<"  E_num="<<std::setw(12)<<E_num
		         <<"  F_num="<<F_num<<"\n";
	}
	if(electrons.root()){ std::ofstream s("results/run_summary.txt");
		s<<std::setprecision(10)<<"run = classical_highdensity_sv/force_test\n"
		 <<"test = analytic_two_gaussian_force\nengine = inq\nboundary = finite (periodicity 0)\n"
		 <<"L_bohr = "<<L<<"  spacing = "<<SPACING<<"\n"
		 <<"sigma_pot = "<<SIGMA_POT<<"  sigma_s = "<<SIGMA_S<<"  delta_fd = "<<DELTA<<"\n"
		 <<"source_norm = "<<src_norm<<"\ncsv = results/force_test.csv\nrun_completed = true\n"; }
	std::cout<<"Done.\n"; return 0;
}
