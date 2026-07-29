// ============================================================================
// classical_highdensity_sv / force_native_analytic / run.cpp
//
// TEST 1 (closed-form check). Validates the INQ-NATIVE analytic Hellmann-Feynman
// force inqkit::dynamics::projectile_force_analytic_z (density-gradient form,
//   F = -int V_proj . grad n dr,  V_proj = poisson(gaussian_density))
// against BOTH:
//   (a) the existing finite-difference operator projectile_force_z, and
//   (b) the closed-form two-Gaussian Coulomb force (compared in python).
//
// Setup (identical geometry to force_test): a FIXED Gaussian "density" n at z=0
// (width SIGMA_S, unit norm) stands in for the electron density; the projectile
// Gaussian (width SIGMA_POT) is swept along z. For each center z_c:
//   F_ours = projectile_force_analytic_z(n, cell, {0,0,z_c}, SIGMA_POT)
//   F_fd   = projectile_force_z(poisson(n), {0,0,z_c}, SIGMA_POT, DELTA)
// Written to results/force_native_analytic.csv as z_c,F_ours,F_fd.
//
// FINITE box (periodicity 0) => free-space Poisson => closed form applies.
// Env: FT_L(85) FT_SPACING(0.5) FT_SIGMA_POT(0.35355) FT_SIGMA_S(0.5)
//      FT_DELTA(0.05) FT_ZMIN(1.0) FT_ZMAX(14.0) FT_ZSTEP(0.5)
// No inq/ or inq-study/ edit -- wrapper-only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density
#include <inqkit/dynamics/projectile_force.hpp>              // projectile_force_z, _analytic_z

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

	std::cout<<std::setprecision(10)<<"\n=== force_native_analytic  FINITE(periodicity 0)  L="<<L
	         <<" dx="<<SPACING<<"  sigma_pot="<<SIGMA_POT<<" sigma_s="<<SIGMA_S
	         <<" delta="<<DELTA<<" ===\n";

	auto cell=systems::cell::orthorhombic(35.0_b,35.0_b,L*1.0_b).finite();
	auto ions=systems::ions(cell);
	auto electrons=systems::electrons(ions,
		options::electrons{}.spacing(SPACING*1.0_b).extra_electrons(2),
		input::kpoints::gamma());
	ground_state::initial_guess(ions, electrons);
	auto basis=electrons.density().basis();

	// Fixed source Gaussian "density" at origin; phi_drag = poisson(n_source) for FD.
	auto n_source=inqkit::jellium::gaussian_density(basis, {0.0,0.0,0.0}, SIGMA_S);
	const double src_norm=operations::integral(n_source);
	auto phi_drag=solvers::poisson::solve(n_source);
	std::cout<<"  source norm (should be 1) = "<<src_norm<<"\n";

	std::filesystem::create_directories("results");
	std::ofstream csv("results/force_native_analytic.csv");
	csv<<std::setprecision(12)<<"z_c,F_ours,F_fd\n";
	for(double zc=ZMIN; zc<=ZMAX+1e-9; zc+=ZSTEP){
		inq::vector3<double> center{0.0,0.0,zc};
		const double F_ours=inqkit::dynamics::projectile_force_analytic_z(n_source, cell, center, SIGMA_POT);
		const double F_fd  =inqkit::dynamics::projectile_force_z(phi_drag, center, SIGMA_POT, DELTA);
		csv<<zc<<","<<F_ours<<","<<F_fd<<"\n";
		std::cout<<"  z_c="<<std::setw(6)<<zc<<"  F_ours="<<std::setw(14)<<F_ours
		         <<"  F_fd="<<F_fd<<"\n";
	}
	if(electrons.root()){ std::ofstream s("results/run_summary.txt");
		s<<std::setprecision(10)<<"run = classical_highdensity_sv/force_native_analytic\n"
		 <<"test = native_analytic_vs_fd_vs_closed_form\nengine = inq\nboundary = finite (periodicity 0)\n"
		 <<"L_bohr = "<<L<<"  spacing = "<<SPACING<<"\n"
		 <<"sigma_pot = "<<SIGMA_POT<<"  sigma_s = "<<SIGMA_S<<"  delta_fd = "<<DELTA<<"\n"
		 <<"source_norm = "<<src_norm<<"\ncsv = results/force_native_analytic.csv\nrun_completed = true\n"; }
	std::cout<<"Done.\n"; return 0;
}
