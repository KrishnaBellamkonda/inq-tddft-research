// ============================================================================
// classical_highdensity_sv / force_vs_native / run.cpp
//
// TEST 2 (DECISIVE). Proves that our INQ-native analytic Hellmann-Feynman force
// operator inqkit::dynamics::projectile_force_analytic (F = -int V_proj . grad n)
// equals INQ's OWN native local ionic force on a real ghost-UPF ion carrying the
// SAME Gaussian local potential V_loc(r) = erf(r/(sqrt2.sigma_pot))/r.
//
// A ghost ion (H symbol, ghost_sigma0p354.upf: z_valence=0, 0 projectors, local
// potential = +erf/r) contributes to the KS Hamiltonian a PURE local potential and
// zero electrons. INQ's ground-state force on it (result.forces[ghost]) is therefore
// exactly the local Hellmann-Feynman term -int V_loc . grad n. We build a nontrivial
// density from a neutral He atom at the origin, place the ghost at z_c, run the GS,
// and compare F_INQ (native) to F_ours (our analytic operator on electrons.density()).
//
// PASS = F_ours.z ~ F_INQ.z to <~1% for several z_c.
//
// Env: FV_L(35) FV_SPACING(0.5) FV_SIGMA_POT(0.35355)
//      FV_ZC list is hard-coded (3.0, 3.5, 4.0, 5.0).
//      FV_GHOST_UPF (path to ghost_sigma0p354.upf) FV_HE_UPF (optional He pseudo)
// No inq/ or inq-study/ edit -- wrapper-only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density
#include <inqkit/dynamics/projectile_force.hpp>              // projectile_force_analytic

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace inq;
using namespace inq::magnitude;
static double env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static std::string env_s(const char* k,const std::string& d){const char*v=std::getenv(k);return v?std::string(v):d;}

int main(){
	const double L=env_d("FV_L",35), SPACING=env_d("FV_SPACING",0.5);
	const double SIGMA_POT=env_d("FV_SIGMA_POT",0.35355);
	const std::string GHOST_UPF=env_s("FV_GHOST_UPF","ghost_sigma0p354.upf");
	const std::string HE_UPF   =env_s("FV_HE_UPF","");   // empty => let inq resolve "He"

	const std::vector<double> zcs = {3.0, 3.5, 4.0, 5.0};

	std::cout<<std::setprecision(10)<<"\n=== force_vs_native  FINITE(periodicity 0)  L="<<L
	         <<" dx="<<SPACING<<"  sigma_pot="<<SIGMA_POT<<"  ghost="<<GHOST_UPF<<" ===\n";

	std::filesystem::create_directories("results");
	std::ofstream csv("results/force_vs_native.csv");
	csv<<std::setprecision(12)<<"z_c,F_INQ_x,F_INQ_y,F_INQ_z,F_ours_x,F_ours_y,F_ours_z,ratio_z\n";

	for(double zc : zcs){
		// FINITE box; He at origin (real density), ghost at {0,0,zc}.
		auto cell=systems::cell::orthorhombic(L*1.0_b,L*1.0_b,L*1.0_b).finite();
		auto ions=systems::ions(cell);

		auto he = HE_UPF.empty() ? ionic::species("He") : ionic::species("He").pseudo_file(HE_UPF);
		ions.insert(he, {0.0_b, 0.0_b, 0.0_b});
		// ghost projectile ion (z_valence=0, 0 projectors, local +erf/r)
		auto ghost = ionic::species("H").pseudo_file(GHOST_UPF);
		ions.insert(ghost, {0.0_b, 0.0_b, zc*1.0_b});
		const int GHOST_IDX = 1;

		auto electrons=systems::electrons(ions,
			options::electrons{}.spacing(SPACING*1.0_b),
			input::kpoints::gamma());
		ground_state::initial_guess(ions, electrons);
		auto result = ground_state::calculate(ions, electrons,
			options::theory{}.lda(),
			options::ground_state{}.calculate_forces());

		inq::vector3<double> F_INQ{ result.forces[GHOST_IDX][0],
		                            result.forces[GHOST_IDX][1],
		                            result.forces[GHOST_IDX][2] };

		auto F_ours = inqkit::dynamics::projectile_force_analytic(
			electrons.density(), cell, {0.0,0.0,zc}, SIGMA_POT);

		const double ratio = F_INQ[2]!=0.0 ? F_ours[2]/F_INQ[2] : 0.0;

		std::cout<<"  z_c="<<zc<<"\n"
		         <<"    F_INQ  = ("<<F_INQ[0]<<", "<<F_INQ[1]<<", "<<F_INQ[2]<<")\n"
		         <<"    F_ours = ("<<F_ours[0]<<", "<<F_ours[1]<<", "<<F_ours[2]<<")\n"
		         <<"    ratio_z (ours/INQ) = "<<ratio<<"\n";
		csv<<zc<<","<<F_INQ[0]<<","<<F_INQ[1]<<","<<F_INQ[2]<<","
		   <<F_ours[0]<<","<<F_ours[1]<<","<<F_ours[2]<<","<<ratio<<"\n";
	}

	{ std::ofstream s("results/run_summary.txt");
	  s<<std::setprecision(10)<<"run = classical_highdensity_sv/force_vs_native\n"
	   <<"test = analytic_force_vs_inq_native_ghost_upf\nengine = inq\n"
	   <<"boundary = finite (periodicity 0)\nL_bohr = "<<L<<"  spacing = "<<SPACING<<"\n"
	   <<"sigma_pot = "<<SIGMA_POT<<"  ghost_upf = "<<GHOST_UPF<<"\n"
	   <<"csv = results/force_vs_native.csv\nrun_completed = true\n"; }
	std::cout<<"Done.\n"; return 0;
}
