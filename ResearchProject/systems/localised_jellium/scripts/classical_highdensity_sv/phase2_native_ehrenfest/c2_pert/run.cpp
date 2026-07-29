// ============================================================================
// Phase 2 Test C — RUN C2: PERTURBATION projectile (our hand-rolled scheme).
//
// IDENTICAL physical setup to C1 (native): FINITE box, same He atom fixed at the
// centre, same dx/dt/N_STEPS/launch_z/v0/mass/sigma_pot. The ONLY difference is
// the projectile representation:
//   - NO ghost ion in the KS system.
//   - A rigid Gaussian charge realised as a moving_gaussian_projectile_perturbation
//     (adds +poisson(n_proj), the same +erf/r the ghost UPF carries).
//   - Advanced each RT step by our velocity-Verlet inqkit::dynamics::Projectile,
//     using the ANALYTIC Hellmann-Feynman force
//     inqkit::dynamics::projectile_force_analytic_z(density, cell, center, sigma_pot)
//     (F = -int V_proj . grad n) — the exact local HF force, NOT the FD variant.
//
// This is our whole hand-rolled Ehrenfest scheme (force + velocity-Verlet + our
// callback intra-step ordering). Comparing native.csv (C1) vs pert.csv (C2)
// isolates the one remaining difference: intra-step ordering (O(dt)).
//
// Records each step to pert.csv (same columns as native.csv):
//   step,time,z,vz,E_elec,E_total,E_kin,E_hartree,E_external,E_nonlocal,E_xc,E_ion
//
// Env: same TC_* as C1 (shared by the dispatcher).
// No inq/ or inq-study/ edit — wrapper-only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/dynamics/projectile_force.hpp>
#include <inqkit/dynamics/moving_gaussian_projectile_perturbation.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
static double env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static int    env_i(const char* k,int d){const char*v=std::getenv(k);return v?std::atoi(v):d;}
static std::string env_s(const char* k,const std::string& d){const char*v=std::getenv(k);return v?std::string(v):d;}

int main(){
	const double L=env_d("TC_L",35), SPACING=env_d("TC_SPACING",0.4);
	const double SIGMA_POT=env_d("TC_SIGMA_POT",0.35355);
	const double LAUNCH_Z=env_d("TC_LAUNCH_Z",-8.0), V0=env_d("TC_V0",1.2);
	const double MASS_AU=env_d("TC_MASS_AU",1.0);        // projectile mass in atomic units (m_e)
	const int    N_STEPS=env_i("TC_N_STEPS",400);
	const double DT=env_d("TC_DT",0.02);
	const std::string HE_UPF=env_s("TC_HE_UPF","");

	std::cout<<std::setprecision(12)
	         <<"\n=== C2 PERTURBATION  FINITE  L="<<L<<" dx="<<SPACING
	         <<"  launch_z="<<LAUNCH_Z<<" v0="<<V0<<" mass_au="<<MASS_AU
	         <<" N="<<N_STEPS<<" dt="<<DT<<"  sigma_pot="<<SIGMA_POT<<" ===\n";

	auto cell=systems::cell::orthorhombic(L*1.0_b,L*1.0_b,L*1.0_b).finite();
	auto ions=systems::ions(cell);

	// SAME neutral He at centre; NO projectile ion.
	auto he = HE_UPF.empty() ? ionic::species("He") : ionic::species("He").pseudo_file(HE_UPF);
	ions.insert(he, {0.0_b, 0.0_b, 0.0_b});

	auto electrons=systems::electrons(ions,
		options::electrons{}.spacing(SPACING*1.0_b),
		input::kpoints::gamma());
	ground_state::initial_guess(ions, electrons);
	auto gs = ground_state::calculate(ions, electrons, options::theory{}.lda(),
		options::ground_state{}.energy_tolerance(1e-8_Ha));
	std::cout<<"  GS energy = "<<gs.energy.total()<<" Ha\n";

	// Live projectile (mass 1 a.u., charge -1) + moving Gaussian perturbation.
	inqkit::dynamics::Projectile proj(MASS_AU, -1.0,
		inqkit::detail::Vec3{0.0,0.0,LAUNCH_Z},
		inqkit::detail::Vec3{0.0,0.0,V0});
	inqkit::dynamics::moving_gaussian_projectile_perturbation proj_pert(proj, SIGMA_POT);

	std::filesystem::create_directories("results");
	std::ofstream csv;
	if(electrons.root()){
		csv.open("results/pert.csv"); csv<<std::setprecision(12)
		 <<"step,time,z,vz,E_elec,E_total,E_kin,E_hartree,E_external,E_nonlocal,E_xc,E_ion\n";
	}

	auto func=[&](auto data){
		auto en = data.energy();
		// record the CURRENT center (R_n, the center used for this step's density)
		auto Rn = proj.R();
		const double zc = Rn.z, vz_before = proj.V().z;
		const double E_elec = en.kinetic()+en.hartree()+en.external()+en.non_local()+en.xc();
		// ANALYTIC HF force at R_n on the post-step electron density.
		inq::vector3<double> center{Rn.x, Rn.y, Rn.z};
		double Fz = inqkit::dynamics::projectile_force_analytic_z(
			electrons.density(), ions.cell(), center, SIGMA_POT);
		proj.advance(inqkit::detail::Vec3{0.0,0.0,Fz}, DT);   // V->V_n, R->R_{n+1}
		if(electrons.root()){
			// vz recorded = the integer-step velocity V_n after completion kick.
			csv<<data.iter()<<","<<data.time()<<","<<zc<<","<<proj.V().z<<","
			   <<E_elec<<","<<en.total()<<","<<en.kinetic()<<","<<en.hartree()<<","
			   <<en.external()<<","<<en.non_local()<<","<<en.xc()<<","<<en.ion()<<"\n";
			(void)vz_before;
		}
	};

	auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
	real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts, proj_pert);

	if(electrons.root()){
		csv.close();
		std::ofstream s("results/run_summary.txt");
		s<<std::setprecision(12)<<"run = phase2_native_ehrenfest/c2_pert\n"
		 <<"scheme = perturbation projectile (analytic HF force + inqkit velocity-Verlet, callback ordering)\n"
		 <<"engine = inq\nboundary = finite (periodicity 0)\n"
		 <<"L_bohr = "<<L<<"  spacing = "<<SPACING<<"  sigma_pot = "<<SIGMA_POT<<"\n"
		 <<"launch_z = "<<LAUNCH_Z<<"  v0 = "<<V0<<"  mass_au = "<<MASS_AU<<"\n"
		 <<"n_steps = "<<N_STEPS<<"  dt = "<<DT<<"\n"
		 <<"z_final = "<<proj.R().z<<"  vz_final = "<<proj.V().z<<"\n"
		 <<"csv = results/pert.csv\nrun_completed = true\n";
	}
	std::cout<<"  done  z_final="<<proj.R().z<<" vz_final="<<proj.V().z<<"\n";
	return 0;
}
