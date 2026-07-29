// ============================================================================
// Phase 2 Test C — RUN C1: INQ NATIVE Ehrenfest ion dynamics.
//
// A mass-1 ghost-UPF ion (H symbol, ghost_sigma0p354.upf: z_valence=0, 0
// projectors, local V_loc = +erf(r/(sqrt2.sigma_pot))/r) is launched at v0 along
// +z from one side of a FINITE box, past a fixed neutral He atom at the centre.
// Ion motion is driven by INQ's OWN native Ehrenfest integrator
// (options::real_time{}.ehrenfest() => ionic::propagator::molecular_dynamics,
// velocity-Verlet with a = F_localHF / species.mass(), advanced inside ETRS).
//
// The mass is set host-side via .mass(1.0/1822.8885) so species.mass() == 1.0
// atomic unit (= m_e). The ghost carries no electrons (z_valence=0) so its only
// coupling to the KS system is its local +erf/r potential — exactly the same
// V_proj the perturbation twin (C2) uses.
//
// Records each step: ions.positions()[ghost].z, ions.velocities()[ghost].z, and
// the electronic energy breakdown, to native.csv:
//   step,time,z,vz,E_elec,E_total,E_kin,E_hartree,E_external,E_nonlocal,E_xc,E_ion
// E_elec := kinetic+hartree+external+nonlocal+xc (electronic only, cross-run comparable).
//
// GATE: the ghost ion MUST move (z changes). If ion_dynamics=EHRENFEST does not
// move a z_valence=0 ghost, that is the decisive finding — reported by z(t) flat.
//
// Env: TC_L(35) TC_SPACING(0.4) TC_SIGMA_POT(0.35355) TC_LAUNCH_Z(-8.0)
//      TC_V0(1.2) TC_MASS_AU(1.0) TC_N_STEPS(400) TC_DT(0.02)
//      TC_GHOST_UPF(ghost_sigma0p354.upf) TC_HE_UPF(optional)
// No inq/ or inq-study/ edit — wrapper-only.
// ============================================================================
#include <inq/inq.hpp>

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
	const std::string GHOST_UPF=env_s("TC_GHOST_UPF","ghost_sigma0p354.upf");
	const std::string HE_UPF   =env_s("TC_HE_UPF","");

	std::cout<<std::setprecision(12)
	         <<"\n=== C1 NATIVE Ehrenfest  FINITE  L="<<L<<" dx="<<SPACING
	         <<"  launch_z="<<LAUNCH_Z<<" v0="<<V0<<" mass_au="<<MASS_AU
	         <<" N="<<N_STEPS<<" dt="<<DT<<"  ghost="<<GHOST_UPF<<" ===\n";

	auto cell=systems::cell::orthorhombic(L*1.0_b,L*1.0_b,L*1.0_b).finite();
	auto ions=systems::ions(cell);

	// Neutral He at centre = fixed density source the projectile passes through.
	auto he = HE_UPF.empty() ? ionic::species("He") : ionic::species("He").pseudo_file(HE_UPF);
	ions.insert(he, {0.0_b, 0.0_b, 0.0_b});               // index 0

	// Ghost projectile ion: mass 1.0 a.u. (= m_e) via .mass(1/1822.8885 amu).
	auto ghost = ionic::species("H").pseudo_file(GHOST_UPF).mass(MASS_AU/1822.8885);
	ions.insert(ghost, {0.0_b, 0.0_b, LAUNCH_Z*1.0_b});   // index 1
	const int GH = 1;

	// Initial velocity: +z toward/past He.
	ions.velocities()[GH] = vector3<double>{0.0, 0.0, V0};

	auto electrons=systems::electrons(ions,
		options::electrons{}.spacing(SPACING*1.0_b),
		input::kpoints::gamma());
	ground_state::initial_guess(ions, electrons);
	auto gs = ground_state::calculate(ions, electrons, options::theory{}.lda(),
		options::ground_state{}.energy_tolerance(1e-8_Ha));
	std::cout<<"  GS energy = "<<gs.energy.total()<<" Ha\n";

	std::filesystem::create_directories("results");
	std::ofstream csv;
	if(electrons.root()){
		csv.open("results/native.csv"); csv<<std::setprecision(12)
		 <<"step,time,z,vz,E_elec,E_total,E_kin,E_hartree,E_external,E_nonlocal,E_xc,E_ion\n";
	}

	auto func=[&](auto data){
		auto pos = data.positions();
		auto vel = data.velocities();
		auto en  = data.energy();
		const double zc = pos[GH][2], vz = vel[GH][2];
		const double E_elec = en.kinetic()+en.hartree()+en.external()+en.non_local()+en.xc();
		if(electrons.root()){
			csv<<data.iter()<<","<<data.time()<<","<<zc<<","<<vz<<","
			   <<E_elec<<","<<en.total()<<","<<en.kinetic()<<","<<en.hartree()<<","
			   <<en.external()<<","<<en.non_local()<<","<<en.xc()<<","<<en.ion()<<"\n";
		}
	};

	auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime).ehrenfest();
	real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts);

	const double z_final = ions.positions()[GH][2], vz_final = ions.velocities()[GH][2];
	if(electrons.root()){
		csv.close();
		std::ofstream s("results/run_summary.txt");
		s<<std::setprecision(12)<<"run = phase2_native_ehrenfest/c1_native\n"
		 <<"scheme = INQ native Ehrenfest (ionic::propagator::molecular_dynamics, velocity-Verlet in ETRS)\n"
		 <<"engine = inq\nboundary = finite (periodicity 0)\n"
		 <<"L_bohr = "<<L<<"  spacing = "<<SPACING<<"  sigma_pot = "<<SIGMA_POT<<"\n"
		 <<"launch_z = "<<LAUNCH_Z<<"  v0 = "<<V0<<"  mass_au = "<<MASS_AU<<"\n"
		 <<"n_steps = "<<N_STEPS<<"  dt = "<<DT<<"\n"
		 <<"z_final = "<<z_final<<"  vz_final = "<<vz_final<<"\n"
		 <<"ghost_upf = "<<GHOST_UPF<<"\ncsv = results/native.csv\nrun_completed = true\n";
	}
	std::cout<<"  done  z_final="<<z_final<<" vz_final="<<vz_final<<"\n";
	return 0;
}
