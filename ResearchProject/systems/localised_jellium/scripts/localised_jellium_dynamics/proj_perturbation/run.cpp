// ============================================================================
// localised_jellium_dynamics / proj_perturbation / run.cpp
//
// Classical projectile as a STATIONARY GAUSSIAN CHARGE PERTURBATION (no UPF ghost,
// no ion, no r_cut). The projectile's Poisson potential v_proj = +poisson(n_proj)
// is added to the KS potential alongside the background well, composed via
// inq::perturbations::sum(background, projectile). This removes the pseudopotential
// aliasing entirely so d(E_H+E_ext) − U_proj_bg can be read cleanly.
//
// Loads the bare jellium slab GS, applies (background + projectile), does N_STEPS
// RT steps, tabulates E_total/kin/H/xc/ext, and computes the CLEAN projectile↔
// background diagnostic U_proj_bg = −∫ n_proj·φ₊ (== the r_cut-invariant "ideal").
//
// Env: LJ_LX(50) LJ_LY(50) LJ_LZ(120) LJ_HALF(12.5) LJ_N(82) LJ_EDGE_W(0)
//      LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_SIGMA(0.5) LJ_LAUNCH_Z(-24.5)
//      LJ_N_STEPS(2) LJ_DT(0.01) LJ_GS_DIR(REQUIRED) LJ_OUT(proj_pert).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/gaussian_projectile_perturbation.hpp>

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
static std::string env_s(const char* k,const std::string&d){const char*v=std::getenv(k);return v?std::string(v):d;}

int main(){
    const double HA=27.211386;
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA_WP=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.5);
    const int N_STEPS=env_i("LJ_N_STEPS",2); const double DT=env_d("LJ_DT",0.01);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","proj_pert");
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);
    const double N0=double(N)/(LX*LY*(2.0*HALF));
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}

    std::cout<<std::setprecision(12)<<"\n=== proj_perturbation OUT="<<OUT<<" z="<<LAUNCH_Z<<" spacing="<<SPACING<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);                       // NO projectile ion
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(GS_DIR);

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);
    inqkit::jellium::gaussian_projectile_perturbation proj_pert({0.0,0.0,LAUNCH_Z}, SIGMA_POT);
    auto pert = perturbations::sum(bg_pert, proj_pert);

    // CLEAN U_proj_bg = −∫ n_proj·φ₊   (projectile −1 charge ↔ +background; == "ideal", r_cut-free)
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);
    auto nproj   = proj_pert.charge_density(basis);
    const double nproj_norm = operations::integral(nproj);
    const double U_proj_bg  = -operations::integral_product(nproj, phiplus);
    std::cout<<"  n_proj_norm = "<<nproj_norm<<"   U_proj_bg = "<<U_proj_bg*HA<<" eV (clean, no r_cut)\n";

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables.csv",sel); obs.write_header();
    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){obs.append(ctx);});
    auto step_fn=[&](auto const&data){rt.step(data);};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,pert);

    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)
         <<"run = proj_perturbation (Gaussian-charge classical projectile)\nengine = inq-study\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  N = "<<N<<"  spacing = "<<SPACING<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  sigma_wp = "<<SIGMA_WP<<"  sigma_pot = "<<SIGMA_POT<<"\n"
         <<"n_proj_norm = "<<nproj_norm<<"\n"
         <<"U_proj_bg_ha = "<<U_proj_bg<<"   U_proj_bg_ev = "<<U_proj_bg*HA<<"   (clean ideal, r_cut-free)\n"
         <<"gs_dir = "<<GS_DIR<<"\nrun_completed = true\n";}
    std::cout<<"  done\n"; return 0;
}
