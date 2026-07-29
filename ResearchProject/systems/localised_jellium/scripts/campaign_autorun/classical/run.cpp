// ============================================================================
// localised_jellium / scripts/campaign_autorun / classical / run.cpp
//
// ENV-PARAMETRISED stationary matched-ghost E_total(0) run for H5 (and H0-style).
// Loads a matching GS (LJ_GS_DIR), inserts a stationary classical Gaussian ghost
// (electron_gaussian_wpsigma0p5.upf, z_valence 0) at z=LJ_LAUNCH_Z, re-applies
// the background, propagates a few steps; E_total(0) from step 0. The omitted
// int v_ghost*n_+ is re-added in analysis (H5).
//
// Env: LJ_LX LJ_LY LJ_LZ LJ_HALF LJ_N LJ_EDGE_W LJ_PERIODICITY(3|2) LJ_SPACING
//      LJ_LAUNCH_Z LJ_GS_DIR(REQUIRED) LJ_OUT(REQUIRED) LJ_N_STEPS(2) LJ_DT(0.01).
// Build against INQ_SOURCE=inq-study; runtime shares from inq/install/share.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include <chrono>
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
static const char* PROJ_PSEUDO=
  "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main(){
    auto t0=std::chrono::steady_clock::now();
    constexpr double M_PROJ=1.0/1822.8885;
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",3);
    const double SPACING=env_d("LJ_SPACING",0.5), LAUNCH_Z=env_d("LJ_LAUNCH_Z",-16.5);
    const int N_STEPS=env_i("LJ_N_STEPS",2); const double DT=env_d("LJ_DT",0.01);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","cl");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    std::cout<<"\n=== autorun classical OUT="<<OUT<<" per="<<PER<<" Lz="<<LZ<<" z="<<LAUNCH_Z<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto sp=ionic::species("H").pseudo_file(env_s("LJ_PROJ_UPF",PROJ_PSEUDO)).mass(M_PROJ);
    ions.insert(sp,{0.0*1.0_b,0.0*1.0_b,LAUNCH_Z*1.0_b});
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(GS_DIR);
    ions.velocities()[0]=vector3<double>{0.0,0.0,0.0};

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true; sel.density_l2=false;
    // full energy decomposition (external term makes the ghost's V_ext change explicit)
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=sel.energy_ion_kinetic=true;
    sel.energy_exact_exchange=sel.energy_nvxc=sel.energy_eigenvalues=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables.csv",sel); obs.write_header();
    std::ofstream nlog(OUT+"/raw/observables/electron_number.csv"); nlog<<"step,time_au,N_total\n";
    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){obs.append(ctx);});
    auto step_fn=[&](auto const&data){rt.step(data); if(data.root()) nlog<<data.iter()<<","<<(data.iter()*DT)<<","<<data.num_electrons()<<"\n";};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert);

    double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = campaign_autorun/classical/"<<env_s("LJ_OUT","cl")<<"\nengine = inq-study\n"
         <<"projectile = ghost wpsigma0p5 (z_valence 0), stationary\n"
         <<"proj_upf = "<<env_s("LJ_PROJ_UPF",PROJ_PSEUDO)<<"\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  launch_z = "<<LAUNCH_Z<<"  r_from_face = "<<(-LAUNCH_Z-HALF)<<"\n"
         <<"ghost_background_term_omitted = true\n"
         <<"gs_dir = "<<GS_DIR<<"  dt = "<<DT<<"  n_steps = "<<N_STEPS<<"\nwall_time_s = "<<wall<<"\nrun_completed = true\n";}
    std::cout<<"  done wall="<<wall<<"s\n"; return 0;
}
