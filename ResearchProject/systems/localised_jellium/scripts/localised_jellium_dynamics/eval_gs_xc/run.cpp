// ============================================================================
// localised_jellium_dynamics / eval_gs_xc / run.cpp
//
// Single-point energy evaluation of the BARE jellium-slab ground state (NO
// projectile at all) to read its E_xc directly, for the WP-added-XC test
// (campaign localised-jellium-dynamics-analysis). Loads gs_p2_lz120, applies the
// background, does ONE RT step, and streams the full energy decomposition. The
// step-0 E_xc is the bare-GS E_xc (density unchanged from the load).
//
// Env: LJ_LX(50) LJ_LY(50) LJ_LZ(120) LJ_HALF(12.5) LJ_N(82) LJ_EDGE_W(0)
//      LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_GS_DIR(REQUIRED) LJ_OUT(gs_eval).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

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
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","gs_eval");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    std::cout<<"\n=== eval_gs_xc (BARE jellium GS, no projectile) ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);                       // NO projectile
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(GS_DIR);

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables.csv",sel); obs.write_header();
    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){obs.append(ctx);});
    auto step_fn=[&](auto const&data){rt.step(data);};
    auto opts=options::real_time{}.num_steps(1).dt(0.01*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert);

    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = eval_gs_xc (bare jellium GS)\nengine = inq-study\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  N = "<<N<<"\ngs_dir = "<<GS_DIR<<"\nrun_completed = true\n";}
    std::cout<<"  done\n"; return 0;
}
