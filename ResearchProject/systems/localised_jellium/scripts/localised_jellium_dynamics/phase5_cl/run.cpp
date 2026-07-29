// ============================================================================
// localised_jellium_dynamics / phase5_cl / run.cpp
//
// Phase 5 (screening): RT-propagate a stationary classical ghost near the slab
// and save per-step total densities for the screening GIFs. From the SAME GS as
// the WP partner. Saves frames/total/density_tNNNN.vti = electrons.density()
// every SAVE_EVERY steps (the ghost adds NO density — this is the bath only, so
// induced_CL(t) = n(t) − n(0) with n(0) = the bare GS). ETRS propagator.
//
// Env: as phase5_wp but no LJ_K0; adds LJ_PROJ_UPF (default full ghost UPF).
// Build vs INQ_SOURCE=inq-study.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
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
#include <sstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
static double env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static int    env_i(const char* k,int d){const char*v=std::getenv(k);return v?std::atoi(v):d;}
static std::string env_s(const char* k,const std::string&d){const char*v=std::getenv(k);return v?std::string(v):d;}
static std::string tag4(int n){std::ostringstream o;o<<std::setw(6)<<std::setfill('0')<<n;return o.str();}
static const char* PROJ_PSEUDO_DEFAULT=
  "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main(){
    auto t0=std::chrono::steady_clock::now();
    constexpr double M_PROJ=1.0/1822.8885;
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.5);
    const int N_STEPS=env_i("LJ_N_STEPS",500); const double DT=env_d("LJ_DT",0.01);
    const int SAVE_EVERY=env_i("LJ_SAVE_EVERY",25);
    const std::string PROJ_UPF=env_s("LJ_PROJ_UPF",PROJ_PSEUDO_DEFAULT);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","cl");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    std::cout<<"\n=== phase5_cl OUT="<<OUT<<" z="<<LAUNCH_Z<<" steps="<<N_STEPS<<" save="<<SAVE_EVERY<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto sp=ionic::species("H").pseudo_file(PROJ_UPF).mass(M_PROJ);
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
    std::filesystem::create_directories(OUT+"/frames/total");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables.csv",sel); obs.write_header();

    auto save_opt=[&](const std::string&nm){return inqkit::io::RealField3DLayout{.field_name=nm,
        .include_meta=false,.emit_raw=false,.emit_vti=true,
        .vti_format=inqkit::io::VTIWriteOptions::Format::binary};};
    auto save_frame=[&](int step){
        auto n_tot=inqkit::fields::density::total(electrons);
        using W=inqkit::io::RealField3DWriter;
        { W wr(OUT+"/frames/total",save_opt("density"),{.overwrite=true}); wr.write(n_tot,"density_t"+tag4(step)); }
    };

    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){
        obs.append(ctx);
        if(ctx.step % SAVE_EVERY == 0) save_frame(ctx.step);
    });
    auto step_fn=[&](auto const&data){rt.step(data);};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert);

    double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = localised_jellium_dynamics/phase5_cl/"<<env_s("LJ_OUT","cl")<<"\nengine = inq-study\n"
         <<"mode = classical ghost (at rest)\nperiodicity = "<<PER<<"  Lz = "<<LZ<<"  launch_z = "<<LAUNCH_Z
         <<"  r_from_face = "<<(-LAUNCH_Z-HALF)<<"\nproj_upf = "<<PROJ_UPF<<"\n"
         <<"n_steps = "<<N_STEPS<<"  dt = "<<DT<<"  save_every = "<<SAVE_EVERY<<"\n"
         <<"gs_dir = "<<GS_DIR<<"\nwall_time_s = "<<wall<<"\nrun_completed = true\n";}
    std::cout<<"  done wall="<<wall<<"s\n"; return 0;
}
