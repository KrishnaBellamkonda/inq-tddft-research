// ============================================================================
// localised_jellium / scripts/campaign_autorun / wp / run.cpp
//
// ENV-PARAMETRISED stationary/moving WP E_total(0) run for H4 (and H0-style).
// Loads a matching GS (LJ_GS_DIR), injects a WP at z=LJ_LAUNCH_Z with k0=LJ_K0,
// re-applies the background well, propagates a few steps; E_total(0) from step 0.
//
// Env: LJ_LX(50) LJ_LY(50) LJ_LZ(120) LJ_HALF(12.5) LJ_N(82) LJ_EDGE_W(0)
//      LJ_PERIODICITY(3|2) LJ_SPACING(0.5) LJ_SIGMA(0.5) LJ_LAUNCH_Z LJ_K0(0)
//      LJ_GS_DIR(REQUIRED) LJ_OUT(REQUIRED) LJ_N_STEPS(2) LJ_DT(0.01).
// Build against INQ_SOURCE=inq-study; runtime shares from inq/install/share.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
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

int main(){
    auto t0=std::chrono::steady_clock::now();
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",3);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-16.5), K0=env_d("LJ_K0",0.0);
    const int N_STEPS=env_i("LJ_N_STEPS",2); const double DT=env_d("LJ_DT",0.01);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","wp");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    std::cout<<"\n=== autorun wp OUT="<<OUT<<" per="<<PER<<" Lz="<<LZ<<" z="<<LAUNCH_Z<<" k0="<<K0<<" sig="<<SIGMA<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(GS_DIR);

    auto wp=inqkit::WavePacket{}.center(0.0,0.0,LAUNCH_Z).sigma(SIGMA).k0(0.0,0.0,K0)
                .orthogonalise_against_occupied(electrons);
    auto report=wp.inject_into_last_extra_state(electrons,1.0);
    std::cout<<"  WP idx="<<report.state_index<<" norm_after="<<report.norm_after<<"\n";

    // Screening / WP-potential test (LJ_SAVE_DENSITY=1): write the t=0 densities
    // BEFORE propagation so the injected WP source charge is captured exactly.
    //   density_wp   = |psi_WP|^2               (the WP's own charge, ~norm 1)
    //   density_total= electrons.density()      (slab + WP)
    //   density_bath = total - |psi_WP|^2       (slab only; screening baseline)
    // n_WP feeds the Coulomb-potential comparison vs the classical Gaussian ghost.
    if(env_i("LJ_SAVE_DENSITY",0)==1){
        auto n_wp   = inqkit::fields::density::orbital(electrons,report.state_index);
        auto n_tot  = inqkit::fields::density::total(electrons);
        auto n_bath = inqkit::fields::density::total_excluding_orbital(n_tot,n_wp,1.0);
        using W=inqkit::io::RealField3DWriter;
        auto opt=[&](const std::string&nm){return inqkit::io::RealField3DLayout{.field_name=nm,
            .include_meta=false,.emit_raw=false,.emit_vti=true,
            .vti_format=inqkit::io::VTIWriteOptions::Format::binary};};
        std::filesystem::create_directories(OUT+"/density_wp");
        { W wr(OUT+"/density_wp",   opt("density"),{.overwrite=true}); wr.write(n_wp,  "density_wp"); }
        std::filesystem::create_directories(OUT+"/density_total");
        { W wr(OUT+"/density_total",opt("density"),{.overwrite=true}); wr.write(n_tot, "density_total"); }
        std::filesystem::create_directories(OUT+"/density_bath");
        { W wr(OUT+"/density_bath", opt("density"),{.overwrite=true}); wr.write(n_bath,"density_bath"); }
        std::cout<<"  saved t=0 densities (wp/total/bath)\n";
    }

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true; sel.density_l2=false;
    // full energy decomposition (kinetic/Hartree/xc rise on WP insertion; stream the rest too)
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=sel.energy_ion_kinetic=true;
    sel.energy_exact_exchange=sel.energy_nvxc=sel.energy_eigenvalues=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables.csv",sel); obs.write_header();
    std::ofstream nlog(OUT+"/raw/observables/electron_number.csv"); nlog<<"step,time_au,N_total\n";
    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){obs.append(ctx);});
    auto step_fn=[&](auto const&data){rt.step(data); if(data.root()) nlog<<data.iter()<<","<<(data.iter()*DT)<<","<<data.num_electrons()<<"\n";};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert);

    const double zp=3.0/(4.0*SIGMA*SIGMA);
    double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = campaign_autorun/wp/"<<env_s("LJ_OUT","wp")<<"\nengine = inq-study\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  sigma_WP = "<<SIGMA<<"  k0 = "<<K0<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  r_from_face = "<<(-LAUNCH_Z-HALF)<<"\n"
         <<"wp_norm_after = "<<report.norm_after<<"  zero_point_ha = "<<zp<<"\n"
         <<"gs_dir = "<<GS_DIR<<"  dt = "<<DT<<"  n_steps = "<<N_STEPS<<"\nwall_time_s = "<<wall<<"\nrun_completed = true\n";}
    std::cout<<"  done wall="<<wall<<"s\n"; return 0;
}
