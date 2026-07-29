// ============================================================================
// localised_jellium_dynamics / phase5_wp / run.cpp
//
// Phase 5 (screening): RT-propagate an injected wavepacket AT REST near the slab
// and save per-step densities for the screening GIFs. From the SAME GS as the
// classical partner. Saves every SAVE_EVERY steps:
//   frames/total/density_tNNNN.vti  = electrons.density()      (slab + WP)
//   frames/wp/density_tNNNN.vti     = |psi_WP|^2               (WP orbital)
// Bath = total − wp; induced_bath(t) = bath(t) − bath(0) computed in the notebook.
// ETRS propagator (CN renormalises the WP each step — reference_inq_propagator).
//
// Env: LJ_LX(50) LJ_LY(50) LJ_LZ(120) LJ_HALF(12.5) LJ_N(82) LJ_EDGE_W(0)
//      LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_SIGMA(0.5) LJ_LAUNCH_Z(-24.5,r=12)
//      LJ_K0(0) LJ_GS_DIR(REQUIRED) LJ_OUT(REQUIRED) LJ_N_STEPS(500) LJ_DT(0.01)
//      LJ_SAVE_EVERY(25). Build vs INQ_SOURCE=inq-study.
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
#include <inqkit/jellium/interaction_energies.hpp>   // WP pairwise P/S/B decomposition

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
static double read_state_d(const std::string& path,const char* key,double def){
    std::ifstream f(path); std::string line; std::string k=std::string(key)+"=";
    while(std::getline(f,line)){ auto p=line.find(k);
        if(p!=std::string::npos) return std::atof(line.substr(p+k.size()).c_str()); }
    return def;
}

int main(){
    auto t0=std::chrono::steady_clock::now();
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.5), K0=env_d("LJ_K0",0.0);
    const int N_STEPS=env_i("LJ_N_STEPS",500); const double DT=env_d("LJ_DT",0.01);
    const int SAVE_EVERY=env_i("LJ_SAVE_EVERY",25);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","wp");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    // ---- checkpoint / resume (rule final-timestep-checkpoint.md) ----
    const bool RESUME=env_i("LJ_RESUME",0)!=0;
    const std::string CKPT=OUT+"/checkpoint", RT_STATE=OUT+"/rt_state.txt";
    int START=0;
    if(RESUME){
        START=(int)read_state_d(RT_STATE,"last_step",-1);
        if(START<0){std::cerr<<"FATAL: LJ_RESUME=1 but no readable "<<RT_STATE<<"\n";return 2;}
        if(START>=N_STEPS){std::cout<<"Already at/after target ("<<START<<">="<<N_STEPS<<"); nothing to do.\n";return 0;}
    }
    const std::string SEG=(START>0)?(".from"+std::to_string(START)):std::string("");

    std::cout<<"\n=== phase5_wp OUT="<<OUT<<" z="<<LAUNCH_Z<<" k0="<<K0<<" START="<<START
             <<" -> "<<N_STEPS<<(RESUME?"  [RESUME]":"")<<" save="<<SAVE_EVERY<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(RESUME?CKPT:GS_DIR);                     // resume: ckpt already holds the propagated WP

    int wp_idx; double wp_norm_after=1.0;
    if(RESUME){                                             // do NOT re-inject — WP is in the checkpoint
        wp_idx=(int)read_state_d(RT_STATE,"wp_idx",-1);
        std::cout<<"  [RESUME] loaded ckpt at step "<<START<<"  wp_idx="<<wp_idx<<"\n";
    } else {
        auto wp=inqkit::WavePacket{}.center(0.0,0.0,LAUNCH_Z).sigma(SIGMA).k0(0.0,0.0,K0)
                    .orthogonalise_against_occupied(electrons);
        auto report=wp.inject_into_last_extra_state(electrons,1.0);
        wp_idx=report.state_index; wp_norm_after=report.norm_after;
        std::cout<<"  WP idx="<<wp_idx<<" norm_after="<<wp_norm_after<<"\n";
    }

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    std::filesystem::create_directories(OUT+"/raw/observables");
    std::filesystem::create_directories(OUT+"/frames/total");
    std::filesystem::create_directories(OUT+"/frames/wp");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables"+SEG+".csv",sel); obs.write_header();

    // pairwise interaction decomposition (P=WP, S=slab, B=background); E_BB once
    auto ie_basis   = electrons.density().basis();
    auto ie_nplus   = bg_pert.background_density(ie_basis);
    auto ie_phiplus = solvers::poisson::solve(ie_nplus);
    const double E_BB = inqkit::jellium::background_self_energy(ie_nplus, ie_phiplus);
    std::ofstream ix;   // Hartree; e_ss+e_ps+e_pp must == E_hartree, e_sb+e_pb == E_external
    if(electrons.root()){ ix.open(OUT+"/raw/observables/interactions"+SEG+".csv");
        ix<<std::setprecision(12)<<"step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
            "e_hartree_check,e_external_check,norm_wp,norm_total\n"; }

    auto save_opt=[&](const std::string&nm){return inqkit::io::RealField3DLayout{.field_name=nm,
        .include_meta=false,.emit_raw=false,.emit_vti=true,
        .vti_format=inqkit::io::VTIWriteOptions::Format::binary};};
    auto save_frame=[&](int step){
        auto n_wp  = inqkit::fields::density::orbital(electrons,wp_idx);
        auto n_tot = inqkit::fields::density::total(electrons);
        using W=inqkit::io::RealField3DWriter;
        { W wr(OUT+"/frames/total",save_opt("density"),{.overwrite=true}); wr.write(n_tot,"density_t"+tag4(step)); }
        { W wr(OUT+"/frames/wp",   save_opt("density"),{.overwrite=true}); wr.write(n_wp, "density_t"+tag4(step)); }
    };

    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){
        obs.append(ctx);
        if(ctx.step % SAVE_EVERY == 0) save_frame(ctx.step);
        auto n_wp = inqkit::jellium::orbital_density_field(electrons, wp_idx);
        auto ct   = inqkit::jellium::compute_coulomb_wp(electrons.density(), n_wp, ie_phiplus);
        if(electrons.root())
            ix<<ctx.step<<","<<ctx.time_au<<","<<ct.e_ss<<","<<ct.e_pp<<","<<ct.e_ps<<","
              <<ct.e_sb<<","<<ct.e_pb<<","<<E_BB<<","<<ct.e_hartree_check<<","
              <<ct.e_external_check<<","<<ct.norm_wp<<","<<ct.norm_total<<"\n";
    });
    auto step_fn=[&](auto const&data){rt.step(data);};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert,START);

    // FINAL checkpoint (rule final-timestep-checkpoint.md): the checkpoint holds the
    // propagated WP so LJ_RESUME=1 + a larger LJ_N_STEPS extends without recomputing.
    electrons.save(CKPT);
    if(electrons.root()){ std::ofstream st(RT_STATE);
        st<<std::setprecision(12)<<"last_step="<<N_STEPS<<"\ntime_au="<<(N_STEPS*DT)
          <<"\ndt="<<DT<<"\nwp_idx="<<wp_idx<<"\n"; }

    if(electrons.root()) ix.close();
    double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = localised_jellium_dynamics/phase5_wp/"<<env_s("LJ_OUT","wp")<<"\nengine = inq-study\n"
         <<"mode = wavepacket (k0="<<K0<<")\nperiodicity = "<<PER<<"  Lz = "<<LZ<<"  sigma_WP = "<<SIGMA<<"  k0 = "<<K0<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  r_from_face = "<<(-LAUNCH_Z-HALF)<<"  wp_idx = "<<wp_idx<<"\n"
         <<"wp_norm_after = "<<wp_norm_after<<"\n"
         <<"start_step = "<<START<<"  n_steps = "<<N_STEPS<<"  dt = "<<DT<<"  save_every = "<<SAVE_EVERY<<"\n"
         <<"gs_dir = "<<GS_DIR<<"\nwall_time_s = "<<wall<<"\nrun_completed = true\n";}
    std::cout<<"  done wall="<<wall<<"s\n"; return 0;
}
