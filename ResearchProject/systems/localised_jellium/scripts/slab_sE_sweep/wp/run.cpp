// ============================================================================
// slab_sE_sweep / wp / run.cpp   (slab wavepacket twin, WITH CAP + σ_r(t))
//
// Identical to slab_bulk_limit/wp/run.cpp EXCEPT it also writes the per-step WP
// real-space moments (wp_real_space_stats.csv) so the transverse/z spread σ_r(t)
// can be time-averaged (campaign slab_sE_sweep, S(E) over 5 energies × 5 L_z).
//
// Clone of localised_jellium_dynamics/phase5_wp/run.cpp (inject a σ_WP Gaussian
// wavepacket as an extra KS orbital; full pairwise P/S/B decomposition with the
// closure checks; checkpoint/resume) with:
//   * a MOVING packet: k0 = m*v0 passed via LJ_K0 (phase5_wp defaulted k0=0),
//   * a TWO-SIDED sin^2 CAP at the z-boundaries, geometry IDENTICAL to the
//     classical twin, so the transmitted WP is absorbed (no reflection).
//
// E_absorbed (WP) = |E_total(t) - E_total(0)| read AFTER the CAP absorbs the
// transmitted WP + excited density. INQ reports a normalized scalar under a CAP
// (reference_inq_reports_normalized_energy) -> cross-check against the extensive
// bath decomposition Δ(e_ss) in interactions.csv. ETRS propagator (CAP absorption;
// CN would renormalise the WP each step). Build vs INQ_SOURCE=inq-study.
//
// Env: LJ_LX(50) LJ_LY(50) LJ_LZ(70) LJ_HALF(5) LJ_N(32) LJ_EDGE_W(1)
//      LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_SIGMA(2.0) LJ_LAUNCH_Z(-15)
//      LJ_K0(2.71) LJ_N_STEPS(2500) LJ_DT(0.04) LJ_SAVE_EVERY(50)
//      LJ_CAP_ETA(-1.0) LJ_CAP_L(10) LJ_GS_DIR(REQUIRED) LJ_OUT(wp).
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
#include <inqkit/observables/orbital_kinetic_stats.hpp> // all-orbital extensive (norm-corrected) kinetic
#include <inqkit/observables/wp_real_space_stats.hpp>    // per-step WP σ_r(t): ⟨r⟩, σ_x/σ_y/σ_z (+circular)

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
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",70),HALF=env_d("LJ_HALF",5);
    const int N=env_i("LJ_N",32); const double EDGE_W=env_d("LJ_EDGE_W",1); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA=env_d("LJ_SIGMA",2.0);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-15.0), K0=env_d("LJ_K0",2.71);
    const int N_STEPS=env_i("LJ_N_STEPS",2500); const double DT=env_d("LJ_DT",0.04);
    const int SAVE_EVERY=env_i("LJ_SAVE_EVERY",50);
    const int CKPT_EVERY=env_i("LJ_CKPT_EVERY",200);   // interior RT checkpoints (rule checkpoint-dont-block)
    const double CAP_ETA=env_d("LJ_CAP_ETA",-1.0), CAP_L=env_d("LJ_CAP_L",10.0);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","wp");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    // ---- two-sided sin^2 CAP ----
    const bool   CAP_ON        = (CAP_ETA != 0.0);
    const double CAP_WIDTH_FRAC = CAP_L / LZ;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC/2.0;
    const double z_cap_in       = LZ/2.0 - CAP_L;

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

    std::cout<<"\n=== slab_bulk_limit/wp OUT="<<OUT<<" z="<<LAUNCH_Z<<" k0="<<K0<<" sigma="<<SIGMA
             <<" CAP="<<(CAP_ON?"ON":"OFF")<<" eta="<<CAP_ETA<<" L="<<CAP_L<<"/face"
             <<" START="<<START<<" -> "<<N_STEPS<<(RESUME?"  [RESUME]":"")<<" save="<<SAVE_EVERY<<" ===\n";
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

    // ALL-ORBITAL norm-division fix (rule: INQ reports per-particle-normalized
    // kinetic; under the CAP a draining orbital pins the reported kinetic at the
    // remnant's mean -> inflates energy_total). kin_bare_total = Σ occ_i T_i is the
    // extensive kinetic (no norm division). E_corr = E_total - kinetic_INQ + kin_bare.
    inqkit::observables::OrbitalKineticStats extkin(
        OUT+"/raw/observables/orbital_kinetic_stats"+SEG+".csv",
        inqkit::observables::OrbitalKineticStatsConfig{.write_every=1});

    // WP real-space moments: ⟨r⟩(t) + per-axis spread σ_x,σ_y,σ_z (naive + circular).
    // σ_r(t) = sqrt(σ_x²+σ_y²) (transverse radial width) is time-averaged offline;
    // σ_z tracked separately. wp_idx is valid on resume (read from rt_state).
    inqkit::observables::WPRealSpaceStats wprs(
        OUT+"/raw/observables/wp_real_space_stats"+SEG+".csv", wp_idx,
        inqkit::observables::WPRealSpaceStatsConfig{.write_every=1});

    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){
        obs.append(ctx);
        if(SAVE_EVERY>0 && ctx.step % SAVE_EVERY == 0) save_frame(ctx.step);
        auto n_wp = inqkit::jellium::orbital_density_field(electrons, wp_idx);
        auto ct   = inqkit::jellium::compute_coulomb_wp(electrons.density(), n_wp, ie_phiplus);
        if(electrons.root())
            ix<<ctx.step<<","<<ctx.time_au<<","<<ct.e_ss<<","<<ct.e_pp<<","<<ct.e_ps<<","
              <<ct.e_sb<<","<<ct.e_pb<<","<<E_BB<<","<<ct.e_hartree_check<<","
              <<ct.e_external_check<<","<<ct.norm_wp<<","<<ct.norm_total<<"\n";
    });
    auto step_fn=[&](auto const&data){
        rt.step(data); extkin.maybe_accumulate(data); wprs.maybe_accumulate(data);
        const int it=(int)data.iter();                     // interior RT checkpoint (rule checkpoint-dont-block)
        if(it>START && it%CKPT_EVERY==0 && it<N_STEPS){
            electrons.save(CKPT);                          // collective; overwrites the checkpoint
            if(electrons.root()){ std::ofstream st(RT_STATE,std::ios::trunc);
                st<<std::setprecision(12)<<"last_step="<<it<<"\ntime_au="<<(it*DT)
                  <<"\ndt="<<DT<<"\nwp_idx="<<wp_idx<<"\n";
                std::cout<<"  [ckpt] step "<<it<<" (t="<<it*DT<<")\n"; }
        }
    };
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    if(CAP_ON){
        perturbations::absorbing cap_lo(CAP_ETA*1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
        perturbations::absorbing cap_hi(CAP_ETA*1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
        real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,
            perturbations::sum(bg_pert, perturbations::sum(cap_lo,cap_hi)),START);
    } else {
        real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert,START);
    }

    // FINAL checkpoint (rule final-timestep-checkpoint.md): the checkpoint holds the
    // propagated WP so LJ_RESUME=1 + a larger LJ_N_STEPS extends without recomputing.
    electrons.save(CKPT);
    if(electrons.root()){ std::ofstream st(RT_STATE);
        st<<std::setprecision(12)<<"last_step="<<N_STEPS<<"\ntime_au="<<(N_STEPS*DT)
          <<"\ndt="<<DT<<"\nwp_idx="<<wp_idx<<"\n"; }

    if(electrons.root()) ix.close();
    double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = slab_bulk_limit/wp/"<<env_s("LJ_OUT","wp")<<"\nengine = inq-study\n"
         <<"mode = wavepacket (k0="<<K0<<")\nperiodicity = "<<PER<<"  Lx = "<<LX<<"  Ly = "<<LY<<"  Lz = "<<LZ<<"\n"
         <<"half_width = "<<HALF<<"  N = "<<N<<"  n0 = "<<N0<<"  edge_width = "<<EDGE_W<<"  spacing = "<<SPACING<<"\n"
         <<"sigma_WP = "<<SIGMA<<"  k0 = "<<K0<<"  launch_z = "<<LAUNCH_Z<<"  r_from_face = "<<(-LAUNCH_Z-HALF)<<"\n"
         <<"wp_idx = "<<wp_idx<<"  wp_norm_after = "<<wp_norm_after<<"\n"
         <<"cap_on = "<<(CAP_ON?"yes":"no")<<"  cap_eta_ha = "<<CAP_ETA<<"  cap_width_bohr = "<<CAP_L<<" /face\n"
         <<"cap_band_hi_bohr = ["<<z_cap_in<<","<<LZ/2.0<<"]  cap_band_lo_bohr = ["<<-LZ/2.0<<","<<-z_cap_in<<"]\n"
         <<"start_step = "<<START<<"  n_steps = "<<N_STEPS<<"  dt = "<<DT<<"  save_every = "<<SAVE_EVERY<<"\n"
         <<"gs_dir = "<<GS_DIR<<"\nwall_time_s = "<<wall<<"\nrun_completed = true\n";}
    std::cout<<"  done wall="<<wall<<"s\n"; return 0;
}
