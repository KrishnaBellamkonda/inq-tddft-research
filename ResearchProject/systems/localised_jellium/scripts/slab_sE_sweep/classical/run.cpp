// ============================================================================
// slab_bulk_limit / classical / run.cpp   (slab classical twin, WITH CAP)
//
// Clone of localised_jellium_dynamics/proj_dyn/run.cpp (moving Gaussian charge
// advanced by EHRENFEST dynamics; full pairwise energy decomposition + U_proj_bg
// + checkpoint/resume) with a TWO-SIDED sin^2 CAP added at the z-boundaries so the
// slab-face->CAP geometry is IDENTICAL to the WP twin (campaign slab_bulk_limit).
//
// E_absorbed (classical) = [E_total(t) - E_total(0)] + U_proj_bg, read in the
// plateau window (projectile clear of the slab, bath not yet at the CAP). The
// conserved quantity E_electronic + proj_ke + U_proj_bg is flat in the EARLY
// (pre-CAP) window -> the dt-accuracy gate for the smoke test.
//
// Env (adds LJ_CAP_ETA, LJ_CAP_L to the proj_dyn set):
//   LJ_LX(50) LJ_LY(50) LJ_LZ(70) LJ_HALF(5) LJ_N(32) LJ_EDGE_W(1)
//   LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_SIGMA(2.0) LJ_LAUNCH_Z(-15)
//   LJ_K0(2.71) LJ_MASS(1.0) LJ_DELTA(0.1) LJ_N_STEPS(2500) LJ_DT(0.04)
//   LJ_CAP_ETA(-1.0) LJ_CAP_L(10) LJ_SAVE_EVERY(0)
//   LJ_GS_DIR(REQUIRED) LJ_OUT(classical). Build vs INQ_SOURCE=inq-study.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density
#include <inqkit/jellium/interaction_energies.hpp>           // pairwise P/S/B decomposition
#include <inqkit/observables/orbital_kinetic_stats.hpp>       // all-orbital extensive (norm-corrected) kinetic
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/dynamics/projectile_force.hpp>
#include <inqkit/dynamics/moving_gaussian_projectile_perturbation.hpp>

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
static double read_state_d(const std::string& path,const char* key,double def){
    std::ifstream f(path); std::string line; std::string k=std::string(key)+"=";
    while(std::getline(f,line)){ auto p=line.find(k);
        if(p!=std::string::npos) return std::atof(line.substr(p+k.size()).c_str()); }
    return def;
}
static std::string tag6(int n){std::ostringstream o;o<<std::setw(6)<<std::setfill('0')<<n;return o.str();}

int main(){
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",70),HALF=env_d("LJ_HALF",5);
    const int N=env_i("LJ_N",32); const double EDGE_W=env_d("LJ_EDGE_W",1); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA_WP=env_d("LJ_SIGMA",2.0);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-15.0);
    const double K0=env_d("LJ_K0",2.71), MASS=env_d("LJ_MASS",1.0), DELTA=env_d("LJ_DELTA",0.1);
    const int N_STEPS=env_i("LJ_N_STEPS",2500); const double DT=env_d("LJ_DT",0.04);
    const int SAVE_EVERY=env_i("LJ_SAVE_EVERY",0);   // 0 = no density frames
    const int CKPT_EVERY=env_i("LJ_CKPT_EVERY",200); // interior RT checkpoints (rule checkpoint-dont-block)
    const double CAP_ETA=env_d("LJ_CAP_ETA",-1.0), CAP_L=env_d("LJ_CAP_L",10.0);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","classical");
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);
    const double N0=double(N)/(LX*LY*(2.0*HALF));
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}

    // ---- two-sided sin^2 CAP (bands at +/-CAP_MID_FRAC, width CAP_WIDTH_FRAC) ----
    const bool   CAP_ON        = (CAP_ETA != 0.0);
    const double CAP_WIDTH_FRAC = CAP_L / LZ;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC/2.0;
    const double z_cap_in       = LZ/2.0 - CAP_L;           // inner face of the +z CAP

    // ---- checkpoint / resume (rule final-timestep-checkpoint.md) ----
    const bool RESUME=env_i("LJ_RESUME",0)!=0;
    const std::string CKPT=OUT+"/checkpoint", RT_STATE=OUT+"/rt_state.txt";
    int START=0; double R0z=LAUNCH_Z, V0z=K0/MASS;
    if(RESUME){
        START=(int)read_state_d(RT_STATE,"last_step",-1);
        if(START<0){std::cerr<<"FATAL: LJ_RESUME=1 but no readable "<<RT_STATE<<"\n";return 2;}
        if(START>=N_STEPS){std::cout<<"Already at/after target ("<<START<<">="<<N_STEPS<<"); nothing to do.\n";return 0;}
        R0z=read_state_d(RT_STATE,"proj_z",LAUNCH_Z); V0z=read_state_d(RT_STATE,"proj_vz",K0/MASS);
    }
    const std::string SEG=(START>0)?(".from"+std::to_string(START)):std::string("");

    std::cout<<std::setprecision(12)<<"\n=== slab_bulk_limit/classical OUT="<<OUT<<" z0="<<LAUNCH_Z<<" k0="<<K0
             <<" CAP="<<(CAP_ON?"ON":"OFF")<<" eta="<<CAP_ETA<<" L="<<CAP_L<<"/face band=[+/-"<<z_cap_in<<","<<LZ/2.0<<"]"
             <<" START="<<START<<" -> N_STEPS="<<N_STEPS<<(RESUME?"  [RESUME]":"")<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);                       // NO projectile ion
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(RESUME?CKPT:GS_DIR);                     // resume from RT ckpt, else GS

    // background well + live projectile (Ehrenfest) + moving perturbation
    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    inqkit::dynamics::Projectile proj(MASS, -1.0,
        inqkit::detail::Vec3{0.0,0.0,R0z},
        inqkit::detail::Vec3{0.0,0.0,V0z});                 // resume state or v0=k0/m toward slab
    inqkit::dynamics::moving_gaussian_projectile_perturbation proj_pert(proj, SIGMA_POT);
    auto pert = perturbations::sum(bg_pert, proj_pert);

    // constant background fields for the force / U_proj_bg
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);          // φ₊

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables"+SEG+".csv",sel); obs.write_header();

    std::ofstream pj;   // energies in HARTREE (native), matching observables.csv
    if(electrons.root()){ pj.open(OUT+"/raw/observables/projectile"+SEG+".csv");
        pj<<std::setprecision(12)<<"step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal\n"; }

    // pairwise interaction-energy decomposition (P/S/B); E_BB constant computed once
    const double E_BB = inqkit::jellium::background_self_energy(nplus, phiplus);
    std::ofstream ix;   // Hartree units; sums close to INQ scalars (see interaction_energies.hpp)
    if(electrons.root()){ ix.open(OUT+"/raw/observables/interactions"+SEG+".csv");
        ix<<std::setprecision(12)<<"step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,norm_slab,norm_proj\n"; }

    // density frames (classical: slab electrons responding to the moving projectile)
    if(SAVE_EVERY>0) std::filesystem::create_directories(OUT+"/frames/total");
    auto save_frame=[&](int step){
        if(SAVE_EVERY<=0) return;
        auto n_tot = inqkit::fields::density::total(electrons);
        inqkit::io::RealField3DLayout lay{.field_name="density",.include_meta=false,
            .emit_raw=false,.emit_vti=true,.vti_format=inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter wr(OUT+"/frames/total",lay,{.overwrite=true});
        wr.write(n_tot,"density_t"+tag6(step));
    };

    // ALL-ORBITAL norm-division fix: under the CAP the bath orbitals also drain
    // (norm_total 64 -> 63.996 here), so INQ's per-particle-normalized kinetic is
    // slightly inflated on the classical half too. kin_bare_total = Σ occ_i T_i is
    // the extensive kinetic; E_corr = E_total - kinetic_INQ + kin_bare (§skill).
    inqkit::observables::OrbitalKineticStats extkin(
        OUT+"/raw/observables/orbital_kinetic_stats"+SEG+".csv",
        inqkit::observables::OrbitalKineticStatsConfig{.write_every=1});

    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const& ctx){
        obs.append(ctx);
        if(SAVE_EVERY>0 && ctx.step % SAVE_EVERY == 0) save_frame(ctx.step);
        // Ehrenfest projectile: HF force from the (electrons − background) field.
        auto Rn = proj.R();
        inq::vector3<double> center{Rn.x,Rn.y,Rn.z};
        auto phi_e = solvers::poisson::solve(electrons.density());          // φ_e = poisson(n_e)
        const double Fz = inqkit::dynamics::projectile_force_z(phi_e,   center, SIGMA_POT, DELTA)
                        - inqkit::dynamics::projectile_force_z(phiplus, center, SIGMA_POT, DELTA);
        // U_proj_bg(t) at the current center (ideal, r_cut-free)
        auto nproj = inqkit::jellium::gaussian_density(basis, center, SIGMA_POT);
        const double Uprojbg = -operations::integral_product(nproj, phiplus);
        // pairwise decomposition at the CURRENT step (n_slab = slab electrons; n_P = n_proj)
        auto ct = inqkit::jellium::compute_coulomb(electrons.density(), nproj, phiplus);
        proj.advance(inqkit::detail::Vec3{0.0,0.0,Fz}, DT);                 // V→V_n, R→R_{n+1}
        if(electrons.root()){
            pj<<ctx.step<<","<<ctx.time_au<<","<<Rn.z<<","<<proj.V().z<<","
              <<proj.ke()<<","<<Uprojbg<<"\n";
            ix<<ctx.step<<","<<ctx.time_au<<","<<ct.e_ss<<","<<ct.e_pp<<","<<ct.e_ps<<","
              <<ct.e_sb<<","<<ct.e_pb<<","<<E_BB<<","<<ct.norm_slab<<","<<ct.norm_p<<"\n";
        }
    });

    auto step_fn=[&](auto const&data){
        rt.step(data); extkin.maybe_accumulate(data);
        const int it=(int)data.iter();                     // interior RT checkpoint (rule checkpoint-dont-block)
        if(it>START && it%CKPT_EVERY==0 && it<N_STEPS){
            electrons.save(CKPT);                          // collective; overwrites the checkpoint
            if(electrons.root()){ std::ofstream st(RT_STATE,std::ios::trunc);
                st<<std::setprecision(12)<<"last_step="<<it<<"\ntime_au="<<(it*DT)<<"\ndt="<<DT
                  <<"\nproj_z="<<proj.R().z<<"\nproj_vz="<<proj.V().z
                  <<"\nproj_mass="<<MASS<<"\nproj_charge="<<-1.0<<"\n";
                std::cout<<"  [ckpt] step "<<it<<" (t="<<it*DT<<")\n"; }
        }
    };
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    if(CAP_ON){
        perturbations::absorbing cap_lo(CAP_ETA*1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
        perturbations::absorbing cap_hi(CAP_ETA*1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
        real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,
            perturbations::sum(pert, perturbations::sum(cap_lo,cap_hi)),START);
    } else {
        real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,pert,START);
    }

    // FINAL checkpoint (rule final-timestep-checkpoint.md)
    electrons.save(CKPT);
    if(electrons.root()){ std::ofstream st(RT_STATE);
        st<<std::setprecision(12)
          <<"last_step="<<N_STEPS<<"\ntime_au="<<(N_STEPS*DT)<<"\ndt="<<DT
          <<"\nproj_z="<<proj.R().z<<"\nproj_vz="<<proj.V().z
          <<"\nproj_mass="<<MASS<<"\nproj_charge="<<-1.0<<"\n"; }

    if(electrons.root()){ pj.close(); ix.close(); std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)
         <<"run = slab_bulk_limit/classical (moving Gaussian-charge projectile, Ehrenfest, CAP)\nengine = inq-study\n"
         <<"representation = perturbation\nprojectile = moving Gaussian charge (velocity-Verlet)\n"
         <<"periodicity = "<<PER<<"  Lx = "<<LX<<"  Ly = "<<LY<<"  Lz = "<<LZ<<"  half_width = "<<HALF<<"\n"
         <<"N = "<<N<<"  n0 = "<<N0<<"  edge_width = "<<EDGE_W<<"  spacing = "<<SPACING<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  sigma_wp = "<<SIGMA_WP<<"  sigma_pot = "<<SIGMA_POT<<"\n"
         <<"k0 = "<<K0<<"  mass = "<<MASS<<"  delta_fd = "<<DELTA<<"\n"
         <<"cap_on = "<<(CAP_ON?"yes":"no")<<"  cap_eta_ha = "<<CAP_ETA<<"  cap_width_bohr = "<<CAP_L<<" /face\n"
         <<"cap_band_hi_bohr = ["<<z_cap_in<<","<<LZ/2.0<<"]  cap_band_lo_bohr = ["<<-LZ/2.0<<","<<-z_cap_in<<"]\n"
         <<"start_step = "<<START<<"  n_steps = "<<N_STEPS<<"  dt = "<<DT
         <<"  proj_z_final = "<<proj.R().z<<"  proj_vz_final = "<<proj.V().z<<"\n"
         <<"gs_dir = "<<GS_DIR<<"\nrun_completed = true\n";}
    std::cout<<"  done  proj_z_final="<<proj.R().z<<" proj_vz_final="<<proj.V().z<<"\n"; return 0;
}
