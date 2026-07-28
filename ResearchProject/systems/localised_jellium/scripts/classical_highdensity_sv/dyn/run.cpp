// ============================================================================
// classical_slab_stopping / run.cpp
//   Classical stopping baseline for the LOCALISED jellium slab (twin of the WP
//   qsp_phase5 S(E) sweep). Cloned from localised_jellium_dynamics/proj_dyn.
//
// MOVING classical projectile = a rigid Gaussian charge represented as a moving
// perturbation (no ghost UPF, no r_cut aliasing). Two drive modes:
//   LJ_CONST_V=0 (Phase 1, Ehrenfest): each RT step compute the Hellmann-Feynman
//     force from the (electrons − background) field, velocity-Verlet advance —
//     a LIGHT electron decelerates; S(v0) = initial-drag slope.
//   LJ_CONST_V=1 (Phase 2, prescribed constant velocity): zero force ⇒ a=0 ⇒
//     R=R0+V0·t, V≡V0 — external drive; S = ΔE_deposited / L_slab over the
//     transit. Size N_STEPS to stop before the Gaussian wraps the periodic box.
//
// Emits: full energy ledger (observables.csv) + a projectile trajectory
// (projectile.csv: step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal)
// that twin-run-analysis merges. Conserved quantity (correctness gate):
// E_electronic + energy_proj_ke + U_proj_bg must be flat.
//
// Env: LJ_LX(50) LJ_LY(50) LJ_LZ(120) LJ_HALF(12.5) LJ_N(82) LJ_EDGE_W(0)
//      LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_SIGMA(0.5) LJ_LAUNCH_Z(-24.5)
//      LJ_K0(1.0) LJ_MASS(1.0) LJ_DELTA(0.1) LJ_N_STEPS(50) LJ_DT(0.05)
//      LJ_GS_DIR(REQUIRED) LJ_OUT(classical_slab_stopping) LJ_CONST_V(0)
//      LJ_SAVE_EVERY(0) LJ_RESUME(0).
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
    const double HA=27.211386;
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA_WP=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.5);
    const double K0=env_d("LJ_K0",1.0), MASS=env_d("LJ_MASS",1.0), DELTA=env_d("LJ_DELTA",0.1);
    const int N_STEPS=env_i("LJ_N_STEPS",50); const double DT=env_d("LJ_DT",0.05);
    const int SAVE_EVERY=env_i("LJ_SAVE_EVERY",0);   // 0 = no density frames
    const bool CONST_V=env_i("LJ_CONST_V",0)!=0;     // 0=Ehrenfest (force-driven), 1=prescribed const velocity
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","classical_slab_stopping");
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);
    const double N0=double(N)/(LX*LY*(2.0*HALF));
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}

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

    std::cout<<std::setprecision(12)<<"\n=== classical_slab_stopping OUT="<<OUT<<" drive="<<(CONST_V?"const_v":"ehrenfest")
             <<" z0="<<LAUNCH_Z<<" k0="<<K0
             <<" START="<<START<<" -> N_STEPS="<<N_STEPS<<(RESUME?"  [RESUME]":"")<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);                       // NO projectile ion
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(24).temperature(0.00862*1.0_eV),input::kpoints::gamma());
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

    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const& ctx){
        obs.append(ctx);
        if(SAVE_EVERY>0 && ctx.step % SAVE_EVERY == 0) save_frame(ctx.step);
        // Ehrenfest projectile: HF force from the (electrons − background) field.
        // projectile_force_z is linear in the field, so F(φ_e−φ₊)=F(φ_e)−F(φ₊)
        // (no field subtraction needed; φ_e re-solved each step, φ₊ constant).
        auto Rn = proj.R();
        inq::vector3<double> center{Rn.x,Rn.y,Rn.z};
        // Ehrenfest: HF force from (electrons − background). CONST_V: zero force ⇒
        // velocity-Verlet with a=0 ⇒ R=R0+V0·t, V≡V0 (prescribed constant velocity).
        double Fz = 0.0;
        if(not CONST_V){
            auto phi_e = solvers::poisson::solve(electrons.density());       // φ_e = poisson(n_e)
            Fz = inqkit::dynamics::projectile_force_z(phi_e,   center, SIGMA_POT, DELTA)
               - inqkit::dynamics::projectile_force_z(phiplus, center, SIGMA_POT, DELTA);
        }
        // U_proj_bg(t) at the current center (ideal, r_cut-free)
        auto nproj = inqkit::jellium::gaussian_density(basis, center, SIGMA_POT);
        const double Uprojbg = -operations::integral_product(nproj, phiplus);
        // pairwise decomposition at the CURRENT step (n_slab = slab electrons; n_P = n_proj)
        auto ct = inqkit::jellium::compute_coulomb(electrons.density(), nproj, phiplus);
        proj.advance(inqkit::detail::Vec3{0.0,0.0,Fz}, DT);                 // V→V_n, R→R_{n+1}
        if(electrons.root()){
            pj<<ctx.step<<","<<ctx.time_au<<","<<Rn.z<<","<<proj.V().z<<","
              <<proj.ke()<<","<<Uprojbg<<"\n";   // Hartree (engine converts uniformly)
            ix<<ctx.step<<","<<ctx.time_au<<","<<ct.e_ss<<","<<ct.e_pp<<","<<ct.e_ps<<","
              <<ct.e_sb<<","<<ct.e_pb<<","<<E_BB<<","<<ct.norm_slab<<","<<ct.norm_p<<"\n";
        }
    });

    auto step_fn=[&](auto const&data){rt.step(data);};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,pert,START);

    // FINAL checkpoint (rule final-timestep-checkpoint.md): save RT state so this
    // run can be EXTENDED to more steps via LJ_RESUME=1 + a larger LJ_N_STEPS.
    electrons.save(CKPT);
    if(electrons.root()){ std::ofstream st(RT_STATE);
        st<<std::setprecision(12)
          <<"last_step="<<N_STEPS<<"\ntime_au="<<(N_STEPS*DT)<<"\ndt="<<DT
          <<"\nproj_z="<<proj.R().z<<"\nproj_vz="<<proj.V().z
          <<"\nproj_mass="<<MASS<<"\nproj_charge="<<-1.0<<"\n"; }

    if(electrons.root()){ pj.close(); ix.close(); std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)
         <<"run = proj_dyn (moving Gaussian-charge classical projectile, Ehrenfest)\nengine = inq-study\n"
         <<"representation = perturbation\nprojectile = moving Gaussian charge (velocity-Verlet)\n"
         <<"drive = "<<(CONST_V?"const_velocity":"ehrenfest")<<"\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  N = "<<N<<"  spacing = "<<SPACING<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  sigma_wp = "<<SIGMA_WP<<"  sigma_pot = "<<SIGMA_POT<<"\n"
         <<"k0 = "<<K0<<"  mass = "<<MASS<<"  delta_fd = "<<DELTA<<"\n"
         <<"start_step = "<<START<<"  n_steps = "<<N_STEPS<<"  dt = "<<DT
         <<"  proj_z_final = "<<proj.R().z<<"  proj_vz_final = "<<proj.V().z<<"\n"
         <<"gs_dir = "<<GS_DIR<<"\nrun_completed = true\n";}
    std::cout<<"  done  proj_z_final="<<proj.R().z<<" proj_vz_final="<<proj.V().z<<"\n"; return 0;
}
