// ============================================================================
// classical_highdensity_sv / pilot_direct / run.cpp   (Phase 3 pilot — RUN A')
//
// DIRECT-POTENTIAL A/B variant of the Phase-3 pilot (RUN A). IDENTICAL in every
// physical parameter to scripts/classical_highdensity_sv/pilot/run.cpp except the
// projectile representation, which is the DIRECT free-space erf/r potential (no
// charge density, no Poisson solve, no periodic neutralizing background) instead
// of the Poisson-based moving_gaussian_projectile_perturbation.
//
// EXACTLY TWO changes vs pilot/run.cpp:
//   (a) proj_pert = inqkit::dynamics::moving_gaussian_projectile_potential
//       (adds gaussian_potential = erf(|r-R|/(sqrt2 sigma))/|r-R| directly) instead
//       of moving_gaussian_projectile_perturbation (poisson(gaussian_density)).
//   (b) Ehrenfest force: projectile_force_direct_z(...) (same HF integrand
//       −∫V_proj·∇n but V_proj = the direct erf/r field) instead of
//       projectile_force_analytic_z(...), for BOTH the electrons and the
//       background term:  F = force_direct(n_e) − force_direct(n_+).
//
// Motivation: the Poisson representation carries the periodic G=0 neutralizing
// background + transverse images, whose offset lurches abruptly as the clipped
// charge crosses the far box face (proj_z ~ +12.5..+18) — the "exit transient"
// spike (~+40 eV/step) in ΔE_total. The direct erf/r potential has no charge in
// the cell to clip, so the transient should be GONE. The force gradient is
// insensitive to that offset, so the trajectory should be nearly identical.
//
// Ehrenfest, mass 1, charge -1, sigma_pot=0.35355, launch_z=-30, K0=2 (v=2).
// The light electron DECELERATES (expected); we do NOT gate on velocity drift.
// GATE: abort on NaN/complex energy; expect TRANSIT (proj_z crosses +12.5 far
// slab face with v>0) and E_electronic PLATEAU after exit.
//
// Full energy ledger: observables.csv + interactions.csv + projectile.csv, and
// the conservation column E_elec + KE_proj + U_proj_bg (must be flat).
// SAVE ~300 density frames (LJ_SAVE_EVERY) for the density-evolution GIF.
// Final checkpoint + resume (final-timestep-checkpoint.md).
//
// Env (defaults = this pilot): LJ_LX(35) LJ_LY(35) LJ_LZ(85) LJ_HALF(12.5)
//   LJ_N(100) LJ_EDGE_W(1.0) LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_SIGMA(0.5)
//   LJ_LAUNCH_Z(-30) LJ_K0(2.0) LJ_MASS(1.0) LJ_N_STEPS(1600) LJ_DT(0.04)
//   LJ_SAVE_EVERY(5) LJ_GS_DIR(REQUIRED) LJ_OUT(pilot_direct) LJ_RESUME(0).
// No inq/ or inq-study/ edit — wrapper-only.
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
#include <inqkit/dynamics/projectile_force.hpp>              // analytic + FD + DIRECT force
#include <inqkit/dynamics/moving_gaussian_projectile_potential.hpp>   // (a) DIRECT erf/r perturbation

#include <cmath>
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
    const double LX=env_d("LJ_LX",35),LY=env_d("LJ_LY",35),LZ=env_d("LJ_LZ",85),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",100); const double EDGE_W=env_d("LJ_EDGE_W",1.0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA_WP=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-30.0);
    const double K0=env_d("LJ_K0",2.0), MASS=env_d("LJ_MASS",1.0);
    const int N_STEPS=env_i("LJ_N_STEPS",1600); const double DT=env_d("LJ_DT",0.04);
    const int SAVE_EVERY=env_i("LJ_SAVE_EVERY",5);   // ~1600/5 = 320 frames
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","pilot_direct");
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

    std::cout<<std::setprecision(12)<<"\n=== pilot_direct (DIRECT erf/r potential Ehrenfest, DIRECT force) OUT="<<OUT
             <<" z0="<<LAUNCH_Z<<" k0="<<K0<<" mass="<<MASS
             <<" START="<<START<<" -> N_STEPS="<<N_STEPS<<(RESUME?"  [RESUME]":"")<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);                       // NO projectile ion (perturbation)
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(24).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(RESUME?CKPT:GS_DIR);                     // resume from RT ckpt, else GS

    // background well + live projectile (Ehrenfest) + moving DIRECT-potential perturbation
    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    inqkit::dynamics::Projectile proj(MASS, -1.0,
        inqkit::detail::Vec3{0.0,0.0,R0z},
        inqkit::detail::Vec3{0.0,0.0,V0z});                 // resume state or v0=k0/m toward slab
    // (a) DIRECT free-space erf/r perturbation (no charge/Poisson/background).
    inqkit::dynamics::moving_gaussian_projectile_potential proj_pert(proj, SIGMA_POT);
    auto pert = perturbations::sum(bg_pert, proj_pert);

    // constant background fields (n_+ for the direct force & U_proj_bg)
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);       // n_+ (constant field)
    auto phiplus = solvers::poisson::solve(nplus);          // phi_+  (for U_proj_bg & pairwise)

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables"+SEG+".csv",sel); obs.write_header();

    std::ofstream pj;   // energies in HARTREE (native), matching observables.csv
    if(electrons.root()){ pj.open(OUT+"/raw/observables/projectile"+SEG+".csv");
        pj<<std::setprecision(12)<<"step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal,force_z\n"; }

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
        auto Rn = proj.R();
        inq::vector3<double> center{Rn.x,Rn.y,Rn.z};
        // Ehrenfest: DIRECT erf/r HF force from (electrons - background). The direct
        // force uses V_proj = gaussian_potential (erf/r), consistent with the DIRECT
        // perturbation; it is linear in the density field, so the drag = F(n_e) - F(n_+),
        // both evaluated at the SAME projectile center.  (b) direct vs analytic.
        const double Fz = inqkit::dynamics::projectile_force_direct_z(electrons.density(), cell, center, SIGMA_POT)
                        - inqkit::dynamics::projectile_force_direct_z(nplus,             cell, center, SIGMA_POT);
        // U_proj_bg(t) at the current center (ideal, r_cut-free)
        auto nproj = inqkit::jellium::gaussian_density(basis, center, SIGMA_POT);
        const double Uprojbg = -operations::integral_product(nproj, phiplus);
        // pairwise decomposition at the CURRENT step (n_slab = slab electrons; n_P = n_proj)
        auto ct = inqkit::jellium::compute_coulomb(electrons.density(), nproj, phiplus);
        proj.advance(inqkit::detail::Vec3{0.0,0.0,Fz}, DT);                 // V->V_n, R->R_{n+1}
        if(electrons.root()){
            pj<<ctx.step<<","<<ctx.time_au<<","<<Rn.z<<","<<proj.V().z<<","
              <<proj.ke()<<","<<Uprojbg<<","<<Fz<<"\n";   // Hartree
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
         <<"run = classical_highdensity_sv/pilot_direct (DIRECT erf/r potential Ehrenfest, direct force)\nengine = inq-study\n"
         <<"representation = direct_potential\nprojectile = moving Gaussian direct erf/r potential (velocity-Verlet)\n"
         <<"force = projectile_force_direct (erf/r density-gradient HF, no charge/Poisson/background)\n"
         <<"drive = ehrenfest\n"
         <<"periodicity = "<<PER<<"  Lx = "<<LX<<"  Ly = "<<LY<<"  Lz = "<<LZ<<"  N = "<<N<<"  spacing = "<<SPACING<<"\n"
         <<"slab_half = "<<HALF<<"  edge_w = "<<EDGE_W<<"  n0 = "<<N0<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  sigma_wp = "<<SIGMA_WP<<"  sigma_pot = "<<SIGMA_POT<<"\n"
         <<"k0 = "<<K0<<"  mass = "<<MASS<<"  v0 = "<<(K0/MASS)<<"\n"
         <<"start_step = "<<START<<"  n_steps = "<<N_STEPS<<"  dt = "<<DT
         <<"  proj_z_final = "<<proj.R().z<<"  proj_vz_final = "<<proj.V().z<<"\n"
         <<"gs_dir = "<<GS_DIR<<"\nrun_completed = true\n";}
    std::cout<<"  done  proj_z_final="<<proj.R().z<<" proj_vz_final="<<proj.V().z<<"\n"; return 0;
}
