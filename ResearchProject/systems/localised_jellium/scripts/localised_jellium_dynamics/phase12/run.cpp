// ============================================================================
// localised_jellium_dynamics / phase12 / run.cpp
//
// Classical insertion run carrying the E_proj_bg diagnostic columns (campaign
// localised-jellium-dynamics-analysis, Phase 1 ledger + Phase 2 r_cut sweep).
// Loads the matching GS, inserts a stationary classical Gaussian ghost at
// z=LJ_LAUNCH_Z with UPF LJ_PROJ_UPF (truncated variants for Phase 2), computes
// projectile_background_energy(...) and streams energy_proj_bg_{ideal,impl} as
// two observables.csv columns (DIAGNOSTIC — not in energy_total). Static ghost →
// the term is a per-run constant (set once via obs.set_proj_bg). 2-step t=0
// insertion, identical geometry to the A1 ledger so dE_WP is reused.
//
// Env: LJ_LX(50) LJ_LY(50) LJ_LZ(120) LJ_HALF(12.5) LJ_N(82) LJ_EDGE_W(0)
//      LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_SIGMA(0.5) LJ_LAUNCH_Z(REQUIRED)
//      LJ_PROJ_UPF(default full) LJ_GS_DIR(REQUIRED) LJ_OUT(REQUIRED)
//      LJ_N_STEPS(2) LJ_DT(0.01). Build vs INQ_SOURCE=inq-study.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>

#include <chrono>
#include <cmath>
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
static const char* PROJ_PSEUDO_DEFAULT=
  "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main(){
    auto t0=std::chrono::steady_clock::now();
    constexpr double M_PROJ=1.0/1822.8885;
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA_WP=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.5);   // default r=12
    const int N_STEPS=env_i("LJ_N_STEPS",2); const double DT=env_d("LJ_DT",0.01);
    const std::string PROJ_UPF=env_s("LJ_PROJ_UPF",PROJ_PSEUDO_DEFAULT);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","cl");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}
    if(!std::filesystem::exists(PROJ_UPF)){std::cerr<<"FATAL: UPF missing: "<<PROJ_UPF<<"\n";return 2;}
    const double N0=double(N)/(LX*LY*(2.0*HALF));
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);

    std::cout<<"\n=== phase12 classical OUT="<<OUT<<" per="<<PER<<" z="<<LAUNCH_Z<<" upf="<<PROJ_UPF<<" ===\n";
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

    // E_proj_bg (per-run constant, static ghost). ideal = ∫ n_proj·v_bg (r_cut-invariant),
    // impl = −∫ n₊·v_ion (r_cut-dependent). DIAGNOSTIC — not added to any energy.
    auto epb=inqkit::jellium::projectile_background_energy(
        bg_pert, electrons, ions, {0.0,0.0,LAUNCH_Z}, SIGMA_POT);
    const double HA=27.211386;
    std::cout<<"  E_proj_bg ideal="<<epb.ideal*HA<<" eV  impl="<<epb.impl*HA
             <<" eV  n_proj_norm="<<epb.n_proj_norm<<"\n";

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true; sel.density_l2=false;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=sel.energy_ion_kinetic=true;
    sel.energy_exact_exchange=sel.energy_nvxc=sel.energy_eigenvalues=true;
    sel.energy_proj_bg_ideal=sel.energy_proj_bg_impl=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables.csv",sel);
    obs.set_proj_bg(epb.ideal, epb.impl);
    obs.write_header();
    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){obs.append(ctx);});
    auto step_fn=[&](auto const&data){rt.step(data);};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert);

    double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = localised_jellium_dynamics/phase12/"<<env_s("LJ_OUT","cl")<<"\nengine = inq-study\n"
         <<"projectile = ghost sigma_wp="<<SIGMA_WP<<" (sigma_pot="<<SIGMA_POT<<"), stationary\n"
         <<"proj_upf = "<<PROJ_UPF<<"\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  launch_z = "<<LAUNCH_Z<<"  r_from_face = "<<(-LAUNCH_Z-HALF)<<"\n"
         <<"e_proj_bg_ideal_ha = "<<epb.ideal<<"  e_proj_bg_impl_ha = "<<epb.impl<<"\n"
         <<"e_proj_bg_ideal_ev = "<<epb.ideal*HA<<"  e_proj_bg_impl_ev = "<<epb.impl*HA<<"\n"
         <<"n_proj_norm = "<<epb.n_proj_norm<<"   (guard: must be ~1.0)\n"
         <<"gs_dir = "<<GS_DIR<<"  dt = "<<DT<<"  n_steps = "<<N_STEPS<<"\nwall_time_s = "<<wall<<"\nrun_completed = true\n";}
    std::cout<<"  done wall="<<wall<<"s\n"; return 0;
}
