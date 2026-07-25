// ============================================================================
// localised_jellium / scripts/ke_check / classical / run.cpp
//
// KE-bookkeeping experiment (campaign localised_jellium_parameter_study_2):
// a MOVING classical Gaussian ghost, to test whether INQ records the
// projectile's kinetic energy as an ionic energy. Copy of campaign_autorun/
// classical/run.cpp with two added knobs:
//   LJ_VZ   projectile velocity along z (Bohr/atomictime), default 0
//   LJ_EHRENFEST  1 = Ehrenfest ion dynamics (ion actually moves), default 0
// Everything else (geometry, background, GS, observables) is identical, so the
// energy columns are directly comparable to the campaign insertion runs.
// Build against INQ_SOURCE=inq-study.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
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
static const char* PROJ_PSEUDO=
  "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main(){
    constexpr double M_PROJ=1.0/1822.8885;   // amu -> exactly 1 electron mass in atomic units
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), LAUNCH_Z=env_d("LJ_LAUNCH_Z",-16.5);
    const double VZ=env_d("LJ_VZ",0.0); const int EHRENFEST=env_i("LJ_EHRENFEST",0);
    const int N_STEPS=env_i("LJ_N_STEPS",4); const double DT=env_d("LJ_DT",0.01);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","cl");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    std::cout<<"\n=== ke_check classical OUT="<<OUT<<" per="<<PER<<" vz="<<VZ<<" ehrenfest="<<EHRENFEST<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto sp=ionic::species("H").pseudo_file(env_s("LJ_PROJ_UPF",PROJ_PSEUDO)).mass(M_PROJ);
    ions.insert(sp,{0.0*1.0_b,0.0*1.0_b,LAUNCH_Z*1.0_b});
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(GS_DIR);
    ions.velocities()[0]=vector3<double>{0.0,0.0,VZ};

    // analytic ionic KE = 1/2 m v^2 with m = species.mass() (atomic units = 1 m_e here)
    const double ke_ion_analytic=0.5*sp.mass()*VZ*VZ;
    std::cout<<"  species.mass()="<<sp.mass()<<" (atomic units)  KE_ion=1/2 m v^2 = "
             <<ke_ion_analytic<<" Ha = "<<ke_ion_analytic*27.211386<<" eV\n";

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true; sel.density_l2=false;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=sel.energy_ion_kinetic=true;
    sel.energy_exact_exchange=sel.energy_nvxc=sel.energy_eigenvalues=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables.csv",sel); obs.write_header();
    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){obs.append(ctx);});
    auto step_fn=[&](auto const&data){rt.step(data);};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    if(EHRENFEST) opts=opts.ehrenfest();
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert);

    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = ke_check/classical/"<<env_s("LJ_OUT","cl")<<"\nengine = inq-study\n"
         <<"vz = "<<VZ<<"  ehrenfest = "<<EHRENFEST<<"  species_mass_au = "<<sp.mass()<<"\n"
         <<"ke_ion_analytic_ha = "<<ke_ion_analytic<<"  ke_ion_analytic_ev = "<<ke_ion_analytic*27.211386<<"\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  launch_z = "<<LAUNCH_Z<<"\n"
         <<"gs_dir = "<<GS_DIR<<"  dt = "<<DT<<"  n_steps = "<<N_STEPS<<"\nrun_completed = true\n";}
    std::cout<<"  done\n"; return 0;
}
