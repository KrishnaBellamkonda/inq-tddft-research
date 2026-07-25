// ============================================================================
// localised_jellium_dynamics / proj_ghost / run.cpp
//
// Classical projectile as a GHOST-UPF PSEUDOPOTENTIAL (the "pseudopotential"
// representation), the twin of proj_perturbation's Gaussian-charge case. Inserts
// electron_gaussian_wpsigma0p5.upf (z_valence 0, local erf(r/0.5)/r, charge std
// sigma_WP/sqrt2 = 0.354) as a STATIONARY ghost ion at LJ_LAUNCH_Z, re-applies the
// background well, propagates a few short steps, and tabulates the FULL energy
// ledger (total/kin/hartree/xc/external/nonlocal/ion) so the twin-run-analysis
// engine can read the pseudopotential-representation residual (expected ~7.4 eV;
// the ~14 eV shortfall vs the 21.7 eV self-Hartree is the KNOWN ghost-UPF tail
// aliasing — reference_ghost_upf_tail_aliasing).
//
// U_proj_bg is the CLEAN IDEAL term -int n_proj*phi_+ (r_cut-free), NEVER the impl
// term. Config MATCHES the existing WP twin (p2 open-z, GS h2/gs_p2_lz120,
// launch_z -24.5, N 82, spacing 0.5) so the pair is valid for check_twin.py.
//
// Env: LJ_LX(50) LJ_LY(50) LJ_LZ(120) LJ_HALF(12.5) LJ_N(82) LJ_EDGE_W(0)
//      LJ_PERIODICITY(2) LJ_SPACING(0.5) LJ_SIGMA(0.5) LJ_LAUNCH_Z(-24.5)
//      LJ_N_STEPS(2) LJ_DT(0.01) LJ_GS_DIR(REQUIRED) LJ_OUT(proj_ghost).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density

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

static const char* PROJ_PSEUDO =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
    "shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main(){
    const double HA=27.211386;
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA_WP=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.5);
    const int N_STEPS=env_i("LJ_N_STEPS",2); const double DT=env_d("LJ_DT",0.01);
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","proj_ghost");
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);
    const double N0=double(N)/(LX*LY*(2.0*HALF));
    constexpr double M_PROJ=1.0/1822.8885;   // electron mass in amu (irrelevant at v=0)
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}

    std::cout<<std::setprecision(12)<<"\n=== proj_ghost OUT="<<OUT<<" z="<<LAUNCH_Z<<" spacing="<<SPACING<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto sp=ionic::species("H").pseudo_file(PROJ_PSEUDO).mass(M_PROJ);   // GHOST projectile ion
    ions.insert(sp, {0.0*1.0_b, 0.0*1.0_b, LAUNCH_Z*1.0_b});
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(GS_DIR);
    ions.velocities()[0]=vector3<double>{0.0,0.0,0.0};   // stationary

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);   // ghost ion supplies proj potential

    // CLEAN U_proj_bg = -int n_proj*phi_+   (ideal, r_cut-free; NEVER the impl term for a ghost)
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);
    auto nproj   = inqkit::jellium::gaussian_density(basis, {0.0,0.0,LAUNCH_Z}, SIGMA_POT);
    const double nproj_norm = operations::integral(nproj);
    const double U_proj_bg  = -operations::integral_product(nproj, phiplus);
    std::cout<<"  n_proj_norm = "<<nproj_norm<<"   U_proj_bg = "<<U_proj_bg*HA<<" eV (clean ideal)\n";

    std::filesystem::create_directories(OUT+"/raw/observables");
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables.csv",sel); obs.write_header();
    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const&ctx){obs.append(ctx);});
    auto step_fn=[&](auto const&data){rt.step(data);};
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,bg_pert);

    if(electrons.root()){std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)
         <<"run = proj_ghost (ghost-UPF pseudopotential classical projectile)\nengine = inq-study\n"
         <<"representation = pseudopotential\nprojectile = ghost UPF wpsigma0p5 (z_valence 0), stationary\n"
         <<"periodicity = "<<PER<<"  Lz = "<<LZ<<"  N = "<<N<<"  spacing = "<<SPACING<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  sigma_wp = "<<SIGMA_WP<<"  sigma_pot = "<<SIGMA_POT<<"\n"
         <<"n_proj_norm = "<<nproj_norm<<"\n"
         <<"U_proj_bg_ha = "<<U_proj_bg<<"   U_proj_bg_ev = "<<U_proj_bg*HA<<"   (clean ideal, r_cut-free)\n"
         <<"gs_dir = "<<GS_DIR<<"\nrun_completed = true\n";}
    std::cout<<"  done\n"; return 0;
}
