// ============================================================================
// localised_jellium_dynamics / eval_projpot / run.cpp
//
// Single-point diagnostic (NO dynamics) to understand why the electrostatic
// residual reads ~7.4 eV instead of the clean 21.5 eV self-Hartree. Computes,
// in INQ's OWN p2 open-z G=0 convention, every projectile potential term:
//   ideal = ∫ n_proj·v_bg              (true Gaussian charge; r_cut-invariant)
//   impl  = −∫ n₊·v_ion                 (as-implemented pseudopotential)
//   gap   = ideal − impl               (pure pseudopotential representation error)
//   e_proj      = ∫ n_e·v_ion          (electron-projectile, as in E_ext)   [GS only]
//   e_proj_ideal= ∫ n_e·V_proj_ideal   (electron vs analytic Gaussian pot)  [GS only]
//   slabmbg_dv  = ∫ (n_e−n₊)·(v_ion−V_proj_ideal)  = residual − self-Hartree [GS only]
//
// Sweep LJ_SPACING (grid) and LJ_PROJ_UPF (r_cut) to test grid/cutoff sensitivity.
// If LJ_GS_DIR empty → basis-only (ideal/impl/gap; no electron-density terms).
// Env mirrors phase12; no inq/ or inq-study edit — wrapper only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>

#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
static double env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static int    env_i(const char* k,int d){const char*v=std::getenv(k);return v?std::atoi(v):d;}
static std::string env_s(const char* k,const std::string&d){const char*v=std::getenv(k);return v?std::string(v):d;}
static const char* UPF_DEFAULT=
  "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main(){
    const double HA=27.211386;
    const double LX=env_d("LJ_LX",50),LY=env_d("LJ_LY",50),LZ=env_d("LJ_LZ",120),HALF=env_d("LJ_HALF",12.5);
    const int N=env_i("LJ_N",82); const double EDGE_W=env_d("LJ_EDGE_W",0); const int PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",0.5), SIGMA_WP=env_d("LJ_SIGMA",0.5);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.5);
    const std::string PROJ_UPF=env_s("LJ_PROJ_UPF",UPF_DEFAULT);
    const std::string GS_DIR=env_s("LJ_GS_DIR","");
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);
    const double N0=double(N)/(LX*LY*(2.0*HALF));
    if(!std::filesystem::exists(PROJ_UPF)){std::cerr<<"FATAL: UPF missing: "<<PROJ_UPF<<"\n";return 2;}

    std::cout<<std::setprecision(10);
    std::cout<<"\n=== eval_projpot  spacing="<<SPACING<<"  upf="<<PROJ_UPF<<"  z="<<LAUNCH_Z<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);
    auto sp=ionic::species("H").pseudo_file(PROJ_UPF).mass(1.0/1822.8885);
    ions.insert(sp,{0.0*1.0_b,0.0*1.0_b,LAUNCH_Z*1.0_b});
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(20).temperature(0.00862*1.0_eV),input::kpoints::gamma());

    bool loaded=false;
    if(!GS_DIR.empty() && std::filesystem::exists(GS_DIR)){ electrons.load(GS_DIR); loaded=true; }

    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    auto basis   = electrons.density().basis();
    auto & comm  = electrons.states_comm();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);               // φ₊ ; v_bg = −φ₊
    auto nproj   = inqkit::jellium::gaussian_density(basis,{0.0,0.0,LAUNCH_Z},SIGMA_POT);
    auto vproj   = solvers::poisson::solve(nproj);               // V_proj_ideal = poisson(n_proj)
    auto & atomic_pot = electrons.atomic_pot();
    auto v_ion = operations::add(
        solvers::poisson::solve(atomic_pot.ionic_density(comm,basis,ions)),
        atomic_pot.local_potential(comm,basis,ions));

    const double ideal    = -operations::integral_product(nproj,phiplus);   // ∫n_proj·v_bg
    const double impl     = -operations::integral_product(nplus,v_ion);     // −∫n₊·v_ion
    const double bg_vion  =  operations::integral_product(nplus,v_ion);     // ∫n₊·v_ion (=−impl)
    const double bg_videal=  operations::integral_product(nplus,vproj);     // ∫n₊·V_proj_ideal
    const double nprojnorm=  operations::integral(nproj);

    std::cout<<"  n_proj_norm = "<<nprojnorm<<"\n";
    std::cout<<"  ideal            ∫n_proj·v_bg      = "<<ideal*HA<<" eV\n";
    std::cout<<"  impl            −∫n₊·v_ion         = "<<impl*HA<<" eV\n";
    std::cout<<"  gap = ideal−impl (pseudopot error) = "<<(ideal-impl)*HA<<" eV\n";
    std::cout<<"  ∫n₊·v_ion                          = "<<bg_vion*HA<<" eV\n";
    std::cout<<"  ∫n₊·V_proj_ideal                   = "<<bg_videal*HA<<" eV\n";

    if(loaded){
        auto n_e = electrons.density();
        const double e_proj       = operations::integral_product(n_e,v_ion);
        const double e_proj_ideal = operations::integral_product(n_e,vproj);
        // residual − self_Hartree = ∫(n_e−n₊)·(v_ion−V_proj_ideal)
        const double slabmbg_dv = (e_proj - e_proj_ideal) - (bg_vion - bg_videal);
        std::cout<<"  e_proj           ∫n_e·v_ion        = "<<e_proj*HA<<" eV\n";
        std::cout<<"  e_proj_ideal     ∫n_e·V_proj_ideal  = "<<e_proj_ideal*HA<<" eV\n";
        std::cout<<"  ∫(n_e−n₊)·(v_ion−V_proj_ideal)     = "<<slabmbg_dv*HA<<" eV   <== residual − self_Hartree\n";
    } else {
        std::cout<<"  (no GS loaded — electron-density terms skipped)\n";
    }
    std::cout<<"  DONE\n";
    return 0;
}
