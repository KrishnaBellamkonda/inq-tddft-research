// ============================================================================
// systems/localised_jellium/scripts/ng_mass_ladder/classical/run.cpp
//
// CLASSICAL half of the Nazarov-Gross mass ladder — the M -> infinity anchor.
// Plan: docs/plans/nazarov-gross-slab-mass-ladder.md
//
// The projectile is a rigid Gaussian CHARGE advanced by Ehrenfest velocity-Verlet
// and applied as a moving PERTURBATION -- NOT a UPF/ion (user decision
// 2026-08-05). sigma_pot is therefore a runtime argument, so the sigma sweep of
// plan step 10 costs nothing extra on this side, and there is no r_cut aliasing.
//
// WHY IT IS THE ANCHOR. A rigid classical cloud has NO hbar: it cannot disperse,
// so its width is frozen at sigma_pot for the whole transit and its stopping
// CANNOT depend on mass. Every quantum rung's departure from this run is
// therefore unambiguously a quantum effect -- which is exactly the argument
// against the "it's just a Gaussian form factor" objection (plan section 7.7).
//
// TWO CLASSICAL RUNS ARE INTENDED, differing only in NG_MASS:
//   * NG_MASS = 1e6  -> effectively infinite: constant velocity, the NG M=inf rung
//   * NG_MASS = 1.0  -> matched-mass control: recoils like the WP does, so the
//                       WP-minus-this difference isolates quantum-ness at fixed M
// Mass here affects ONLY the Ehrenfest recoil, never the charge distribution.
//
// Conserved quantity (correctness gate): E_electronic + proj_ke + U_proj_bg flat.
//
// Env: NG_MASS(1e6) NG_V(Cfg::V0_AU) NG_SIGMA_WP(Cfg::WP_SIGMA_BOHR)
//      NG_SPACING(0.50) NG_GS_DIR(REQUIRED) NG_OUT(cl_inf) NG_N_STEPS
//      NG_DT(0.02) NG_LAUNCH_Z(-25) NG_SAVE_EVERY(50) NG_CKPT_EVERY(200)
//      NG_RESUME(0) NG_CAP(1) NG_DELTA(0.1)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>   // gaussian_density
#include <inqkit/jellium/interaction_energies.hpp>           // pairwise P/S/B
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/dynamics/projectile_force.hpp>
#include <inqkit/dynamics/moving_gaussian_projectile_perturbation.hpp>

#include "../../../shared/configs/slab_n206_L30x30x120_rs2p5.hpp"

#include <algorithm>
#include <chrono>
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
using Cfg = localised_jellium::config::SlabN206_L30x30x120_rs2p5;

static double      env_d(const char* k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static int         env_i(const char* k,int d){const char*v=std::getenv(k);return v?std::atoi(v):d;}
static std::string env_s(const char* k,const std::string&d){const char*v=std::getenv(k);return v?std::string(v):d;}
static double read_state_d(const std::string& path,const char* key,double def){
    std::ifstream f(path); std::string line; std::string k=std::string(key)+"=";
    while(std::getline(f,line)) if(line.rfind(k,0)==0) return std::atof(line.c_str()+k.size());
    return def;
}
static std::string tag6(int n){std::ostringstream o;o<<std::setw(6)<<std::setfill('0')<<n;return o.str();}

int main(){
    auto t0 = std::chrono::steady_clock::now();
    const double HA = 27.211386245988;

    const double MASS      = env_d("NG_MASS", 1.0e6);            // 1e6 => M -> infinity
    const double VEL       = env_d("NG_V", Cfg::V0_AU);
    const double SIGMA_WP  = env_d("NG_SIGMA_WP", Cfg::WP_SIGMA_BOHR);
    const double SIGMA_POT = SIGMA_WP / std::sqrt(2.0);           // rule sigma-wp-convention.md
    const double SPACING   = env_d("NG_SPACING", Cfg::SPACING_BOHR);
    const std::string TAG  = env_s("NG_OUT", "cl_inf");
    const std::string OUT  = "results/" + TAG;
    const double DT        = env_d("NG_DT", 0.08 * 1.0 * Cfg::SPACING_BOHR * Cfg::SPACING_BOHR);
    const int    N_STEPS   = env_i("NG_N_STEPS", 2560);
    // DISK CADENCES (user instruction 2026-08-05). Measured on this grid:
    // a density VTI is 6.9 MB and an RT checkpoint 1.7 GB (123 states x 864k
    // points x 16 B). Interior checkpoints ROLL INTO ONE DIRECTORY, so
    // N_STEPS/3 gives 2 interior writes + the mandatory FINAL one at 1.7 GB
    // total, not 3 x 1.7 GB. ~24 density frames is ample for the GIF builders,
    // which sample ~30.
    const int    SAVE_EVERY= env_i("NG_SAVE_EVERY", std::max(1, N_STEPS/24));
    const int    CKPT_EVERY= env_i("NG_CKPT_EVERY", std::max(1, N_STEPS/3));
    const double LAUNCH_Z  = env_d("NG_LAUNCH_Z", Cfg::WP_CZ_BOHR);
    const double DELTA     = env_d("NG_DELTA", 0.1);              // finite-difference force step
    const bool   RESUME    = env_i("NG_RESUME",0)!=0;
    const bool   USE_CAP   = env_i("NG_CAP",1)!=0;
    const double E_EV      = 0.5*MASS*VEL*VEL*HA;

    const std::string GS_DIR = env_s("NG_GS_DIR","");
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: NG_GS_DIR missing: '"<<GS_DIR<<"'\n";return 2;}

    // dt guard: the BATH orbitals have m = 1, so the ceiling is 0.08*1*h^2
    // regardless of the (classical, off-grid) projectile mass.
    const double dt_ceiling = 0.08 * 1.0 * SPACING * SPACING;
    if(DT > dt_ceiling*1.001){
        std::cerr<<"FATAL: dt = "<<DT<<" exceeds the bath ceiling "<<dt_ceiling
                 <<" = 0.08*h^2 (h="<<SPACING<<").\n"; return 2; }

    const std::string CKPT=OUT+"/checkpoint", RT_STATE=OUT+"/rt_state.txt";
    int START=0; double R0z=LAUNCH_Z, V0z=VEL;
    if(RESUME){
        START=(int)read_state_d(RT_STATE,"last_step",-1);
        if(START<0){std::cerr<<"FATAL: NG_RESUME=1 but no readable "<<RT_STATE<<"\n";return 2;}
        if(START>=N_STEPS){std::cout<<"Already at/after target ("<<START<<">="<<N_STEPS<<"); nothing to do.\n";return 0;}
        R0z=read_state_d(RT_STATE,"proj_z",LAUNCH_Z); V0z=read_state_d(RT_STATE,"proj_vz",VEL);
    }
    const std::string SEG=(START>0)?(".from"+std::to_string(START)):std::string("");

    std::cout<<std::setprecision(12)
             <<"\n=== ng_mass_ladder classical (out="<<OUT<<") ===\n"
             <<"  mass="<<MASS<<" v="<<VEL<<" ("<<VEL/Cfg::KF_AU<<" v_F, "
             <<VEL/Cfg::V_BRAGG_AU<<" of Bragg peak) E="<<E_EV<<" eV\n"
             <<"  sigma_WP="<<SIGMA_WP<<" -> sigma_pot="<<SIGMA_POT<<"  h="<<SPACING
             <<"  dt="<<DT<<" (ceiling "<<dt_ceiling<<")\n"
             <<"  START="<<START<<" -> N_STEPS="<<N_STEPS<<(RESUME?"  [RESUME]":"")
             <<"  launch_z="<<LAUNCH_Z<<"  cap="<<(USE_CAP?"on":"off")<<"\n";

    auto cell0=systems::cell::orthorhombic(Cfg::LX_BOHR*1.0_b,Cfg::LY_BOHR*1.0_b,Cfg::LZ_BOHR*1.0_b);
    auto cell=(Cfg::PERIODICITY==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);                       // NO projectile ion
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(Cfg::N_ELECTRONS).extra_states(Cfg::EXTRA_STATES)
        .temperature(Cfg::TEMPERATURE_EV*1.0_eV),input::kpoints::gamma());
    electrons.load(RESUME?CKPT:GS_DIR);

    // background well + live projectile (Ehrenfest) + moving perturbation + CAPs
    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=Cfg::N0;
    bg.half_width=Cfg::SLAB_HALF_WIDTH; bg.slab_axis=Cfg::SLAB_AXIS;
    bg.center={0.0,0.0,Cfg::SLAB_CENTER_BOHR}; bg.edge_width=Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    inqkit::dynamics::Projectile proj(MASS, -1.0,
        inqkit::detail::Vec3{0.0,0.0,R0z},
        inqkit::detail::Vec3{0.0,0.0,V0z});
    inqkit::dynamics::moving_gaussian_projectile_perturbation proj_pert(proj, SIGMA_POT);

    perturbations::absorbing cap_lo(Cfg::CAP_ETA_HA*1.0_Ha, -Cfg::CAP_MID_FRAC, Cfg::CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(Cfg::CAP_ETA_HA*1.0_Ha,  Cfg::CAP_MID_FRAC, Cfg::CAP_WIDTH_FRAC);

    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);
    const double E_BB = inqkit::jellium::background_self_energy(nplus, phiplus);

    std::filesystem::create_directories(OUT+"/raw/observables");
    std::filesystem::create_directories(OUT+"/raw/vti/density_total");
    // FULL energy decomposition, matching the WP half column-for-column so the
    // two ledgers can be differenced directly. energy_external is non-zero here
    // (the projectile IS an external potential) and identically zero in the WP
    // half — which is exactly why interactions.csv, not this file, carries the
    // comparable numbers (rule decomposed-interaction-energies.md).
    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=true;
    sel.energy_ion_kinetic=sel.energy_exact_exchange=true;
    sel.energy_nvxc=sel.energy_eigenvalues=true;
    sel.cod_x=sel.cod_y=sel.cod_z=true;
    inqkit::io::ObservablesWriter obs(OUT+"/raw/observables/observables"+SEG+".csv",sel); obs.write_header();

    std::ofstream pj, ix, nlog;
    if(electrons.root()){
        pj.open(OUT+"/raw/observables/projectile"+SEG+".csv");
        pj<<std::setprecision(12)<<"step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal\n";
        ix.open(OUT+"/raw/observables/interactions"+SEG+".csv");
        ix<<std::setprecision(12)<<"step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,norm_slab,norm_proj\n";
        nlog.open(OUT+"/raw/observables/electron_number"+SEG+".csv");
        nlog<<std::setprecision(12)<<"step,time_au,N_total\n";
    }

    inqkit::io::RealField3DLayout lay{.field_name="density",.include_meta=false,
        .emit_raw=false,.emit_vti=true,.vti_format=inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter total_wr(OUT+"/raw/vti/density_total", lay, {.overwrite = !RESUME});

    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const& ctx){
        obs.append(ctx);
        if(SAVE_EVERY>0 && ctx.step % SAVE_EVERY == 0)
            total_wr.write(inqkit::fields::density::total(*ctx.electrons), ctx.time_au, ctx.step);

        // Ehrenfest projectile: HF force from the (electrons - background) field.
        // projectile_force_z is linear in the field, so F(phi_e - phi+) =
        // F(phi_e) - F(phi+); phi_e is re-solved each step, phi+ is constant.
        auto Rn = proj.R();
        inq::vector3<double> center{Rn.x,Rn.y,Rn.z};
        auto phi_e = solvers::poisson::solve(ctx.electrons->density());
        const double Fz = inqkit::dynamics::projectile_force_z(phi_e,   center, SIGMA_POT, DELTA)
                        - inqkit::dynamics::projectile_force_z(phiplus, center, SIGMA_POT, DELTA);
        auto nproj = inqkit::jellium::gaussian_density(basis, center, SIGMA_POT);
        const double Uprojbg = -operations::integral_product(nproj, phiplus);
        auto ct = inqkit::jellium::compute_coulomb(ctx.electrons->density(), nproj, phiplus);
        proj.advance(inqkit::detail::Vec3{0.0,0.0,Fz}, DT);       // V->V_n, R->R_{n+1}
        if(ctx.electrons->root()){
            pj<<ctx.step<<","<<ctx.time_au<<","<<Rn.z<<","<<proj.V().z<<","
              <<proj.ke()<<","<<Uprojbg<<"\n";
            ix<<ctx.step<<","<<ctx.time_au<<","<<ct.e_ss<<","<<ct.e_pp<<","<<ct.e_ps<<","
              <<ct.e_sb<<","<<ct.e_pb<<","<<E_BB<<","<<ct.norm_slab<<","<<ct.norm_p<<"\n";
            nlog<<ctx.step<<","<<ctx.time_au<<","<<ct.norm_slab<<"\n";
        }
    });

    auto step_fn=[&](auto const& data){
        rt.step(data);
        const int it = data.iter();
        if(it > START && it % CKPT_EVERY == 0 && it < N_STEPS){
            electrons.save(CKPT);
            if(data.root()){
                std::ofstream st(RT_STATE, std::ios::trunc);
                st<<std::setprecision(12)
                  <<"last_step="<<it<<"\ntime_au="<<(it*DT)<<"\ndt="<<DT
                  <<"\nproj_z="<<proj.R().z<<"\nproj_vz="<<proj.V().z
                  <<"\nproj_mass="<<MASS<<"\nproj_charge="<<-1.0<<"\n";
                std::cout<<"  [ckpt] step "<<it<<" (t="<<it*DT<<")\n";
            }
        }
    };

    auto pert_cap = perturbations::sum(bg_pert, perturbations::sum(proj_pert,
                        perturbations::sum(cap_lo, cap_hi)));
    auto pert_raw = perturbations::sum(bg_pert, proj_pert);
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    if(USE_CAP) real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,pert_cap,START);
    else        real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,pert_raw,START);

    // FINAL checkpoint (rule final-timestep-checkpoint.md)
    electrons.save(CKPT);
    if(electrons.root()){
        std::ofstream st(RT_STATE, std::ios::trunc);
        st<<std::setprecision(12)
          <<"last_step="<<N_STEPS<<"\ntime_au="<<(N_STEPS*DT)<<"\ndt="<<DT
          <<"\nproj_z="<<proj.R().z<<"\nproj_vz="<<proj.V().z
          <<"\nproj_mass="<<MASS<<"\nproj_charge="<<-1.0<<"\n";
    }

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if(electrons.root()){
        pj.close(); ix.close(); nlog.close();
        std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)
         <<"run = localised_jellium/ng_mass_ladder/"<<TAG<<"\n"
         <<"plan = docs/plans/nazarov-gross-slab-mass-ladder.md\n"
         <<"engine = inq-study\nxc = LDA (no SIC)\n"
         <<"representation = gaussian_charge_perturbation\n"
         <<"projectile = moving Gaussian charge (velocity-Verlet Ehrenfest), NOT a UPF\n"
         <<"sigma_wp = "<<SIGMA_WP<<"  sigma_pot = "<<SIGMA_POT<<"\n"
         <<"mass = "<<MASS<<"  velocity = "<<VEL<<"  E = "<<E_EV<<" eV\n"
         <<"v_over_vF = "<<VEL/Cfg::KF_AU<<"  v_over_vBragg = "<<VEL/Cfg::V_BRAGG_AU<<"\n"
         <<"cap = "<<(USE_CAP?"on":"off")<<"  cap_eta_ha = "<<Cfg::CAP_ETA_HA
         <<"  cap_region_bohr = +/-["<<Cfg::CAP_INNER_FACE_BOHR<<","<<Cfg::LZ_BOHR/2.0<<"]\n"
         <<"cell_bohr = "<<Cfg::LX_BOHR<<"x"<<Cfg::LY_BOHR<<"x"<<Cfg::LZ_BOHR
         <<"  periodicity = "<<Cfg::PERIODICITY<<"  spacing = "<<SPACING<<"\n"
         <<"background = slab half_width "<<Cfg::SLAB_HALF_WIDTH<<" edge "<<Cfg::EDGE_WIDTH_BOHR<<"\n"
         <<"n0_a0m3 = "<<Cfg::N0<<"  r_s = "<<Cfg::RS_BOHR<<"  kf_au = "<<Cfg::KF_AU<<"\n"
         <<"n_electrons = "<<Cfg::N_ELECTRONS<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  delta_fd = "<<DELTA<<"\n"
         <<"dt = "<<DT<<"  dt_ceiling = "<<dt_ceiling<<"\n"
         <<"start_step = "<<START<<"  n_steps = "<<N_STEPS
         <<"  segment_suffix = "<<(SEG.empty()?"(none)":SEG)<<"\n"
         <<"proj_z_final = "<<proj.R().z<<"  proj_vz_final = "<<proj.V().z<<"\n"
         <<"e_bb_ha = "<<E_BB<<"\n"
         <<"gs_dir = "<<GS_DIR<<"\n"
         <<"wall_time_s = "<<wall<<" (this segment)\nrun_completed = true\n";
    }
    std::cout<<"  done  proj_z_final="<<proj.R().z<<" proj_vz_final="<<proj.V().z
             <<"  wall="<<wall<<"s\n";
    return 0;
}
