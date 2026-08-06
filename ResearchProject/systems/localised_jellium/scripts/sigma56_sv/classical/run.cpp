// ============================================================================
// localised_jellium / scripts/sigma56_sv / classical / run.cpp
//
// CLASSICAL half of the sigma_WP = 5 and 6 Bohr twin campaign.
// Plan: docs/plans/sigma56-sv-twin.md
//
// Derived from scripts/classical_highdensity_sv/dyn_direct/run.cpp (the DIRECT
// erf/r representation, user-selected 2026-08-02) with four changes:
//   (1) TWO ABSORBING BANDS at the z faces, so this half carries the SAME
//       absorber as its wavepacket twin and E_absorbed/L_slab is the same
//       estimator on both. LJ_CAP_ETA=0 turns them into an exact no-op, which is
//       how the CAP-free control at v = 3.0 is produced (no second binary).
//   (2) Cfg-backed defaults for the L_z = 105 box (launch -27.5, dx 0.40).
//   (3) Density frames written to the CANONICAL raw/vti/density_total/, with
//       frames/total left as a symlink. The ancestor wrote only frames/total/,
//       and the notebook builder looks in raw/vti/ — that mismatch silently
//       produced 8 classical notebooks with NO density GIFs in the slab_ks_wrap
//       campaign. Fixed at the source here.
//   (4) Interior + retained numbered checkpoints (the ancestor had only a final
//       one), per .claude/rules/checkpoint-dont-block.md — these runs are hours
//       long and a kill must cost at most one interval.
//
// ---------------------------------------------------------------------------
// WHY THE DIRECT erf/r POTENTIAL AND NOT THE POISSON PERTURBATION
// The projectile is added as V_proj = erf(|r-R|/(sqrt2 sigma_pot))/|r-R| directly
// (moving_gaussian_projectile_potential), NOT as poisson(gaussian_density).
// Consequences, all wanted here:
//   * there is NO projectile charge inside the cell to clip, so the projectile
//     may sit ANYWHERE — including outside the box — and the in-box potential is
//     simply the tail of its free-space field;
//   * no periodic G=0 neutralising sheet, hence none of the abrupt ledger lurch
//     (~+40 eV/step) the Poisson form shows as the clipped charge crosses the far
//     face — the "exit transient" that motivated this variant.
// The slab's own terms stay Poisson-exact; only the projectile is direct.
//
// ---------------------------------------------------------------------------
// THE WIDTH MAPPING (.claude/rules/sigma-wp-convention.md)
// LJ_SIGMA is sigma_WP — the LABEL, shared with the wavepacket twin. The Gaussian
// actually used here is sigma_pot = sigma_WP/sqrt(2), derived below, because that
// is the width at which the classical charge cloud equals the t=0 DENSITY of a
// sigma_WP wavepacket. sigma_WP = 5 -> sigma_pot = 3.53553; 6 -> 4.24264.
// NEVER pass sigma_pot in as LJ_SIGMA: it would be wrong by sqrt(2) and would
// silently mis-scale every pairwise energy.
//
// ---------------------------------------------------------------------------
// WHY THIS HALF NOW CARRIES A CAP (user decision 2026-08-02)
// The original classical benchmark was deliberately CAP-free: periodicity(2)
// makes z open for the ELECTROSTATICS, the moving charge leaves the box, energy
// is conserved and E_absorbed = plateau - E_GS is exact. Its wavepacket twin
// CANNOT be CAP-free (a KS orbital wraps on the FFT grid), so comparing the two
// compared a conserved quantity against an absorbed one — the documented reason
// the WP deposit curve sits 3-5x below the classical one. Putting the same
// absorber on both halves makes them the same measurement. The cost is measured,
// not assumed: one CAP-free control per sigma at v = 3.0 (LJ_CAP_ETA=0) gives the
// difference directly. Slab electrons are bound and the bands start at |z| = 40,
// so only genuinely emitted electrons should reach them.
//
// ---------------------------------------------------------------------------
// perturbations::absorbing takes FRACTIONAL cell coordinates, NOT Bohr (it
// compares point_op.rvector()[2], which lies in [-0.5,0.5)). For 12.5 Bohr per
// face on L_z = 105: width 0.119047619048, mid 0.440476190476 (= 46.25 Bohr),
// bands z in [+40,+52.5] and [-52.5,-40]. eta < 0 ABSORBS.
// ENGINE: inq-study — the CAP needs the complexified scalar potential; it does
// not compile against stock inq.
//
// ---------------------------------------------------------------------------
// LIGHT PROJECTILE (.claude/rules/light-projectile-stopping.md). This is a mass-1
// electron under free Ehrenfest dynamics: it DECELERATES strongly and that is
// physics, not a fault. Nothing here gates on velocity drift.
//
// Full energy ledger: observables.csv + interactions.csv (P/S/B decomposition,
// .claude/rules/decomposed-interaction-energies.md) + projectile.csv, every step.
//
// Env: LJ_LX(35) LJ_LY(35) LJ_LZ(105) LJ_HALF(12.5) LJ_N(100) LJ_EDGE_W(1.0)
//      LJ_PERIODICITY(2) LJ_SPACING(0.40) LJ_SIGMA(6.0) LJ_LAUNCH_Z(-27.5)
//      LJ_K0(2.0) LJ_MASS(1.0) LJ_N_STEPS(4360) LJ_DT(0.04) LJ_SAVE_EVERY(14)
//      LJ_CKPT_EVERY(0=auto N/5) LJ_CAP_ETA(-1.0) LJ_CAP_L(12.5)
//      LJ_GS_DIR(REQUIRED) LJ_OUT(REQUIRED) LJ_RESUME(0)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/projectile_background_energy.hpp>
#include <inqkit/jellium/gaussian_potential.hpp>             // direct erf/r V_proj
#include <inqkit/jellium/interaction_energies.hpp>           // compute_coulomb_direct
#include <inqkit/jellium/analytics.hpp>
#include <inqkit/dynamics/projectile.hpp>
#include <inqkit/dynamics/projectile_force.hpp>
#include <inqkit/dynamics/moving_gaussian_projectile_potential.hpp>

#include "../../../shared/configs/slab_n100_L35x35x105.hpp"

#include <algorithm>
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
using Cfg = localised_jellium::config::SlabN100_L35x35x105;

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
    // ---- parameters ---------------------------------------------------------
    const double LX=env_d("LJ_LX",Cfg::LX_BOHR), LY=env_d("LJ_LY",Cfg::LY_BOHR);
    const double LZ=env_d("LJ_LZ",Cfg::LZ_BOHR), HALF=env_d("LJ_HALF",Cfg::SLAB_HALF_WIDTH);
    const int    N=env_i("LJ_N",Cfg::N_ELECTRONS);
    const double EDGE_W=env_d("LJ_EDGE_W",Cfg::EDGE_WIDTH_BOHR);
    const int    PER=env_i("LJ_PERIODICITY",2);
    const double SPACING=env_d("LJ_SPACING",Cfg::SPACING_BOHR);
    // sigma_WP — the LABEL shared with the WP twin. sigma_pot is derived, never passed.
    const double SIGMA_WP=env_d("LJ_SIGMA",Cfg::WP_SIGMA_BOHR);
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",Cfg::LAUNCH_Z_BOHR);
    const double K0=env_d("LJ_K0",2.0), MASS=env_d("LJ_MASS",Cfg::PROJ_MASS);
    const int    N_STEPS=env_i("LJ_N_STEPS",4360);
    const double DT=env_d("LJ_DT",0.04);
    const int    SAVE_EVERY=env_i("LJ_SAVE_EVERY",14);
    const double CAP_ETA=env_d("LJ_CAP_ETA",Cfg::CAP_ETA_HA);   // 0 => CAP-free control
    const double CAP_L  =env_d("LJ_CAP_L",  Cfg::CAP_L_BOHR);
    const std::string GS_DIR=env_s("LJ_GS_DIR","");
    const std::string OUT="results/"+env_s("LJ_OUT","classical");
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    int CKPT_EVERY=env_i("LJ_CKPT_EVERY",0);
    if(CKPT_EVERY<=0) CKPT_EVERY=std::max(1,N_STEPS/5);          // >= 4 retained

    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){
        std::cerr<<"FATAL: LJ_GS_DIR missing or unset: '"<<GS_DIR<<"'\n"; return 2; }

    // CAP geometry — FRACTIONAL, derived from the ACTUAL L_z in use so it follows
    // the box automatically (do not hard-code 0.119...; LJ_LZ is an env knob).
    const double CAP_WIDTH_FRAC = CAP_L/LZ;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC/2.0;
    const double z_cap_in       = LZ/2.0 - CAP_L;
    const bool   CAP_ON         = (CAP_ETA != 0.0);
    const double rs             = inqkit::jellium::rs_from_n0(N0);

    // ---- checkpoint / resume (rule final-timestep-checkpoint.md) ------------
    const bool RESUME=env_i("LJ_RESUME",0)!=0;
    const std::string CKPT=OUT+"/checkpoint", RT_STATE=OUT+"/rt_state.txt";
    int START=0; double R0z=LAUNCH_Z, V0z=K0/MASS;
    if(RESUME){
        START=(int)read_state_d(RT_STATE,"last_step",-1);
        if(START<0){std::cerr<<"FATAL: LJ_RESUME=1 but no readable "<<RT_STATE<<"\n";return 2;}
        if(START>=N_STEPS){std::cout<<"Already at/after target ("<<START<<">="<<N_STEPS<<"); nothing to do.\n";return 0;}
        // The projectile is NOT part of `electrons`, so its R and V must be
        // restored by hand or the resumed run silently restarts it at the launch.
        R0z=read_state_d(RT_STATE,"proj_z",LAUNCH_Z); V0z=read_state_d(RT_STATE,"proj_vz",K0/MASS);
    }
    const std::string SEG=(START>0)?(".from"+std::to_string(START)):std::string("");

    std::cout<<std::setprecision(12)
             <<"\n=== sigma56_sv CLASSICAL (direct erf/r, Ehrenfest)  OUT="<<OUT<<" ===\n"
             <<"  cell      = "<<LX<<" x "<<LY<<" x "<<LZ<<" Bohr, periodicity("<<PER<<"), dx="<<SPACING<<"\n"
             <<"  slab      = "<<2*HALF<<" Bohr (half "<<HALF<<", edge "<<EDGE_W<<"), N="<<N
             <<", n0="<<N0<<", r_s="<<rs<<"\n"
             <<"  projectile= DIRECT erf/r, sigma_WP="<<SIGMA_WP<<" -> sigma_pot="<<SIGMA_POT
             <<"  mass="<<MASS<<"  charge="<<Cfg::PROJ_CHARGE<<"\n"
             <<"  launch    = z="<<LAUNCH_Z<<"  v0=k0="<<K0
             <<"  (standoff "<<(-HALF-LAUNCH_Z)<<" Bohr from the face at "<<-HALF<<")\n"
             <<"  CAP       = "<<(CAP_ON?"ON":"OFF (control)")<<"  eta="<<CAP_ETA
             <<" Ha  width="<<CAP_L<<" Bohr/face  bands +/-["<<z_cap_in<<","<<LZ/2.0<<"]\n"
             <<"  CAP frac  : mid="<<CAP_MID_FRAC<<"  width="<<CAP_WIDTH_FRAC<<"\n"
             <<"  steps     = "<<START<<" -> "<<N_STEPS<<"  dt="<<DT
             <<"  (t_final="<<N_STEPS*DT<<" a.u.)"<<(RESUME?"  [RESUME]":"")<<"\n"
             <<"  cadence   = density/"<<SAVE_EVERY<<"  ckpt/"<<CKPT_EVERY<<"  stats/1\n\n";

    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);                       // NO projectile ion (perturbation)
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(Cfg::EXTRA_STATES).temperature(Cfg::TEMPERATURE_EV*1.0_eV),
        input::kpoints::gamma());
    electrons.load(RESUME?CKPT:GS_DIR);                  // resume from RT ckpt, else GS

    // background well + live projectile (Ehrenfest) + moving DIRECT-potential
    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    inqkit::dynamics::Projectile proj(MASS, Cfg::PROJ_CHARGE,
        inqkit::detail::Vec3{0.0,0.0,R0z},
        inqkit::detail::Vec3{0.0,0.0,V0z});
    inqkit::dynamics::moving_gaussian_projectile_potential proj_pert(proj, SIGMA_POT);

    // The absorbing bands are summed in UNCONDITIONALLY: eta = 0 makes them an
    // exact identity, so the CAP-free control needs no separate code path (and
    // therefore cannot drift away from the production binary).
    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
    auto pert = perturbations::sum(perturbations::sum(bg_pert, proj_pert),
                                   perturbations::sum(cap_lo, cap_hi));

    // constant background fields (n_+ for the direct force & the pairwise ledger)
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);
    const double E_BB = inqkit::jellium::background_self_energy(nplus, phiplus);

    const std::string OBS = OUT+"/raw/observables";
    const std::string VTI = OUT+"/raw/vti";
    std::filesystem::create_directories(OBS);

    inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
    sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
    sel.energy_external=sel.energy_nonlocal=sel.energy_ion=true;
    inqkit::io::ObservablesWriter obs(OBS+"/observables"+SEG+".csv",sel); obs.write_header();

    std::ofstream pj;   // energies in HARTREE (native), matching observables.csv
    if(electrons.root()){ pj.open(OBS+"/projectile"+SEG+".csv");
        pj<<std::setprecision(12)
          <<"step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal,force_z\n"; }

    std::ofstream ix;
    if(electrons.root()){ ix.open(OBS+"/interactions"+SEG+".csv");
        ix<<std::setprecision(12)
          <<"step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,norm_slab,norm_proj\n"; }

    // ---- density frames -----------------------------------------------------
    // CANONICAL path raw/vti/density_total (what make_density_gif_battery and the
    // notebook builders read). frames/total is kept as a SYMLINK for any consumer
    // still hard-coding the ancestor's layout; the guard refuses to clobber a real
    // directory, so re-running over old data is safe.
    if(SAVE_EVERY>0){
        std::filesystem::create_directories(VTI+"/density_total");
        const std::filesystem::path legacy = OUT+"/frames/total";
        if(!std::filesystem::exists(legacy)){
            std::filesystem::create_directories(OUT+"/frames");
            std::error_code ec;
            std::filesystem::create_directory_symlink("../raw/vti/density_total", legacy, ec);
            if(ec) std::cout<<"  [warn] could not create frames/total symlink: "<<ec.message()<<"\n";
        }
    }
    inqkit::io::RealField3DLayout vti_layout{
        .field_name="density", .include_meta=false, .emit_raw=false,
        .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter total_wr(VTI+"/density_total", vti_layout,
                                           {.overwrite=(START==0)});

    auto write_rt_state=[&](int step){
        if(!electrons.root()) return;
        std::ofstream st(RT_STATE);
        st<<std::setprecision(12)
          <<"last_step="<<step<<"\ntime_au="<<(step*DT)<<"\ndt="<<DT
          <<"\nproj_z="<<proj.R().z<<"\nproj_vz="<<proj.V().z
          <<"\nproj_mass="<<MASS<<"\nproj_charge="<<Cfg::PROJ_CHARGE<<"\n";
    };

    // t = 0 frame, so the GIF starts at the true initial condition
    if(START==0 && SAVE_EVERY>0)
        total_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);

    inqkit::RealTimeSession rt(ions,electrons,1);
    rt.add([&](inqkit::StepContext const& ctx){
        obs.append(ctx);
        auto Rn = proj.R();
        inq::vector3<double> center{Rn.x,Rn.y,Rn.z};
        // Ehrenfest: DIRECT erf/r HF force from (electrons - background). Linear in
        // the density, so the drag is F(n_e) - F(n_+) at the SAME centre.
        const double Fz = inqkit::dynamics::projectile_force_direct_z(electrons.density(), cell, center, SIGMA_POT)
                        - inqkit::dynamics::projectile_force_direct_z(nplus,             cell, center, SIGMA_POT);
        // DIRECT-potential pairwise ledger: build V_proj ONCE and use it for both
        // the projectile pairwise terms and U_proj_bg. Slab terms stay Poisson-exact.
        auto vproj = inqkit::jellium::gaussian_potential(basis, center, SIGMA_POT);
        auto ct = inqkit::jellium::compute_coulomb_direct(electrons.density(), vproj, nplus, phiplus, SIGMA_POT);
        const double Uprojbg = ct.e_pb;
        proj.advance(inqkit::detail::Vec3{0.0,0.0,Fz}, DT);       // V->V_n, R->R_{n+1}
        if(electrons.root()){
            pj<<ctx.step<<","<<ctx.time_au<<","<<Rn.z<<","<<proj.V().z<<","
              <<proj.ke()<<","<<Uprojbg<<","<<Fz<<"\n";
            // norm_slab is the CAP witness on this half: with the absorber on it
            // decays as emitted electrons leave, with eta = 0 it must stay at N.
            ix<<ctx.step<<","<<ctx.time_au<<","<<ct.e_ss<<","<<ct.e_pp<<","<<ct.e_ps<<","
              <<ct.e_sb<<","<<ct.e_pb<<","<<E_BB<<","<<ct.norm_slab<<","<<ct.norm_p<<"\n";
        }
        if(SAVE_EVERY>0 && ctx.step % SAVE_EVERY == 0)
            total_wr.write(inqkit::fields::density::total(electrons), ctx.time_au, ctx.step);
    });

    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,
        [&](auto const& data){
            rt.step(data);
            // Interior checkpoints: rolling `checkpoint` (what LJ_RESUME loads)
            // plus a RETAINED numbered snapshot, so >= 4 resume points survive.
            if(data.iter()>0 && data.iter()%CKPT_EVERY==0){
                electrons.save(CKPT);
                write_rt_state(data.iter());
                electrons.save(OUT+"/ckpt_step"+tag6(data.iter()));
            }
        },
        options::theory{}.lda(),opts,pert,START);

    // FINAL checkpoint (rule final-timestep-checkpoint.md)
    electrons.save(CKPT);
    write_rt_state(N_STEPS);

    const double n_final = operations::integral(electrons.density());

    if(electrons.root()){
        pj.close(); ix.close();
        std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)
         <<"run = localised_jellium/sigma56_sv/classical/"<<env_s("LJ_OUT","classical")<<"\n"
         <<"plan = docs/plans/sigma56-sv-twin.md\n"
         <<"engine = inq-study\nxc = LDA\n"
         <<"representation = direct_potential\n"
         <<"projectile = moving Gaussian direct erf/r potential (velocity-Verlet Ehrenfest)\n"
         <<"force = projectile_force_direct (erf/r density-gradient HF)\n"
         <<"ledger = compute_coulomb_direct (P/S/B pairwise)\n"
         <<"periodicity = "<<PER<<"  Lx = "<<LX<<"  Ly = "<<LY<<"  Lz = "<<LZ
         <<"  N = "<<N<<"  spacing = "<<SPACING<<"\n"
         <<"slab_half = "<<HALF<<"  edge_w = "<<EDGE_W<<"  n0 = "<<N0<<"  r_s = "<<rs<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  sigma_wp = "<<SIGMA_WP<<"  sigma_pot = "<<SIGMA_POT<<"\n"
         <<"k0 = "<<K0<<"  mass = "<<MASS<<"  v0 = "<<(K0/MASS)<<"\n"
         <<"cap = "<<(CAP_ON?"on":"off")<<"  cap_eta_ha = "<<CAP_ETA
         <<"  cap_width_bohr = "<<CAP_L<<" per face\n"
         <<"cap_mid_frac = "<<CAP_MID_FRAC<<"  cap_width_frac = "<<CAP_WIDTH_FRAC
         <<"  cap_z_inner = "<<z_cap_in<<"\n"
         <<"start_step = "<<START<<"  n_steps = "<<N_STEPS<<"  dt = "<<DT
         <<"  t_final_au = "<<(N_STEPS*DT)<<"\n"
         <<"save_every = "<<SAVE_EVERY<<"  ckpt_every = "<<CKPT_EVERY<<"\n"
         <<"proj_z_final = "<<proj.R().z<<"  proj_vz_final = "<<proj.V().z<<"\n"
         <<"electron_count_final = "<<n_final<<"  (initial "<<N<<")\n"
         <<"gs_dir = "<<GS_DIR<<"\nrun_completed = true\n";
    }
    std::cout<<"  done  proj_z_final="<<proj.R().z<<"  proj_vz_final="<<proj.V().z
             <<"  integral_n_final="<<n_final<<" (initial "<<N<<")\n";
    return 0;
}
