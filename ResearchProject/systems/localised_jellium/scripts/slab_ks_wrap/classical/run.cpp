// ============================================================================
// localised_jellium / scripts/slab_ks_wrap / classical / run.cpp
//
// WRAPPED CLASSICAL TWIN of the CAP-free wrap-around slab KS-stopping study.
// Plan: docs/plans/slab-ks-orbital-stopping-wrap.md
//
// Forked from scripts/classical_highdensity_sv/dyn/run.cpp. The projectile is
// still a moving Gaussian CHARGE (a perturbation, not a ghost UPF) advanced by
// velocity-Verlet Ehrenfest from its own Hellmann-Feynman force. What changes:
//
//   1. IT WRAPS (LJ_WRAP_Z=1). When the centre passes +L_z/2 it is re-entered
//      at -L_z/2, so it makes the SAME ~14 slab crossings the wavepacket makes.
//   2. sigma_WP = 2.0 (sigma_pot = sigma_WP/sqrt2 = 1.41421), matching the WP
//      half. The published benchmark curve is sigma_WP = 0.5 and cannot serve
//      as the reference for a sigma = 2 study.
//   3. Interior checkpoints every LJ_CKPT_EVERY steps, not just a final one
//      (.claude/rules/checkpoint-dont-block.md — these runs are hours long).
//
// ---------------------------------------------------------------------------
// WHY THE WRAP NEEDS THE MINIMUM-IMAGE GAUSSIAN
// inqkit::jellium::gaussian_density builds the blob from a PLAIN Cartesian
// displacement, so a projectile sitting on the +z face loses the half of its
// charge that falls outside the grid: the integral drops to Phi((L/2-b)/sigma)
// and the force is wrong for as long as it straddles. The wavepacket has no
// such problem — the wavefunction basis is a 3-D FFT, periodic in all three
// directions, so it wraps exactly. Using the clipped kernel here would make the
// two twins differ precisely at the boundary this study introduces on purpose.
// LJ_WRAP_Z therefore ALSO switches the perturbation (and the ledger's n_proj)
// to gaussian_density_minimum_image. Pinned by
// inq-stack/tests/include/inqkit/jellium/test_gaussian_minimum_image_engine.cpp.
//
// Note that wrapping the CHARGE does not change the Poisson boundary condition:
// under periodicity(2) a straddling blob is solved as two lumps at opposite ends
// of a z-open box — which is exactly what the solver also does with the
// straddling wavepacket density, so the twins stay matched.
//
// ---------------------------------------------------------------------------
// EXPECTED BEHAVIOUR — NOT A BUG (.claude/rules/light-projectile-stopping.md)
// A mass-1 electron at v = 2.0 carries KE = 54 eV and the benchmark measured
// 27 eV deposited per 25-Bohr slab crossing, so it STOPS after about two
// crossings and then sits. At v = 3.5 (KE = 167 eV, 12.7 eV per crossing) it
// survives the whole run. S is therefore extracted as the INITIAL DRAG over the
// early near-constant-velocity window, never as a full-run regression.
//
// ---------------------------------------------------------------------------
// OUTPUT. projectile.csv gains a proj_z_unwrapped column: proj_z is the wrapped
// coordinate (it jumps by one L_z at each wrap) and proj_z_unwrapped is the
// continuous path, so post-processing never has to guess where a wrap happened.
// Conserved quantity (correctness gate): E_electronic + energy_proj_ke +
// U_proj_bg must be flat.
//
// Env: LJ_LX(35) LJ_LY(35) LJ_LZ(85) LJ_HALF(12.5) LJ_N(100) LJ_EDGE_W(1.0)
//      LJ_PERIODICITY(2) LJ_SPACING(0.40) LJ_SIGMA(2.0) LJ_LAUNCH_Z(-24)
//      LJ_K0(2.0) LJ_MASS(1.0) LJ_DELTA(0.1) LJ_N_STEPS(4529) LJ_DT(0.04)
//      LJ_WRAP_Z(1) LJ_CKPT_EVERY(0=auto N/5) LJ_EXTRA_STATES(24)
//      LJ_GS_DIR(REQUIRED) LJ_OUT(REQUIRED) LJ_CONST_V(0) LJ_SAVE_EVERY(15)
//      LJ_RESUME(0)
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

#include <algorithm>
#include <vector>
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
    const double SPACING=env_d("LJ_SPACING",0.40), SIGMA_WP=env_d("LJ_SIGMA",2.0);
    const double LAUNCH_Z=env_d("LJ_LAUNCH_Z",-24.0);
    const double K0=env_d("LJ_K0",2.0), MASS=env_d("LJ_MASS",1.0), DELTA=env_d("LJ_DELTA",0.1);
    const int N_STEPS=env_i("LJ_N_STEPS",4529); const double DT=env_d("LJ_DT",0.04);
    const int SAVE_EVERY=env_i("LJ_SAVE_EVERY",15);
    const int EXTRA_STATES=env_i("LJ_EXTRA_STATES",24);
    // Periodic wrap: re-enter the projectile at the opposite face so it makes the
    // same repeated slab crossings the wavepacket makes on the FFT grid.
    const bool WRAP_Z=env_i("LJ_WRAP_Z",1)!=0;
    int CKPT_EVERY=env_i("LJ_CKPT_EVERY",0);
    if(CKPT_EVERY<=0) CKPT_EVERY=std::max(1,N_STEPS/3);      // -> 3 retained snapshots
    const bool CONST_V=env_i("LJ_CONST_V",0)!=0;     // 0=Ehrenfest (force-driven), 1=prescribed const velocity
    const std::string GS_DIR=env_s("LJ_GS_DIR",""), OUT="results/"+env_s("LJ_OUT","classical_slab_stopping");
    const double SIGMA_POT=SIGMA_WP/std::sqrt(2.0);
    const double N0=double(N)/(LX*LY*(2.0*HALF));
    if(GS_DIR.empty()||!std::filesystem::exists(GS_DIR)){std::cerr<<"FATAL: GS missing: "<<GS_DIR<<"\n";return 2;}

    // ---- checkpoint / resume (rule final-timestep-checkpoint.md) ----
    const bool RESUME=env_i("LJ_RESUME",0)!=0;
    const std::string CKPT=OUT+"/checkpoint", RT_STATE=OUT+"/rt_state.txt";
    int START=0; double R0z=LAUNCH_Z, V0z=K0/MASS, Z_UNWRAP0=LAUNCH_Z; int NWRAP0=0;
    if(RESUME){
        START=(int)read_state_d(RT_STATE,"last_step",-1);
        if(START<0){std::cerr<<"FATAL: LJ_RESUME=1 but no readable "<<RT_STATE<<"\n";return 2;}
        if(START>=N_STEPS){std::cout<<"Already at/after target ("<<START<<">="<<N_STEPS<<"); nothing to do.\n";return 0;}
        R0z=read_state_d(RT_STATE,"proj_z",LAUNCH_Z); V0z=read_state_d(RT_STATE,"proj_vz",K0/MASS);
        Z_UNWRAP0=read_state_d(RT_STATE,"proj_z_unwrapped",R0z);
        NWRAP0=(int)read_state_d(RT_STATE,"n_wraps",0);
    }
    const std::string SEG=(START>0)?(".from"+std::to_string(START)):std::string("");

    std::cout<<std::setprecision(12)<<"\n=== classical_slab_stopping OUT="<<OUT<<" drive="<<(CONST_V?"const_v":"ehrenfest")
             <<" z0="<<LAUNCH_Z<<" k0="<<K0
             <<" START="<<START<<" -> N_STEPS="<<N_STEPS<<(RESUME?"  [RESUME]":"")<<" ===\n";
    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);                       // NO projectile ion
    auto electrons=systems::electrons(ions,options::electrons{}.spacing(SPACING*1.0_b)
        .extra_electrons(N).extra_states(EXTRA_STATES).temperature(0.00862*1.0_eV),input::kpoints::gamma());
    electrons.load(RESUME?CKPT:GS_DIR);                     // resume from RT ckpt, else GS

    // background well + live projectile (Ehrenfest) + moving perturbation
    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF; bg.slab_axis=2;
    bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    inqkit::dynamics::Projectile proj(MASS, -1.0,
        inqkit::detail::Vec3{0.0,0.0,R0z},
        inqkit::detail::Vec3{0.0,0.0,V0z});                 // resume state or v0=k0/m toward slab
    // WRAP_Z also selects the minimum-image charge kernel: a wrapping blob that
    // is clipped at the face is not the twin of a wavepacket that wraps exactly.
    inqkit::dynamics::moving_gaussian_projectile_perturbation proj_pert(proj, SIGMA_POT, WRAP_Z);
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
        pj<<std::setprecision(12)
          <<"step,time_au,proj_z,proj_vz,energy_proj_ke,energy_proj_bg_ideal,proj_z_unwrapped,n_wraps\n"; }

    // pairwise interaction-energy decomposition (P/S/B); E_BB constant computed once
    const double E_BB = inqkit::jellium::background_self_energy(nplus, phiplus);
    std::ofstream ix;   // Hartree units; sums close to INQ scalars (see interaction_energies.hpp)
    if(electrons.root()){ ix.open(OUT+"/raw/observables/interactions"+SEG+".csv");
        ix<<std::setprecision(12)<<"step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,norm_slab,norm_proj\n"; }

    // density frames (classical: slab electrons responding to the moving projectile)
    // FRAME LAYOUT. The WP binary and inqview's GIF battery both expect
    // raw/vti/density_total/; this binary inherited a bespoke frames/total/ from
    // its ancestor, so the builder silently found NO frames and produced
    // classical notebooks with zero density GIFs (2026-08-01). Written to the
    // canonical path now, with frames/total kept as a symlink so any existing
    // post-processing that hard-codes the old path keeps working.
    const std::string FRAMES = OUT+"/raw/vti/density_total";
    if(SAVE_EVERY>0){
        std::filesystem::create_directories(FRAMES);
        std::error_code ec;
        if(!std::filesystem::exists(OUT+"/frames/total")){
            std::filesystem::create_directories(OUT+"/frames");
            std::filesystem::create_directory_symlink("../raw/vti/density_total",
                                                      OUT+"/frames/total", ec);
        }
    }
    auto save_frame=[&](int step){
        if(SAVE_EVERY<=0) return;
        auto n_tot = inqkit::fields::density::total(electrons);
        inqkit::io::RealField3DLayout lay{.field_name="density",.include_meta=false,
            .emit_raw=false,.emit_vti=true,.vti_format=inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter wr(FRAMES,lay,{.overwrite=true});
        wr.write(n_tot,"density_t"+tag6(step));
    };

    // Continuous-path bookkeeping across wraps.
    double z_unwrapped = Z_UNWRAP0, z_prev_wrapped = R0z;
    int    n_wraps     = NWRAP0;

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
        auto nproj = WRAP_Z
            ? inqkit::jellium::gaussian_density_minimum_image(basis, center, SIGMA_POT)
            : inqkit::jellium::gaussian_density(basis, center, SIGMA_POT);
        const double Uprojbg = -operations::integral_product(nproj, phiplus);
        // pairwise decomposition at the CURRENT step (n_slab = slab electrons; n_P = n_proj)
        auto ct = inqkit::jellium::compute_coulomb(electrons.density(), nproj, phiplus);
        proj.advance(inqkit::detail::Vec3{0.0,0.0,Fz}, DT);                 // V→V_n, R→R_{n+1}

        // PERIODIC WRAP. Relabelling the position across a face is not a
        // dynamical event: the velocity and the accumulated path are untouched,
        // only the reported coordinate moves by one lattice vector. z_unwrapped
        // is the continuous path the stopping fit uses; proj_z is what the
        // Gaussian is actually centred on.
        if(WRAP_Z){
            const double z_before = proj.R().z;
            if(proj.wrap_into_cell(inqkit::detail::Vec3{0.0,0.0,LZ})){
                ++n_wraps;
                std::cout<<"  [wrap] step "<<ctx.step<<"  z "<<z_before
                         <<" -> "<<proj.R().z<<"  (wrap "<<n_wraps<<")\n"<<std::flush;
            }
        }
        z_unwrapped += proj.R().z - z_prev_wrapped;
        if(WRAP_Z && (proj.R().z - z_prev_wrapped) < -0.5*LZ) z_unwrapped += LZ;
        z_prev_wrapped = proj.R().z;

        if(electrons.root()){
            pj<<ctx.step<<","<<ctx.time_au<<","<<Rn.z<<","<<proj.V().z<<","
              <<proj.ke()<<","<<Uprojbg<<","<<z_unwrapped<<","<<n_wraps<<"\n";
            ix<<ctx.step<<","<<ctx.time_au<<","<<ct.e_ss<<","<<ct.e_pp<<","<<ct.e_ps<<","
              <<ct.e_sb<<","<<ct.e_pb<<","<<E_BB<<","<<ct.norm_slab<<","<<ct.norm_p<<"\n";
        }
    });

    // rt_state carries the FULL dynamical state, including what `electrons` does
    // NOT hold: the projectile position/velocity and the wrap bookkeeping.
    // At most MAX_CKPT retained numbered snapshots (user instruction 2026-07-31).
    // Zero-padded names => lexicographic order IS step order, oldest first. The
    // rolling `checkpoint` that LJ_RESUME loads is separate and never pruned.
    const int MAX_CKPT=env_i("LJ_MAX_CKPT",3);
    auto prune_ckpts=[&](){
        if(!electrons.root()) return;
        std::vector<std::filesystem::path> snaps;
        for(auto const& e : std::filesystem::directory_iterator(OUT))
            if(e.is_directory() && e.path().filename().string().rfind("ckpt_step",0)==0)
                snaps.push_back(e.path());
        std::sort(snaps.begin(),snaps.end());
        while((int)snaps.size()>MAX_CKPT){
            std::error_code ec;
            std::filesystem::remove_all(snaps.front(),ec);
            std::cout<<"  [ckpt] pruned "<<snaps.front().filename().string()
                     <<" (keeping newest "<<MAX_CKPT<<")\n"<<std::flush;
            snaps.erase(snaps.begin());
        }
    };

    auto write_rt_state=[&](int last){
        if(!electrons.root()) return;
        std::ofstream st(RT_STATE);
        st<<std::setprecision(12)
          <<"last_step="<<last<<"\ntime_au="<<(last*DT)<<"\ndt="<<DT
          <<"\nproj_z="<<proj.R().z<<"\nproj_vz="<<proj.V().z
          <<"\nproj_z_unwrapped="<<z_unwrapped<<"\nn_wraps="<<n_wraps
          <<"\nproj_mass="<<MASS<<"\nproj_charge="<<-1.0<<"\n";
    };

    auto step_fn=[&](auto const&data){
        rt.step(data);
        // Interior checkpoints: rolling `checkpoint` (what LJ_RESUME loads) plus a
        // RETAINED numbered snapshot, so a killed job loses at most one interval.
        if(data.iter()>0 && data.iter()%CKPT_EVERY==0){
            electrons.save(CKPT);
            write_rt_state(data.iter());
            const std::string snap=OUT+"/ckpt_step"+tag6(data.iter());
            electrons.save(snap);
            prune_ckpts();
            std::cout<<"  [ckpt] step "<<data.iter()<<" -> "<<snap<<"\n"<<std::flush;
        }
    };
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    real_time::propagate(ions,electrons,step_fn,options::theory{}.lda(),opts,pert,START);

    // FINAL checkpoint (rule final-timestep-checkpoint.md): save RT state so this
    // run can be EXTENDED to more steps via LJ_RESUME=1 + a larger LJ_N_STEPS.
    // The LAST TIMESTEP is saved TWICE on purpose: the rolling `checkpoint` that
    // LJ_RESUME=1 loads, AND a STEP-STAMPED ckpt_step<N_STEPS> so the final state
    // is identifiable by step number rather than being an anonymous directory.
    // Pruning runs after it is written and it sorts last, so it is never pruned.
    electrons.save(CKPT);
    {
        const std::string final_snap=OUT+"/ckpt_step"+tag6(N_STEPS);
        electrons.save(final_snap);
        prune_ckpts();
        std::cout<<"  [ckpt] FINAL step "<<N_STEPS<<" -> "<<final_snap
                 <<"  (t = "<<(N_STEPS*DT)<<" a.u.)\n"<<std::flush;
    }
    write_rt_state(N_STEPS);

    if(electrons.root()){ pj.close(); ix.close(); std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)
         <<"run = localised_jellium/slab_ks_wrap/classical/"<<env_s("LJ_OUT","classical")<<"\n"
         <<"campaign = slab-ks-orbital-stopping-wrap (CAP-free, wrap-around)\n"
         <<"plan = docs/plans/slab-ks-orbital-stopping-wrap.md\n"
         <<"engine = inq-study\n"
         <<"wrap_around = "<<(WRAP_Z?"yes":"no")
         <<"  (minimum-image Gaussian charge; proj_z wrapped, proj_z_unwrapped continuous)\n"
         <<"n_wraps = "<<n_wraps<<"  proj_z_unwrapped_final = "<<z_unwrapped<<"\n"
         <<"slab_half_width = "<<HALF<<"  edge_width = "<<EDGE_W<<"\n"
         <<"n0_a0m3 = "<<N0<<"\n"
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
