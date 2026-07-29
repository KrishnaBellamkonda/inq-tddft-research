// ============================================================================
// classical_highdensity_sv / pilot_native / run.cpp   (Phase 3 pilot — RUN B)
//
// CONTRAST pilot: single-transit INQ NATIVE Ehrenfest of a REAL ghost-UPF ion
// through the SAME high-density localised jellium slab as RUN A (pilot/). This
// is the real-system faithfulness check for the perturbation projectile:
// identical slab, identical projectile potential (V_proj = clean +erf/r ghost
// UPF, same sigma_pot), but the ion is moved by INQ's OWN native Ehrenfest
// integrator (options::real_time{}.ehrenfest() => velocity-Verlet inside ETRS,
// a = F_localHF / species.mass()) instead of our perturbation velocity-Verlet.
//
// The slab jellium is a BACKGROUND PERTURBATION (as in the GS build and RUN A);
// the ghost projectile is a real mass-1 ion (z_valence=0, only its +erf/r local
// potential couples to the KS system). GS is recomputed HERE with the ghost ion
// present at launch_z AND the slab background perturbation (matching RUN A's GS
// recipe) — the perturbation GS checkpoint can't be reused (it has no ghost ion).
//
// GATE / decisive fact: the ghost ion MUST move (z changes). If native Ehrenfest
// cannot move a z_valence=0 ion, z(t) stays flat and THAT is the finding — the
// run records it clearly and exits 0; RUN A + its notebook still ship.
//
// Records per step -> native.csv:
//   step,time,z,vz,E_elec,E_total,E_kin,E_hartree,E_external,E_nonlocal,E_xc,E_ion
// E_elec := kinetic+hartree+external+nonlocal+xc (electronic only). Density frames
// (frames/total/density_t*.vti) for the density-evolution GIF.
//
// Env: TC_LX(35) TC_LY(35) TC_LZ(85) TC_HALF(12.5) TC_N(100) TC_EDGE_W(1.0)
//   TC_PERIODICITY(2) TC_SPACING(0.5) TC_SIGMA_POT(0.35355) TC_LAUNCH_Z(-30)
//   TC_V0(2.0) TC_MASS_AU(1.0) TC_N_STEPS(1600) TC_DT(0.04) TC_SAVE_EVERY(5)
//   TC_GHOST_UPF(REQUIRED, ghost_sigma0p354.upf) TC_OUT(pilot_native)
// No inq/ or inq-study/ edit — wrapper-only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

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
static std::string env_s(const char* k,const std::string& d){const char*v=std::getenv(k);return v?std::string(v):d;}
static std::string tag6(int n){std::ostringstream o;o<<std::setw(6)<<std::setfill('0')<<n;return o.str();}

int main(){
    const double LX=env_d("TC_LX",35),LY=env_d("TC_LY",35),LZ=env_d("TC_LZ",85),HALF=env_d("TC_HALF",12.5);
    const int    N=env_i("TC_N",100); const double EDGE_W=env_d("TC_EDGE_W",1.0);
    const int    PER=env_i("TC_PERIODICITY",2);
    const double SPACING=env_d("TC_SPACING",0.5);
    const double SIGMA_POT=env_d("TC_SIGMA_POT",0.35355);
    const double LAUNCH_Z=env_d("TC_LAUNCH_Z",-30.0), V0=env_d("TC_V0",2.0);
    const double MASS_AU=env_d("TC_MASS_AU",1.0);        // projectile mass in atomic units (m_e)
    const int    N_STEPS=env_i("TC_N_STEPS",1600);
    const double DT=env_d("TC_DT",0.04);
    const int    SAVE_EVERY=env_i("TC_SAVE_EVERY",5);
    const std::string GHOST_UPF=env_s("TC_GHOST_UPF","ghost_sigma0p354.upf");
    const std::string OUT="results/"+env_s("TC_OUT","pilot_native");
    const double N0=double(N)/(LX*LY*(2.0*HALF));

    if(!std::filesystem::exists(GHOST_UPF)){std::cerr<<"FATAL: ghost UPF missing: "<<GHOST_UPF<<"\n";return 2;}

    std::cout<<std::setprecision(12)
             <<"\n=== B: NATIVE Ehrenfest ghost-UPF ion  slab_n100 35x35x85  dx="<<SPACING
             <<"  launch_z="<<LAUNCH_Z<<" v0="<<V0<<" mass_au="<<MASS_AU
             <<" N="<<N_STEPS<<" dt="<<DT<<"  ghost="<<GHOST_UPF<<" ===\n";

    auto cell0=systems::cell::orthorhombic(LX*1.0_b,LY*1.0_b,LZ*1.0_b);
    auto cell=(PER==2)?cell0.periodicity(2):cell0.periodic();
    auto ions=systems::ions(cell);

    // Ghost projectile ion: mass 1.0 a.u. (= m_e) via .mass(1/1822.8885 amu),
    // z_valence=0 (no electrons), local V_loc = +erf(r/(sqrt2.sigma_pot))/r.
    auto ghost = ionic::species("H").pseudo_file(GHOST_UPF).mass(MASS_AU/1822.8885);
    ions.insert(ghost, {0.0_b, 0.0_b, LAUNCH_Z*1.0_b});   // index 0
    const int GH = 0;

    auto electrons=systems::electrons(ions,
        options::electrons{}.spacing(SPACING*1.0_b)
            .extra_electrons(N).extra_states(24).temperature(0.00862*1.0_eV),
        input::kpoints::gamma());

    // slab jellium background as a perturbation (as in RUN A + GS build)
    inqkit::jellium::localised_background_params bg;
    bg.shape=inqkit::jellium::background_shape::slab; bg.n0=N0; bg.half_width=HALF;
    bg.slab_axis=2; bg.center={0.0,0.0,0.0}; bg.edge_width=EDGE_W;
    inqkit::jellium::localised_background_perturbation bgpert(bg);

    // GS WITH ghost ion present + slab background (own GS; can't reuse pert GS).
    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(ions, electrons, options::theory{}.lda(),
        options::ground_state{}.energy_tolerance(1e-4_Ha).max_steps(300)
            .broyden_mixing().mixing_ndim(8).mixing(0.1),
        bgpert);
    std::cout<<"  GS energy = "<<gs.energy.total()<<" Ha\n";

    // Initial velocity: +z toward/through slab.
    ions.velocities()[GH] = vector3<double>{0.0, 0.0, V0};

    std::filesystem::create_directories(OUT);
    std::filesystem::create_directories(OUT+"/frames/total");
    std::ofstream csv;
    if(electrons.root()){
        csv.open(OUT+"/native.csv"); csv<<std::setprecision(12)
         <<"step,time,z,vz,E_elec,E_total,E_kin,E_hartree,E_external,E_nonlocal,E_xc,E_ion\n";
    }
    auto save_frame=[&](int step){
        if(SAVE_EVERY<=0) return;
        auto n_tot = inqkit::fields::density::total(electrons);
        inqkit::io::RealField3DLayout lay{.field_name="density",.include_meta=false,
            .emit_raw=false,.emit_vti=true,.vti_format=inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter wr(OUT+"/frames/total",lay,{.overwrite=true});
        wr.write(n_tot,"density_t"+tag6(step));
    };

    auto func=[&](auto data){
        auto pos = data.positions();
        auto vel = data.velocities();
        auto en  = data.energy();
        const double zc = pos[GH][2], vz = vel[GH][2];
        const double E_elec = en.kinetic()+en.hartree()+en.external()+en.non_local()+en.xc();
        if(electrons.root()){
            csv<<data.iter()<<","<<data.time()<<","<<zc<<","<<vz<<","
               <<E_elec<<","<<en.total()<<","<<en.kinetic()<<","<<en.hartree()<<","
               <<en.external()<<","<<en.non_local()<<","<<en.xc()<<","<<en.ion()<<"\n";
        }
        if(SAVE_EVERY>0 && data.iter() % SAVE_EVERY == 0) save_frame(data.iter());
    };

    // native Ehrenfest + slab background perturbation
    auto opts=options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime).ehrenfest();
    real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts, bgpert);

    const double z_final = ions.positions()[GH][2], vz_final = ions.velocities()[GH][2];
    const bool moved = std::abs(z_final - LAUNCH_Z) > 1e-3;
    if(electrons.root()){
        csv.close();
        std::ofstream s(OUT+"/run_summary.txt");
        s<<std::setprecision(12)<<"run = classical_highdensity_sv/pilot_native (native Ehrenfest ghost-UPF ion)\n"
         <<"scheme = INQ native Ehrenfest (ionic::propagator::molecular_dynamics, velocity-Verlet in ETRS)\n"
         <<"engine = inq\nperiodicity = "<<PER<<"\n"
         <<"cell_bohr = "<<LX<<"x"<<LY<<"x"<<LZ<<"  slab_half = "<<HALF<<"  edge_w = "<<EDGE_W<<"\n"
         <<"N = "<<N<<"  n0 = "<<N0<<"  spacing = "<<SPACING<<"  sigma_pot = "<<SIGMA_POT<<"\n"
         <<"launch_z = "<<LAUNCH_Z<<"  v0 = "<<V0<<"  mass_au = "<<MASS_AU<<"\n"
         <<"n_steps = "<<N_STEPS<<"  dt = "<<DT<<"\n"
         <<"z_final = "<<z_final<<"  vz_final = "<<vz_final<<"  moved = "<<(moved?"true":"false")<<"\n"
         <<"ghost_upf = "<<GHOST_UPF<<"\ncsv = native.csv\nrun_completed = true\n";
    }
    std::cout<<"  done  z_final="<<z_final<<" vz_final="<<vz_final<<"  moved="<<(moved?"YES":"NO (ghost did not move)")<<"\n";
    return 0;
}
