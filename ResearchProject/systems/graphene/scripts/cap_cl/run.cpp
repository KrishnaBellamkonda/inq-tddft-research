// ============================================================================
// graphene CAP scattering — CLASSICAL projectile (Ehrenfest point electron).
//
// Companion to run.cpp (WP mode). A classical electron = a moving ion with the
// Gaussian-smeared electron pseudopotential (sigma=1.47 Bohr = WP width) and
// mass = m_e, launched at the WP launch point with v = k0 (+z). Ensemble member
// drawn from Gaussian position+momentum (std::mt19937) matching the WP, per
// Yao & Schliefe. Same two-sided sin^2 CAP (eta=-0.5, L=20) on the z-ends.
//
// Build vs inq-study: INQ_SOURCE=.../inq-study inq-run
// Env: GR_E_EV(100), GR_CX, GR_CY (impact pt), GR_CAP(1), GR_OUTDIR, GR_DT(0.02),
//      GR_NSTEPS(0=auto), GR_SEED(0=central/no-jitter; >0 = random draw), GR_TAG.
// All CAP results PROVISIONAL until inq-study Task #7.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/eigenvalue_dump.hpp>
#include <inqkit/observables/minimum_observable_set.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>
#include "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared/configs/graphene_gs.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace fs=std::filesystem; namespace cfg=graphene_cfg; namespace obs_=inqkit::observables;
static double env_d(const char*k,double d){const char*v=std::getenv(k);return v?std::atof(v):d;}
static int    env_i(const char*k,int d){const char*v=std::getenv(k);return v?std::atoi(v):d;}
static std::string env_s(const char*k,const std::string&d){const char*v=std::getenv(k);return v?std::string(v):d;}

int main(){
  const double HA_TO_EV=27.211386245988;
  const double E_eV=env_d("GR_E_EV",100.0);
  const double cx=env_d("GR_CX",0.0), cy=env_d("GR_CY",0.0);
  const bool capon=env_i("GR_CAP",1)!=0;
  const double dt=env_d("GR_DT",0.02);
  const int seed=env_i("GR_SEED",0);
  const std::string outdir=env_s("GR_OUTDIR","results"), tag=env_s("GR_TAG","cl");

  const double sigma=1.47, k0=std::sqrt(2.0*E_eV/HA_TO_EV);
  const double Lz=cfg::LZ_BOHR, Lcap=20.0, Lhalf=Lcap/2.0, z_in=(Lz-Lcap)/2.0;
  const double z0=-z_in+5.0*sigma, eta=-0.5;
  const double mid_frac=(z_in+Lhalf/2.0)/Lz, width_frac=Lhalf/Lz;
  const double tau=((z_in-z0)+Lcap+8.0)/k0+4.0;
  int N_STEPS=env_i("GR_NSTEPS",0); if(N_STEPS<=0) N_STEPS=std::max(1,(int)std::llround(tau/dt));
  const int WRITE_EVERY=std::max(1,N_STEPS/60), OBS_EVERY=std::max(1,N_STEPS/200), SNAP_EVERY=std::max(1,N_STEPS/20);
  const double sigma_k=1.0/(std::sqrt(2.0)*sigma);

  // ensemble jitter (Gaussian pos+mom matching the WP); seed=0 => central member
  double dx=0,dy=0,dz=0,dvx=0,dvy=0,dvz=0;
  if(seed>0){ std::mt19937 rng(seed); std::normal_distribution<double> gp(0.0,sigma), gk(0.0,sigma_k);
    dx=gp(rng); dy=gp(rng); dz=gp(rng); dvx=gk(rng); dvy=gk(rng); dvz=gk(rng); }
  const double px=cx+dx, py=cy+dy, pz=z0+dz;
  const double vx=dvx, vy=dvy, vz=k0+dvz;

  fs::create_directories(outdir+"/raw/observables/eigenvalues");
  fs::create_directories(outdir+"/raw/vti"); fs::create_directories(outdir+"/raw/screens/total");
  std::printf("\n=== graphene CLASSICAL [%s] E=%.1f k0=%.4f CAP=%d seed=%d ===\n",tag.c_str(),E_eV,k0,(int)capon,seed);
  std::printf("  proj r=(%.3f,%.3f,%.3f) v=(%.3f,%.3f,%.3f) N=%d\n",px,py,pz,vx,vy,vz,N_STEPS);

  auto cell=systems::cell::orthorhombic(cfg::LX_BOHR*1.0_b,cfg::LY_BOHR*1.0_b,cfg::LZ_BOHR*1.0_b).periodic();
  auto ions=systems::ions::parse(cfg::GEOMETRY_XYZ,cell);
  const int n_carbon=ions.size();
  // insert classical-electron projectile (Gaussian UPF, mass m_e).
  // The UPF carries z_valence=-1 (consistent with its +1/r repulsive local
  // potential) so the projectile is a proper -1 charge in the ion-ion Ewald
  // sum (interaction.hpp:329 uses valence_charge()) and thus FEELS the carbon
  // nuclei (attraction), cancelling its electron-cloud repulsion at long range.
  // With z_valence=0 it felt only the electrons -> spurious vacuum deceleration.
  auto sp=ionic::species("H").pseudo_file(
            "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared/pseudopotentials/electron_gaussian_sigma1p47_zm1.upf")
          .mass(1.0/1822.8885);
  ions.insert(sp,{px*1.0_b,py*1.0_b,pz*1.0_b});
  const int proj=ions.size()-1;  // projectile index (last)

  // .extra_electrons(+1) compensates the projectile's z_valence=-1 so the
  // QUANTUM electron count stays at the graphene value (96); the cell then
  // carries the physical net -1 charge (projectile) + a +1 uniform background.
  auto electrons=systems::electrons(ions,options::electrons{}.cutoff(cfg::CUTOFF_HA*1.0_Ha)
                       .extra_electrons(1.0)
                       .extra_states(cfg::EXTRA_STATES).temperature(cfg::TEMPERATURE_EV*1.0_eV));
  if(!fs::exists(cfg::GS_CHECKPOINT_DIR)){std::fprintf(stderr,"FATAL: no GS\n");return 2;}
  electrons.load(cfg::GS_CHECKPOINT_DIR);
  ions.velocities()[proj]=vector3<double>{vx,vy,vz};
  std::printf("  carbons=%d proj_idx=%d proj_mass_au=%.4f (expect 1) v_set=(%.3f,%.3f,%.3f)\n",
              n_carbon,proj,ions.species(proj).mass(),vx,vy,vz);

  { std::ofstream s(outdir+"/run_summary.txt"); s<<"run_completed = false\ntag = "<<tag<<"\n"; }
  obs_::write_manifest(outdir+"/observables_manifest.json",obs_::RunType::jellium_classical,WRITE_EVERY,N_STEPS);

  using RLay=inqkit::io::RealField3DLayout; using ROpt=inqkit::io::RealField3DWriteOptions;
  const auto VBIN=inqkit::io::VTIWriteOptions::Format::binary;
  auto mkR=[&](const std::string&d,const std::string&n){return inqkit::io::RealField3DWriter(d,
      RLay{.field_name=n,.include_meta=false,.emit_raw=false,.emit_vti=true,.vti_format=VBIN},ROpt{.overwrite=true});};
  auto total_wr=mkR(outdir+"/raw/vti/density_rt_total","density");
  auto system_wr=mkR(outdir+"/raw/vti/density_rt_system","density");
  { auto g=mkR(outdir+"/raw/vti/density_gs_system","density"); g.write(inqkit::fields::density::total(electrons),"density_gs_system"); }

  inqkit::io::ObservableSelection sel; sel.step=sel.time_au=true;
  sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true; sel.density_l2=true;
  sel.current_x=sel.current_y=sel.current_z=true; sel.dipole_x=sel.dipole_y=sel.dipole_z=true;
  inqkit::io::ObservablesWriter obs_writer(outdir+"/raw/observables/observables.csv",sel); obs_writer.write_header();
  obs_::DensityDelta density_delta(outdir+"/raw/vti/density_delta",outdir+"/raw/vti/density_delta_coarse",{.compute_l2=true,.coarse_bin_bohr=3.0});
  std::ofstream etrack(outdir+"/raw/observables/electron_track.csv"); etrack<<std::setprecision(12);
  etrack<<"step,time_au,x,y,z,vx,vy,vz,fx,fy,fz\n";

  using inqkit::screens::LeedPatternAccumulator; using inqkit::screens::PlaneScreen;
  const std::array<double,8> sz={-16,-12,-8,-4,4,8,12,16}; std::array<LeedPatternAccumulator,8> acc;
  for(int k=0;k<8;++k){char lb[16];std::snprintf(lb,sizeof lb,"scr_%+03d",(int)sz[k]);acc[k]=LeedPatternAccumulator(PlaneScreen{sz[k],std::string(lb)});}

  auto cap=perturbations::absorbing(eta*1.0_Ha,mid_frac,width_frac)+perturbations::absorbing(eta*1.0_Ha,-mid_frac,width_frac);

  { auto t0=inqkit::fields::density::total(electrons); total_wr.write(t0,0.0,0); system_wr.write(t0,0.0,0);
    auto const&p=ions.positions()[proj]; auto const&v=ions.velocities()[proj];
    etrack<<0<<",0,"<<p[0]<<','<<p[1]<<','<<p[2]<<','<<v[0]<<','<<v[1]<<','<<v[2]<<",0,0,0\n"; }

  auto t_start=std::chrono::steady_clock::now();
  auto cb=[&](auto const&data){
    const int step=data.iter(); const double t=data.time();
    { auto const&p=ions.positions()[proj]; auto const&v=ions.velocities()[proj];
      etrack<<step<<','<<t<<','<<p[0]<<','<<p[1]<<','<<p[2]<<','<<v[0]<<','<<v[1]<<','<<v[2]<<",0,0,0\n"; }
    if(step%OBS_EVERY==0){
      inqkit::StepContext ctx; ctx.step=step; ctx.time_au=t; ctx.ions=&ions; ctx.electrons=&electrons;
      ctx.energy_total=data.energy().total(); ctx.energy_kinetic=data.energy().kinetic();
      ctx.energy_hartree=data.energy().hartree(); ctx.energy_xc=data.energy().xc();
      try{auto c=data.current();ctx.current={c[0],c[1],c[2]};}catch(...){}
      try{auto dp=data.dipole();ctx.dipole={dp[0],dp[1],dp[2]};}catch(...){}
      auto tt=inqkit::fields::density::total(electrons); ctx.density_l2=density_delta.snapshot(tt,t,step);
      obs_writer.append(ctx);
    }
    if(step%WRITE_EVERY==0){ auto tt=inqkit::fields::density::total(electrons); total_wr.write(tt,t,step); system_wr.write(tt,t,step); }
    for(auto&a:acc) a.accumulate(electrons,dt);
  };

  auto theory=options::theory{}.lda();
  // .ehrenfest() is REQUIRED: INQ freezes ions by default; without it the
  // projectile moves ballistically (no deceleration => zero stopping).
  auto rt=options::real_time{}.num_steps(N_STEPS).dt(dt*1.0_atomictime).ehrenfest();
  if(capon) real_time::propagate(ions,electrons,cb,theory,rt,cap);
  else      real_time::propagate(ions,electrons,cb,theory,rt);
  etrack.close();
  for(int k=0;k<8;++k){char fn[256];std::snprintf(fn,sizeof fn,"%s/raw/screens/total/scr%+03d.dat",outdir.c_str(),(int)sz[k]);acc[k].save(std::string(fn));}
  obs_::dump_eigenvalues(electrons,outdir+"/raw/observables/eigenvalues/eigenvalues.csv");

  double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t_start).count();
  auto const&pf=ions.positions()[proj]; auto const&vf=ions.velocities()[proj];
  const double KE_f=0.5*(vf[0]*vf[0]+vf[1]*vf[1]+vf[2]*vf[2])*HA_TO_EV;
  { std::ofstream s(outdir+"/run_summary.txt"); s.precision(12);
    s<<"run_completed = true\nrun_type = graphene-classical-cap\ntag = "<<tag<<"\nE_eV = "<<E_eV
     <<"\nk0 = "<<k0<<"\ncap = "<<(int)capon<<"\neta_Ha = "<<eta<<"\nL_cap = "<<Lcap<<"\nseed = "<<seed
     <<"\nsigma = "<<sigma<<"\nproj_r0 = "<<px<<' '<<py<<' '<<pz<<"\nproj_v0 = "<<vx<<' '<<vy<<' '<<vz
     <<"\nproj_rf = "<<pf[0]<<' '<<pf[1]<<' '<<pf[2]<<"\nproj_vf = "<<vf[0]<<' '<<vf[1]<<' '<<vf[2]
     <<"\nKE_final_eV = "<<KE_f<<"\nKE_loss_eV = "<<(E_eV-KE_f)<<"\nn_carbon = "<<n_carbon<<"\nproj_idx = "<<proj
     <<"\ndt = "<<dt<<"\nN_STEPS = "<<N_STEPS<<"\nwrite_every = "<<WRITE_EVERY
     <<"\ncell_bohr = "<<cfg::LX_BOHR<<' '<<cfg::LY_BOHR<<' '<<cfg::LZ_BOHR
     <<"\npropagator = etrs (Ehrenfest)\nwall_s = "<<wall<<"\nPROVISIONAL = inq-study Task#7\n"; }
  std::printf("  [%s] KE_final=%.2f eV KE_loss=%.2f eV wall=%.1fs -> %s\n",tag.c_str(),KE_f,E_eV-KE_f,wall,outdir.c_str());
  return 0;
}
