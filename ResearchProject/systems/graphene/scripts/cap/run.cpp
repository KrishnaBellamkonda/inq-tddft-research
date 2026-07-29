// ============================================================================
// graphene CAP scattering — env-driven build-once binary (WP projectile).
//
// Feasibility replica of Yao & Schleife. Plan: docs/plans/graphene-cap.md.
// Loads the converged graphene GS (24 C, 3x2), injects a Gaussian electron WP,
// propagates (ETRS) with an OPTIONAL two-sided sin^2 CAP on the z-ends, and
// emits the coronene-style observable set + survival fraction + LEED screens.
//
// Build against inq-study (the CAP needs the complexified scalar potential):
//   INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study inq-run
//
// Env (all optional, sensible defaults):
//   GR_E_EV     projectile kinetic energy (default 100)
//   GR_CX,GR_CY lateral impact point Bohr (default 0,0 = centroid/atom;
//               channeling/hollow set per geometry)
//   GR_CAP      1 = two-sided CAP on (default), 0 = no-CAP baseline
//   GR_OUTDIR   output dir (default results)
//   GR_DT       timestep a.u. (default 0.02)
//   GR_NSTEPS   override N_steps (default 0 => auto from travel)
//   GR_TAG      label string for run_summary
//
// CAP geometry (locked): two-sided sin^2, total L=20 Bohr (10/end), W=-0.5 Ha,
//   free region |z|<20 (z_in=20). z-cell 60 Bohr.
// WP: sigma_r=1.47 Bohr, launch z0 = -20 + 5*sigma = -12.65 Bohr, +z.
// All CAP results PROVISIONAL until inq-study engine regression (Task #7).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/absorbers/mask_absorber.hpp>       // inner_region_norm[_twosided]
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/eigenvalue_dump.hpp>
#include <inqkit/observables/minimum_observable_set.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include "/local/data/public/skcb2/tddft/ResearchProject/systems/graphene/shared/configs/graphene_gs.hpp"

#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace fs   = std::filesystem;
namespace cfg  = graphene_cfg;
namespace abs_ = inqkit::absorbers;
namespace obs_ = inqkit::observables;

static double env_d(const char *k, double d){ const char*v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char *k, int d)   { const char*v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char*k,const std::string&d){ const char*v=std::getenv(k); return v?std::string(v):d; }

int main() {
  const double HA_TO_EV = 27.211386245988;
  const double E_eV  = env_d("GR_E_EV", 100.0);
  const double cx    = env_d("GR_CX", 0.0);
  const double cy    = env_d("GR_CY", 0.0);
  const bool   capon = env_i("GR_CAP", 1) != 0;
  const double dt    = env_d("GR_DT", 0.02);
  const std::string outdir = env_s("GR_OUTDIR", "results");
  const std::string tag    = env_s("GR_TAG", "wp");

  // --- geometry / kinematics (locked) ---
  const double sigma = 1.47;                          // WP sigma_r (d=1.1 Ang)
  const double k0    = std::sqrt(2.0 * E_eV / HA_TO_EV);
  const double Lz    = cfg::LZ_BOHR;                  // 60
  const double Lcap  = 20.0;                          // total CAP width
  const double Lhalf = Lcap / 2.0;                    // 10 per end
  const double z_in  = (Lz - Lcap) / 2.0;             // 20  (free-region half-width)
  const double z0    = -z_in + 5.0 * sigma;           // launch 5 sigma inside free edge = -12.65
  const double eta   = -0.5;                          // W (CAP depth, Ha)
  const double mid_frac   = (z_in + Lhalf/2.0) / Lz;  // 25/60 = 0.4167
  const double width_frac = Lhalf / Lz;               // 10/60 = 0.1667
  // travel: z0 -> +z_in -> through far CAP, plus absorption tail.
  const double tau = ((z_in - z0) + Lcap + 8.0) / k0 + 4.0;
  int N_STEPS = env_i("GR_NSTEPS", 0);
  if (N_STEPS <= 0) N_STEPS = std::max(1, (int)std::llround(tau / dt));
  const int WRITE_EVERY = std::max(1, N_STEPS / 60);
  const int OBS_EVERY   = std::max(1, N_STEPS / 200);
  const int SNAP_EVERY  = std::max(1, N_STEPS / 20);

  fs::create_directories(outdir + "/raw/observables/eigenvalues");
  fs::create_directories(outdir + "/raw/vti");
  fs::create_directories(outdir + "/raw/screens/total");
  fs::create_directories(outdir + "/raw/screens/instantaneous");

  std::printf("\n=== graphene CAP [%s] E=%.1f eV k0=%.4f CAP=%d ===\n", tag.c_str(), E_eV, k0, (int)capon);
  std::printf("  z0=%.3f z_in=%.1f mid_frac=%.4f width_frac=%.4f dt=%.3f N=%d WE=%d cx=%.2f cy=%.2f\n",
              z0, z_in, mid_frac, width_frac, dt, N_STEPS, WRITE_EVERY, cx, cy);

  // --- cell + ions + GS ---
  auto cell = systems::cell::orthorhombic(cfg::LX_BOHR*1.0_b, cfg::LY_BOHR*1.0_b, cfg::LZ_BOHR*1.0_b).periodic();
  auto ions = systems::ions::parse(cfg::GEOMETRY_XYZ, cell);
  auto electrons = systems::electrons(
      ions, options::electrons{}.cutoff(cfg::CUTOFF_HA*1.0_Ha)
                                .extra_states(cfg::EXTRA_STATES)
                                .temperature(cfg::TEMPERATURE_EV*1.0_eV));
  if (!fs::exists(cfg::GS_CHECKPOINT_DIR)) { std::fprintf(stderr,"FATAL: no GS checkpoint\n"); return 2; }
  electrons.load(cfg::GS_CHECKPOINT_DIR);
  const int n_states = electrons.states().num_states();
  const int n_elec   = electrons.states().num_electrons();
  std::printf("  loaded GS: states=%d electrons=%d\n", n_states, n_elec);

  // stub summary
  { std::ofstream s(outdir+"/run_summary.txt");
    s<<"run_completed = false\ntag = "<<tag<<"\nE_eV = "<<E_eV<<"\ncap = "<<(int)capon<<"\n"; }

  // --- WP injection (orthogonalised against the occupied graphene subspace) ---
  auto rep = inqkit::WavePacket{}
                 .center(cx, cy, z0).sigma(sigma).k0(0.0, 0.0, k0)
                 .orthogonalise_against_occupied(electrons)
                 .inject_into_last_extra_state(electrons, 1.0);
  const long wp_idx = rep.state_index;
  const double N0 = abs_::inner_region_norm(electrons, 2, +1e12, wp_idx);
  std::printf("  WP injected: idx=%ld norm_after=%.5f max_overlap=%.2e N0=%.5f\n",
              wp_idx, rep.norm_after, rep.max_overlap, N0);

  obs_::write_manifest(outdir+"/observables_manifest.json", obs_::RunType::coronene, WRITE_EVERY, N_STEPS);

  // --- VTI writers ---
  using RLay = inqkit::io::RealField3DLayout; using ROpt = inqkit::io::RealField3DWriteOptions;
  const auto VBIN = inqkit::io::VTIWriteOptions::Format::binary;
  auto mkR=[&](const std::string&dir,const std::string&name){
    return inqkit::io::RealField3DWriter(dir, RLay{.field_name=name,.include_meta=false,.emit_raw=false,.emit_vti=true,.vti_format=VBIN}, ROpt{.overwrite=true}); };
  auto total_wr  = mkR(outdir+"/raw/vti/density_rt_total","density");
  auto system_wr = mkR(outdir+"/raw/vti/density_rt_system","density");
  auto wp_wr     = mkR(outdir+"/raw/vti/density_rt_wp","density");
  // GS system density once
  { auto g=mkR(outdir+"/raw/vti/density_gs_system","density");
    g.write(inqkit::fields::density::total_excluding_orbital(electrons, wp_idx, 1.0), "density_gs_system"); }
  // WP complex wavefunction (enables WP current density in post-processing)
  inqkit::io::ComplexField3DWriter wp_wf_wr(outdir+"/raw/vti/wavefunction_wp_rt",
      {.field_name="wavefunction",.include_meta=false,.emit_raw=false,.emit_vti=true,.vti_format=VBIN}, {.overwrite=true});

  // --- observables CSV ---
  inqkit::io::ObservableSelection sel;
  sel.step=sel.time_au=true; sel.energy_total=sel.energy_kinetic=sel.energy_hartree=sel.energy_xc=true;
  sel.density_l2=true; sel.current_x=sel.current_y=sel.current_z=true; sel.dipole_x=sel.dipole_y=sel.dipole_z=true;
  inqkit::io::ObservablesWriter obs_writer(outdir+"/raw/observables/observables.csv", sel);
  obs_writer.write_header();
  obs_::WPMomentumStats wp_mom(outdir+"/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every=OBS_EVERY});
  obs_::MomentumDistribution mom_dist(outdir+"/raw/observables/momentum_distribution.csv", wp_idx, Lz,
                                      {.n_bins=64,.k_max_bohr_inv=0.0,.write_every=OBS_EVERY});
  inqkit::observables::OrbitalOverlapMatrix overlap_obs(electrons, wp_idx, outdir+"/raw/observables/overlap");
  obs_::DensityDelta density_delta(outdir+"/raw/vti/density_delta", outdir+"/raw/vti/density_delta_coarse",
                                   {.compute_l2=true,.coarse_bin_bohr=3.0});
  std::ofstream inner_csv(outdir+"/raw/observables/inner_norm_vs_time.csv");
  inner_csv << "step,time_au,total_wp_norm,survival_inner_over_N0,absorbed_fraction\n";

  // --- LEED screens: 8 (4 per side) at z = +/-4,8,12,16 ---
  using inqkit::screens::LeedPatternAccumulator; using inqkit::screens::PlaneScreen;
  const std::array<double,8> screen_z = {-16,-12,-8,-4,4,8,12,16};
  std::array<LeedPatternAccumulator,8> acc;
  for (int k=0;k<8;++k){ char lb[16]; std::snprintf(lb,sizeof lb,"scr_%+03d",(int)screen_z[k]);
    acc[k]=LeedPatternAccumulator(PlaneScreen{screen_z[k], std::string(lb)}); }

  // --- two-sided CAP (built only if capon) ---
  auto cap = perturbations::absorbing(eta*1.0_Ha,  mid_frac, width_frac)
           + perturbations::absorbing(eta*1.0_Ha, -mid_frac, width_frac);

  // t=0 density frames
  { auto sys0=inqkit::fields::density::total_excluding_orbital(electrons,wp_idx,1.0);
    auto wp0 =inqkit::fields::density::orbital(electrons,wp_idx);
    auto tot0=inqkit::fields::density::total(electrons);
    system_wr.write(sys0,0.0,0); wp_wr.write(wp0,0.0,0); total_wr.write(tot0,0.0,0); }

  auto t_start=std::chrono::steady_clock::now();
  auto step_cb=[&](auto const&data){
    const int step=data.iter(); const double t=data.time();
    if (step % OBS_EVERY == 0){
      inqkit::StepContext ctx; ctx.step=step; ctx.time_au=t; ctx.ions=&ions; ctx.electrons=&electrons;
      ctx.energy_total=data.energy().total(); ctx.energy_kinetic=data.energy().kinetic();
      ctx.energy_hartree=data.energy().hartree(); ctx.energy_xc=data.energy().xc();
      try{auto c=data.current(); ctx.current={c[0],c[1],c[2]};}catch(...){}
      try{auto dp=data.dipole(); ctx.dipole={dp[0],dp[1],dp[2]};}catch(...){}
      { auto bath=inqkit::fields::density::total_excluding_orbital(electrons,wp_idx,1.0);
        ctx.density_l2=density_delta.snapshot(bath,t,step); }
      obs_writer.append(ctx);
      overlap_obs.snapshot_wp_only(electrons,t,step);
      double tot=abs_::inner_region_norm(electrons,2,+1e12,wp_idx);
      double inn=abs_::inner_region_norm_twosided(electrons,2,z_in,wp_idx)/N0;
      inner_csv<<step<<','<<t<<','<<tot<<','<<inn<<','<<(1.0-tot/N0)<<'\n';
    }
    if (step % WRITE_EVERY == 0){
      auto sys=inqkit::fields::density::total_excluding_orbital(electrons,wp_idx,1.0);
      auto wp =inqkit::fields::density::orbital(electrons,wp_idx);
      auto tot=inqkit::fields::density::total(electrons);
      system_wr.write(sys,t,step); wp_wr.write(wp,t,step); total_wr.write(tot,t,step);
      auto psi=inqkit::fields::orbital::wavefunction(electrons,wp_idx);
      char wf[48]; std::snprintf(wf,sizeof wf,"wavefunction_t%06d",step); wp_wf_wr.write(psi,std::string(wf));
    }
    for (auto&a:acc) a.accumulate(electrons,dt);
    if (step % SNAP_EVERY == 0){
      for (int k=0;k<8;++k){ auto sl=acc[k].screen().extract(electrons);
        char fn[256]; std::snprintf(fn,sizeof fn,"%s/raw/screens/instantaneous/scr%+03d_t%06d.dat",
                                    outdir.c_str(),(int)screen_z[k],step);
        acc[k].screen().save(sl,t,std::string(fn)); } }
    mom_dist.maybe_accumulate(data);
    wp_mom.maybe_accumulate(data);
  };

  auto theory = options::theory{}.lda();
  auto rt = options::real_time{}.num_steps(N_STEPS).dt(dt*1.0_atomictime);  // ETRS default
  if (capon) real_time::propagate(ions, electrons, step_cb, theory, rt, cap);
  else       real_time::propagate(ions, electrons, step_cb, theory, rt);

  inner_csv.close();
  for (int k=0;k<8;++k){ char fn[256]; std::snprintf(fn,sizeof fn,"%s/raw/screens/total/scr%+03d.dat",outdir.c_str(),(int)screen_z[k]);
    acc[k].save(std::string(fn)); }
  obs_::dump_eigenvalues(electrons, outdir+"/raw/observables/eigenvalues/eigenvalues.csv");

  const double inner_tau=abs_::inner_region_norm_twosided(electrons,2,z_in,wp_idx);
  const double total_tau=abs_::inner_region_norm(electrons,2,+1e12,wp_idx);
  const double eps=inner_tau/N0;
  double wall=std::chrono::duration<double>(std::chrono::steady_clock::now()-t_start).count();
  { std::ofstream s(outdir+"/run_summary.txt"); s.precision(12);
    s<<"run_completed = true\nrun_type = graphene-wp-cap\ntag = "<<tag<<"\nE_eV = "<<E_eV
     <<"\nk0 = "<<k0<<"\ncap = "<<(int)capon<<"\neta_Ha = "<<eta<<"\nL_cap = "<<Lcap
     <<"\nsigma = "<<sigma<<"\nz0 = "<<z0<<"\nz_in = "<<z_in<<"\ncx = "<<cx<<"\ncy = "<<cy
     <<"\nepsilon_survival = "<<eps<<"\nabsorbed_fraction = "<<(1.0-total_tau/N0)
     <<"\nN0 = "<<N0<<"\nwp_norm_tau = "<<total_tau<<"\nnorm_after_inject = "<<rep.norm_after
     <<"\nmax_overlap = "<<rep.max_overlap<<"\nn_electrons = "<<n_elec<<"\nwp_idx = "<<wp_idx
     <<"\ndt = "<<dt<<"\nN_STEPS = "<<N_STEPS<<"\nwrite_every = "<<WRITE_EVERY
     <<"\ncell_bohr = "<<cfg::LX_BOHR<<' '<<cfg::LY_BOHR<<' '<<cfg::LZ_BOHR
     <<"\npropagator = etrs\nwall_s = "<<wall<<"\nPROVISIONAL = inq-study Task#7\n"; }
  std::printf("  [%s] epsilon=%.5f absorbed=%.4f wall=%.1fs -> %s\n", tag.c_str(), eps, 1.0-total_tau/N0, wall, outdir.c_str());
  return 0;
}
