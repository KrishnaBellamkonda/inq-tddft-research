// ============================================================================
// run_classical_n162_L50_sv_sigma*/run.cpp — σ-convergence S(v) sweep.
//
// Full-observable classical-electron projectile (erf-smoothed Gaussian, m_e)
// through the cubic 50^3 Bohr periodic jellium bath (N=162, dx=0.40), launched
// +z on-axis from z=-20 Bohr. One build serves the whole (σ, v) sweep; σ via
// SV_PSEUDO, velocity via PROJ_V0, length via SV_N_STEPS, output via SV_OUT_ROOT.
//
// Observables = ADR-0006 jellium_classical minimum set + center-of-density,
// with cadences DECOUPLED so the heavy density VTI stays at ~6 frames while the
// scalar energy/track stays dense for Method-A stopping extraction:
//   - electron_track.csv          every step          (Method A/B path & v)
//   - observables.csv (E,J,μ,l2,COD) every SV_OBS_EVERY (dense, cheap)
//   - density_total/system VTI    every SV_VTI_EVERY  (≈6 frames; memory cap)
//   - state_energies.csv          every 5·SV_OBS_EVERY
//   - overlap_full                t=0 and t=end only
//   - gs eigenvalues/occupations + gs density VTI (once, from checkpoint)
//   - observables_manifest.json   (RunType::jellium_classical)
//
// Runtime env (one build serves the sweep):
//   SV_PSEUDO     projectile UPF path                 (default sigma0p15)
//   PROJ_V0       initial +z velocity (a.u.)          (default 1.0)
//   SV_N_STEPS    number of propagation steps         (default 700)
//   SV_VTI_EVERY  density-VTI cadence (6 frames)      (default N_STEPS/5)
//   SV_OBS_EVERY  scalar-observable cadence           (default max(1,N_STEPS/80))
//   SV_OUT_ROOT   output root dir                     (default results/run)
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/center_of_density.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/observables/minimum_observable_set.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>

#include "../shared/configs/sv_ladder_L50_sigma0p5.hpp"
#include "../shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::SV_Ladder_L50_sigma0p5;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }
static std::string iso_now(){ auto t=std::time(nullptr); auto tm=*std::localtime(&t);
    char b[64]; std::strftime(b,sizeof(b),"%Y-%m-%dT%H:%M:%S",&tm); return std::string(b); }

int main() {
    auto t0wall = std::chrono::steady_clock::now();

    const double      V0        = env_d("PROJ_V0", Cfg::PROJ_VEL_Z_DEFAULT);
    const int         N_STEPS   = env_i("SV_N_STEPS", 700);
    const int         VTI_EVERY = env_i("SV_VTI_EVERY", std::max(1, N_STEPS/5));
    const int         OBS_EVERY = env_i("SV_OBS_EVERY", std::max(1, N_STEPS/80));
    const std::string PSEUDO    = env_s("SV_PSEUDO", Cfg::PROJ_PSEUDO_PATH);
    const std::string OUT       = env_s("SV_OUT_ROOT", "results/run");
    const std::string RAW       = OUT + "/raw";
    const std::string OBSDIR    = RAW + "/observables";

    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p40";

    std::cout << "\n=== sv_sweep classical ===\n"
              << "  v0=" << V0 << "  N_STEPS=" << N_STEPS
              << "  dt=" << Cfg::DT_AU << "\n"
              << "  VTI_EVERY=" << VTI_EVERY << " (~" << (N_STEPS/VTI_EVERY + 1)
              << " frames)  OBS_EVERY=" << OBS_EVERY << "\n"
              << "  psp=" << PSEUDO << "\n  out=" << OUT
              << "\n  launch_z=" << Cfg::PROJ_LAUNCH_Z << " Bohr\n";

    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing\n"; return 2; }

    std::filesystem::create_directories(OBSDIR);
    std::filesystem::create_directories(OBSDIR + "/overlap_full");
    std::filesystem::create_directories(RAW + "/vti/density_total");
    std::filesystem::create_directories(RAW + "/vti/density_system");
    std::filesystem::create_directories(RAW + "/vti/density_gs_system");

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto sp_e = ionic::species(Cfg::PROJ_SPECIES_SYMBOL).pseudo_file(PSEUDO).mass(Cfg::PROJ_MASS_AMU);
    ions.insert(sp_e, {Cfg::PROJ_LAUNCH_X * 1.0_b, Cfg::PROJ_LAUNCH_Y * 1.0_b, Cfg::PROJ_LAUNCH_Z * 1.0_b});

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);
    jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OBSDIR + "/eigenvalues");

    ions.velocities()[0] = vector3<double>{0.0, 0.0, V0};
    std::cout << "  mass_au=" << ions.species(0).mass() << " (expect 1.0)  KE0="
              << ions.kinetic_energy() << " Ha\n";

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();

    // ----- manifest (ADR 0006) -----
    inqkit::observables::write_manifest(
        OUT + "/observables_manifest.json",
        inqkit::observables::RunType::jellium_classical, OBS_EVERY, N_STEPS);

    // stub run_summary
    {
        std::ofstream s(OUT + "/run_summary.txt");
        s << "run_name        = sv_sweep_classical\n"
          << "date_started    = " << iso_now() << "\n"
          << "projectile      = classical electron, erf Gaussian, m_e, ehrenfest\n"
          << "rs              = 5.69 (N=162, L=50, dx=0.40)\n"
          << "psp             = " << PSEUDO << "\n"
          << "v0_au           = " << V0 << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "n_steps         = " << N_STEPS << "\n"
          << "vti_every       = " << VTI_EVERY << "\n"
          << "obs_every       = " << OBS_EVERY << "\n"
          << "run_completed   = false\n";
    }

    inqkit::io::RealField3DLayout vti_layout{
        .field_name="density", .include_meta=false, .emit_raw=false,
        .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary};

    // t=0 GS density VTI
    auto gs_density = inqkit::fields::density::total(electrons);
    {
        inqkit::io::RealField3DWriter gs_sys(RAW + "/vti/density_gs_system", vti_layout, {.overwrite=true});
        gs_sys.write(gs_density, "density_gs_system");
    }

    // overlap (start + end only)
    const int n_ref = n_states - 1;
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(electrons, n_ref, OBSDIR + "/overlap_full");

    inqkit::io::RealField3DWriter total_wr (RAW + "/vti/density_total",  vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(RAW + "/vti/density_system", vti_layout, {.overwrite=true});
    system_wr.write(gs_density, 0.0, 0);
    total_wr.write (gs_density, 0.0, 0);

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x  = sel.dipole_y  = sel.dipole_z  = true;
    sel.density_l2 = true;
    sel.cod_x = sel.cod_y = sel.cod_z = true;
    inqkit::io::ObservablesWriter obs_writer(OBSDIR + "/observables.csv", sel);
    obs_writer.write_header();

    inqkit::observables::StateEnergyWriter state_energy_wr(OBSDIR + "/state_energies.csv", true);

    // density_delta used for the L2 metric ONLY (no VTI series — memory cap).
    inqkit::observables::DensityDelta density_delta(
        RAW + "/vti/_ddelta_unused", RAW + "/vti/_ddelta_unused_coarse",
        {.emit_raw_vti=false, .emit_coarse_vti=false, .compute_l2=true, .coarse_bin_bohr=3.0});
    density_delta.set_reference(gs_density);

    std::ofstream trk(OBSDIR + "/electron_track.csv");
    trk << std::setprecision(16) << "step,time_au,x,y,z,vx,vy,vz,fx,fy,fz\n";
    {
        auto const&p=ions.positions()[0]; auto const&v=ions.velocities()[0];
        trk<<0<<",0,"<<p[0]<<","<<p[1]<<","<<p[2]<<","<<v[0]<<","<<v[1]<<","<<v[2]<<",0,0,0\n";
    }
    overlap_full_obs.snapshot(electrons, 0.0, 0);

    // ---- callbacks ----
    // scalar observables (dense): energy/current/dipole + L2 + COD
    inqkit::RealTimeSession rt_obs(ions, electrons, OBS_EVERY);
    rt_obs.add([&](inqkit::StepContext const &ctx){
        auto f = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(f, ctx.time_au, ctx.step);  // l2 only
        auto cod = inqkit::observables::center_of_density(f);
        inqkit::StepContext o = ctx;
        o.density_l2 = l2;
        o.wp_center  = vector3<double>{cod.center_bohr.x, cod.center_bohr.y, cod.center_bohr.z};
        obs_writer.append(o);
    });

    // density VTI (≈6 frames)
    inqkit::RealTimeSession rt_vti(ions, electrons, VTI_EVERY);
    rt_vti.add([&](inqkit::StepContext const &ctx){
        auto f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(f, ctx.time_au, ctx.step);
        total_wr.write (f, ctx.time_au, ctx.step);
    });

    // electron track (every step)
    inqkit::RealTimeSession rt_track(ions, electrons, 1);
    rt_track.add([&](inqkit::StepContext const &ctx){
        auto const&p=ctx.ions->positions()[0]; auto const&v=ctx.ions->velocities()[0];
        trk<<ctx.step<<","<<ctx.time_au<<","<<p[0]<<","<<p[1]<<","<<p[2]<<","
           <<v[0]<<","<<v[1]<<","<<v[2]<<",0,0,0\n";
    });

    real_time::propagate(
        ions, electrons,
        [&](auto const &data){
            rt_obs.step(data);
            rt_vti.step(data);
            rt_track.step(data);
            if (data.iter() % (5*OBS_EVERY) == 0) state_energy_wr.snapshot(data);
        },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .ehrenfest()
            .observables_current()
            .observables_dipole());

    overlap_full_obs.snapshot(electrons, Cfg::DT_AU * N_STEPS, N_STEPS);  // end

    trk.flush();
    auto t1wall = std::chrono::steady_clock::now();
    double wall = std::chrono::duration<double>(t1wall - t0wall).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt", std::ios::app);
        s << std::setprecision(16)
          << "date_finished   = " << iso_now() << "\n"
          << "wall_time_s     = " << wall << "\n"
          << "final_z         = " << ions.positions()[0][2] << "\n"
          << "final_vz        = " << ions.velocities()[0][2] << "\n"
          << "run_completed   = true\n";
    }
    std::cout << "  done. final_z=" << ions.positions()[0][2]
              << " final_vz=" << ions.velocities()[0][2]
              << " wall=" << wall << "s\n";
    return 0;
}
