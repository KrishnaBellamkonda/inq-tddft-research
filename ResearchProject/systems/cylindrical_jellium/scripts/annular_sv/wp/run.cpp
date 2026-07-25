// ============================================================================
// cylindrical_jellium / scripts/annular_sv / wp / run.cpp
//
// Phase-4 QUANTUM rung: an electron WAVEPACKET (sigma_WP, drift k0) gliding
// on-axis down the bore of the PERIODIC annular jellium tube — the quantum
// counterpart of the classical Gaussian-electron glide. The WP is injected as an
// extra electron state (net −1 projectile, INQ's G=0 background compensates), the
// annulus background re-applied. Compared against the matched classical ghost
// (the v=PROJ_V0 production run at the same r_s).
//
// S from ΔE_system(t) regression vs the WP centroid path (wp_real_space_stats);
// quantum spreading from wp_momentum_stats / momentum_distribution.
//
// Env (orchestrator/runner supplies all; defaults = r_s=6, v=0.30, sigma_WP=0.5):
//   geometry: CJ_LXY(40) CJ_LZ(48) CJ_RIN(5) CJ_ROUT(13) CJ_N(24) CJ_EDGE_W(1.0)
//             CJ_SPACING(0.5)
//   run:      CJ_GS_DIR(REQUIRED) PROJ_V0(0.30) CJ_WP_SIGMA(0.5)
//             SV_N_STEPS(2000) SV_WRITE_EVERY(=N_STEPS/300) SV_OUT_SUBDIR(wp_run)
//             CJ_LAUNCH_Z(=-(L_z/2)+1)
// Build against INQ_SOURCE=inq-study; runtime shares from inq/install/share.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
    auto t0 = std::chrono::steady_clock::now();

    const double LXY = env_d("CJ_LXY", 40.0), LZ = env_d("CJ_LZ", 48.0);
    const double RIN = env_d("CJ_RIN", 5.0), ROUT = env_d("CJ_ROUT", 13.0);
    const int    N = env_i("CJ_N", 24);
    const double EDGE_W = env_d("CJ_EDGE_W", 1.0);
    const double SPACING = env_d("CJ_SPACING", 0.5);
    const double V0 = env_d("PROJ_V0", 0.30);            // WP drift momentum k0 (m_e=1)
    const double WP_SIGMA = env_d("CJ_WP_SIGMA", 0.5);   // envelope width sigma_WP
    const int    N_STEPS = env_i("SV_N_STEPS", 2000);
    const int    WRITE_EVERY = env_i("SV_WRITE_EVERY", std::max(1, (int)std::lround(N_STEPS/300.0)));
    const double DT_AU = 0.020;
    const double LAUNCH_Z = env_d("CJ_LAUNCH_Z", -(LZ/2.0) + 1.0);
    const std::string OUT = "results/" + env_s("SV_OUT_SUBDIR", "wp_run");
    const std::string GS_DIR = env_s("CJ_GS_DIR", "");
    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2;
    }
    const double V_ann = M_PI * (ROUT*ROUT - RIN*RIN) * LZ;
    const double N0 = double(N) / V_ann;

    std::cout << "\n=== annular_sv WP OUT=" << OUT << " ===\n"
              << "  tube R_in=" << RIN << " R_out=" << ROUT << " L_z=" << LZ << " (periodic)\n"
              << "  N=" << N << " n0=" << N0 << " sigma_WP=" << WP_SIGMA << " k0=" << V0
              << " launch_z=" << LAUNCH_Z << "\n"
              << "  N_STEPS=" << N_STEPS << " dt=" << DT_AU << " write_every=" << WRITE_EVERY << "\n";

    auto cell = systems::cell::orthorhombic(LXY*1.0_b, LXY*1.0_b, LZ*1.0_b).periodic();
    auto ions = systems::ions(cell);   // no classical projectile ion — the WP is the projectile
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING*1.0_b)
            .extra_electrons(N)
            .extra_states(20)
            .temperature(0.00862*1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);
    const int n_states = electrons.states().num_states();

    // ----- WP injection (the projectile electron) -----------------------
    auto wp = inqkit::WavePacket{}
                  .center(0.0, 0.0, LAUNCH_Z).sigma(WP_SIGMA)
                  .k0(0.0, 0.0, V0).orthogonalise_against_occupied(electrons);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << report.norm_after
              << " max_overlap=" << report.max_overlap << "\n";

    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::annulus;
    bg.n0 = N0; bg.half_width = ROUT; bg.inner_radius = RIN; bg.slab_axis = 2;
    bg.center = {0.0, 0.0, 0.0}; bg.edge_width = EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    // ----- output skeleton + writers ------------------------------------
    for (auto sub : {"density_system","density_delta","density_delta_coarse","density_wp"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
    std::filesystem::create_directories(OUT + "/raw/observables");

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter wp_density_wr(OUT + "/raw/vti/density_wp", vti_layout, {.overwrite=true});
    { auto s0 = inqkit::fields::density::total(electrons); system_wr.write(s0, 0.0, 0); }

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
    obs_writer.write_header();
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});
    inqkit::observables::MomentumDistribution momentum_dist(
        OUT + "/raw/observables/momentum_distribution.csv", wp_idx, LZ,
        {.n_bins=64, .k_max_bohr_inv=0.0, .write_every=WRITE_EVERY});
    inqkit::observables::WPMomentumStats wp_momentum_stats(
        OUT + "/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every=WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every=WRITE_EVERY});
    std::ofstream nlog(OUT + "/raw/observables/electron_number.csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext c = ctx; c.density_l2 = l2; obs_writer.append(c);
        wp_density_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, ctx.step);
    });

    auto step_fn = [&](auto const& data) {
        rt_obs.step(data);
        momentum_dist.maybe_accumulate(data);
        wp_momentum_stats.maybe_accumulate(data);
        wp_real_space_stats.maybe_accumulate(data);
        if (data.root()) nlog << data.iter() << "," << (data.iter()*DT_AU) << "," << data.num_electrons() << "\n";
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU*1.0_atomictime)
                       .observables_current().observables_dipole();
    real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, bg_pert);

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = annular_sv/wp/" << env_s("SV_OUT_SUBDIR","wp_run") << "\nengine = inq-study\n"
          << "projectile = electron WAVEPACKET sigma_WP " << WP_SIGMA << " k0 " << V0 << " (quantum)\n"
          << "geometry = annular_tube  R_in=" << RIN << " R_out=" << ROUT << " L_z=" << LZ << " (periodic)\n"
          << "cell_bohr = " << LXY << " x " << LXY << " x " << LZ << "  spacing = " << SPACING << "\n"
          << "n_electrons_bath = " << N << "  n0 = " << N0 << "  n_states = " << n_states
          << "  wp_state_index = " << wp_idx << "\n"
          << "wp_norm_after = " << report.norm_after << "  wp_max_overlap = " << report.max_overlap << "\n"
          << "launch_z = " << LAUNCH_Z << "  k0 = " << V0 << "  drift_ke_ha = " << (0.5*V0*V0) << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
          << "gs_dir = " << GS_DIR << "\nwall_time_s = " << wall << "\nrun_completed = true\n";
    }
    std::cout << "  done. wall=" << wall << "s\n";
    return 0;
}
