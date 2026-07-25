// ============================================================================
// localised_jellium / scripts/nazarov_gross / wp / run.cpp
//
// Fixed-VELOCITY mass-sweep wavepacket run for the Nazarov-Gross campaign
// (docs/campaigns/nazarov_gross_comparison/nazarov-gross-comparison.md).
// Clone of scripts/fullsuite_wp/run.cpp (the p3 baseline, CAP OFF) with:
//   * NG_MASS     projectile mass in m_e (default 1.0) -> per-state inverse
//                 mass on the WP slot via the inq-study fork (the SAME call
//                 validated in muon_mass_fork: effmass_pair/quantum/run.cpp:81)
//   * NG_V        projectile VELOCITY a.u. (default Cfg::WP_K0 = 2.711, i.e.
//                 the p3 electron velocity). k0 = NG_MASS * NG_V per the
//                 fixed-velocity sweep convention; E = m v^2 / 2.
//   * NG_SPACING  grid spacing Bohr (default 0.35)
//   * NG_GS_DIR   GS checkpoint (must match NG_SPACING)
//   * NG_OUT / NG_N_STEPS(880) / NG_LAUNCH_Z(-23) / NG_DT(0.02) /
//     NG_WRITE_EVERY(10) / NG_WF_EVERY(10)
// No CAP (p3 is the closed-box baseline). Bath = mass-1 electrons.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include "../../../shared/configs/slab_n234_L50.hpp"
#include "../../../../jellium/shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN234L50;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

// Checkpoint state parsers (cloned from sigma1_massonly/wp/run.cpp:55-73).
static int read_last_step(const std::string& path) {
    std::ifstream f(path); if (!f) return -1;
    std::string line; int last = -1;
    while (std::getline(f, line)) {
        auto p = line.find("last_step=");
        if (p == 0) last = std::atoi(line.c_str() + 10);
    }
    return last;
}
static int read_kv_int(const std::string& path, const std::string& key, int dflt) {
    std::ifstream f(path); if (!f) return dflt;
    std::string line;
    while (std::getline(f, line)) {
        auto p = line.find(key + "=");
        if (p == 0) return std::atoi(line.c_str() + key.size() + 1);
    }
    return dflt;
}

int main() {
    auto t0 = std::chrono::steady_clock::now();

    const double MASS        = env_d("NG_MASS", 1.0);
    const double VEL         = env_d("NG_V", Cfg::WP_K0);       // fixed-velocity invariant
    const double KZ          = MASS * VEL;                       // k0 = m*v
    const double INV_MASS    = 1.0 / MASS;
    const double SPACING     = env_d("NG_SPACING", 0.35);
    const std::string OUT    = "results/" + env_s("NG_OUT", "null_m1");
    const double DT_AU       = env_d("NG_DT", 0.02);
    const int    N_STEPS     = env_i("NG_N_STEPS", 880);
    const int    WRITE_EVERY = env_i("NG_WRITE_EVERY", 10);
    const int    WF_EVERY    = env_i("NG_WF_EVERY", 10);
    const double LAUNCH_Z    = env_d("NG_LAUNCH_Z", -23.0);
    const int    CKPT_EVERY  = env_i("NG_CKPT_EVERY", 200);     // interior RT checkpoints
    const bool   RESUME      = env_i("NG_RESUME", 0) != 0;
    const double E_EV        = 0.5 * MASS * VEL * VEL * 27.211386245988;
    const std::string RT_CKPT_DIR  = OUT + "/rt_ckpt";
    const std::string RT_STATE_TXT = RT_CKPT_DIR + "/rt_state.txt";

    const std::string GS_DIR = env_s("NG_GS_DIR",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "shared_gs/slab_n234_L50_h0p35");
    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

    // ----- resume bookkeeping (pattern: sigma1_massonly, bit-faithful ETRS
    //       restart: jellium has no nuclei, background is static) -------------
    int START = 0;
    if (RESUME) {
        if (!std::filesystem::exists(RT_CKPT_DIR)) { std::cerr << "FATAL: resume but no RT ckpt: " << RT_CKPT_DIR << "\n"; return 2; }
        START = read_last_step(RT_STATE_TXT);
        if (START < 0) { std::cerr << "FATAL: resume but unreadable " << RT_STATE_TXT << "\n"; return 2; }
        if (START >= N_STEPS) { std::cout << "Already at/after target (" << START << ">=" << N_STEPS << "); nothing to do.\n"; return 0; }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << "\n=== nazarov_gross wp (out=" << OUT << ") ===\n"
              << "  mass=" << MASS << " v=" << VEL << " k0=" << KZ << " E=" << E_EV << " eV\n"
              << "  h=" << SPACING << " N_STEPS=" << N_STEPS << " dt=" << DT_AU
              << " launch_z=" << LAUNCH_Z << "\n";

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    int wp_idx_pre = -1;
    if (RESUME) {
        electrons.load(RT_CKPT_DIR);                     // already holds the propagated WP
        wp_idx_pre = read_kv_int(RT_STATE_TXT, "wp_idx", electrons.states().num_states() - 1);
        std::cout << "  [RESUME] loaded RT ckpt at step " << START << " from " << RT_CKPT_DIR << "\n";
    } else {
        electrons.load(GS_DIR);
        std::cout << "  Loaded GS (inq-study) from " << GS_DIR << "\n";
        jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OUT + "/raw/observables/eigenvalues");
    }
    const int n_states = electrons.states().num_states();

    // ----- output skeleton ----------------------------------------------
    for (auto sub : {"density_total","density_system","density_gs_system","density_wp",
                     "wavefunction_wp","density_delta","density_delta_coarse"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
    std::filesystem::create_directories(OUT + "/raw/observables/overlap");
    std::filesystem::create_directories(OUT + "/raw/observables/overlap_full");

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};

    if (!RESUME) { inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
      gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system"); }

    // ----- WP injection + per-state mass (the inq-study fork) ------------
    int wp_idx = wp_idx_pre;
    double wp_norm_after = 1.0;
    if (!RESUME) {
        auto wp = inqkit::WavePacket{}
                      .center(0.0, 0.0, LAUNCH_Z).sigma(Cfg::WP_SIGMA_BOHR)
                      .k0(0.0, 0.0, KZ).orthogonalise_against_occupied(electrons);
        auto report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index;
        wp_norm_after = report.norm_after;
        std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << report.norm_after
                  << " max_overlap=" << report.max_overlap << "\n";
    }
    electrons.inverse_mass()[0][wp_idx] = INV_MASS;      // the mass fork (re-applied on resume)
    std::cout << "  inverse_mass[" << wp_idx << "]=" << INV_MASS << "\n";

    // ----- background well ------------------------------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::slab;
    bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
    bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    // ----- density writers ----------------------------------------------
    inqkit::io::RealField3DWriter total_wr (OUT + "/raw/vti/density_total",  vti_layout, {.overwrite = !RESUME});
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite = !RESUME});
    if (!RESUME) { auto s0 = inqkit::fields::density::total(electrons); total_wr.write(s0,0.0,0); system_wr.write(s0,0.0,0); }
    inqkit::io::RealField3DWriter wp_density_wr(OUT + "/raw/vti/density_wp", vti_layout, {.overwrite = !RESUME});
    inqkit::io::ComplexField3DWriter wp_wf_wr(
        OUT + "/raw/vti/wavefunction_wp",
        {.field_name="wavefunction", .include_meta=false, .emit_raw=false,
         .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite = !RESUME});

    // ----- scalar/observable writers ------------------------------------
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.density_l2 = true;
    // Segment-suffixed CSVs on resume (post-processing concatenates segments;
    // density_delta L2 baseline resets per segment — noted in run_summary).
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables" + SEG + ".csv", sel);
    obs_writer.write_header();
    inqkit::observables::StateEnergyWriter state_energy_wr(OUT + "/raw/observables/state_energies" + SEG + ".csv", true);
    inqkit::observables::OccupationsWriter occupations_wr(OUT + "/raw/observables/occupations_vs_time" + SEG + ".csv");
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});

    // ----- WP-specific observables --------------------------------------
    inqkit::observables::OrbitalOverlapMatrix overlap_obs(electrons, wp_idx, OUT + "/raw/observables/overlap");
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(electrons, n_states - 1, OUT + "/raw/observables/overlap_full");
    if (!RESUME) {
        overlap_full_obs.snapshot(electrons, 0.0, 0);
        overlap_obs.snapshot_wp_only(electrons, 0.0, 0);
    }
    inqkit::observables::MomentumDistribution momentum_dist(
        OUT + "/raw/observables/momentum_distribution" + SEG + ".csv", wp_idx, Cfg::L_BOHR,
        {.n_bins=64, .k_max_bohr_inv=0.0, .write_every=WRITE_EVERY});
    inqkit::observables::WPMomentumStats wp_momentum_stats(
        OUT + "/raw/observables/wp_momentum_stats" + SEG + ".csv", wp_idx, {.write_every=WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        OUT + "/raw/observables/wp_real_space_stats" + SEG + ".csv", wp_idx, {.write_every=WRITE_EVERY});

    std::ofstream nlog(OUT + "/raw/observables/electron_number" + SEG + ".csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";
    std::filesystem::create_directories(RT_CKPT_DIR);

    // ----- real-time density/observable session -------------------------
    inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        total_wr.write (sys_f, ctx.time_au, ctx.step);
        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext c = ctx; c.density_l2 = l2; obs_writer.append(c);
        if (ctx.step % WF_EVERY == 0) {
            wp_density_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, ctx.step);
            char nm[64]; std::snprintf(nm, sizeof(nm), "wavefunction_t%06d", ctx.step);
            wp_wf_wr.write(inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx), std::string(nm));
        }
        if (ctx.step % 10 == 0) overlap_obs.snapshot_wp_only(*ctx.electrons, ctx.time_au, ctx.step);
    });

    auto step_fn = [&](auto const& data) {
        rt_obs.step(data);
        const int it = data.iter();
        if (it % (5 * WRITE_EVERY) == 0) { state_energy_wr.snapshot(data); occupations_wr.snapshot(data); }
        momentum_dist.maybe_accumulate(data);
        wp_momentum_stats.maybe_accumulate(data);
        wp_real_space_stats.maybe_accumulate(data);
        if (data.root()) nlog << it << "," << (it*DT_AU) << "," << data.num_electrons() << "\n";
        if (it > START && it % CKPT_EVERY == 0 && it < N_STEPS) {
            electrons.save(RT_CKPT_DIR);      // collective; overwrites the checkpoint
            if (data.root()) {
                std::ofstream st(RT_STATE_TXT, std::ios::trunc);
                st << "last_step=" << it << "\ntime_au=" << (it*DT_AU) << "\nwp_idx=" << wp_idx
                   << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS
                   << "\ndt_au=" << DT_AU << "\n";
                std::cout << "  [ckpt] saved at step " << it << " (t=" << it*DT_AU << ")\n";
            }
        }
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime)
                       .observables_current().observables_dipole();
    real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, bg_pert, START);

    // final checkpoint (so a follow-on resume is a clean no-op)
    electrons.save(RT_CKPT_DIR);
    if (electrons.root()) {
        std::ofstream st(RT_STATE_TXT, std::ios::trunc);
        st << "last_step=" << N_STEPS << "\ntime_au=" << N_STEPS*DT_AU << "\nwp_idx=" << wp_idx
           << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS << "\ndt_au=" << DT_AU << "\n";
    }

    overlap_full_obs.snapshot(electrons, DT_AU * N_STEPS, N_STEPS);

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        // run_summary.txt is ALWAYS the final-state file (orchestrator keys its
        // done() on it); segments record start_step + segment CSV suffix inside.
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/nazarov_gross/" << env_s("NG_OUT","null_m1") << "\n"
          << "engine = inq-study (mass fork)\n"
          << "projectile = wavepacket sigma " << Cfg::WP_SIGMA_BOHR << " mass " << MASS
          << " velocity " << VEL << " k0 " << KZ << " E " << E_EV << " eV\n"
          << "inverse_mass = " << INV_MASS << "\n"
          << "cap = off\n"
          << "cell_bohr = " << Cfg::L_BOHR << "^3  spacing = " << SPACING << "\n"
          << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "  wp_state_index = " << wp_idx << "\n"
          << "wp_norm_after = " << wp_norm_after << "  launch_z = " << LAUNCH_Z << "\n"
          << "gs_dir = " << GS_DIR << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
          << "ckpt_every = " << CKPT_EVERY << "  rt_ckpt_dir = " << RT_CKPT_DIR << "\n"
          << "resume = " << (RESUME?"true":"false") << "  start_step = " << START
          << "  segment_suffix = " << (SEG.empty()?"(none)":SEG) << "\n"
          << "wall_time_s = " << wall << " (this segment)\nrun_completed = true\n";
    }
    std::cout << "  done. wall=" << wall << "s\n";
    return 0;
}
