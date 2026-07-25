// ============================================================================
// localised_jellium / scripts/wp_cap_energy_plateau / wp / run.cpp
//
// WP projectile through the localised-jellium slab (N=102, 25×25×140, w=0.5),
// run WITHOUT (WP_CAP_ETA=0) and WITH (WP_CAP_ETA=-0.7) a TWO-SIDED sin² CAP,
// for the energy-plateau diagnostic
// (docs/campaigns/localised_jellium/wp_cap_energy_plateau.md):
//   * no-CAP  = closed periodic box → energy_total CONSERVED, plateau = all
//               deposited energy retained.
//   * CAP     = escaping flux drained at ±60 → energy_total DROPS; the plateau
//               gap between the two runs = energy radiated to the boundaries.
//
// Records ALL KS energy components each step (energies.csv), the momentum
// distribution EVERY step, the projectile wavefunction every 10 steps, and
// density frames (for the density GIF). Checkpointed (every 200) + resumable.
// Engine: inq-study (the CAP's imaginary potential is functional only there;
// stock inq will not compile a real-time absorbing run).
//
// Env: WP_CAP_ETA(0=no CAP; -0.7=CAP) WP_OUT WP_GS_DIR WP_N_STEPS(5000)
//      WP_DT(0.02) WP_WF_EVERY(10) WP_DENS_EVERY(20) WP_CKPT_EVERY(200)
//      WP_RESUME(0) WP_COMPILE_PROBE=1
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include "../../../shared/configs/slab_n102_L25x25x140_w0p5.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace obs_ = inqkit::observables;
using Cfg = localised_jellium::config::SlabN102_L25x25x140_w0p5;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static int read_last_step(const std::string& path) {
    std::ifstream f(path); if (!f) return -1;
    std::string line; int last = -1;
    while (std::getline(f, line)) { auto p = line.find("last_step="); if (p == 0) last = std::atoi(line.c_str() + 10); }
    return last;
}
static int read_kv_int(const std::string& path, const std::string& key, int dflt) {
    std::ifstream f(path); if (!f) return dflt;
    std::string line;
    while (std::getline(f, line)) { auto p = line.find(key + "="); if (p == 0) return std::atoi(line.c_str() + key.size() + 1); }
    return dflt;
}

int main() {
    if (std::getenv("WP_COMPILE_PROBE")) { std::cout << "compile probe ok\n"; return 2; }
    auto t0 = std::chrono::steady_clock::now();

    const double CAP_ETA     = env_d("WP_CAP_ETA", 0.0);
    const bool   CAP_ON      = (CAP_ETA != 0.0);
    const std::string OUT    = "results/" + env_s("WP_OUT", (CAP_ON? "cap":"nocap"));
    const double DT_AU       = env_d("WP_DT", 0.02);
    const int    N_STEPS     = env_i("WP_N_STEPS", 5000);        // 100 a.u.
    const int    WF_EVERY    = env_i("WP_WF_EVERY", 10);
    const int    DENS_EVERY  = env_i("WP_DENS_EVERY", 20);
    const int    CKPT_EVERY  = env_i("WP_CKPT_EVERY", 200);
    const bool   RESUME      = env_i("WP_RESUME", 0) != 0;
    const std::string GS_DIR = env_s("WP_GS_DIR", "");
    if (GS_DIR.empty()) { std::cerr << "FATAL: WP_GS_DIR not set\n"; return 2; }
    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

    const std::string RT_CKPT_DIR  = OUT + "/rt_ckpt";
    const std::string RT_STATE_TXT = RT_CKPT_DIR + "/rt_state.txt";

    int START = 0;
    if (RESUME) {
        if (!std::filesystem::exists(RT_CKPT_DIR)) { std::cerr << "FATAL: resume but no RT ckpt\n"; return 2; }
        START = read_last_step(RT_STATE_TXT);
        if (START < 0) { std::cerr << "FATAL: unreadable " << RT_STATE_TXT << "\n"; return 2; }
        if (START >= N_STEPS) { std::cout << "Already at target (" << START << ">=" << N_STEPS << ").\n"; return 0; }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(6)
              << "\n=== wp_cap_energy_plateau WP (out=" << OUT << ") ===\n"
              << "  CAP " << (CAP_ON?"ON":"OFF") << " eta=" << CAP_ETA
              << "  sigma=" << Cfg::WP_SIGMA_BOHR << " E=" << Cfg::WP_EKIN_EV << " eV k0=" << Cfg::WP_K0 << "\n"
              << "  launch_z=" << Cfg::WP_CZ_BOHR << " N_STEPS=" << N_STEPS << " dt=" << DT_AU
              << " (start=" << START << ")\n";

    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR*1.0_b, Cfg::LY_BOHR*1.0_b, Cfg::LZ_BOHR*1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR*1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV*1.0_eV),
        input::kpoints::gamma());

    int wp_idx = -1;
    double wp_norm_after = 1.0;
    if (RESUME) {
        electrons.load(RT_CKPT_DIR);
        wp_idx = read_kv_int(RT_STATE_TXT, "wp_idx", electrons.states().num_states() - 1);
        std::cout << "  [RESUME] loaded RT ckpt at step " << START << "\n";
    } else {
        electrons.load(GS_DIR);
        auto wp = inqkit::WavePacket{}
                      .center(0.0, 0.0, Cfg::WP_CZ_BOHR).sigma(Cfg::WP_SIGMA_BOHR)
                      .k0(0.0, 0.0, Cfg::WP_KZ).orthogonalise_against_occupied(electrons);
        auto report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index;
        wp_norm_after = report.norm_after;
        std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << report.norm_after
                  << " max_overlap=" << report.max_overlap << "\n";
    }
    const int n_states = electrons.states().num_states();

    // ----- background well + optional two-sided CAP -------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::slab;
    bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
    bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);
    perturbations::absorbing cap_lo(CAP_ETA*1.0_Ha, -Cfg::CAP_MID_FRAC, Cfg::CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(CAP_ETA*1.0_Ha, +Cfg::CAP_MID_FRAC, Cfg::CAP_WIDTH_FRAC);
    auto pert = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

    // ----- output skeleton --------------------------------------------------
    for (auto sub : {"density_total","density_wp","wavefunction_wp","density_delta","density_delta_coarse"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
    std::filesystem::create_directories(OUT + "/raw/observables");
    std::filesystem::create_directories(RT_CKPT_DIR);

    inqkit::io::RealField3DLayout vti_layout{
        .field_name="density", .include_meta=false, .emit_raw=false,
        .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter total_wr(OUT + "/raw/vti/density_total", vti_layout, {.overwrite=!RESUME});
    inqkit::io::RealField3DWriter wp_density_wr(OUT + "/raw/vti/density_wp", vti_layout, {.overwrite=!RESUME});
    inqkit::io::ComplexField3DWriter wp_wf_wr(
        OUT + "/raw/vti/wavefunction_wp",
        {.field_name="wavefunction", .include_meta=false, .emit_raw=false,
         .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary}, {.overwrite=!RESUME});
    if (!RESUME) {
        auto s0 = inqkit::fields::density::total(electrons);
        total_wr.write(s0, 0.0, 0);
        wp_density_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    }
    obs_::DensityDelta density_delta(OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});

    // ----- ALL KS energy components, every step -----------------------------
    std::ofstream en(OUT + "/raw/observables/energies" + SEG + ".csv");
    en << std::setprecision(12)
       << "step,time_au,total,kinetic,hartree,external,non_local,xc,exact_exchange,ion,ion_kinetic,density_l2,N_total\n";

    // ----- momentum EVERY step; WP real-space stats at WF cadence ------------
    obs_::MomentumDistribution mom_dist(OUT + "/raw/observables/momentum_distribution" + SEG + ".csv",
        wp_idx, Cfg::LZ_BOHR, {.n_bins=64, .k_max_bohr_inv=0.0, .write_every=1});
    obs_::WPMomentumStats  wp_mom(OUT + "/raw/observables/wp_momentum_stats" + SEG + ".csv", wp_idx, {.write_every=1});
    obs_::WPRealSpaceStats wp_rs (OUT + "/raw/observables/wp_real_space_stats" + SEG + ".csv", wp_idx, {.write_every=WF_EVERY});

    auto step_fn = [&](auto const& data) {
        const int it = data.iter();
        auto e = data.energy();
        auto sys_f = inqkit::fields::density::total(electrons);
        const double l2 = density_delta.snapshot(sys_f, it*DT_AU, it);
        if (data.root()) {
            en << it << ',' << (it*DT_AU) << ',' << e.total() << ',' << e.kinetic() << ',' << e.hartree()
               << ',' << e.external() << ',' << e.non_local() << ',' << e.xc() << ',' << e.exact_exchange()
               << ',' << e.ion() << ',' << e.ion_kinetic() << ',' << l2 << ',' << data.num_electrons() << '\n';
        }
        mom_dist.maybe_accumulate(data);
        wp_mom.maybe_accumulate(data);
        wp_rs.maybe_accumulate(data);
        if (it % DENS_EVERY == 0) {
            total_wr.write(sys_f, it*DT_AU, it);
            wp_density_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), it*DT_AU, it);
        }
        if (it % WF_EVERY == 0) {
            char nm[64]; std::snprintf(nm, sizeof(nm), "wavefunction_t%06d", it);
            wp_wf_wr.write(inqkit::fields::orbital::wavefunction(electrons, wp_idx), std::string(nm));
        }
        if (it > START && it % CKPT_EVERY == 0 && it < N_STEPS) {
            electrons.save(RT_CKPT_DIR);
            if (data.root()) {
                std::ofstream st(RT_STATE_TXT, std::ios::trunc);
                st << "last_step=" << it << "\ntime_au=" << (it*DT_AU) << "\nwp_idx=" << wp_idx
                   << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS << "\ndt_au=" << DT_AU << "\n";
                std::cout << "  [ckpt] step " << it << " (t=" << it*DT_AU << ")\n";
            }
        }
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU*1.0_atomictime)
                       .observables_current().observables_dipole();
    real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, pert, START);
    en.close();

    // final checkpoint (resume becomes a clean no-op)
    electrons.save(RT_CKPT_DIR);
    if (electrons.root()) {
        std::ofstream st(RT_STATE_TXT, std::ios::trunc);
        st << "last_step=" << N_STEPS << "\ntime_au=" << N_STEPS*DT_AU << "\nwp_idx=" << wp_idx
           << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS << "\ndt_au=" << DT_AU << "\n";
    }

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/wp_cap_energy_plateau/" << env_s("WP_OUT",(CAP_ON?"cap":"nocap")) << "\n"
          << "engine = inq-study\ntheory = lda\n"
          << "cap = " << (CAP_ON?"on":"off") << "  cap_eta_ha = " << CAP_ETA
          << "  cap_L_bohr = " << Cfg::CAP_L_BOHR << "  cap_mid_frac = +/-" << Cfg::CAP_MID_FRAC << "\n"
          << "projectile = wavepacket sigma " << Cfg::WP_SIGMA_BOHR << " E " << Cfg::WP_EKIN_EV
          << " eV k0 " << Cfg::WP_K0 << " mass 1\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR
          << "  spacing = " << Cfg::SPACING_BOHR << "\n"
          << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " edge_width " << Cfg::EDGE_WIDTH_BOHR << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "  wp_state_index = " << wp_idx << "\n"
          << "wp_norm_after = " << wp_norm_after << "  launch_z = " << Cfg::WP_CZ_BOHR << "\n"
          << "gs_dir = " << GS_DIR << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  wf_every = " << WF_EVERY
          << "  dens_every = " << DENS_EVERY << "  mom_every = 1\n"
          << "resume = " << (RESUME?"true":"false") << "  start_step = " << START
          << "  segment_suffix = " << (SEG.empty()?"(none)":SEG) << "\n"
          << "wall_time_s = " << wall << " (this segment)\nrun_completed = true\n";
    }
    std::cout << "  done. wall=" << wall << "s -> " << OUT << "\n";
    return 0;
}
