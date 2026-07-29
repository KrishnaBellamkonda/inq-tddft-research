// ============================================================================
// run_wp_n162_L50_E100_sigma8_wf/run.cpp — σ=8 WP at E=100 eV, L=50 bath.
// 2026-05-31 σ-sweep rerun: adds density_wp + wavefunction_wp saving (every
// density frame) AND bath-only system density via total_excluding_orbital, so
// the saved density_system is the TARGET RESPONSE with the WP removed (162 e),
// directly comparable to the classical run and usable for momentum/loss-fn.
//
// Cfg: jellium::config::Electron_Proj_E100_L50_sigma8_WP_dx0p40 (WF_WRITE_EVERY
// = WRITE_EVERY). GS: checkpoints/gs_L50_cubic_N162_dx0p40.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/center_of_density.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/shells.hpp>

#include "../shared/configs/electron_proj_E100_L50_cubic.hpp"
#include "../shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::Electron_Proj_E100_L50_sigma8_WP_dx0p40;

static std::string iso_now() {
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

int main() {
    auto t_wallclock_start = std::chrono::steady_clock::now();

    const std::string RUN_NAME = "run_wp_n162_L50_E100_sigma8_wf";
    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p40";

    std::cout << "\n=== " << RUN_NAME << " ===\n"
              << "  cell = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
              << "  WP: sigma=" << Cfg::WP_SIGMA_BOHR << "  k0=" << Cfg::WP_K0
              << "  E=" << Cfg::WP_EKIN_EV << " eV  launch z=" << Cfg::WP_CZ_BOHR << "\n"
              << "  dt=" << Cfg::DT_AU << " N_steps=" << Cfg::N_STEPS
              << " WRITE_EVERY=" << Cfg::WRITE_EVERY
              << " WF_WRITE_EVERY=" << Cfg::WF_WRITE_EVERY << "\n";

    if (!std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: GS checkpoint missing.\n"; return 2;
    }

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);
    jellium::eigenvalues::copy_from_checkpoint(
        GS_DIR, "results/raw/observables/eigenvalues");

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    const int n_occupied  = n_electrons / 2;
    std::cout << "  num_states=" << n_states << " num_electrons=" << n_electrons
              << " n_occupied=" << n_occupied << "\n";

    std::filesystem::create_directories("results/raw/observables");
    std::filesystem::create_directories("results/raw/observables/overlap");
    std::filesystem::create_directories("results/raw/observables/overlap_full");
    std::filesystem::create_directories("results/raw/vti/density_total");
    std::filesystem::create_directories("results/raw/vti/density_system");
    std::filesystem::create_directories("results/raw/vti/density_delta");
    std::filesystem::create_directories("results/raw/vti/density_delta_coarse");
    std::filesystem::create_directories("results/raw/vti/density_gs_system");
    std::filesystem::create_directories("results/raw/vti/density_wp");
    std::filesystem::create_directories("results/raw/vti/wavefunction_wp");

    {
        std::ofstream s("results/run_summary.txt");
        s << "RUN SUMMARY (stub - written at start)\n=====================================\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "date_started    = " << iso_now() << "\n"
          << "cell_bohr       = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "n_electrons     = " << n_electrons << "\n"
          << "n_states        = " << n_states << "\n"
          << "wp_ekin_ev      = " << Cfg::WP_EKIN_EV << "\n"
          << "wp_sigma_bohr   = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "n_steps         = " << Cfg::N_STEPS << "\n"
          << "run_completed   = false\n";
    }

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true,
        .vti_format = inqkit::io::VTIWriteOptions::Format::binary};

    {
        inqkit::io::RealField3DWriter gs_sys(
            "results/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
        gs_sys.write(inqkit::fields::density::total(electrons), "density_gs_system");
    }

    auto wp = inqkit::WavePacket{}
                  .center(Cfg::WP_CX_BOHR, Cfg::WP_CY_BOHR, Cfg::WP_CZ_BOHR)
                  .sigma(Cfg::WP_SIGMA_BOHR)
                  .k0(Cfg::WP_KX, Cfg::WP_KY, Cfg::WP_KZ)
                  .orthogonalise_against_occupied(electrons);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx
              << " norm_after=" << report.norm_after
              << " max_overlap=" << report.max_overlap << "\n";

    {
        std::ofstream f("results/raw/observables/wp_config.txt");
        f << "wp_center_bohr = " << Cfg::WP_CX_BOHR << " " << Cfg::WP_CY_BOHR
                                 << " " << Cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr  = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv = " << Cfg::WP_KX << " " << Cfg::WP_KY << " "
                                 << Cfg::WP_KZ << "\n"
          << "wp_energy_ev   = " << Cfg::WP_EKIN_EV << "\n"
          << "wp_state_index = " << wp_idx << "\n";
    }
    {
        std::ofstream f("results/raw/observables/wp_injection_report.txt");
        f << std::setprecision(16);
        f << "norm_before      = " << report.norm_before << "\n"
          << "norm_after       = " << report.norm_after << "\n"
          << "max_overlap      = " << report.max_overlap << "\n";
    }

    inqkit::io::RealField3DWriter total_wr(
        "results/raw/vti/density_total", vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(
        "results/raw/vti/density_system", vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter wp_density_wr(
        "results/raw/vti/density_wp", vti_layout, {.overwrite=true});
    inqkit::io::ComplexField3DWriter wp_wf_wr(
        "results/raw/vti/wavefunction_wp",
        {.field_name = "wavefunction", .include_meta = false, .emit_raw = false,
         .emit_vti = true,
         .vti_format = inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite=true});

    // t=0: system = bath (full - wp); total = full DFT density; wp = wp alone.
    {
        auto full0 = inqkit::fields::density::total(electrons);
        auto wp0   = inqkit::fields::density::orbital(electrons, wp_idx);
        auto bath0 = inqkit::fields::density::total_excluding_orbital(
            full0, wp0, 1.0);
        system_wr.write(bath0, 0.0, 0);
        wp_density_wr.write(wp0, 0.0, 0);
        total_wr.write(full0, 0.0, 0);
    }

    inqkit::observables::OrbitalOverlapMatrix overlap_obs(
        electrons, wp_idx, "results/raw/observables/overlap");
    auto shell_table = inqkit::jellium::shells::enumerate_for_n_states(n_states);
    auto proxies     = inqkit::jellium::shells::pick_proxies(shell_table, 2);
    inqkit::observables::OrbitalOverlapMatrix overlap_proxy_obs(
        electrons, wp_idx, "results/raw/observables/overlap_proxies");
    inqkit::jellium::shells::write_shells_csv(
        shell_table, proxies, "results/raw/observables/overlap_proxies");
    const int PROXY_SNAPSHOT_STRIDE = 5 * Cfg::WRITE_EVERY;

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x  = sel.dipole_y  = sel.dipole_z  = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(
        "results/raw/observables/observables.csv", sel);
    obs_writer.write_header();

    inqkit::observables::StateEnergyWriter state_energy_wr(
        "results/raw/observables/state_energies.csv", true);
    inqkit::observables::OccupationsWriter occupations_wr(
        "results/raw/observables/occupations_vs_time.csv");
    inqkit::observables::MomentumDistribution momentum_dist(
        "results/raw/observables/momentum_distribution.csv", wp_idx, Cfg::L_BOHR,
        {.n_bins = 64, .k_max_bohr_inv = 0.0, .write_every = 10 * Cfg::WRITE_EVERY});
    inqkit::observables::WPMomentumStats wp_momentum_stats(
        "results/raw/observables/wp_momentum_stats.csv", wp_idx,
        {.write_every = Cfg::WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        "results/raw/observables/wp_real_space_stats.csv", wp_idx,
        {.write_every = Cfg::WRITE_EVERY});
    inqkit::observables::DensityDelta density_delta(
        "results/raw/vti/density_delta", "results/raw/vti/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true, .compute_l2 = true,
         .coarse_bin_bohr = 3.0});
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(
        electrons, n_states - 1, "results/raw/observables/overlap_full");
    overlap_full_obs.snapshot(electrons, 0.0, 0);
    overlap_obs.snapshot_wp_only(electrons, 0.0, 0);
    overlap_proxy_obs.snapshot_proxies(electrons, proxies, 0.0, 0);

    inqkit::RealTimeSession rt_obs(ions, electrons, Cfg::WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const &ctx) {
        // Full DFT density (WP-included) and WP orbital density every frame.
        auto full_f = inqkit::fields::density::total(*ctx.electrons);
        auto wp_f   = inqkit::fields::density::orbital(*ctx.electrons, wp_idx);
        // Bath = full - wp: the target response, WP removed (integrates to N_e).
        auto bath_f = inqkit::fields::density::total_excluding_orbital(
            full_f, wp_f, 1.0);

        system_wr.write(bath_f, ctx.time_au, ctx.step);   // bath only (WP removed)
        total_wr.write( full_f, ctx.time_au, ctx.step);    // true total (WP incl.)
        wp_density_wr.write(wp_f, ctx.time_au, ctx.step);  // WP density alone

        const double l2 = density_delta.snapshot(bath_f, ctx.time_au, ctx.step);
        inqkit::StepContext ctx_out = ctx;
        ctx_out.density_l2 = l2;
        obs_writer.append(ctx_out);

        // Wavefunction (complex) for momentum / loss-function — high cadence.
        if (ctx.step % Cfg::WF_WRITE_EVERY == 0) {
            auto wp_wf = inqkit::fields::orbital::wavefunction(
                *ctx.electrons, wp_idx);
            char wf_name[64];
            std::snprintf(wf_name, sizeof(wf_name), "wavefunction_t%06d", ctx.step);
            wp_wf_wr.write(wp_wf, std::string(wf_name));
        }
        if (ctx.step % 10 == 0)
            overlap_obs.snapshot_wp_only(*ctx.electrons, ctx.time_au, ctx.step);
        if (ctx.step > 0 && ctx.step % PROXY_SNAPSHOT_STRIDE == 0)
            overlap_proxy_obs.snapshot_proxies(*ctx.electrons, proxies,
                                               ctx.time_au, ctx.step);
    });

    real_time::propagate(
        ions, electrons,
        [&](auto const &data) {
            rt_obs.step(data);
            if (data.iter() % (5 * Cfg::WRITE_EVERY) == 0) {
                state_energy_wr.snapshot(data);
                occupations_wr.snapshot(data);
            }
            momentum_dist.maybe_accumulate(data);
            wp_momentum_stats.maybe_accumulate(data);
            wp_real_space_stats.maybe_accumulate(data);
        },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(Cfg::N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .observables_current()
            .observables_dipole());

    overlap_full_obs.snapshot(electrons, Cfg::DT_AU * Cfg::N_STEPS, Cfg::N_STEPS);

    auto t_wallclock_end = std::chrono::steady_clock::now();
    double wall_seconds =
        std::chrono::duration<double>(t_wallclock_end - t_wallclock_start).count();
    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "RUN SUMMARY\n===========\n\n1. Run identity\n---------------\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "run_type        = wave-packet projectile, jellium TDDFT (ALDA)\n"
          << "date_finished   = " << iso_now() << "\n"
          << "wall_time_s     = " << wall_seconds << "\n\n"
          << "3. System configuration\n-----------------------\n"
          << "cell_bohr       = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
          << "n_electrons     = " << n_electrons << "\n"
          << "n_occupied      = " << n_occupied << "\n"
          << "extra_states    = " << Cfg::EXTRA_STATES << "\n"
          << "wp_state_index  = " << wp_idx << "\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "xc_functional   = LDA (ALDA in TDDFT)\n\n"
          << "5. Wavepacket configuration and injection\n"
          << "-----------------------------------------\n"
          << "wp_enabled      = yes\n"
          << "wp_center_bohr  = " << Cfg::WP_CX_BOHR << " " << Cfg::WP_CY_BOHR
                                  << " " << Cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr   = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv  = " << Cfg::WP_KX << " " << Cfg::WP_KY << " "
                                  << Cfg::WP_KZ << "\n"
          << "wp_energy_ev    = " << Cfg::WP_EKIN_EV << "\n"
          << "norm_after      = " << report.norm_after << "\n"
          << "max_overlap     = " << report.max_overlap << "\n\n"
          << "6. Real-time configuration\n--------------------------\n"
          << "rt_num_steps    = " << Cfg::N_STEPS << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "total_time_au   = " << (Cfg::DT_AU * Cfg::N_STEPS) << "\n"
          << "write_every     = " << Cfg::WRITE_EVERY << "\n"
          << "wf_write_every  = " << Cfg::WF_WRITE_EVERY << "\n\n"
          << "9. End-of-run diagnostics\n-------------------------\n"
          << "run_completed   = true\n";
    }
    std::cout << "Done. Wall time " << wall_seconds << " s.\n";
    return 0;
}
