// ============================================================================
// run_wp_n162_L50_E50_v2/run.cpp — Gaussian wave-packet projectile in a cubic
// 50^3 Bohr periodic jellium bath at N=162, dx=0.40, E_kin=100 eV.
//
// Run-3 (WP half) of the 2026-05-21 Emilio meeting campaign — fastest path to an
// emailable result; validates the whole infra pipeline end-to-end (see
// docs/plans/jellium-meeting-2026-05-21.md §"Per-family motivations"
// "Run-1 — 25 eV WP"). No matched classical sibling (Run-9 is the
// post-Knudsen post-meeting companion).
//
// Cfg: jellium::config::Backlog_L50_cubic_WP_E50.
// WP is a Gaussian at launch_z = -L/2 + 4σ = -5 Bohr (updated 2026-05-17
// to the universal boundary rule), k_0=2.711 bohr^-1 in +z, sigma_r=5.0
// Bohr, KE=100 eV. N_STEPS=462, WRITE_EVERY=2 (~231 frames).
//
// PHYSICS: v=2.71 sits between v_F=0.55 (Lindhard) and v=10.5 (Bethe).
// Plasmon at this geometry is at ~3.53 eV (m=1, from
// run_plasmon_n162_L50_E15). The WP's bandwidth sigma_k=0.20 spans
// k=2.51..2.91, well above the Lindhard particle-hole continuum top edge
// but well below the high-k cutoff. Expect SOME plasmon coupling
// (Bohm-Gross dispersion permits it at this k) and some single-particle
// e-h excitation.
//
// GRID: dx=0.40 gives k_Nyquist=π/0.40=7.85 bohr^-1, comfortable
// headroom for k_max=k_0+3sigma_k=3.31. The dx=0.30 grid used by the
// E=1500 WP attempt deadlocked at init; at dx=0.40 the orbital memory
// is roughly half, so init should fit in 24 GB GPU.
//
// CALLBACK DESIGN: density::total per WRITE_EVERY uses the working
// bulk-copy patch (inqkit/fields/density.hpp), so density_total VTI
// series and density_delta tracking are ENABLED for this run (contrast
// the E=1500 WP attempt which stripped them due to dx=0.30 memory
// pressure). density::orbital (per-element loop) remains stripped —
// the WP-density VTI series and cod_z observable are reconstructed in
// postprocess from density_total - density_bath_t0 if needed.
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

#include "../shared/configs/electron_proj_backlog_L50_cubic.hpp"
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
using Cfg = jellium::config::Backlog_L50_cubic_WP_E50;

static std::string iso_now() {
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

int main() {
    auto t_wallclock_start = std::chrono::steady_clock::now();

    const std::string RUN_NAME = "run_wp_n162_L50_E50_v2";
    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p40";

    std::cout << "\n=== " << RUN_NAME << " ===\n"
              << "  cell = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
              << "  N_electrons = " << Cfg::N_ELECTRONS << "\n"
              << "  spacing = " << Cfg::SPACING_BOHR << " Bohr\n"
              << "  WP: sigma=" << Cfg::WP_SIGMA_BOHR
              << "  k0=" << Cfg::WP_K0
              << "  E=" << Cfg::WP_EKIN_EV << " eV"
              << "  launch (" << Cfg::WP_CX_BOHR << ", " << Cfg::WP_CY_BOHR
              << ", " << Cfg::WP_CZ_BOHR << ") Bohr\n"
              << "  dt=" << Cfg::DT_AU << " a.u., N_steps="
              << Cfg::N_STEPS << ", t_total="
              << (Cfg::DT_AU * Cfg::N_STEPS) << " a.u.\n"
              << "  checkpoint = " << GS_DIR << "\n\n";

    if (!std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: GS checkpoint missing — run save_gs first.\n";
        return 2;
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
    std::cout << "  Loaded GS from " << GS_DIR << "\n";

    jellium::eigenvalues::copy_from_checkpoint(
        GS_DIR, "results/raw/observables/eigenvalues");

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    const int n_occupied  = n_electrons / 2;
    std::cout << "  num_states=" << n_states
              << "  num_electrons=" << n_electrons
              << "  n_occupied=" << n_occupied << "\n";

    // ----- Output skeleton -----------------------------------------------
    std::filesystem::create_directories("results/raw/observables");
    std::filesystem::create_directories("results/raw/observables/overlap");
    std::filesystem::create_directories("results/raw/observables/overlap_full");
    std::filesystem::create_directories("results/raw/vti/density_total");
    std::filesystem::create_directories("results/raw/vti/density_system");
    std::filesystem::create_directories("results/raw/vti/density_delta");
    std::filesystem::create_directories("results/raw/vti/density_delta_coarse");
    std::filesystem::create_directories("results/raw/vti/density_gs_system");

    // Stub run_summary so a crash leaves metadata.
    {
        std::ofstream s("results/run_summary.txt");
        s << "RUN SUMMARY (stub - written at start)\n"
          << "=====================================\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "date_started    = " << iso_now() << "\n"
          << "checkpoint_dir  = " << GS_DIR << "\n"
          << "cell_bohr       = " << Cfg::L_BOHR
                                  << "^3 (cubic, periodic)\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "n_electrons     = " << n_electrons << "\n"
          << "n_states        = " << n_states << "\n"
          << "wp_ekin_ev      = " << Cfg::WP_EKIN_EV << "\n"
          << "wp_sigma_bohr   = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "n_steps         = " << Cfg::N_STEPS << "\n"
          << "run_completed   = false\n";
    }

    // ----- VTI layout used by all density writers ------------------------
    inqkit::io::RealField3DLayout vti_layout{
        .field_name   = "density",
        .include_meta = false,
        .emit_raw     = false,
        .emit_vti     = true,
        .vti_format   = inqkit::io::VTIWriteOptions::Format::binary,
    };

    // ----- t=0 GS density VTI -------------------------------------------
    {
        inqkit::io::RealField3DWriter gs_sys(
            "results/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
        gs_sys.write(inqkit::fields::density::total(electrons),
                     "density_gs_system");
    }

    // ----- WP injection --------------------------------------------------
    auto wp = inqkit::WavePacket{}
                  .center(Cfg::WP_CX_BOHR, Cfg::WP_CY_BOHR, Cfg::WP_CZ_BOHR)
                  .sigma(Cfg::WP_SIGMA_BOHR)
                  .k0(Cfg::WP_KX, Cfg::WP_KY, Cfg::WP_KZ)
                  .orthogonalise_against_occupied(electrons);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx
              << "  norm_after=" << report.norm_after
              << "  max_overlap=" << report.max_overlap << "\n";

    // ----- WP artefacts --------------------------------------------------
    {
        std::ofstream f("results/raw/observables/wp_config.txt");
        f << "wp_center_bohr = " << Cfg::WP_CX_BOHR << " "
                                 << Cfg::WP_CY_BOHR << " "
                                 << Cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr  = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv = " << Cfg::WP_KX << " " << Cfg::WP_KY << " "
                                 << Cfg::WP_KZ << "\n"
          << "wp_energy_ev   = " << Cfg::WP_EKIN_EV << "\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "orthogonalised = " << (report.orthogonalised ? "yes" : "no") << "\n";
    }
    {
        std::ofstream f("results/raw/observables/wp_injection_report.txt");
        f << std::setprecision(16);
        f << "norm_before      = " << report.norm_before << "\n"
          << "norm_after       = " << report.norm_after << "\n"
          << "max_overlap      = " << report.max_overlap << "\n"
          << "passed_tolerance = " << (report.passed_tolerance ? "yes" : "no") << "\n";
    }
    // Initial WP density / wavefunction VTI writes SKIPPED. Same reason as
    // E=1500 WP run: density::orbital uses a per-element GPU->host loop
    // (~minutes per call at this grid size).

    // ----- Real-time outputs --------------------------------------------
    // overlap_obs (WP-only): 1 wavefunction extraction per call, used
    // every 10 propagator steps to track WP decoherence into the GS basis.
    inqkit::observables::OrbitalOverlapMatrix overlap_obs(
        electrons, wp_idx, "results/raw/observables/overlap");

    // overlap_proxy_obs: GS-basis-projected occupations via 2 fixed
    // proxies per shell. Snapshot every PROXY_SNAPSHOT_STRIDE steps.
    auto shell_table = inqkit::jellium::shells::enumerate_for_n_states(n_states);
    auto proxies     = inqkit::jellium::shells::pick_proxies(shell_table, 2);
    inqkit::observables::OrbitalOverlapMatrix overlap_proxy_obs(
        electrons, wp_idx, "results/raw/observables/overlap_proxies");
    inqkit::jellium::shells::write_shells_csv(
        shell_table, proxies, "results/raw/observables/overlap_proxies");
    std::cout << "  Shell table: " << shell_table.size() << " shells, "
              << proxies.size() << " proxies (2 per shell)\n";

    const int PROXY_SNAPSHOT_STRIDE = 5 * Cfg::WRITE_EVERY;

    inqkit::io::RealField3DWriter total_wr(
        "results/raw/vti/density_total", vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(
        "results/raw/vti/density_system", vti_layout, {.overwrite=true});

    // t=0 system density frame (we DO write the bath density at t=0 via
    // bulk-copy density::total). The WP-only density at t=0 is skipped.
    {
        auto sys0 = inqkit::fields::density::total(electrons);
        system_wr.write(sys0, 0.0, 0);
        total_wr.write( sys0, 0.0, 0);
    }

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
        "results/raw/observables/state_energies.csv", /*emit_variance=*/true);

    inqkit::observables::OccupationsWriter occupations_wr(
        "results/raw/observables/occupations_vs_time.csv");

    inqkit::observables::MomentumDistribution momentum_dist(
        "results/raw/observables/momentum_distribution.csv",
        wp_idx, Cfg::L_BOHR,
        {.n_bins = 64, .k_max_bohr_inv = 0.0,
         .write_every = 10 * Cfg::WRITE_EVERY});

    // New per-step WP momentum/position statistics (Infra-4/5, validated
    // by Heisenberg known-case tests in Tutorial/wp-{momentum,real-space}-
    // stats-test). Method-B (Knudsen) stopping power consumes e_kin_ha
    // from wp_momentum_stats.csv via the inqview knudsen_ke phase.
    inqkit::observables::WPMomentumStats wp_momentum_stats(
        "results/raw/observables/wp_momentum_stats.csv",
        wp_idx, {.write_every = Cfg::WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        "results/raw/observables/wp_real_space_stats.csv",
        wp_idx, {.write_every = Cfg::WRITE_EVERY});

    inqkit::observables::DensityDelta density_delta(
        "results/raw/vti/density_delta",
        "results/raw/vti/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0});

    // Full-matrix overlap at t=0 (n_ref × n_evolved). Used by analyse.py's
    // post-run heatmap + effective-GS-occupation analysis. The mid-trajectory
    // snapshot used by the classical run is skipped here because each full
    // snapshot incurs n_ref + 1 orbital extractions (slow under the
    // per-element GPU→host loop). Only t=0 and t=N_STEPS are taken — the
    // pair the user-requested overlap analysis needs.
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(
        electrons, n_states - 1, "results/raw/observables/overlap_full");
    overlap_full_obs.snapshot(electrons, 0.0, 0);
    std::cout << "  Full overlap snapshot at step 0\n";

    // t=0 sanity snapshots:
    overlap_obs.snapshot_wp_only(electrons, 0.0, 0);
    overlap_proxy_obs.snapshot_proxies(electrons, proxies, 0.0, 0);

    // ----- Real-time session callbacks ----------------------------------
    //
    // ENABLED per-step observables (vs the E=1500 WP stripped callback):
    //   - system density VTI via bulk-copy density::total (cheap)
    //   - density_delta tracking (uses density::total, cheap)
    //   - obs_writer CSV row (cheap, INQ viewables)
    //   - WP-only overlap every 10 steps (1 wavefunction extraction)
    //   - proxy overlap every PROXY_SNAPSHOT_STRIDE steps (~18 extractions
    //     per snapshot)
    //
    // STILL STRIPPED (per-element loop, expensive at any grid):
    //   - density::orbital(electrons, wp_idx) per step → no WP-density
    //     VTI series, no cod_z observable (reconstruct in postprocess
    //     from density_total - density_bath_t0).
    inqkit::RealTimeSession rt_obs(ions, electrons, Cfg::WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const &ctx) {
        // System density via bulk-copy patch (cheap).
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        total_wr.write( sys_f, ctx.time_au, ctx.step);

        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext ctx_out = ctx;
        ctx_out.density_l2 = l2;
        obs_writer.append(ctx_out);

        // WP-only overlap every 10 propagator steps.
        if (ctx.step % 10 == 0) {
            overlap_obs.snapshot_wp_only(*ctx.electrons,
                                         ctx.time_au, ctx.step);
        }

        // Proxy overlap every PROXY_SNAPSHOT_STRIDE steps.
        if (ctx.step > 0 && ctx.step % PROXY_SNAPSHOT_STRIDE == 0) {
            overlap_proxy_obs.snapshot_proxies(
                *ctx.electrons, proxies, ctx.time_au, ctx.step);
        }
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

    // Final full-matrix overlap snapshot for analyse.py's overlap-heatmap
    // and effective-GS-projected-occupations bar chart.
    overlap_full_obs.snapshot(electrons,
        Cfg::DT_AU * Cfg::N_STEPS, Cfg::N_STEPS);
    std::cout << "  Full overlap snapshot at step " << Cfg::N_STEPS
              << " (end)\n";

    auto t_wallclock_end = std::chrono::steady_clock::now();
    double wall_seconds =
        std::chrono::duration<double>(t_wallclock_end - t_wallclock_start).count();
    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "RUN SUMMARY\n===========\n\n"
          << "1. Run identity\n---------------\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "run_type        = wave-packet projectile, jellium TDDFT (ALDA)\n"
          << "date_finished   = " << iso_now() << "\n"
          << "wall_time_s     = " << wall_seconds << "\n\n"
          << "3. System configuration\n-----------------------\n"
          << "cell_bohr       = " << Cfg::L_BOHR
                                  << "^3 (cubic, periodic)\n"
          << "n_electrons     = " << n_electrons << "\n"
          << "n_occupied      = " << n_occupied << "\n"
          << "extra_states    = " << Cfg::EXTRA_STATES << "\n"
          << "wp_state_index  = " << wp_idx << "\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "xc_functional   = LDA (ALDA in TDDFT)\n\n"
          << "5. Wavepacket configuration and injection\n-----------------------------------------\n"
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
          << "write_every     = " << Cfg::WRITE_EVERY << "\n\n"
          << "9. End-of-run diagnostics\n-------------------------\n"
          << "run_completed   = true\n";
    }

    std::cout << "Done. Wall time " << wall_seconds << " s.\n";
    return 0;
}
