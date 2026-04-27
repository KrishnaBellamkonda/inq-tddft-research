// ============================================================================
// shared/cpp/run_template.hpp
//
// The propagation body shared by every coronene run.cpp. Each run.cpp is a
// thin wrapper that:
//   1. Includes its variant config header (or just the base);
//   2. Defines RUN_NAME, GS_CHECKPOINT_DIR, plus any parameter overrides;
//   3. Calls coronene::run_template::run_propagation(...).
//
// The function is templated only over a Cfg struct so each run can supply its
// own override values without runtime branching. The outputs land under
// ./results/ per docs/results_folder_structure_spec.md and are produced by
// scripts/coronene_postprocess.py downstream.
//
// Conventions:
//   * Loads electrons from <gs_checkpoint_dir>;
//   * Builds the cell and ions from cfg.LX_BOHR / cfg.LY_BOHR / cfg.LZ_BOHR
//     and shared/geometry/coronene.xyz;
//   * Injects the wave packet into the last extra state with MGS
//     orthogonalisation against the occupied subspace;
//   * Writes the GS density (system) and the GS orbital densities once to
//     results/raw/ground_state/ and results/raw/vti/density_gs_*/;
//   * Writes the initial WP density and complex wavefunction to
//     results/raw/wavepacket/density_wp_initial/ and wavefunction_wp_initial/;
//   * Writes time-stepped density_rt_total / _system / _wp into
//     results/raw/vti/density_rt_*/ at every WRITE_EVERY step;
//   * Writes observables CSV (energy, current, dipole) at every step;
//   * Records WP-only orbital overlap (snapshot_wp_only) every step;
//   * Accumulates 20 LEED screens with one full-time accumulator and
//     one physics window per screen (centroid arrives -> trailing edge
//     passes; backscattering screens use the rebound model) plus the paper
//     window applied uniformly to all 20 screens;
//   * Saves instantaneous LEED snapshots every SCREEN_SNAP_EVERY steps as
//     results/raw/screens/instantaneous/screen_NN_tXXXXXX.dat (flat).
//
// The driver only writes raw artefacts. All plotting / GIFs / VTI of the
// derived 2D fields is done by coronene_postprocess.py from those raws.
// ============================================================================
#pragma once

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include "leed_screen_layout.hpp"
#include "results_paths.hpp"

#include <array>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace coronene::run_template {

// Produce results/raw/density/density_rt_total/<basename>.vti by adding the
// system and WP densities pointwise. Both inputs share the same grid (they
// come from the same electrons object), so a flat add on .values is correct.
inline inqkit::fields::RealField3D
add_real_fields(inqkit::fields::RealField3D const &a,
                inqkit::fields::RealField3D const &b) {
    if (a.values.size() != b.values.size()
        || a.nx != b.nx || a.ny != b.ny || a.nz != b.nz) {
        throw std::runtime_error("add_real_fields: grid mismatch");
    }
    inqkit::fields::RealField3D out = a;  // copy metadata
    for (std::size_t i = 0; i < a.values.size(); ++i) {
        out.values[i] = a.values[i] + b.values[i];
    }
    return out;
}

inline std::string iso_now() {
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

// Per-run config. Cfg is the variant struct (e.g. coronene::config::E30) so
// every override value is statically known. Cfg is duck-typed: it exposes
// LX_BOHR, LY_BOHR, LZ_BOHR, CUTOFF_HA, EXTRA_STATES, WP_CX_BOHR, WP_CY_BOHR,
// WP_CZ_BOHR, WP_SIGMA_BOHR, WP_KX, WP_KY, WP_KZ, WP_EKIN_EV, DT_AU, N_STEPS,
// WRITE_EVERY, SCREEN_SNAP_EVERY, T1_AU, T2_AU, plus the constants from base.
//
// `geometry_xyz_path` is the absolute path to shared/geometry/coronene.xyz
// (passed in by run.cpp because relative-path resolution in inq-run is
// brittle).
//
// `gs_checkpoint_dir` is the absolute path to checkpoints/<sig>/ produced by
// the corresponding save_gs/<sig>/run.cpp.
//
// `run_name` is the short tag (e.g. "run_E30") that appears in
// run_summary.txt and in plot titles produced downstream.
template <class Cfg>
int run_propagation(std::string const &run_name,
                    std::string const &geometry_xyz_path,
                    std::string const &gs_checkpoint_dir) {
    using namespace inq;
    using namespace inq::magnitude;

    auto t_wallclock_start = std::chrono::steady_clock::now();

    std::cout << "\n=== " << run_name << " ===\n"
              << "  cell = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR
              << " x " << Cfg::LZ_BOHR << " Bohr\n"
              << "  geometry  = " << geometry_xyz_path << "\n"
              << "  checkpoint = " << gs_checkpoint_dir << "\n"
              << "  WP sigma = " << Cfg::WP_SIGMA_BOHR
              << " Bohr, |k| = " << Cfg::WP_KZ_MAGNITUDE
              << " Bohr^-1, dt = " << Cfg::DT_AU
              << " a.u., N_steps = " << Cfg::N_STEPS << "\n";

    if (!std::filesystem::exists(gs_checkpoint_dir)) {
        std::cerr << "FATAL: checkpoint directory does not exist: "
                  << gs_checkpoint_dir << "\n"
                  << "Run the corresponding save_gs/<sig>/run.cpp first.\n";
        return 2;
    }

    // ----- Cell + atoms ----------------------------------------------------
    auto cell = systems::cell::orthorhombic(
        Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b
    ).finite();
    auto ions = systems::ions::parse(geometry_xyz_path, cell);
    std::cout << "  Atoms: " << ions.size() << "\n";

    // ----- Electrons (skeleton) + load checkpoint --------------------------
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .cutoff(Cfg::CUTOFF_HA * 1.0_Ha)
            .extra_states(Cfg::EXTRA_STATES));
    electrons.load(gs_checkpoint_dir);
    std::cout << "  Loaded GS from " << gs_checkpoint_dir << "\n";

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    const int n_occupied  = n_electrons / 2;
    std::cout << "  num_states = " << n_states
              << "  num_electrons = " << n_electrons
              << "  n_occupied = " << n_occupied << "\n";

    // ----- Stub run_summary.txt so a crashed run still has metadata --------
    {
        std::ofstream s(coronene::results::run_summary_path());
        s << "RUN SUMMARY (stub - written at start)\n"
          << "=====================================\n"
          << "run_name        = " << run_name << "\n"
          << "date_started    = " << iso_now() << "\n"
          << "geometry_file   = " << geometry_xyz_path << "\n"
          << "checkpoint_dir  = " << gs_checkpoint_dir << "\n"
          << "cell_bohr       = " << Cfg::LX_BOHR << " " << Cfg::LY_BOHR << " "
                                  << Cfg::LZ_BOHR << "\n"
          << "boundary        = finite\n"
          << "xc              = LDA (ALDA in TDDFT)\n"
          << "cutoff_ha       = " << Cfg::CUTOFF_HA << "\n"
          << "extra_states    = " << Cfg::EXTRA_STATES << "\n"
          << "num_states      = " << n_states << "\n"
          << "num_electrons   = " << n_electrons << "\n"
          << "n_occupied      = " << n_occupied << "\n"
          << "wp_sigma_bohr   = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_offset_bohr  = " << Cfg::WP_CZ_BOHR << "\n"
          << "wp_ekin_ev      = " << Cfg::WP_EKIN_EV << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "n_steps         = " << Cfg::N_STEPS << "\n"
          << "write_every     = " << Cfg::WRITE_EVERY << "\n"
          << "run_completed   = false\n";
    }

    // ----- Ground-state artefacts (before WP injection) --------------------
    inqkit::io::RealField3DLayout vti_layout{
        .field_name  = "density",
        .include_meta = false,
        .emit_raw    = false,
        .emit_vti    = true,
        .vti_format  = inqkit::io::VTIWriteOptions::Format::binary,
    };
    {
        // GS total density (system, before any WP)
        inqkit::io::RealField3DWriter gs_sys(
            coronene::results::vti_density_gs_system_dir(), vti_layout,
            {.overwrite = true});
        gs_sys.write(inqkit::fields::density::total(electrons),
                     "density_gs_system");

        // GS orbital densities (one VTI per occupied orbital)
        inqkit::io::RealField3DWriter gs_orb(
            coronene::results::vti_density_gs_orbitals_dir(), vti_layout,
            {.overwrite = true});
        for (int i = 0; i < n_occupied; ++i) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "orbital_%04d", i);
            gs_orb.write(inqkit::fields::density::orbital(electrons, i),
                         std::string(buf));
        }
    }

    // ----- WP injection ----------------------------------------------------
    auto wp = inqkit::WavePacket{}
                  .center(Cfg::WP_CX_BOHR, Cfg::WP_CY_BOHR, Cfg::WP_CZ_BOHR)
                  .sigma(Cfg::WP_SIGMA_BOHR)
                  .k0(Cfg::WP_KX, Cfg::WP_KY, Cfg::WP_KZ)
                  .orthogonalise_against_occupied(electrons);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected: state_index = " << wp_idx
              << "  norm_after = "  << report.norm_after
              << "  max_overlap = " << report.max_overlap << "\n";

    // ----- Wavepacket artefacts -------------------------------------------
    {
        std::ofstream f(coronene::results::wp_config_path());
        f << "wp_center_bohr   = " << Cfg::WP_CX_BOHR << " "
                                   << Cfg::WP_CY_BOHR << " "
                                   << Cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr    = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv   = " << Cfg::WP_KX << " " << Cfg::WP_KY << " "
                                   << Cfg::WP_KZ << "\n"
          << "wp_energy_ev     = " << Cfg::WP_EKIN_EV << "\n"
          << "wp_state_index   = " << wp_idx << "\n"
          << "orthogonalised   = " << (report.orthogonalised ? "yes" : "no") << "\n";
    }
    {
        std::ofstream f(coronene::results::wp_injection_report_path());
        f << std::setprecision(16);
        f << "norm_before      = " << report.norm_before << "\n"
          << "norm_after       = " << report.norm_after << "\n"
          << "max_overlap      = " << report.max_overlap << "\n"
          << "passed_tolerance = " << (report.passed_tolerance ? "yes" : "no") << "\n";
    }
    {
        // Initial WP density (real) and wavefunction (complex)
        inqkit::io::RealField3DWriter wp_d(
            coronene::results::wp_density_initial_dir(), vti_layout,
            {.overwrite = true});
        wp_d.write(inqkit::fields::density::orbital(electrons, wp_idx),
                   "density_wp_initial");

        inqkit::io::ComplexField3DLayout cvti_layout{
            .field_name = "psi_wp",
            .include_meta = false,
            .emit_raw    = false,
            .emit_vti    = true,
            .vti_format  = inqkit::io::VTIWriteOptions::Format::binary,
        };
        inqkit::io::ComplexField3DWriter wp_psi(
            coronene::results::wp_wavefunction_initial_dir(), cvti_layout,
            {.overwrite = true});
        wp_psi.write(inqkit::fields::orbital::wavefunction(electrons, wp_idx),
                     "wavefunction_wp_initial");
    }

    // ----- Real-time outputs setup ----------------------------------------
    inqkit::observables::OrbitalOverlapMatrix overlap_obs(
        electrons, wp_idx, coronene::results::overlap_dir());

    inqkit::io::RealField3DWriter total_wr(
        coronene::results::vti_density_total_dir(), vti_layout, {.overwrite = true});
    inqkit::io::RealField3DWriter system_wr(
        coronene::results::vti_density_system_dir(), vti_layout, {.overwrite = true});
    inqkit::io::RealField3DWriter wp_wr(
        coronene::results::vti_density_wp_dir(), vti_layout, {.overwrite = true});

    // t = 0 frames
    {
        auto sys0   = inqkit::fields::density::total(electrons);
        auto wp0    = inqkit::fields::density::orbital(electrons, wp_idx);
        auto total0 = add_real_fields(sys0, wp0);
        system_wr.write(sys0,   0.0, 0);
        wp_wr.write(    wp0,    0.0, 0);
        total_wr.write( total0, 0.0, 0);
    }

    // Observables CSV
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x  = sel.dipole_y  = sel.dipole_z  = true;
    inqkit::io::ObservablesWriter obs_writer(
        coronene::results::observables_csv_path(), sel);
    obs_writer.write_header();

    // ----- Screens: 20 positions, three accumulators per screen ----------
    //   acc_full[k]          — every step (always-on);
    //   acc_screen_window[k] — physics window: WP centroid arrives -> trailing
    //                          edge passes (forward) / back-scattered trailing
    //                          edge passes (backscattering side);
    //   acc_paper[k]         — single global paper window [T1_AU, T2_AU].
    auto screen_z = coronene::layout::screen_z_positions(Cfg::LZ_BOHR);
    {
        std::cout << "  Screen z-positions (Bohr):";
        for (auto z : screen_z) std::cout << " " << z;
        std::cout << "\n";
    }

    using inqkit::screens::LeedPatternAccumulator;
    using inqkit::screens::PlaneScreen;

    std::array<LeedPatternAccumulator, coronene::layout::N_SCREENS> acc_full;
    std::array<LeedPatternAccumulator, coronene::layout::N_SCREENS> acc_screen_window;
    std::array<LeedPatternAccumulator, coronene::layout::N_SCREENS> acc_paper;

    // Per-screen physics windows derived from compute_screen_window().
    std::array<coronene::layout::ScreenWindow, coronene::layout::N_SCREENS>
        screen_windows;

    for (int k = 0; k < coronene::layout::N_SCREENS; ++k) {
        auto label = coronene::layout::screen_label(k);
        acc_full[k]          = LeedPatternAccumulator(PlaneScreen{screen_z[k], label});
        acc_screen_window[k] = LeedPatternAccumulator(PlaneScreen{screen_z[k], label});
        acc_paper[k]         = LeedPatternAccumulator(PlaneScreen{screen_z[k], label});
        screen_windows[k] = coronene::layout::compute_screen_window(
            screen_z[k], Cfg::WP_CZ_BOHR, Cfg::WP_SIGMA_BOHR, Cfg::WP_K0,
            Cfg::LZ_BOHR, Cfg::DT_AU * Cfg::N_STEPS,
            Cfg::WP_ENVELOPE_SIGMAS);
    }

    // Persist the configuration so the postprocess can reproduce filenames.
    {
        std::ofstream f(coronene::results::screens_config_path());
        f << "screen_index,z_bohr,label,window_kind,t_start_au,t_end_au\n";
        for (int k = 0; k < coronene::layout::N_SCREENS; ++k) {
            f << k << "," << std::fixed << std::setprecision(6)
              << screen_z[k] << "," << coronene::layout::screen_label(k)
              << "," << (screen_windows[k].is_back ? "back" : "forward")
              << "," << screen_windows[k].t_start_au
              << "," << screen_windows[k].t_end_au << "\n";
        }
    }
    {
        std::ofstream f(coronene::results::screens_window_ranges_path());
        f << "kind,screen_index,t_start_au,t_end_au\n";
        for (int k = 0; k < coronene::layout::N_SCREENS; ++k) {
            f << (screen_windows[k].is_back ? "back" : "forward") << ","
              << k << ","
              << std::fixed << std::setprecision(6)
              << screen_windows[k].t_start_au << ","
              << screen_windows[k].t_end_au << "\n";
        }
        f << "paper,-1,"
          << std::fixed << std::setprecision(6)
          << Cfg::T1_AU << "," << Cfg::T2_AU << "\n";
    }

    // ----- Real-time session: density frames + observables + overlap + screens
    inqkit::RealTimeSession rt_dens(ions, electrons, Cfg::WRITE_EVERY);
    rt_dens.add([&](inqkit::StepContext const &ctx) {
        auto sys_f   = inqkit::fields::density::total(*ctx.electrons);
        auto wp_f    = inqkit::fields::density::orbital(*ctx.electrons, wp_idx);
        auto total_f = add_real_fields(sys_f, wp_f);
        system_wr.write(sys_f,   ctx.time_au, ctx.step);
        wp_wr.write(    wp_f,    ctx.time_au, ctx.step);
        total_wr.write( total_f, ctx.time_au, ctx.step);
    });

    inqkit::RealTimeSession rt_obs(ions, electrons, /*write_every=*/1);
    rt_obs.add([&](inqkit::StepContext const &ctx) {
        obs_writer.append(ctx);
        overlap_obs.snapshot_wp_only(*ctx.electrons, ctx.time_au, ctx.step);

        // Always-on accumulators.
        for (auto &a : acc_full) a.accumulate(*ctx.electrons, Cfg::DT_AU);

        // Per-screen physics window (one bracket per screen).
        for (int k = 0; k < coronene::layout::N_SCREENS; ++k) {
            const auto &w = screen_windows[k];
            if (ctx.time_au >= w.t_start_au && ctx.time_au <= w.t_end_au) {
                acc_screen_window[k].accumulate(*ctx.electrons, Cfg::DT_AU);
            }
        }

        // Paper-window accumulator (single global bracket, all screens).
        if (ctx.time_au >= Cfg::T1_AU && ctx.time_au <= Cfg::T2_AU) {
            for (auto &a : acc_paper) a.accumulate(*ctx.electrons, Cfg::DT_AU);
        }

        if (ctx.step % Cfg::SCREEN_SNAP_EVERY == 0) {
            // Flat: results/raw/screens/instantaneous/screen_NN_tXXXXXX.dat
            const std::string base = coronene::results::screens_instantaneous_dir();
            for (int k = 0; k < coronene::layout::N_SCREENS; ++k) {
                auto slice = acc_full[k].screen().extract(*ctx.electrons);
                const std::string filename =
                    base + "/" + coronene::layout::screen_label(k)
                    + "_t" + coronene::layout::zero_pad6(ctx.step) + ".dat";
                acc_full[k].screen().save(slice, ctx.time_au, filename);
            }
        }
    });

    real_time::propagate(
        ions, electrons,
        [&](auto const &data) {
            rt_dens.step(data);
            rt_obs.step(data);
        },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(Cfg::N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .observables_current()
            .observables_dipole());

    // ----- Save accumulated screens ---------------------------------------
    auto screens_total_dir         = coronene::results::screens_total_dir();
    auto screens_time_windowed_dir = coronene::results::screens_time_windowed_dir();
    {
        // total/ — full-time accumulator, one .dat per screen.
        for (int k = 0; k < coronene::layout::N_SCREENS; ++k) {
            acc_full[k].save(screens_total_dir + "/"
                             + coronene::layout::screen_label(k) + ".dat");
        }

        // time_windowed/ — physics window per screen plus the paper window.
        const int total_steps = Cfg::N_STEPS;
        const int step_paper_lo =
            static_cast<int>(std::round(Cfg::T1_AU / Cfg::DT_AU));
        const int step_paper_hi =
            static_cast<int>(std::round(Cfg::T2_AU / Cfg::DT_AU));

        for (int k = 0; k < coronene::layout::N_SCREENS; ++k) {
            const auto &w = screen_windows[k];
            int step_lo = static_cast<int>(std::round(w.t_start_au / Cfg::DT_AU));
            int step_hi = static_cast<int>(std::round(w.t_end_au   / Cfg::DT_AU));
            if (step_lo < 0)             step_lo = 0;
            if (step_hi > total_steps)   step_hi = total_steps;
            const std::string tag = w.is_back ? "back" : "forward";
            std::ostringstream fn;
            fn << screens_time_windowed_dir << "/"
               << coronene::layout::screen_label(k)
               << "_t" << coronene::layout::zero_pad6(step_lo)
               << "_to_t" << coronene::layout::zero_pad6(step_hi)
               << "_" << tag << ".dat";
            acc_screen_window[k].save(fn.str());

            // Paper window: same filename pattern + '_paper' tag, every
            // screen sees the same global [T1_AU, T2_AU] bracket.
            std::ostringstream fnp;
            fnp << screens_time_windowed_dir << "/"
                << coronene::layout::screen_label(k)
                << "_t" << coronene::layout::zero_pad6(step_paper_lo)
                << "_to_t" << coronene::layout::zero_pad6(step_paper_hi)
                << "_paper.dat";
            acc_paper[k].save(fnp.str());
        }
    }

    // ----- Final run_summary.txt ------------------------------------------
    auto t_wallclock_end = std::chrono::steady_clock::now();
    double wall_seconds =
        std::chrono::duration<double>(t_wallclock_end - t_wallclock_start).count();
    if (electrons.root()) {
        std::ofstream s(coronene::results::run_summary_path());
        s << std::setprecision(16);
        s << "RUN SUMMARY\n"
          << "===========\n\n"
          << "1. Run identity\n"
          << "---------------\n"
          << "run_name        = " << run_name << "\n"
          << "run_type        = wave-packet RT-LEED (TDDFT, ALDA)\n"
          << "date_finished   = " << iso_now() << "\n"
          << "wall_time_s     = " << wall_seconds << "\n"
          << "executable      = run.cpp built via inq-run\n"
          << "geometry_file   = " << geometry_xyz_path << "\n"
          << "checkpoint_dir  = " << gs_checkpoint_dir << "\n\n"
          << "3. System configuration\n"
          << "-----------------------\n"
          << "cell_bohr       = " << Cfg::LX_BOHR << " " << Cfg::LY_BOHR << " "
                                  << Cfg::LZ_BOHR << "\n"
          << "boundary        = finite (centred [-L/2,+L/2])\n"
          << "n_ions          = " << ions.size() << "\n"
          << "n_electrons     = " << n_electrons << "\n"
          << "n_occupied      = " << n_occupied << "\n"
          << "extra_states    = " << Cfg::EXTRA_STATES << "\n"
          << "wp_state_index  = " << wp_idx << "\n"
          << "cutoff_ha       = " << Cfg::CUTOFF_HA << "\n"
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
          << "wp_direction    = -z\n"
          << "wp_occupation   = 1.0\n"
          << "orthogonalised  = " << (report.orthogonalised ? "yes" : "no") << "\n"
          << "norm_before     = " << report.norm_before << "\n"
          << "norm_after      = " << report.norm_after << "\n"
          << "max_overlap     = " << report.max_overlap << "\n"
          << "passed_tol      = " << (report.passed_tolerance ? "yes" : "no") << "\n\n"
          << "6. Real-time configuration\n"
          << "--------------------------\n"
          << "rt_num_steps    = " << Cfg::N_STEPS << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "total_time_au   = " << Cfg::DT_AU * Cfg::N_STEPS << "\n"
          << "write_every     = " << Cfg::WRITE_EVERY << "\n"
          << "screen_snap_every = " << Cfg::SCREEN_SNAP_EVERY << "\n\n"
          << "7. Screen configuration\n"
          << "-----------------------\n"
          << "n_screens       = " << coronene::layout::N_SCREENS << "\n"
          << "screen_orientation = z (constant-z planes)\n"
          << "screen_windows  = per-screen physics-derived (forward / back)\n"
          << "paper_window_au = " << Cfg::T1_AU << " " << Cfg::T2_AU << "\n";
        for (int k = 0; k < coronene::layout::N_SCREENS; ++k) {
            s << "screen_z[" << k << "] = " << screen_z[k] << "\n";
        }
        s << "\n"
          << "9. End-of-run diagnostics\n"
          << "-------------------------\n"
          << "run_completed   = true\n"
          << "final_time_au   = " << Cfg::DT_AU * Cfg::N_STEPS << "\n"
          << "vti_format      = binary\n"
          << "raw_emitted     = no\n"
          << "wp_overlap_mode = wp_only (snapshot_wp_only)\n";
    }

    std::cout << "Done. Wall time " << wall_seconds << " s. Output in results/\n";
    return 0;
}

}  // namespace coronene::run_template
