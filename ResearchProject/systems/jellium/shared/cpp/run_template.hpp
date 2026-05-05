// ============================================================================
// shared/cpp/run_template.hpp  (jellium)
//
// The propagation body shared by every jellium run.cpp. Direct port of
// ResearchProject/systems/coronene/shared/cpp/run_template.hpp into namespace
// `jellium::run_template`. Each run.cpp is a thin wrapper that:
//   1. Includes its variant config header (or just base.hpp);
//   2. Calls jellium::run_template::run_propagation<Cfg>(run_name, gs_path).
//
// Differences vs the coronene template:
//   * Cell is `cubic(L_BOHR).periodic()` with empty `systems::ions(cell)`
//     (jellium has no atoms).
//   * `systems::electrons` is constructed with `extra_electrons(N_ELECTRONS)`,
//     `extra_states(EXTRA_STATES)`, `temperature(TEMPERATURE_EV)`,
//     `spacing(SPACING_BOHR)` — there is no `cutoff` for jellium. Γ-only
//     k-points (`input::kpoints::gamma()`).
//   * `run_propagation(...)` takes only `run_name` and `gs_checkpoint_dir`
//     (no geometry path).
//   * The orbital overlap observer calls `snapshot_wp_only()` every 10
//     propagator steps. Records O_{i, wp}(t) = |<psi_i^GS|psi_wp(t)>|^2 for
//     i ∈ [0, n_ref). The full n_ref x n_ref matrix was originally used at
//     WRITE_EVERY cadence but proved prohibitively expensive at high
//     N_states (n_ref ~ 257 ⇒ ~66k GPU reductions per call, ~6 min wall).
//     If the full matrix is needed for diagnostics, take it explicitly
//     (e.g. start/mid/end) outside the propagation hot path.
//   * `run_summary.txt` reports `n_ions = 0` and `geometry_file = (none)`.
//
// All other behaviour (raw/analysis split, GS density + orbital writes, WP
// injection + reports, three RT density writers, observables CSV every step,
// three screen accumulators per screen, screen instantaneous snapshots,
// final run_summary) is identical to coronene.
//
// FOLLOW-UP (tracked in docs/plans/jellium_reorg.md §14):
//   * Per-screen physics windows (`compute_screen_window`) for periodic
//     boundaries — currently placeholder.
//   * Compile-time `compute_n_steps` for jellium loop-back.
// ============================================================================
#pragma once

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
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include "eigenvalues_writer.hpp"
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

namespace jellium::run_template {

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

// Cfg duck-typed: exposes L_BOHR, LX_BOHR, LY_BOHR, LZ_BOHR (all = L_BOHR for
// cubic), N_ELECTRONS, EXTRA_STATES, SPACING_BOHR, TEMPERATURE_EV,
// WP_CX_BOHR, WP_CY_BOHR, WP_CZ_BOHR, WP_SIGMA_BOHR, WP_KX, WP_KY, WP_KZ,
// WP_K0, WP_KZ_MAGNITUDE, WP_EKIN_EV, DT_AU, N_STEPS, WRITE_EVERY,
// SCREEN_SNAP_EVERY, T1_AU, T2_AU, WP_ENVELOPE_SIGMAS, CUTOFF_HA (stub).
template <class Cfg>
int run_propagation(std::string const &run_name,
                    std::string const &gs_checkpoint_dir) {
    using namespace inq;
    using namespace inq::magnitude;

    auto t_wallclock_start = std::chrono::steady_clock::now();

    std::cout << "\n=== " << run_name << " ===\n"
              << "  cell = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
              << "  N_electrons = " << Cfg::N_ELECTRONS << "\n"
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

    // ----- Cell + (empty) ions --------------------------------------------
    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    std::cout << "  Atoms: " << ions.size() << " (jellium — no nuclei)\n";

    // ----- Electrons (skeleton) + load checkpoint -------------------------
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(gs_checkpoint_dir);
    std::cout << "  Loaded GS from " << gs_checkpoint_dir << "\n";

    // Mirror eigenvalues + occupations from the checkpoint into this run's
    // results tree (silent no-op if the checkpoint predates the writer).
    jellium::eigenvalues::copy_from_checkpoint(
        gs_checkpoint_dir, "results/raw/observables/eigenvalues");

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    const int n_occupied  = n_electrons / 2;
    std::cout << "  num_states = " << n_states
              << "  num_electrons = " << n_electrons
              << "  n_occupied = " << n_occupied << "\n";

    // ----- Stub run_summary.txt so a crashed run still has metadata -------
    {
        std::ofstream s(jellium::results::run_summary_path());
        s << "RUN SUMMARY (stub - written at start)\n"
          << "=====================================\n"
          << "run_name        = " << run_name << "\n"
          << "date_started    = " << iso_now() << "\n"
          << "geometry_file   = (none, jellium)\n"
          << "checkpoint_dir  = " << gs_checkpoint_dir << "\n"
          << "cell_bohr       = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
          << "boundary        = periodic\n"
          << "xc              = LDA (ALDA in TDDFT)\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "temperature_ev  = " << Cfg::TEMPERATURE_EV << "\n"
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

    // ----- Ground-state artefacts (before WP injection) -------------------
    inqkit::io::RealField3DLayout vti_layout{
        .field_name  = "density",
        .include_meta = false,
        .emit_raw    = false,
        .emit_vti    = true,
        .vti_format  = inqkit::io::VTIWriteOptions::Format::binary,
    };
    {
        // GS total density (jellium background)
        inqkit::io::RealField3DWriter gs_sys(
            jellium::results::vti_density_gs_system_dir(), vti_layout,
            {.overwrite = true});
        gs_sys.write(inqkit::fields::density::total(electrons),
                     "density_gs_system");

        // GS orbital densities (one VTI per occupied orbital)
        inqkit::io::RealField3DWriter gs_orb(
            jellium::results::vti_density_gs_orbitals_dir(), vti_layout,
            {.overwrite = true});
        for (int i = 0; i < n_occupied; ++i) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "orbital_%04d", i);
            gs_orb.write(inqkit::fields::density::orbital(electrons, i),
                         std::string(buf));
        }
    }

    // ----- WP injection ---------------------------------------------------
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

    // ----- Wavepacket artefacts ------------------------------------------
    {
        std::ofstream f(jellium::results::wp_config_path());
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
        std::ofstream f(jellium::results::wp_injection_report_path());
        f << std::setprecision(16);
        f << "norm_before      = " << report.norm_before << "\n"
          << "norm_after       = " << report.norm_after << "\n"
          << "max_overlap      = " << report.max_overlap << "\n"
          << "passed_tolerance = " << (report.passed_tolerance ? "yes" : "no") << "\n";
    }
    {
        // Initial WP density (real) and wavefunction (complex)
        inqkit::io::RealField3DWriter wp_d(
            jellium::results::wp_density_initial_dir(), vti_layout,
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
            jellium::results::wp_wavefunction_initial_dir(), cvti_layout,
            {.overwrite = true});
        wp_psi.write(inqkit::fields::orbital::wavefunction(electrons, wp_idx),
                     "wavefunction_wp_initial");
    }

    // ----- Real-time outputs setup ---------------------------------------
    // Full GS↔evolved overlap matrix observer (the new observable).
    inqkit::observables::OrbitalOverlapMatrix overlap_obs(
        electrons, wp_idx, jellium::results::overlap_dir());

    inqkit::io::RealField3DWriter total_wr(
        jellium::results::vti_density_total_dir(), vti_layout, {.overwrite = true});
    inqkit::io::RealField3DWriter system_wr(
        jellium::results::vti_density_system_dir(), vti_layout, {.overwrite = true});
    inqkit::io::RealField3DWriter wp_wr(
        jellium::results::vti_density_wp_dir(), vti_layout, {.overwrite = true});

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
    sel.cod_x = sel.cod_y = sel.cod_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(
        jellium::results::observables_csv_path(), sel);
    obs_writer.write_header();

    // Per-state energies <phi_i|H|phi_i>(t) and variance, and one-body
    // momentum distribution n(|k|,t) (total + WP). Both are computed via
    // direct access to the propagator's `viewables` (which exposes the
    // Hamiltonian via the local extension to inq/src/real_time/viewables.hpp).
    inqkit::observables::StateEnergyWriter state_energy_wr(
        jellium::results::state_energies_csv_path(), /*emit_variance=*/true);

    // Occupations dump every step. Cheap audit — also enables the
    // energy_balance postprocess to compute occupation-weighted ΔE_bath
    // without joining with a separate GS-occupations CSV.
    inqkit::observables::OccupationsWriter occupations_wr(
        jellium::results::occupations_csv_path());
    // Momentum binning is dominated by the FFT + per-state binning loop
    // (~60^3 * n_states ops per snapshot). Cadence here is decoupled from
    // and coarser than density VTI cadence — 10x WRITE_EVERY by default,
    // i.e. one snapshot per 20 propagation steps.
    inqkit::observables::MomentumDistribution momentum_dist(
        jellium::results::momentum_distribution_csv_path(),
        wp_idx,
        Cfg::L_BOHR,
        {.n_bins = 64, .k_max_bohr_inv = 0.0,
         .write_every = 10 * Cfg::WRITE_EVERY});

    // Density-fluctuation observable (3-tier: raw delta VTI, coarse VTI,
    // integrated L2). Reference density is captured at t=0 below.
    inqkit::observables::DensityDelta density_delta(
        jellium::results::vti_density_delta_dir(),
        jellium::results::vti_density_delta_coarse_dir(),
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0});
    // density_delta captures its reference on the first snapshot() call —
    // see header for the rationale (avoids stale-density issues when the
    // propagator rebuilds the density on its first iteration).

    // ----- Screens: 20 positions, three accumulators per screen ----------
    auto screen_z = jellium::layout::screen_z_positions(Cfg::L_BOHR);
    {
        std::cout << "  Screen z-positions (Bohr):";
        for (auto z : screen_z) std::cout << " " << z;
        std::cout << "\n";
    }

    using inqkit::screens::LeedPatternAccumulator;
    using inqkit::screens::PlaneScreen;

    std::array<LeedPatternAccumulator, jellium::layout::N_SCREENS> acc_full;
    std::array<LeedPatternAccumulator, jellium::layout::N_SCREENS> acc_screen_window;
    std::array<LeedPatternAccumulator, jellium::layout::N_SCREENS> acc_paper;

    std::array<jellium::layout::ScreenWindow, jellium::layout::N_SCREENS>
        screen_windows;

    for (int k = 0; k < jellium::layout::N_SCREENS; ++k) {
        auto label = jellium::layout::screen_label(k);
        acc_full[k]          = LeedPatternAccumulator(PlaneScreen{screen_z[k], label});
        acc_screen_window[k] = LeedPatternAccumulator(PlaneScreen{screen_z[k], label});
        acc_paper[k]         = LeedPatternAccumulator(PlaneScreen{screen_z[k], label});
        screen_windows[k] = jellium::layout::compute_screen_window(
            screen_z[k], Cfg::WP_CZ_BOHR, Cfg::WP_SIGMA_BOHR, Cfg::WP_K0,
            Cfg::L_BOHR, Cfg::DT_AU * Cfg::N_STEPS,
            Cfg::WP_ENVELOPE_SIGMAS);
    }

    // Persist the configuration so the postprocess can reproduce filenames.
    {
        std::ofstream f(jellium::results::screens_config_path());
        f << "screen_index,z_bohr,label,window_kind,t_start_au,t_end_au\n";
        for (int k = 0; k < jellium::layout::N_SCREENS; ++k) {
            f << k << "," << std::fixed << std::setprecision(6)
              << screen_z[k] << "," << jellium::layout::screen_label(k)
              << "," << (screen_windows[k].is_back ? "back" : "forward")
              << "," << screen_windows[k].t_start_au
              << "," << screen_windows[k].t_end_au << "\n";
        }
    }
    {
        std::ofstream f(jellium::results::screens_window_ranges_path());
        f << "kind,screen_index,t_start_au,t_end_au\n";
        for (int k = 0; k < jellium::layout::N_SCREENS; ++k) {
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

    // ----- Real-time sessions
    //
    // We now drive density VTI writes, CoD/dn observables, screens, and
    // overlap from a single session at WRITE_EVERY cadence. The per-step
    // 60^3 host extraction (fields::density::total/orbital) was the
    // dominant cost; doing it once per WRITE_EVERY rather than every step
    // halves runtime. Side effect: observables.csv has ~N_STEPS/WRITE_EVERY
    // rows instead of N_STEPS, which is still adequate for FFT spectra.
    //
    // Screen accumulators take dt = WRITE_EVERY * Cfg::DT_AU because they
    // now run only every WRITE_EVERY steps (used to be every step).
    const double SCREEN_DT_AU = Cfg::DT_AU * Cfg::WRITE_EVERY;
    inqkit::RealTimeSession rt_obs(ions, electrons, Cfg::WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const &ctx) {
        // Single round of density extraction, reused below.
        // NOTE on naming: `sys_f` is `electrons.density()`, which is the
        // full DFT density Σᵢ fᵢ |ψᵢ|² and ALREADY includes the WP slot
        // (occupation 1.0). The `wp_f` is the WP orbital density alone.
        // The legacy `total = sys + wp` line double-counted the WP; we
        // keep it for VTI compatibility (the system_wr/wp_wr/total_wr
        // contract is consumed by existing postprocess) but use the
        // un-doubled `sys_f` for the density_delta observable.
        auto sys_f   = inqkit::fields::density::total(*ctx.electrons);
        auto wp_f    = inqkit::fields::density::orbital(*ctx.electrons, wp_idx);
        auto total_f = add_real_fields(sys_f, wp_f);

        system_wr.write(sys_f,   ctx.time_au, ctx.step);
        wp_wr.write(    wp_f,    ctx.time_au, ctx.step);
        total_wr.write( total_f, ctx.time_au, ctx.step);

        const double l2 =
            density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        const auto cod = inqkit::observables::center_of_density(wp_f);

        inqkit::StepContext ctx_out = ctx;
        ctx_out.wp_center = {cod.x_bohr, cod.y_bohr, cod.z_bohr};
        ctx_out.density_l2 = l2;
        obs_writer.append(ctx_out);

        // WP-only overlap with the GS basis: snapshot the column
        // O_{i, wp}(t) = |<psi_i^GS | psi_wp(t)>|^2 instead of the full
        // n_ref x n_ref matrix. With n_ref ~ 257 (high-density jellium),
        // the full matrix is ~66k GPU reductions per call (~5-6 min wall);
        // the WP-only column is n_ref+1 reductions (~ms wall).
        // Cadence: every 10 propagator steps (= 5 x WRITE_EVERY) to keep
        // the file series compact.
        if (ctx.step % 10 == 0) {
            overlap_obs.snapshot_wp_only(*ctx.electrons,
                                         ctx.time_au, ctx.step);
        }

        for (auto &a : acc_full) a.accumulate(*ctx.electrons, SCREEN_DT_AU);

        for (int k = 0; k < jellium::layout::N_SCREENS; ++k) {
            const auto &w = screen_windows[k];
            if (ctx.time_au >= w.t_start_au && ctx.time_au <= w.t_end_au) {
                acc_screen_window[k].accumulate(*ctx.electrons, SCREEN_DT_AU);
            }
        }

        if (ctx.time_au >= Cfg::T1_AU && ctx.time_au <= Cfg::T2_AU) {
            for (auto &a : acc_paper) a.accumulate(*ctx.electrons, SCREEN_DT_AU);
        }

        if (ctx.step % Cfg::SCREEN_SNAP_EVERY == 0) {
            const std::string base = jellium::results::screens_instantaneous_dir();
            for (int k = 0; k < jellium::layout::N_SCREENS; ++k) {
                auto slice = acc_full[k].screen().extract(*ctx.electrons);
                const std::string filename =
                    base + "/" + jellium::layout::screen_label(k)
                    + "_t" + jellium::layout::zero_pad6(ctx.step) + ".dat";
                acc_full[k].screen().save(slice, ctx.time_au, filename);
            }
        }
    });

    real_time::propagate(
        ions, electrons,
        [&](auto const &data) {
            rt_obs.step(data);
            // State-energy snapshot does 2 ham(phi) applies + variance,
            // ~equivalent to a few extra propagation steps. Run at 5x
            // WRITE_EVERY (every 10 steps for base) — coarse enough to be
            // cheap, dense enough for a nice bar-plot GIF.
            if (data.iter() % (5 * Cfg::WRITE_EVERY) == 0) {
                state_energy_wr.snapshot(data);
                occupations_wr.snapshot(data);
            }
            momentum_dist.maybe_accumulate(data);
        },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(Cfg::N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .observables_current()
            .observables_dipole());

    // ----- Save accumulated screens --------------------------------------
    auto screens_total_dir         = jellium::results::screens_total_dir();
    auto screens_time_windowed_dir = jellium::results::screens_time_windowed_dir();
    {
        for (int k = 0; k < jellium::layout::N_SCREENS; ++k) {
            acc_full[k].save(screens_total_dir + "/"
                             + jellium::layout::screen_label(k) + ".dat");
        }

        const int total_steps = Cfg::N_STEPS;
        const int step_paper_lo =
            static_cast<int>(std::round(Cfg::T1_AU / Cfg::DT_AU));
        const int step_paper_hi =
            static_cast<int>(std::round(Cfg::T2_AU / Cfg::DT_AU));

        for (int k = 0; k < jellium::layout::N_SCREENS; ++k) {
            const auto &w = screen_windows[k];
            int step_lo = static_cast<int>(std::round(w.t_start_au / Cfg::DT_AU));
            int step_hi = static_cast<int>(std::round(w.t_end_au   / Cfg::DT_AU));
            if (step_lo < 0)             step_lo = 0;
            if (step_hi > total_steps)   step_hi = total_steps;
            const std::string tag = w.is_back ? "back" : "forward";
            std::ostringstream fn;
            fn << screens_time_windowed_dir << "/"
               << jellium::layout::screen_label(k)
               << "_t" << jellium::layout::zero_pad6(step_lo)
               << "_to_t" << jellium::layout::zero_pad6(step_hi)
               << "_" << tag << ".dat";
            acc_screen_window[k].save(fn.str());

            std::ostringstream fnp;
            fnp << screens_time_windowed_dir << "/"
                << jellium::layout::screen_label(k)
                << "_t" << jellium::layout::zero_pad6(step_paper_lo)
                << "_to_t" << jellium::layout::zero_pad6(step_paper_hi)
                << "_paper.dat";
            acc_paper[k].save(fnp.str());
        }
    }

    // ----- Final run_summary.txt -----------------------------------------
    auto t_wallclock_end = std::chrono::steady_clock::now();
    double wall_seconds =
        std::chrono::duration<double>(t_wallclock_end - t_wallclock_start).count();
    if (electrons.root()) {
        std::ofstream s(jellium::results::run_summary_path());
        s << std::setprecision(16);
        s << "RUN SUMMARY\n"
          << "===========\n\n"
          << "1. Run identity\n"
          << "---------------\n"
          << "run_name        = " << run_name << "\n"
          << "run_type        = wave-packet RT-LEED on jellium (TDDFT, ALDA)\n"
          << "date_finished   = " << iso_now() << "\n"
          << "wall_time_s     = " << wall_seconds << "\n"
          << "executable      = run.cpp built via inq-run\n"
          << "geometry_file   = (none, jellium)\n"
          << "checkpoint_dir  = " << gs_checkpoint_dir << "\n\n"
          << "3. System configuration\n"
          << "-----------------------\n"
          << "cell_bohr       = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
          << "boundary        = periodic\n"
          << "n_ions          = 0\n"
          << "n_electrons     = " << n_electrons << "\n"
          << "n_occupied      = " << n_occupied << "\n"
          << "extra_states    = " << Cfg::EXTRA_STATES << "\n"
          << "wp_state_index  = " << wp_idx << "\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "temperature_ev  = " << Cfg::TEMPERATURE_EV << "\n"
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
          << "wp_direction    = +z (jellium convention)\n"
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
          << "n_screens       = " << jellium::layout::N_SCREENS << "\n"
          << "screen_orientation = z (constant-z planes)\n"
          << "screen_windows  = per-screen (placeholder, see leed_screen_layout.hpp)\n"
          << "paper_window_au = " << Cfg::T1_AU << " " << Cfg::T2_AU << "\n";
        for (int k = 0; k < jellium::layout::N_SCREENS; ++k) {
            s << "screen_z[" << k << "] = " << screen_z[k] << "\n";
        }
        s << "\n"
          << "9. End-of-run diagnostics\n"
          << "-------------------------\n"
          << "run_completed   = true\n"
          << "final_time_au   = " << Cfg::DT_AU * Cfg::N_STEPS << "\n"
          << "vti_format      = binary\n"
          << "raw_emitted     = no\n"
          << "wp_overlap_mode = full_matrix (snapshot)\n";
    }

    std::cout << "Done. Wall time " << wall_seconds << " s. Output in results/\n";
    return 0;
}

}  // namespace jellium::run_template
