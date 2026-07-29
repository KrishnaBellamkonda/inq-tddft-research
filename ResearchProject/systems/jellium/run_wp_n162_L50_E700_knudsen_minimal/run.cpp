// ============================================================================
// run_wp_n162_L50_E700_knudsen_minimal/run.cpp
//
// Stripped-down attempt at the Run-8 dx=0.30 CUDA blocker (2026-05-20).
//
// The original run_wp_n162_L50_E700_knudsen failed at propagate step 0
// with "CUDA: illegal memory access" — same blocker that killed the
// legacy E=1500 WP. The dx=0.30 grid puts orbital memory near 24 GB.
//
// This minimal variant removes every heavy callback to see whether the
// bare propagation + per-state moments observables fit:
//   STRIPPED:  full overlap matrix snapshot, per-step overlap (wp_only),
//              proxy overlaps, DensityDelta, density VTI series,
//              MomentumDistribution histogram, t=0 density VTI snapshot.
//   KEPT:      WPMomentumStats (small CSV, GPU reduce), WPRealSpaceStats
//              (small CSV, GPU reduce), observables.csv,
//              state_energies, occupations.
//   WP INJECTION orthogonalisation against occupied: KEPT (correctness).
//
// If this works, the Knudsen velocity sweep is recoverable at the cost
// of post-run reconstruction of overlap / density quantities from a
// recomputed reduced output set.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/shells.hpp>

#include "../shared/configs/knudsen_sweep_L50_cubic.hpp"
#include "../shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::KnudsenSweep_L50_cubic_WP_E700;

static std::string iso_now() {
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

int main() {
    auto t0 = std::chrono::steady_clock::now();
    const std::string RUN_NAME = "run_wp_n162_L50_E700_knudsen_minimal";
    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p30";

    std::cout << "\n=== " << RUN_NAME << " (Run-8 minimal dx=0.30 test) ===\n"
              << "  cell    = " << Cfg::L_BOHR << "^3 Bohr, dx=" << Cfg::SPACING_BOHR << "\n"
              << "  E       = " << Cfg::WP_EKIN_EV << " eV (k0=" << Cfg::WP_K0 << ")\n"
              << "  launch  = (0, 0, " << Cfg::WP_CZ_BOHR << ")\n"
              << "  N_STEPS = " << Cfg::N_STEPS << ", dt=" << Cfg::DT_AU << "\n";

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

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    std::cout << "  loaded GS: " << n_states << " states, "
              << n_electrons << " electrons\n";

    std::filesystem::create_directories("results/raw/observables");

    auto wp = inqkit::WavePacket{}
                  .center(Cfg::WP_CX_BOHR, Cfg::WP_CY_BOHR, Cfg::WP_CZ_BOHR)
                  .sigma(Cfg::WP_SIGMA_BOHR)
                  .k0(Cfg::WP_KX, Cfg::WP_KY, Cfg::WP_KZ)
                  .orthogonalise_against_occupied(electrons);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected at state " << wp_idx
              << "  norm=" << report.norm_after << "\n\n";

    // ----- MINIMAL observable set (no overlap, no VTI density) ---------
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_z = true;
    sel.dipole_z = true;
    sel.density_l2 = false;       // skipped — no DensityDelta
    inqkit::io::ObservablesWriter obs_writer(
        "results/raw/observables/observables.csv", sel);
    obs_writer.write_header();

    inqkit::observables::StateEnergyWriter state_energy_wr(
        "results/raw/observables/state_energies.csv", /*emit_variance=*/true);

    inqkit::observables::OccupationsWriter occupations_wr(
        "results/raw/observables/occupations_vs_time.csv");

    inqkit::observables::WPMomentumStats wp_momentum_stats(
        "results/raw/observables/wp_momentum_stats.csv",
        wp_idx, {.write_every = Cfg::WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        "results/raw/observables/wp_real_space_stats.csv",
        wp_idx, {.write_every = Cfg::WRITE_EVERY});

    // run_summary stub
    {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "run_name        = " << RUN_NAME << "\n"
          << "run_type        = WP jellium TDDFT (Run-8 minimal-callback dx=0.30 test)\n"
          << "date_started    = " << iso_now() << "\n"
          << "cell_bohr       = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "n_electrons     = " << n_electrons << "\n"
          << "n_states        = " << n_states << "\n"
          << "wp_state_index  = " << wp_idx << "\n"
          << "wp_ekin_ev      = " << Cfg::WP_EKIN_EV << "\n"
          << "wp_sigma_bohr   = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv  = " << Cfg::WP_KX << " " << Cfg::WP_KY << " " << Cfg::WP_KZ << "\n"
          << "wp_center_bohr  = " << Cfg::WP_CX_BOHR << " " << Cfg::WP_CY_BOHR
                                  << " " << Cfg::WP_CZ_BOHR << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "n_steps         = " << Cfg::N_STEPS << "\n"
          << "norm_after      = " << report.norm_after << "\n"
          << "run_completed   = false\n";
    }

    inqkit::RealTimeSession rt_obs(ions, electrons, Cfg::WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const &ctx) {
        obs_writer.append(ctx);
    });

    real_time::propagate(
        ions, electrons,
        [&](auto const &data) {
            rt_obs.step(data);
            wp_momentum_stats.maybe_accumulate(data);
            wp_real_space_stats.maybe_accumulate(data);
            if (data.iter() % (5 * Cfg::WRITE_EVERY) == 0) {
                state_energy_wr.snapshot(data);
                occupations_wr.snapshot(data);
            }
        },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(Cfg::N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .observables_current()
            .observables_dipole());

    double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    {
        std::ofstream s("results/run_summary.txt", std::ios::app);
        s << "wall_time_s     = " << wall << "\n"
          << "run_completed   = true\n";
    }
    std::cout << "Done. Wall " << wall << " s\n";
    return 0;
}
