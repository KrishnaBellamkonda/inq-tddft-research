// Stage 2 of 2-GPU MPI Run-8 attempt (2026-05-20).
//
// Loads the WP-injected checkpoint (made by Stage 1) under mpirun -np 2 and
// runs real_time::propagate. INQ's parallelization can split states (101)
// or basis (168 slices) across the 2 ranks, effectively doubling memory
// to ~48 GB. Hypothesis: this fixes the single-GPU memory deadlock at
// dx=0.30.
//
// Observables stripped to a minimal set (no inqkit::WavePacket per-step
// calls, no overlap, no VTI density) because:
//   (a) some inqkit observables are single-rank only
//   (b) we want to isolate "does dx=0.30 propagation work at all on 2 GPUs"
//       from "do the heavy callbacks fit". Method-B Knudsen via
//       WPMomentumStats still works since it does its own MPI all_reduce.

#include <inq/inq.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>

#include "../shared/configs/knudsen_sweep_L50_cubic.hpp"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::KnudsenSweep_L50_cubic_WP_E700;

int main() {
    auto t0 = std::chrono::steady_clock::now();
    const std::string IN_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p30_E700_wp_injected";

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
    electrons.load(IN_DIR);

    std::filesystem::create_directories("results/raw/observables");

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_z = sel.dipole_z = true;
    inqkit::io::ObservablesWriter obs_writer(
        "results/raw/observables/observables.csv", sel);
    obs_writer.write_header();

    const int n_states = electrons.states().num_states();
    const int wp_slot = n_states - 1;
    inqkit::observables::WPMomentumStats wp_momentum_stats(
        "results/raw/observables/wp_momentum_stats.csv",
        wp_slot, {.write_every = Cfg::WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        "results/raw/observables/wp_real_space_stats.csv",
        wp_slot, {.write_every = Cfg::WRITE_EVERY});

    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt");
        s << "run_name        = run_wp_n162_L50_E700_mpi_propagate\n"
          << "run_type        = Run-8 MPI 2-GPU attempt\n"
          << "cell_bohr       = " << Cfg::L_BOHR << "^3\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "wp_state_index  = " << wp_slot << "\n"
          << "wp_ekin_ev      = " << Cfg::WP_EKIN_EV << "\n"
          << "wp_sigma_bohr   = " << Cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv  = " << Cfg::WP_KX << " " << Cfg::WP_KY << " " << Cfg::WP_KZ << "\n"
          << "wp_center_bohr  = " << Cfg::WP_CX_BOHR << " " << Cfg::WP_CY_BOHR
                                  << " " << Cfg::WP_CZ_BOHR << "\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "n_steps         = " << Cfg::N_STEPS << "\n"
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
        },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(Cfg::N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .observables_current()
            .observables_dipole());

    double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt", std::ios::app);
        s << "wall_time_s     = " << wall << "\n"
          << "run_completed   = true\n";
    }
    if (electrons.root()) std::cout << "Done. Wall " << wall << " s\n";
    return 0;
}
