// ============================================================================
// run_free_wp_L50_E100/run.cpp — Run-2 of the 2026-05-21 meeting campaign.
//
// FREE-SPACE wave-packet propagation in an empty 50^3 Bohr finite box,
// non-interacting theory (no Hartree, no XC). Same WP injection geometry
// as Run-3 (the jellium E=100 pair) so the σ_r²(t), σ_p²(t), <z>(t)
// curves can be overlaid:
//
//   sigma   = 5 Bohr               (injector parameter)
//   k0_z    = 2.7111 Bohr^-1       (E = 100 eV)
//   launch  = (0, 0, -5)           (boundary_rule launch_z(5, 50) = -5)
//   stop    = (0, 0, +20)          (boundary_rule stop_z = +20)
//   N_STEPS = 462                  (boundary_rule.n_steps_for at k0)
//   WRITE_EVERY = 2                (~231 frames)
//
// Purpose: validates the inqkit::WavePacket injector against analytic
// free-particle Gaussian spreading
//
//     sigma_r(t) = sigma_r(0) * sqrt(1 + (t / (2 * sigma_r(0)^2))^2)
//
// (atomic units, m_e = 1). With injector sigma=5, density sigma_r(0) =
// sigma/sqrt(2) = 3.5355 Bohr; growth at t=9.22 a.u. is
//   sigma_r(9.22) = 3.5355 * sqrt(1 + (9.22 / 25)^2) = 3.5355 * 1.066
//                 = 3.768 Bohr  (=> 6.6 % expansion).
//
// Pairs with Run-2b (Python Schrödinger toy at the same parameters) and
// Run-3 (jellium pair). All three traces of sigma_r²(t) go on one plot
// at meeting time — see plan §"Per-family motivations" Run-2+2b+3.
//
// Cfg: declared inline because this is a one-off free-space run that
// doesn't share Base_N162_L50_E1p5's jellium-specific defaults.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include "../shared/configs/boundary_rule.hpp"

#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;

static constexpr double HA_TO_EV  = 27.21138625;
static constexpr double FS_TO_AU  = 41.34137333518;

namespace cfg {
    constexpr double L_BOHR        = 50.0;
    constexpr double SPACING_BOHR  = 0.40;
    constexpr double DT_AU         = 0.020;

    constexpr double WP_EKIN_EV    = 100.0;
    constexpr double WP_SIGMA_BOHR = 5.0;
    // k0 from kinetic energy in eV: k = sqrt(2 * E_Ha) = sqrt(2 * E_eV / 27.211)
    constexpr double WP_K0         =
        2.711063340102429;                       // sqrt(2 * 100 / 27.21138625)
    constexpr double WP_KX = 0.0, WP_KY = 0.0, WP_KZ = WP_K0;

    constexpr double WP_CX_BOHR = 0.0;
    constexpr double WP_CY_BOHR = 0.0;
    constexpr double WP_CZ_BOHR =
        jellium::config::boundary::launch_z(WP_SIGMA_BOHR, L_BOHR);   // -5
    constexpr int N_STEPS =
        jellium::config::boundary::n_steps_for(WP_SIGMA_BOHR, L_BOHR,
                                               WP_K0, DT_AU);
    constexpr int WRITE_EVERY =
        jellium::config::boundary::write_every_for(N_STEPS);
}

static std::string iso_now() {
    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

int main() {
    auto t_start = std::chrono::steady_clock::now();
    const std::string RUN_NAME = "run_free_wp_L50_E100";

    std::cout << "\n=== " << RUN_NAME << " (Run-2: free-space WP @ 100 eV) ===\n"
              << "  cell    = " << cfg::L_BOHR << "^3 Bohr (finite, NON-interacting)\n"
              << "  spacing = " << cfg::SPACING_BOHR << " Bohr\n"
              << "  WP      = sigma " << cfg::WP_SIGMA_BOHR
              << ", k0 " << cfg::WP_K0 << ", E " << cfg::WP_EKIN_EV << " eV\n"
              << "  launch  = (" << cfg::WP_CX_BOHR << ", " << cfg::WP_CY_BOHR
              << ", " << cfg::WP_CZ_BOHR << ") Bohr (boundary_rule)\n"
              << "  N_STEPS = " << cfg::N_STEPS
              << ", dt = " << cfg::DT_AU
              << ", t_total = " << (cfg::DT_AU * cfg::N_STEPS) << " a.u.\n"
              << "  WRITE_EVERY = " << cfg::WRITE_EVERY << "\n\n";

    std::filesystem::create_directories("results/raw/observables");

    auto cell = systems::cell::orthorhombic(
                    cfg::L_BOHR * 1.0_b,
                    cfg::L_BOHR * 1.0_b,
                    cfg::L_BOHR * 1.0_b)
                .finite();
    auto ions = systems::ions(cell);

    // Same ghost-state trick as the Heisenberg smoke test: INQ requires
    // num_electrons > 0 even when the WP is the only physically meaningful
    // orbital. The ghost lives at state 0; the WP is injected into state 1.
    const double ec_ha = 0.5 * std::pow(M_PI / cfg::SPACING_BOHR, 2.0);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}.cutoff(ec_ha * 1.0_Ha)
                            .extra_states(1)
                            .extra_electrons(2.0));

    ground_state::initial_guess(ions, electrons);

    auto report = inqkit::WavePacket{}
                      .center(cfg::WP_CX_BOHR, cfg::WP_CY_BOHR, cfg::WP_CZ_BOHR)
                      .sigma(cfg::WP_SIGMA_BOHR)
                      .k0(cfg::WP_KX, cfg::WP_KY, cfg::WP_KZ)
                      .inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx
              << "  norm_after=" << report.norm_after << "\n\n";

    // Write the boundary-rule + WP config to a flat text file so the
    // postprocess phases can recover (launch_z, sigma, k0) for IFW shading.
    {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "RUN SUMMARY\n===========\n\n"
          << "1. Run identity\n---------------\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "run_type        = free-space WP (Run-2)\n"
          << "date_started    = " << iso_now() << "\n\n"
          << "3. System configuration\n-----------------------\n"
          << "cell_bohr       = " << cfg::L_BOHR << "^3 (cubic, finite)\n"
          << "spacing_bohr    = " << cfg::SPACING_BOHR << "\n"
          << "xc_functional   = non-interacting\n\n"
          << "5. Wavepacket configuration and injection\n"
          << "-----------------------------------------\n"
          << "wp_enabled      = yes\n"
          << "wp_center_bohr  = " << cfg::WP_CX_BOHR << " " << cfg::WP_CY_BOHR
                                  << " " << cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr   = " << cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv  = " << cfg::WP_KX << " " << cfg::WP_KY << " "
                                  << cfg::WP_KZ << "\n"
          << "wp_energy_ev    = " << cfg::WP_EKIN_EV << "\n"
          << "norm_after      = " << report.norm_after << "\n\n"
          << "6. Real-time configuration\n--------------------------\n"
          << "rt_num_steps    = " << cfg::N_STEPS << "\n"
          << "dt_au           = " << cfg::DT_AU << "\n"
          << "total_time_au   = " << (cfg::DT_AU * cfg::N_STEPS) << "\n"
          << "write_every     = " << cfg::WRITE_EVERY << "\n";
    }

    // Observables — CSV only (no VTI for this validation run).
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    inqkit::io::ObservablesWriter obs_writer(
        "results/raw/observables/observables.csv", sel);
    obs_writer.write_header();

    inqkit::observables::WPMomentumStats wp_momentum_stats(
        "results/raw/observables/wp_momentum_stats.csv",
        wp_idx, {.write_every = cfg::WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        "results/raw/observables/wp_real_space_stats.csv",
        wp_idx, {.write_every = cfg::WRITE_EVERY});

    inqkit::RealTimeSession rt_obs(ions, electrons, cfg::WRITE_EVERY);
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
        options::theory{}.non_interacting(),
        options::real_time{}
            .num_steps(cfg::N_STEPS)
            .dt(cfg::DT_AU * 1.0_atomictime)
            .observables_dipole());

    double wall_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_start).count();

    // Append the completion marker.
    {
        std::ofstream s("results/run_summary.txt", std::ios::app);
        s << "\n9. End-of-run diagnostics\n-------------------------\n"
          << "wall_time_s     = " << wall_seconds << "\n"
          << "date_finished   = " << iso_now() << "\n"
          << "run_completed   = true\n";
    }

    std::cout << "Done. Wall time " << wall_seconds << " s.\n";
    return 0;
}
