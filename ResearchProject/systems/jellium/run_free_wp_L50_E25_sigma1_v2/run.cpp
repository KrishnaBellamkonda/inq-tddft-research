// ============================================================================
// run_free_wp_L50_E25_sigma1_v2/run.cpp — free-space counterpart of
// run_wp_n162_L50_E25_sigma1_v2.
//
// FREE-SPACE wave-packet propagation in a 50^3 Bohr PERIODIC box,
// non-interacting theory (no Hartree, no XC). Same WP injection geometry
// as the jellium E=25 σ=1 v2 run so density VTIs can be subtracted
// frame-by-frame to isolate the bath response:
//
//     δn_bath(r,t) = δn_jellium(r,t) - δn_free(r,t)
//
// Parameters (matching jellium E25 v2 exactly):
//   sigma   = 1 Bohr
//   k0_z    = 1.35553 Bohr^-1       (E = 25 eV)
//   launch  = (0, 0, -21)           (boundary_rule: -L/2 + 4σ)
//   dt      = 0.010 a.u.            (v2 halved timestep)
//   N_STEPS = 1100                  (self-spread capped, matching jellium run)
//   WRITE_EVERY = 4                 (~275 density frames)
//   WF_WRITE_EVERY = 10             (~110 wavefunction frames)
//
// Grid: dx=0.40 Bohr, periodic cubic — matches the jellium GS grid
// gs_L50_cubic_N162_dx0p40 so VTI files are on the same lattice.
//
// Ghost-state trick: INQ requires n_electrons > 0 even for a single-WP
// free-space run. We inject 2 extra_electrons (1 occupied ghost) and
// place the WP into state 1 (the extra_state). The ghost has uniform
// density ~1.6e-5 e/Bohr³ that cancels in δn = n(t) - n(0).
//
// Observables: full set matching the jellium WP run — density VTIs
// (total, system, delta, delta_coarse, wp, wavefunction_wp),
// observables.csv, state_energies, occupations, momentum_distribution,
// wp_momentum_stats, wp_real_space_stats. Overlap omitted (no bath
// states to project against).
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
#include <inqkit/observables/state_energy_writer.hpp>
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

namespace cfg {
    // Cell and grid — must match jellium GS gs_L50_cubic_N162_dx0p40
    constexpr double L_BOHR       = 50.0;
    constexpr double SPACING_BOHR = 0.40;

    // Timestep — v2 halved from 0.020
    constexpr double DT_AU = 0.010;

    // WP — identical to Electron_Proj_E25_L50_sigma1_v2_WP
    constexpr double WP_EKIN_EV    = 25.0;
    constexpr double WP_SIGMA_BOHR = 1.0;
    // k0 = sqrt(2 * E_eV / 27.21138625)
    constexpr double HA_TO_EV      = 27.21138625;
    constexpr double WP_K0         = 1.355531670051215; // sqrt(2 * 25 / 27.21138625)
    constexpr double WP_KX = 0.0, WP_KY = 0.0, WP_KZ = WP_K0;

    constexpr double WP_CX_BOHR = 0.0;
    constexpr double WP_CY_BOHR = 0.0;
    constexpr double WP_CZ_BOHR =
        jellium::config::boundary::launch_z(WP_SIGMA_BOHR, L_BOHR); // -21

    // Duration — matches jellium v2 run exactly
    constexpr int N_STEPS      = 1100;
    constexpr int WRITE_EVERY  = 4;   // ~275 density frames
    constexpr int WF_WRITE_EVERY = 10; // ~110 wavefunction frames
}

static std::string iso_now() {
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

int main() {
    auto t_wallclock_start = std::chrono::steady_clock::now();

    const std::string RUN_NAME = "run_free_wp_L50_E25_sigma1_v2";

    std::cout << "\n=== " << RUN_NAME << " (free-space WP @ 25 eV) ===\n"
              << "  cell    = " << cfg::L_BOHR << "^3 Bohr (cubic, periodic, NON-interacting)\n"
              << "  spacing = " << cfg::SPACING_BOHR << " Bohr\n"
              << "  WP      = sigma " << cfg::WP_SIGMA_BOHR
              << ", k0 " << cfg::WP_K0 << ", E " << cfg::WP_EKIN_EV << " eV\n"
              << "  launch  = (" << cfg::WP_CX_BOHR << ", " << cfg::WP_CY_BOHR
              << ", " << cfg::WP_CZ_BOHR << ") Bohr (boundary_rule)\n"
              << "  N_STEPS = " << cfg::N_STEPS
              << ", dt = " << cfg::DT_AU
              << ", t_total = " << (cfg::DT_AU * cfg::N_STEPS) << " a.u.\n"
              << "  WRITE_EVERY = " << cfg::WRITE_EVERY
              << ", WF_WRITE_EVERY = " << cfg::WF_WRITE_EVERY << "\n\n";

    // ----- Output skeleton ---------------------------------------------------
    std::filesystem::create_directories("results/raw/observables");
    std::filesystem::create_directories("results/raw/vti/density_total");
    std::filesystem::create_directories("results/raw/vti/density_system");
    std::filesystem::create_directories("results/raw/vti/density_delta");
    std::filesystem::create_directories("results/raw/vti/density_delta_coarse");
    std::filesystem::create_directories("results/raw/vti/density_wp");
    std::filesystem::create_directories("results/raw/vti/wavefunction_wp");

    // Stub run_summary so a crash leaves metadata.
    {
        std::ofstream s("results/run_summary.txt");
        s << "RUN SUMMARY (stub - written at start)\n"
          << "=====================================\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "date_started    = " << iso_now() << "\n"
          << "cell_bohr       = " << cfg::L_BOHR << "^3 (cubic, periodic)\n"
          << "spacing_bohr    = " << cfg::SPACING_BOHR << "\n"
          << "xc_functional   = non-interacting\n"
          << "wp_ekin_ev      = " << cfg::WP_EKIN_EV << "\n"
          << "wp_sigma_bohr   = " << cfg::WP_SIGMA_BOHR << "\n"
          << "dt_au           = " << cfg::DT_AU << "\n"
          << "n_steps         = " << cfg::N_STEPS << "\n"
          << "run_completed   = false\n";
    }

    // ----- Cell and electrons (periodic, matching jellium grid) ---------------
    auto cell = systems::cell::cubic(cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);

    // Ghost-state trick: 2 extra_electrons → 1 occupied ghost, +1 extra for WP.
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(2)
            .extra_states(1),
        input::kpoints::gamma());

    ground_state::initial_guess(ions, electrons);

    const int n_states = electrons.states().num_states();
    std::cout << "  num_states=" << n_states << " (1 ghost + 1 WP target)\n";

    // ----- VTI layout --------------------------------------------------------
    inqkit::io::RealField3DLayout vti_layout{
        .field_name   = "density",
        .include_meta = false,
        .emit_raw     = false,
        .emit_vti     = true,
        .vti_format   = inqkit::io::VTIWriteOptions::Format::binary,
    };

    // ----- t=0 "GS" density (ghost only, before WP injection) ----------------
    {
        std::filesystem::create_directories("results/raw/vti/density_gs_system");
        inqkit::io::RealField3DWriter gs_sys(
            "results/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
        gs_sys.write(inqkit::fields::density::total(electrons),
                     "density_gs_system");
    }

    // ----- WP injection (no orthogonalisation — no bath states) --------------
    auto wp = inqkit::WavePacket{}
                  .center(cfg::WP_CX_BOHR, cfg::WP_CY_BOHR, cfg::WP_CZ_BOHR)
                  .sigma(cfg::WP_SIGMA_BOHR)
                  .k0(cfg::WP_KX, cfg::WP_KY, cfg::WP_KZ);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx
              << "  norm_after=" << report.norm_after << "\n\n";

    // ----- WP artefacts ------------------------------------------------------
    {
        std::ofstream f("results/raw/observables/wp_config.txt");
        f << "wp_center_bohr = " << cfg::WP_CX_BOHR << " "
                                 << cfg::WP_CY_BOHR << " "
                                 << cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr  = " << cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv = " << cfg::WP_KX << " " << cfg::WP_KY << " "
                                 << cfg::WP_KZ << "\n"
          << "wp_energy_ev   = " << cfg::WP_EKIN_EV << "\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "orthogonalised = no\n";
    }
    {
        std::ofstream f("results/raw/observables/wp_injection_report.txt");
        f << std::setprecision(16);
        f << "norm_before      = " << report.norm_before << "\n"
          << "norm_after       = " << report.norm_after << "\n"
          << "max_overlap      = " << report.max_overlap << "\n"
          << "passed_tolerance = " << (report.passed_tolerance ? "yes" : "no") << "\n";
    }

    // ----- Density writers ---------------------------------------------------
    inqkit::io::RealField3DWriter total_wr(
        "results/raw/vti/density_total", vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(
        "results/raw/vti/density_system", vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter wp_density_wr(
        "results/raw/vti/density_wp", vti_layout, {.overwrite=true});
    inqkit::io::ComplexField3DWriter wp_wf_wr(
        "results/raw/vti/wavefunction_wp",
        {.field_name = "wavefunction", .include_meta = false,
         .emit_raw = false, .emit_vti = true,
         .vti_format = inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite=true});

    // t=0 density frames
    {
        auto sys0 = inqkit::fields::density::total(electrons);
        system_wr.write(sys0, 0.0, 0);
        total_wr.write( sys0, 0.0, 0);
    }

    // ----- Observables -------------------------------------------------------
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
        wp_idx, cfg::L_BOHR,
        {.n_bins = 64, .k_max_bohr_inv = 0.0,
         .write_every = 10 * cfg::WRITE_EVERY});

    inqkit::observables::WPMomentumStats wp_momentum_stats(
        "results/raw/observables/wp_momentum_stats.csv",
        wp_idx, {.write_every = cfg::WRITE_EVERY});
    inqkit::observables::WPRealSpaceStats wp_real_space_stats(
        "results/raw/observables/wp_real_space_stats.csv",
        wp_idx, {.write_every = cfg::WRITE_EVERY});

    inqkit::observables::DensityDelta density_delta(
        "results/raw/vti/density_delta",
        "results/raw/vti/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0});

    // ----- Real-time session callbacks ---------------------------------------
    inqkit::RealTimeSession rt_obs(ions, electrons, cfg::WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const &ctx) {
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        total_wr.write( sys_f, ctx.time_au, ctx.step);

        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext ctx_out = ctx;
        ctx_out.density_l2 = l2;
        obs_writer.append(ctx_out);

        if (ctx.step % cfg::WF_WRITE_EVERY == 0) {
            auto wp_dens = inqkit::fields::density::orbital(
                *ctx.electrons, wp_idx);
            wp_density_wr.write(wp_dens, ctx.time_au, ctx.step);

            auto wp_wf = inqkit::fields::orbital::wavefunction(
                *ctx.electrons, wp_idx);
            char wf_name[64];
            std::snprintf(wf_name, sizeof(wf_name),
                          "wavefunction_t%06d", ctx.step);
            wp_wf_wr.write(wp_wf, std::string(wf_name));
        }
    });

    // ----- Propagate (non-interacting: p²/2m only) ---------------------------
    real_time::propagate(
        ions, electrons,
        [&](auto const &data) {
            rt_obs.step(data);
            if (data.iter() % (5 * cfg::WRITE_EVERY) == 0) {
                state_energy_wr.snapshot(data);
                occupations_wr.snapshot(data);
            }
            momentum_dist.maybe_accumulate(data);
            wp_momentum_stats.maybe_accumulate(data);
            wp_real_space_stats.maybe_accumulate(data);
        },
        options::theory{}.non_interacting(),
        options::real_time{}
            .num_steps(cfg::N_STEPS)
            .dt(cfg::DT_AU * 1.0_atomictime)
            .observables_current()
            .observables_dipole());

    // ----- Final run_summary -------------------------------------------------
    auto t_wallclock_end = std::chrono::steady_clock::now();
    double wall_seconds =
        std::chrono::duration<double>(t_wallclock_end - t_wallclock_start).count();
    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "RUN SUMMARY\n===========\n\n"
          << "1. Run identity\n---------------\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "run_type        = free-space WP (non-interacting counterpart of E25 v2)\n"
          << "date_finished   = " << iso_now() << "\n"
          << "wall_time_s     = " << wall_seconds << "\n\n"
          << "3. System configuration\n-----------------------\n"
          << "cell_bohr       = " << cfg::L_BOHR << "^3 (cubic, periodic)\n"
          << "n_electrons     = 2 (ghost + WP)\n"
          << "n_states        = " << n_states << "\n"
          << "spacing_bohr    = " << cfg::SPACING_BOHR << "\n"
          << "xc_functional   = non-interacting\n\n"
          << "5. Wavepacket configuration and injection\n-----------------------------------------\n"
          << "wp_enabled      = yes\n"
          << "wp_center_bohr  = " << cfg::WP_CX_BOHR << " " << cfg::WP_CY_BOHR
                                  << " " << cfg::WP_CZ_BOHR << "\n"
          << "wp_sigma_bohr   = " << cfg::WP_SIGMA_BOHR << "\n"
          << "wp_k0_bohr_inv  = " << cfg::WP_KX << " " << cfg::WP_KY << " "
                                  << cfg::WP_KZ << "\n"
          << "wp_energy_ev    = " << cfg::WP_EKIN_EV << "\n"
          << "norm_after      = " << report.norm_after << "\n"
          << "max_overlap     = " << report.max_overlap << "\n\n"
          << "6. Real-time configuration\n--------------------------\n"
          << "rt_num_steps    = " << cfg::N_STEPS << "\n"
          << "dt_au           = " << cfg::DT_AU << "\n"
          << "total_time_au   = " << (cfg::DT_AU * cfg::N_STEPS) << "\n"
          << "write_every     = " << cfg::WRITE_EVERY << "\n"
          << "wf_write_every  = " << cfg::WF_WRITE_EVERY << "\n\n"
          << "9. End-of-run diagnostics\n-------------------------\n"
          << "run_completed   = true\n";
    }

    std::cout << "Done. Wall time " << wall_seconds << " s.\n";
    return 0;
}
