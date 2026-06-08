// ============================================================================
// run_positive_ion_L50_v0p33 — proton (H⁺) projectile through L=50 jellium.
// Matched-velocity (v = v_F ≈ 0.337 a.u.) charge-conjugate companion of
// run_base_n162_L50_E1p5.
//
// Cfg: jellium::config::Positive_Ion_L50_v0p33.
//
// **Status: SCAFFOLD, NOT YET BUILT.** The shared run_template.hpp does NOT
// support the no-WP + ion-projectile mode out of the box; this file is
// therefore a hand-written run.cpp that follows the QBall-INQ Li pattern
// (`QuantumKickExtension/inq-codebase/Li/run_propagate_v0p0123_extensive/
//  run.cpp`) adapted to:
//   - one H ion (not 54 Li atoms),
//   - L=50 cubic jellium bath at N=162,
//   - finer dt (0.005 a.u.) for ion stability.
//
// See docs/plans/jellium_positive_ion_companion.md for the rationale,
// open API questions (N_e accounting, pseudopotential grid, dt convergence),
// and the verdict report scope.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/inqkit.hpp>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "../shared/configs/positive_ion_L50_v0p33.hpp"

using namespace inq;
using namespace inq::magnitude;
namespace cfg = jellium::config;

int main() {
    using Cfg = cfg::Positive_Ion_L50_v0p33;
    auto env = input::environment{};

    std::cout << "=== run_positive_ion_L50_v0p33 ===\n"
              << "  Cfg: Positive_Ion_L50_v0p33\n"
              << "  Bath: L=" << Cfg::L_BOHR << " Bohr, N=" << Cfg::N_ELECTRONS
              << " (closed shell)\n"
              << "  Projectile: H⁺ at (" << Cfg::ION_LAUNCH_X << ", "
              << Cfg::ION_LAUNCH_Y << ", " << Cfg::ION_LAUNCH_Z << ") Bohr, "
              << "v_z = " << Cfg::ION_VELOCITY_Z << " a.u.\n";

    // ----- Cell -----------------------------------------------------------
    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();

    // ----- Ions: one proton at the launch position -----------------------
    auto ions = systems::ions(cell);
    ions.insert(Cfg::ION_SPECIES,
                {Cfg::ION_LAUNCH_X * 1.0_b,
                 Cfg::ION_LAUNCH_Y * 1.0_b,
                 Cfg::ION_LAUNCH_Z * 1.0_b});

    // ----- Electrons: load N=162 closed-shell jellium GS ----------------
    // **NOTE.** This GS was computed *without* the proton, so the
    // wavefunctions are not self-consistent with the proton's Coulomb
    // attractor at t=0. Option G1 of the plan: accept the injection
    // transient. Option G2 (cleaner) requires a fresh GS save with the
    // proton at rest — TODO before the canonical run.
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV));
    electrons.load("/local/data/public/skcb2/tddft/ResearchProject/systems/"
                   "jellium/checkpoints/gs_L50_cubic_N162_dx1p0");
    std::cout << "  Loaded GS from gs_L50_cubic_N162_dx1p0\n"
              << "  num_states    = " << electrons.states().num_states() << "\n"
              << "  num_electrons = " << electrons.states().num_electrons() << "\n";

    // ----- Set proton velocity (impulsive +z kick) ----------------------
    ions.velocities()[0] = vector3<double>{
        Cfg::ION_VELOCITY_X, Cfg::ION_VELOCITY_Y, Cfg::ION_VELOCITY_Z};
    std::cout << "  Velocity set: " << Cfg::ION_VELOCITY_Z << " a.u. +z\n";

    // ----- Output skeleton (raw/ + analysis/ per spec) -------------------
    std::filesystem::create_directories("results/raw/observables");
    std::filesystem::create_directories("results/raw/vti/density_rt_total");
    std::filesystem::create_directories("results/analysis/observables");
    std::filesystem::create_directories("results/analysis/density");

    // ----- Per-step observables (same set as WP run + proton tracking) --
    inqkit::io::ObservableSelection sel;
    sel.step           = sel.time_au       = true;
    sel.energy_total   = sel.energy_kinetic = true;
    sel.energy_hartree = sel.energy_xc      = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x  = sel.dipole_y  = sel.dipole_z  = true;
    inqkit::io::ObservablesWriter obs_writer(
        "results/raw/observables/observables.csv", sel);
    obs_writer.write_header();

    // TODO: extend ObservablesWriter to add proton_{x,y,z}, proton_v{x,y,z},
    //       and force_z_on_proton columns. Workaround for now:
    //       write a sibling proton_track.csv in a separate session.
    std::ofstream proton_csv("results/raw/observables/proton_track.csv");
    proton_csv << std::setprecision(16);
    proton_csv << "step,time_au,proton_x,proton_y,proton_z,proton_vx,proton_vy,proton_vz,force_x,force_y,force_z\n";

    // ----- Density VTI writer -------------------------------------------
    inqkit::io::RealField3DLayout layout{
        .field_name   = "density",
        .include_meta = false,
        .emit_raw     = false,
        .emit_vti     = true,
        .vti_format   = inqkit::io::VTIWriteOptions::Format::binary,
    };
    inqkit::io::RealField3DWriter total_wr(
        "results/raw/vti/density_rt_total", layout, {.overwrite = true});

    inqkit::observables::DensityDelta delta_obs(
        "results/raw/vti/density_rt_delta",
        "results/raw/vti/density_rt_delta_coarse",
        { .emit_raw_vti = true, .emit_coarse_vti = true,
          .compute_l2 = true, .coarse_bin_bohr = 3.0 });

    // ----- Real-time sessions -------------------------------------------
    inqkit::RealTimeSession rt_dens(ions, electrons, Cfg::WRITE_EVERY);
    rt_dens.add([&](inqkit::StepContext const& ctx) {
        auto rho = inqkit::fields::density::total(*ctx.electrons);
        total_wr.write(rho, ctx.time_au, ctx.step);
        (void)delta_obs.snapshot(rho, ctx.time_au, ctx.step);
    });

    inqkit::RealTimeSession rt_obs(ions, electrons, 1);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        obs_writer.append(ctx);
        // Track proton state every step:
        auto const& pos = ctx.ions->positions()[0];
        auto const& vel = ctx.ions->velocities()[0];
        auto const& f   = ctx.forces ? (*ctx.forces)[0] : vector3<double>{0,0,0};
        proton_csv << ctx.step << "," << ctx.time_au << ","
                   << pos[0] << "," << pos[1] << "," << pos[2] << ","
                   << vel[0] << "," << vel[1] << "," << vel[2] << ","
                   << f[0]   << "," << f[1]   << "," << f[2]   << "\n";
    });

    // ----- t=0 frame ----------------------------------------------------
    {
        auto rho0 = inqkit::fields::density::total(electrons);
        total_wr.write(rho0, 0.0, 0);
        (void)delta_obs.snapshot(rho0, 0.0, 0);
    }

    // ----- Propagation: ETRS, impulsive ions ---------------------------
    real_time::propagate(
        ions, electrons,
        [&](auto const& data) {
            rt_dens.step(data);
            rt_obs.step(data);
        },
        options::theory{}.lda(),                          // bath: same xc as the GS
        options::real_time{}
            .num_steps(Cfg::N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .etrs()                                        // explicit; default
            .impulsive()                                   // ion velocity preserved
            .observables_current()
            .observables_dipole());

    // ----- run_summary.txt ----------------------------------------------
    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "RUN SUMMARY\n===========\n\n";
        s << "1. Run identity\n---------------\n";
        s << "run_name        = run_positive_ion_L50_v0p33\n";
        s << "run_type        = positive-ion (proton) projectile through jellium (TDDFT, ALDA)\n";
        s << "executable      = run.cpp built via inq-run\n\n";
        s << "3. System configuration\n-----------------------\n";
        s << "cell_bohr       = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n";
        s << "boundary        = periodic\n";
        s << "n_ions          = 1\n";
        s << "ion_species     = " << Cfg::ION_SPECIES << "\n";
        s << "n_electrons     = " << Cfg::N_ELECTRONS << "\n";
        s << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n";
        s << "xc_functional   = LDA (ALDA in TDDFT)\n\n";
        s << "5. Projectile configuration\n---------------------------\n";
        s << "wp_enabled      = no\n";
        s << "ion_launch_bohr = " << Cfg::ION_LAUNCH_X << " "
                                  << Cfg::ION_LAUNCH_Y << " "
                                  << Cfg::ION_LAUNCH_Z << "\n";
        s << "ion_velocity_au = " << Cfg::ION_VELOCITY_X << " "
                                  << Cfg::ION_VELOCITY_Y << " "
                                  << Cfg::ION_VELOCITY_Z << "\n";
        s << "atoms_dynamics  = impulsive\n\n";
        s << "6. Real-time configuration\n--------------------------\n";
        s << "rt_num_steps    = " << Cfg::N_STEPS << "\n";
        s << "dt_au           = " << Cfg::DT_AU << "\n";
        s << "total_time_au   = " << (Cfg::DT_AU * Cfg::N_STEPS) << "\n";
        s << "write_every     = " << Cfg::WRITE_EVERY << "\n\n";
        s << "9. End-of-run diagnostics\n-------------------------\n";
        s << "run_completed   = true\n";
    }

    return 0;
}
