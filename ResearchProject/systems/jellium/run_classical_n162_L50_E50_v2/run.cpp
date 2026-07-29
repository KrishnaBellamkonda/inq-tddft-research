// ============================================================================
// run_classical_n162_L50_E50_v2/run.cpp — classical-electron projectile (custom
// UPF + mass override = m_e) through the cubic 50^3 Bohr periodic jellium
// bath at N=162, dx=0.40, E_kin=50 eV. v2: dt=0.01 a.u.
//
// Cfg: jellium::config::Electron_Proj_E50_L50_sigma1_v2_Classical.
// GS: checkpoints/gs_L50_cubic_N162_dx0p40.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/center_of_density.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/jellium/shells.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>

#include "../shared/configs/electron_proj_E50_L50_cubic_sigma1_v2.hpp"
#include "../shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::Electron_Proj_E50_L50_sigma1_v2_Classical;

static std::string iso_now() {
    auto t  = std::time(nullptr);
    auto tm = *std::localtime(&t);
    char buf[64];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%S", &tm);
    return std::string(buf);
}

int main() {
    auto t_wallclock_start = std::chrono::steady_clock::now();

    const std::string RUN_NAME = "run_classical_n162_L50_E50_v2";
    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p40";

    std::cout << "\n=== " << RUN_NAME << " ===\n"
              << "  cell = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
              << "  bath: N_e=" << Cfg::N_ELECTRONS << ", dx="
              << Cfg::SPACING_BOHR << " Bohr\n"
              << "  projectile: classical electron at ("
              << Cfg::PROJ_LAUNCH_X << ", " << Cfg::PROJ_LAUNCH_Y << ", "
              << Cfg::PROJ_LAUNCH_Z << ") Bohr,  v=("
              << Cfg::PROJ_VEL_X << ", " << Cfg::PROJ_VEL_Y << ", "
              << Cfg::PROJ_VEL_Z << ") bohr/atu\n"
              << "  KE = " << (0.5 * 1.0 * Cfg::PROJ_VEL_Z * Cfg::PROJ_VEL_Z)
              << " Ha = " << Cfg::WP_EKIN_EV << " eV (matched to WP)\n"
              << "  pseudopotential = " << Cfg::PROJ_PSEUDO_PATH << "\n"
              << "  mass_amu (override) = " << Cfg::PROJ_MASS_AMU
              << " => mass_au = " << (Cfg::PROJ_MASS_AMU * 1822.8885)
              << " (expect 1.0 = m_e)\n"
              << "  dt=" << Cfg::DT_AU << " a.u., N_steps="
              << Cfg::N_STEPS << ", t_total="
              << (Cfg::DT_AU * Cfg::N_STEPS) << " a.u.\n"
              << "  checkpoint = " << GS_DIR << "\n\n";

    if (!std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: GS checkpoint missing — run save_gs first.\n";
        return 2;
    }

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();

    // ----- Insert one classical-electron ion ----------------------------
    auto ions = systems::ions(cell);
    auto sp_e = ionic::species(Cfg::PROJ_SPECIES_SYMBOL)
                    .pseudo_file(Cfg::PROJ_PSEUDO_PATH)
                    .mass(Cfg::PROJ_MASS_AMU);
    ions.insert(sp_e, {Cfg::PROJ_LAUNCH_X * 1.0_b,
                       Cfg::PROJ_LAUNCH_Y * 1.0_b,
                       Cfg::PROJ_LAUNCH_Z * 1.0_b});
    std::cout << "  Inserted projectile ion. mass_au="
              << ions.species(0).mass() << " (expect 1.0)\n";

    // ----- Load the bath GS (computed without the projectile) -----------
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

    // ----- Set the projectile's initial velocity ------------------------
    ions.velocities()[0] = vector3<double>{Cfg::PROJ_VEL_X,
                                            Cfg::PROJ_VEL_Y,
                                            Cfg::PROJ_VEL_Z};
    std::cout << "  Velocity set on ion 0: ("
              << Cfg::PROJ_VEL_X << ", " << Cfg::PROJ_VEL_Y
              << ", " << Cfg::PROJ_VEL_Z << ") bohr/atu\n";
    std::cout << "  KE_check (1/2 m v^2) = " << ions.kinetic_energy()
              << " Ha (expect " << (Cfg::WP_EKIN_EV / 27.211386) << ")\n";

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    const int n_occupied  = n_electrons / 2;
    std::cout << "  num_states=" << n_states
              << "  num_electrons=" << n_electrons
              << "  n_occupied=" << n_occupied << "\n";

    // ----- Output skeleton ----------------------------------------------
    std::filesystem::create_directories("results/raw/observables");
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
          << "projectile_kind = classical electron (custom UPF, m_e mass)\n"
          << "projectile_KE_eV = " << Cfg::WP_EKIN_EV << "\n"
          << "projectile_v    = " << Cfg::PROJ_VEL_X << " "
                                  << Cfg::PROJ_VEL_Y << " "
                                  << Cfg::PROJ_VEL_Z << " bohr/atu\n"
          << "dt_au           = " << Cfg::DT_AU << "\n"
          << "n_steps         = " << Cfg::N_STEPS << "\n"
          << "run_completed   = false\n";
    }

    // ----- VTI layout ----------------------------------------------------
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

    // ----- Real-time output writers --------------------------------------
    const int n_ref = n_states - 1;
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(
        electrons, n_ref, "results/raw/observables/overlap_full");

    // Proxy snapshot infrastructure.
    auto shell_table = inqkit::jellium::shells::enumerate_for_n_states(n_states);
    auto proxies     = inqkit::jellium::shells::pick_proxies(shell_table, 2);
    inqkit::observables::OrbitalOverlapMatrix overlap_proxy_obs(
        electrons, n_ref, "results/raw/observables/overlap_proxies");
    inqkit::jellium::shells::write_shells_csv(
        shell_table, proxies, "results/raw/observables/overlap_proxies");
    std::cout << "  Shell table: " << shell_table.size() << " shells, "
              << proxies.size() << " proxies (2 per shell)\n";

    inqkit::io::RealField3DWriter total_wr(
        "results/raw/vti/density_total", vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(
        "results/raw/vti/density_system", vti_layout, {.overwrite=true});

    // t=0 density frame
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
        /*wp_state_index=*/-1,                      // -1 = no WP slot
        Cfg::L_BOHR,
        {.n_bins = 64, .k_max_bohr_inv = 0.0,
         .write_every = 10 * Cfg::WRITE_EVERY});
    inqkit::observables::DensityDelta density_delta(
        "results/raw/vti/density_delta",
        "results/raw/vti/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0});

    // electron_track.csv: projectile pos / vel every step.
    std::ofstream electron_csv("results/raw/observables/electron_track.csv");
    electron_csv << std::setprecision(16);
    electron_csv << "step,time_au,x,y,z,vx,vy,vz,fx,fy,fz\n";
    {
        auto const &p = ions.positions()[0];
        auto const &v = ions.velocities()[0];
        electron_csv << 0 << "," << 0.0 << ","
                     << p[0] << "," << p[1] << "," << p[2] << ","
                     << v[0] << "," << v[1] << "," << v[2] << ","
                     << 0.0 << "," << 0.0 << "," << 0.0 << "\n";
    }

    // Full-matrix overlap snapshot at t=0 (bath self-overlap ~ I).
    overlap_full_obs.snapshot(electrons, 0.0, 0);
    std::cout << "  Full overlap snapshot at step 0\n";

    overlap_proxy_obs.snapshot_proxies(electrons, proxies, 0.0, 0);

    const int FULL_OVERLAP_STEP_MID = Cfg::N_STEPS / 2;
    const int FULL_OVERLAP_STEP_END = Cfg::N_STEPS;
    const int PROXY_SNAPSHOT_STRIDE = 5 * Cfg::WRITE_EVERY;

    // ----- Real-time session callbacks -----------------------------------
    inqkit::RealTimeSession rt_dens(ions, electrons, Cfg::WRITE_EVERY);
    rt_dens.add([&](inqkit::StepContext const &ctx) {
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        total_wr.write( sys_f, ctx.time_au, ctx.step);

        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext ctx_out = ctx;
        ctx_out.density_l2 = l2;
        obs_writer.append(ctx_out);

        if (ctx.step == FULL_OVERLAP_STEP_MID) {
            overlap_full_obs.snapshot(*ctx.electrons, ctx.time_au, ctx.step);
            std::cout << "  Full overlap snapshot at step "
                      << ctx.step << " (mid)\n";
        }

        if (ctx.step > 0 && ctx.step % PROXY_SNAPSHOT_STRIDE == 0) {
            overlap_proxy_obs.snapshot_proxies(
                *ctx.electrons, proxies, ctx.time_au, ctx.step);
        }
    });

    // electron_track.csv writer — every step.
    inqkit::RealTimeSession rt_track(ions, electrons, 1);
    rt_track.add([&](inqkit::StepContext const &ctx) {
        auto const &p = ctx.ions->positions()[0];
        auto const &v = ctx.ions->velocities()[0];
        electron_csv << ctx.step << "," << ctx.time_au << ","
                     << p[0] << "," << p[1] << "," << p[2] << ","
                     << v[0] << "," << v[1] << "," << v[2] << ","
                     << 0.0 << "," << 0.0 << "," << 0.0 << "\n";
    });

    real_time::propagate(
        ions, electrons,
        [&](auto const &data) {
            rt_dens.step(data);
            rt_track.step(data);
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
            .ehrenfest()                       // electronic forces decelerate ion
            .observables_current()
            .observables_dipole());

    overlap_proxy_obs.snapshot_proxies(
        electrons, proxies,
        Cfg::DT_AU * Cfg::N_STEPS, FULL_OVERLAP_STEP_END);

    overlap_full_obs.snapshot(electrons,
        Cfg::DT_AU * Cfg::N_STEPS, FULL_OVERLAP_STEP_END);
    std::cout << "  Full overlap snapshot at step "
              << FULL_OVERLAP_STEP_END << " (end)\n";

    auto t_wallclock_end = std::chrono::steady_clock::now();
    double wall_seconds =
        std::chrono::duration<double>(t_wallclock_end - t_wallclock_start).count();
    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "RUN SUMMARY\n===========\n\n"
          << "1. Run identity\n---------------\n"
          << "run_name        = " << RUN_NAME << "\n"
          << "run_type        = classical-electron projectile, jellium TDDFT (ALDA)\n"
          << "date_finished   = " << iso_now() << "\n"
          << "wall_time_s     = " << wall_seconds << "\n\n"
          << "3. System configuration\n-----------------------\n"
          << "cell_bohr       = " << Cfg::L_BOHR
                                  << "^3 (cubic, periodic)\n"
          << "n_electrons     = " << n_electrons << "\n"
          << "n_occupied      = " << n_occupied << "\n"
          << "extra_states    = " << Cfg::EXTRA_STATES << "\n"
          << "spacing_bohr    = " << Cfg::SPACING_BOHR << "\n"
          << "xc_functional   = LDA (ALDA in TDDFT)\n\n"
          << "5. Projectile configuration\n---------------------------\n"
          << "projectile_kind = classical electron (custom UPF + mass override)\n"
          << "pseudopotential = " << Cfg::PROJ_PSEUDO_PATH << "\n"
          << "mass_amu        = " << Cfg::PROJ_MASS_AMU
                                  << " (= 1.0 / 1822.8885)\n"
          << "mass_au         = " << ions.species(0).mass()
                                  << " (expect 1.0 = m_e)\n"
          << "launch_bohr     = " << Cfg::PROJ_LAUNCH_X << " "
                                  << Cfg::PROJ_LAUNCH_Y << " "
                                  << Cfg::PROJ_LAUNCH_Z << "\n"
          << "velocity_atu    = " << Cfg::PROJ_VEL_X << " "
                                  << Cfg::PROJ_VEL_Y << " "
                                  << Cfg::PROJ_VEL_Z << "\n"
          << "KE_eV           = " << Cfg::WP_EKIN_EV << "\n"
          << "ion_dynamics    = ehrenfest (electronic forces decelerate ion)\n\n"
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
