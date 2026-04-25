// ============================================================================
// run_07_paper_replica: exact replica of Tsubonoya 2014 with full observable
// suite from docs/observables_reference.md.
//
// System:
//   Coronene C24H12 (D6h), 18.4 x 18.4 x 31.7 A^3 finite cell.
//   Atoms centred at the origin, z=0 plane (post-fix coordinate convention).
// Theory:
//   ALDA (paper), INQ default norm-conserving pseudopotentials.
//   GS SCF tolerance = 1e-6 Ha.
// Wave-packet:
//   sigma = 0.53 A, |k| = sqrt(2 * 200 eV) Bohr^-1 directed -z, centre at
//   (0, 0, +6.35 A) so it hits the flake at z=0 from above.
// Propagation:
//   dt = 0.020 a.u. (~4.84e-4 fs), 10000 steps -> 200 a.u. = 4.84 fs.
// Observables (every step unless noted):
//   - GS total density, all GS orbital densities (once at GS)
//   - RT target density and RT WP density (every WRITE_EVERY=30 steps)
//   - Energies (total, kinetic, Hartree, xc), currents Jx/Jy/Jz,
//     dipole mu_x/mu_y/mu_z (every step, observables.csv)
//   - Overlap matrix O_ij(t) (every step, results/overlap/)
//   - 20 LEED screens with three accumulators each:
//       full-run avg            -> results/screens/screen_NN.dat
//       paper window [t1, t2]   -> results/screens_leed_window/screen_NN.dat
//       per-step snapshots      -> results/screens_snapshots/step_XXXXXX/...
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/config/tsubonoya_2014_coronene.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace cfg = inqkit::config::tsubonoya_2014;

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

int main() {
    std::cout << "\n=== run_07: Tsubonoya 2014 coronene+WP replica ===\n";
    std::cout << "  cell = " << cfg::LX_BOHR << " x " << cfg::LY_BOHR
              << " x " << cfg::LZ_BOHR << " Bohr\n";
    std::cout << "  WP sigma = " << cfg::WP_SIGMA_BOHR
              << " Bohr, |k| = " << cfg::WP_K0 << " Bohr^-1, dt = "
              << cfg::DT_AU << " a.u., N_steps = " << cfg::N_STEPS << "\n";

    // ----- Cell + atoms ----------------------------------------------------
    auto cell = systems::cell::orthorhombic(
        cfg::LX_BOHR * 1.0_b, cfg::LY_BOHR * 1.0_b, cfg::LZ_BOHR * 1.0_b
    ).finite();
    auto ions = systems::ions::parse("coronene_centred.xyz", cell);
    std::cout << "  Atoms: " << ions.size() << "\n";

    // Defensive cell-bounds check (post-fix lesson from run_06).
    {
        const double half_lx = 0.5 * cfg::LX_BOHR;
        const double half_ly = 0.5 * cfg::LY_BOHR;
        const double half_lz = 0.5 * cfg::LZ_BOHR;
        for (int iatom = 0; iatom < static_cast<int>(ions.size()); ++iatom) {
            auto const & p = ions.positions()[iatom];
            if (std::fabs(p[0]) > half_lx ||
                std::fabs(p[1]) > half_ly ||
                std::fabs(p[2]) > half_lz) {
                std::cerr << "FATAL: atom " << iatom << " outside [-L/2, +L/2]: "
                          << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
                return 2;
            }
        }
    }

    // ----- Electrons + GS SCF ---------------------------------------------
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .cutoff(cfg::CUTOFF_HA * 1.0_Ha)
            .extra_states(cfg::EXTRA_STATES)
    );

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),                       // ALDA
        options::ground_state{}
            .energy_tolerance(cfg::SCF_TOL_HA * 1.0_Ha)
            .max_steps(cfg::SCF_MAX_STEPS)
            .broyden_mixing()
            .mixing_ndim(cfg::SCF_MIX_NDIM)
            .mixing(cfg::SCF_MIX_ALPHA)
    );
    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    const int n_occupied  = n_electrons / 2;          // spin-unpolarised, doubly occupied
    std::cout << "  num_states = " << n_states
              << "  num_electrons = " << n_electrons
              << "  n_occupied = " << n_occupied << "\n";

    std::filesystem::create_directories("results");

    // ----- GS total density and per-orbital densities ---------------------
    {
        inqkit::io::RealField3DWriter gs_wr("results/density_gs",
            {.field_name = "density", .include_meta = true},
            {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
    }
    std::filesystem::create_directories("results/density_gs_orbitals");
    for (int i = 0; i < n_states; ++i) {
        const auto out_dir = std::string("results/density_gs_orbitals/orbital_") + zero_pad(i, 4);
        std::filesystem::create_directories(out_dir);
        inqkit::io::RealField3DWriter orb_wr(out_dir,
            {.field_name = "density", .include_meta = true},
            {.overwrite = true});
        orb_wr.write(inqkit::fields::density::orbital(electrons, i), 0.0, 0);
    }
    std::cout << "  Wrote GS density and " << n_states << " orbital densities\n";

    // ----- WP injection ---------------------------------------------------
    auto wp = inqkit::WavePacket{}
        .center(cfg::WP_CX_BOHR, cfg::WP_CY_BOHR, cfg::WP_CZ_BOHR)
        .sigma(cfg::WP_SIGMA_BOHR)
        .k0(cfg::WP_KX, cfg::WP_KY, cfg::WP_KZ)
        .orthogonalise_against_occupied(electrons);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    std::cout << "  WP injected: state_index = " << wp_idx
              << "  norm_after = "  << report.norm_after
              << "  max_overlap = " << report.max_overlap << "\n";

    // ----- Overlap matrix observer ---------------------------------------
    inqkit::observables::OrbitalOverlapMatrix overlap_obs(
        electrons, wp_idx, "results/overlap");

    // ----- RT density writers (target + WP only; total is target+WP) -----
    inqkit::io::RealField3DWriter target_wr("results/density_rt_target",
        {.field_name = "density", .include_meta = true},
        {.overwrite = true});
    inqkit::io::RealField3DWriter wp_wr("results/density_rt_wp",
        {.field_name = "density", .include_meta = true},
        {.overwrite = true});

    // t=0 frames
    target_wr.write(inqkit::fields::density::total(electrons),  0.0, 0);
    wp_wr.write(    inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);

    // ----- Observables CSV (every step) ----------------------------------
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x  = sel.dipole_y  = sel.dipole_z  = true;
    inqkit::io::ObservablesWriter obs_writer("results/observables.csv", sel);
    obs_writer.write_header();

    // ----- 20 LEED screens, each with 2 accumulators -------------------
    // Layout: boundaries inset 1 Bohr from cell faces; 18 interior screens
    // with a deterministic jitter (no two on the same plane).
    std::array<double, cfg::N_SCREENS> screen_z;
    {
        const double half_lz = 0.5 * cfg::LZ_BOHR;
        screen_z[0]                  = -half_lz + 1.0;          // near -Lz/2
        screen_z[cfg::N_SCREENS - 1] = +half_lz - 1.0;          // near +Lz/2
        const int n_interior = cfg::N_SCREENS - 2;
        const double z_lo = -half_lz + 2.5;
        const double z_hi = +half_lz - 2.5;
        const double dz   = (z_hi - z_lo) / static_cast<double>(n_interior - 1);
        for (int k = 0; k < n_interior; ++k) {
            // small deterministic jitter so screens never coincide with grid planes
            const double jitter = ((k % 2) == 0 ? +0.07 : -0.13);
            screen_z[k + 1] = z_lo + k * dz + jitter;
        }
        std::cout << "  Screen z-positions (Bohr):";
        for (auto z : screen_z) std::cout << " " << z;
        std::cout << "\n";
    }

    std::array<inqkit::screens::LeedPatternAccumulator, cfg::N_SCREENS> acc_full;
    std::array<inqkit::screens::LeedPatternAccumulator, cfg::N_SCREENS> acc_window;
    for (int k = 0; k < cfg::N_SCREENS; ++k) {
        const std::string name = "screen_" + zero_pad(k, 2);
        acc_full[k]   = inqkit::screens::LeedPatternAccumulator(
            inqkit::screens::PlaneScreen{screen_z[k], name});
        acc_window[k] = inqkit::screens::LeedPatternAccumulator(
            inqkit::screens::PlaneScreen{screen_z[k], name});
    }

    // ----- Real-time propagation ----------------------------------------
    // Session A: density writes every WRITE_EVERY steps.
    inqkit::RealTimeSession rt_dens(ions, electrons, cfg::WRITE_EVERY);
    rt_dens.add([&](inqkit::StepContext const& ctx) {
        target_wr.write(inqkit::fields::density::total(*ctx.electrons),
                        ctx.time_au, ctx.step);
        wp_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx),
                    ctx.time_au, ctx.step);
    });

    // Session B: observables, screens, overlap every step.
    inqkit::RealTimeSession rt_obs(ions, electrons, 1);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        obs_writer.append(ctx);
        overlap_obs.snapshot(*ctx.electrons, ctx.time_au, ctx.step);

        // Always-on full-run accumulator
        for (auto & a : acc_full)   a.accumulate(*ctx.electrons, cfg::DT_AU);
        // Paper window accumulator: gate on [T1_AU, T2_AU]
        if (ctx.time_au >= cfg::T1_AU && ctx.time_au <= cfg::T2_AU) {
            for (auto & a : acc_window) a.accumulate(*ctx.electrons, cfg::DT_AU);
        }

        // Per-step snapshots every SCREEN_SNAP_EVERY steps
        if (ctx.step % cfg::SCREEN_SNAP_EVERY == 0) {
            const std::string step_dir =
                "results/screens_snapshots/step_" + zero_pad(ctx.step, 6);
            std::filesystem::create_directories(step_dir);
            for (int k = 0; k < cfg::N_SCREENS; ++k) {
                auto slice = acc_full[k].screen().extract(*ctx.electrons);
                acc_full[k].screen().save(slice, ctx.time_au,
                    step_dir + "/screen_" + zero_pad(k, 2) + ".dat");
            }
        }
    });

    real_time::propagate(ions, electrons,
        [&](auto const & data) { rt_dens.step(data); rt_obs.step(data); },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(cfg::N_STEPS)
            .dt(cfg::DT_AU * 1.0_atomictime)
            .observables_current()
            .observables_dipole());

    // ----- Final: save time-averaged screens -----------------------------
    std::filesystem::create_directories("results/screens");
    std::filesystem::create_directories("results/screens_leed_window");
    for (int k = 0; k < cfg::N_SCREENS; ++k) {
        acc_full[k].save(  "results/screens/screen_"             + zero_pad(k, 2) + ".dat");
        acc_window[k].save("results/screens_leed_window/screen_" + zero_pad(k, 2) + ".dat");
    }

    // ----- Run summary file ----------------------------------------------
    if (electrons.root()) {
        std::ofstream summary("results/run_summary.txt");
        summary << std::setprecision(16);
        summary << "run = run_07_paper_replica (Tsubonoya 2014)\n";
        summary << "system = coronene_C24H12\n";
        summary << "geometry_file = coronene_centred.xyz\n";
        summary << "cell_bohr = " << cfg::LX_BOHR << ' ' << cfg::LY_BOHR << ' ' << cfg::LZ_BOHR << "\n";
        summary << "boundary = finite\n";
        summary << "xc = ALDA\n";
        summary << "cutoff_ha = " << cfg::CUTOFF_HA << "\n";
        summary << "extra_states = " << cfg::EXTRA_STATES << "\n";
        summary << "scf_tol_ha = " << cfg::SCF_TOL_HA << "\n";
        summary << "ground_state_energy_ha = " << gs.energy.total() << "\n";
        summary << "num_states = " << n_states << "\n";
        summary << "num_electrons = " << n_electrons << "\n";
        summary << "n_occupied = " << n_occupied << "\n";
        summary << "wp_state_index = " << wp_idx << "\n";
        summary << "wp_norm_after = " << report.norm_after << "\n";
        summary << "wp_max_overlap = " << report.max_overlap << "\n";
        summary << "wp_sigma_bohr = " << cfg::WP_SIGMA_BOHR << "\n";
        summary << "wp_k0 = " << cfg::WP_K0 << "\n";
        summary << "wp_offset_bohr = " << cfg::WP_OFFSET_BOHR << "\n";
        summary << "wp_ekin_ev = " << cfg::WP_EKIN_EV << "\n";
        summary << "dt_au = " << cfg::DT_AU << "\n";
        summary << "n_steps = " << cfg::N_STEPS << "\n";
        summary << "write_every = " << cfg::WRITE_EVERY << "\n";
        summary << "screen_snap_every = " << cfg::SCREEN_SNAP_EVERY << "\n";
        summary << "t1_au = " << cfg::T1_AU << "\n";
        summary << "t2_au = " << cfg::T2_AU << "\n";
        for (int k = 0; k < cfg::N_SCREENS; ++k)
            summary << "screen_z[" << k << "] = " << screen_z[k] << "\n";
    }

    std::cout << "Done. Output in results/\n";
    return 0;
}
