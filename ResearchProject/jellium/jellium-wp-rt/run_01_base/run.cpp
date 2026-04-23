// ============================================================================
// jellium-wp-rt / run_01_base: N=38 closed-shell jellium, 200 eV, σ=0.53 Å, +z
//
// Physical: L=40.0 bohr=21.18 Å; N_gs=38, r_s=7.38 a₀, n=5.94e-4 e/bohr³
// WP: σ=0.53 Å=1.00 bohr, k₀=3.834 bohr⁻¹, v=3.834 bohr/a.u.=83.5 Å/fs
// Loop-back: T_sim=8.34 a.u. / T_loop=10.43 a.u. = 80% → OK
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>
#include <array>
#include <cmath>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
using namespace inq;
using namespace inq::magnitude;

static constexpr double ANG_TO_BOHR = 1.8897259886;
static constexpr double HA_TO_EV    = 27.21138625;
static constexpr double L_BOHR      = 40.0;

static constexpr double WP_SIGMA_ANG  = 0.53;
static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;
static constexpr double WP_EKIN_EV    = 200.0;
static constexpr double WP_EKIN_HA    = WP_EKIN_EV / HA_TO_EV;
static const     double WP_K0         = std::sqrt(2.0 * WP_EKIN_HA);

static const double WP_CX = L_BOHR / 2.0;
static const double WP_CY = L_BOHR / 2.0;
static const double WP_CZ = 5.0 * WP_SIGMA_BOHR;

static constexpr int    N_STEPS          = 417;
static constexpr double DT_AU            = 0.02;
static constexpr int    WRITE_EVERY      = 2;
static constexpr int    SCREEN_SNAP_EVERY = 3;
static constexpr int    N_SCREENS        = 20;

// 20 screen z-positions: boundaries at 0.5 and 39.5 bohr, 18 interior with jitter
static constexpr double SCREEN_Z[N_SCREENS] = {
    0.5,
    2.53, 4.66, 6.78, 8.87, 10.95, 12.97, 15.03, 17.06, 19.09,
    21.07, 23.11, 25.08, 27.12, 29.04, 31.03, 33.01, 34.97, 36.95,
    39.5
};

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

static void add_field_inplace(inqkit::fields::RealField3D& a,
                               inqkit::fields::RealField3D const& b) {
    for (std::size_t i = 0; i < a.values.size(); i++) a.values[i] += b.values[i];
}

int main() {
    std::cout << "\n=== jellium run_01_base: N=38, 200 eV, σ=0.53 Å, +z ===\n";

    auto cell      = systems::cell::cubic(L_BOHR * 1.0_b).periodic();
    auto ions      = systems::ions(cell);
    auto electrons = systems::electrons(ions,
        options::electrons{}
            .spacing(0.50 * 1.0_b)
            .extra_electrons(38)
            .extra_states(3)
            .temperature(0.00862 * 1.0_eV),
        input::kpoints::gamma());

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(1e-4_Ha).max_steps(300)
            .broyden_mixing().mixing_ndim(8).mixing(0.1));
    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    // Save GS total density (jellium, no WP)
    {
        inqkit::io::RealField3DWriter gs_wr("results/density_gs",
            {.field_name = "density", .include_meta = true}, {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
    }

    // Inject WP
    auto wp = inqkit::WavePacket{}
        .center(WP_CX, WP_CY, WP_CZ)
        .sigma(WP_SIGMA_BOHR)
        .k0(0.0, 0.0, WP_K0)
        .orthogonalise_against_occupied(electrons);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    int wp_idx = report.state_index;
    std::cout << "  state_index = " << wp_idx
              << "  norm_after = "  << report.norm_after
              << "  max_overlap = " << report.max_overlap << "\n";

    // Save GS orbital densities (indices 0..wp_idx-1)
    std::filesystem::create_directories("results/density_gs_orbitals");
    for (int i = 0; i < wp_idx; ++i) {
        inqkit::io::RealField3DWriter orb_wr(
            "results/density_gs_orbitals/orbital_" + zero_pad(i, 4),
            {.field_name = "density", .include_meta = true}, {.overwrite = true});
        orb_wr.write(inqkit::fields::density::orbital(electrons, i), 0.0, i);
    }

    // Overlap matrix observer (constructed after injection — GS refs are orbitals 0..wp_idx-1)
    inqkit::observables::OrbitalOverlapMatrix overlap_obs(electrons, wp_idx, "results/overlap");

    // Three RT density writers
    inqkit::io::RealField3DWriter total_wr("results/density_rt_total",
        {.field_name = "density", .include_meta = true}, {.overwrite = true});
    inqkit::io::RealField3DWriter jell_wr("results/density_rt_jellium",
        {.field_name = "density", .include_meta = true}, {.overwrite = true});
    inqkit::io::RealField3DWriter wp_wr("results/density_rt_wp",
        {.field_name = "density", .include_meta = true}, {.overwrite = true});

    // Write t=0 frames for all three series
    {
        auto rho_jell = inqkit::fields::density::total(electrons);
        auto rho_wp_  = inqkit::fields::density::orbital(electrons, wp_idx);
        auto rho_tot  = rho_jell;
        add_field_inplace(rho_tot, rho_wp_);
        total_wr.write(rho_tot,  0.0, 0);
        jell_wr.write(rho_jell,  0.0, 0);
        wp_wr.write(rho_wp_,     0.0, 0);
    }

    // Observables writer
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x  = sel.dipole_y  = sel.dipole_z  = true;
    inqkit::io::ObservablesWriter obs_writer("results/observables.csv", sel);
    obs_writer.write_header();

    // 20 LEED screen accumulators
    std::array<inqkit::screens::LeedPatternAccumulator, N_SCREENS> accums;
    for (int k = 0; k < N_SCREENS; ++k)
        accums[k] = inqkit::screens::LeedPatternAccumulator(
            inqkit::screens::PlaneScreen{SCREEN_Z[k], "screen_" + zero_pad(k, 2)});

    // Session A: density writes at WRITE_EVERY intervals
    inqkit::RealTimeSession rt(ions, electrons, WRITE_EVERY);
    rt.add([&](inqkit::StepContext const& ctx) {
        auto rho_jell = inqkit::fields::density::total(*ctx.electrons);
        auto rho_wp_  = inqkit::fields::density::orbital(*ctx.electrons, wp_idx);
        auto rho_tot  = rho_jell;
        add_field_inplace(rho_tot, rho_wp_);
        total_wr.write(rho_tot,  ctx.time_au, ctx.step);
        jell_wr.write(rho_jell,  ctx.time_au, ctx.step);
        wp_wr.write(rho_wp_,     ctx.time_au, ctx.step);
    });

    // Session B: observables, screens, overlap — every step
    inqkit::RealTimeSession rt_obs(ions, electrons, 1);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        obs_writer.append(ctx);
        for (auto& acc : accums) acc.accumulate(*ctx.electrons, DT_AU);
        overlap_obs.snapshot(*ctx.electrons, ctx.time_au, ctx.step);
        if (ctx.step % SCREEN_SNAP_EVERY == 0) {
            std::string step_dir = "results/screens_snapshots/step_" + zero_pad(ctx.step, 6);
            std::filesystem::create_directories(step_dir);
            for (int k = 0; k < N_SCREENS; ++k) {
                auto slice = accums[k].screen().extract(*ctx.electrons);
                accums[k].screen().save(slice, ctx.time_au,
                    step_dir + "/screen_" + zero_pad(k, 2) + ".dat");
            }
        }
    });

    real_time::propagate(ions, electrons,
        [&](auto const& data) { rt.step(data); rt_obs.step(data); },
        options::theory{}.lda(),
        options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime)
            .observables_current().observables_dipole());

    // Save final time-averaged LEED patterns
    std::filesystem::create_directories("results/screens");
    for (int k = 0; k < N_SCREENS; ++k)
        accums[k].save("results/screens/screen_" + zero_pad(k, 2) + ".dat");

    std::cout << "Done. Output in results/\n";
    return 0;
}
