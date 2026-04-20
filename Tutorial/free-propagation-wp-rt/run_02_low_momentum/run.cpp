// run_02_low_momentum: free WP, 50 eV, sigma=0.53A, -z
#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>
#include <cmath>
#include <filesystem>
#include <iostream>
using namespace inq;
using namespace inq::magnitude;

static constexpr double ANG_TO_BOHR = 1.8897259886;
static constexpr double HA_TO_EV    = 27.21138625;
static constexpr double LX_BOHR = 34.771;
static constexpr double LY_BOHR = 34.771;
static constexpr double LZ_BOHR = 89.856;
static constexpr double WP_SIGMA_ANG  = 0.53;
static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;
static constexpr double WP_EKIN_EV    = 50.0;
static constexpr double WP_EKIN_HA    = WP_EKIN_EV / HA_TO_EV;
static const     double WP_K0         = std::sqrt(2.0 * WP_EKIN_HA);
static const double WP_CX = LX_BOHR / 2.0;
static const double WP_CY = LY_BOHR / 2.0;
static const double WP_CZ = LZ_BOHR - 5.0 * WP_SIGMA_BOHR;
static constexpr int    N_STEPS     = 10000;
static constexpr double DT_AU       = 0.02;
static constexpr int    WRITE_EVERY = 100;

static void add_field_inplace(inqkit::fields::RealField3D & a,
                               inqkit::fields::RealField3D const& b) {
    for (std::size_t i = 0; i < a.values.size(); i++) a.values[i] += b.values[i];
}

int main() {
    std::cout << "\n=== run_02_low_momentum: free WP, 50 eV, sigma=0.53A, -z ===\n";
    auto cell = systems::cell::orthorhombic(LX_BOHR*1.0_b, LY_BOHR*1.0_b, LZ_BOHR*1.0_b).finite();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(ions, options::electrons{}.cutoff(40.0_Ha).extra_states(1));
    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons, options::theory{}.lda(),
        options::ground_state{}.energy_tolerance(1e-4_Ha).max_steps(10));

    auto wp = inqkit::WavePacket{}.center(WP_CX, WP_CY, WP_CZ).sigma(WP_SIGMA_BOHR).k0(0.0, 0.0, -WP_K0);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    std::cout << "  state_index = " << report.state_index << "  norm_after = " << report.norm_after << "\n";

    inqkit::io::RealField3DWriter density_writer("results/density_rt",
        {.field_name = "density", .include_meta = true}, {.overwrite = true});
    auto rho_t0 = inqkit::fields::density::total(electrons);
    auto rho_wp = inqkit::fields::density::orbital(electrons, report.state_index);
    add_field_inplace(rho_t0, rho_wp);
    density_writer.write(rho_t0, 0.0, 0);

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x  = sel.dipole_y  = sel.dipole_z  = true;
    inqkit::io::ObservablesWriter obs_writer("results/observables.csv", sel);
    obs_writer.write_header();

    inqkit::screens::LeedPatternAccumulator sc1(inqkit::screens::PlaneScreen{LZ_BOHR*0.25, "screen_z25"});
    inqkit::screens::LeedPatternAccumulator sc2(inqkit::screens::PlaneScreen{LZ_BOHR*0.50, "screen_z50"});
    inqkit::screens::LeedPatternAccumulator sc3(inqkit::screens::PlaneScreen{LZ_BOHR*0.75, "screen_z75"});

    inqkit::RealTimeSession rt(ions, electrons, WRITE_EVERY);
    rt.add([&](inqkit::StepContext const& ctx) {
        density_writer.write(inqkit::fields::density::total(*ctx.electrons), ctx.time_au, ctx.step);
    });
    inqkit::RealTimeSession rt_obs(ions, electrons, 1);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        obs_writer.append(ctx);
        sc1.accumulate(*ctx.electrons, DT_AU);
        sc2.accumulate(*ctx.electrons, DT_AU);
        sc3.accumulate(*ctx.electrons, DT_AU);
    });

    real_time::propagate(ions, electrons,
        [&](auto const& data) { rt.step(data); rt_obs.step(data); },
        options::theory{}.lda(),
        options::real_time{}.num_steps(N_STEPS).dt(DT_AU*1.0_atomictime)
            .observables_current().observables_dipole());

    std::filesystem::create_directories("results/screens");
    sc1.save("results/screens/screen_z25.dat");
    sc2.save("results/screens/screen_z50.dat");
    sc3.save("results/screens/screen_z75.dat");
    std::cout << "Done.\n";
    return 0;
}
