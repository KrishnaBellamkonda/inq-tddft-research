// ============================================================================
// jellium-wp-rt / run_01_base: WP + N=40 jellium, 200 eV, sigma=0.53A, +z
//
// System: cubic L=40 bohr periodic cell, 40 electrons (no ions = jellium),
//         r_s = 7.26 a0.  extra_states(3) for Broyden SCF stability.
//         Temperature smearing kT = 0.00862 eV.
//
// WP: sigma=0.53A, E_kin=200 eV, center=(L/2, L/2, 5*sigma), direction +z.
//     Orthogonalised against occupied KS states.
//
// TDDFT: dt=0.02 a.u., 10000 steps ≈ 4.83 fs
//   density every 100 steps, observables every step, 3 screens (z=10,20,30)
// ============================================================================

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
static constexpr double L_BOHR      = 40.0;

static constexpr double WP_SIGMA_ANG  = 0.53;
static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;
static constexpr double WP_EKIN_EV    = 200.0;
static constexpr double WP_EKIN_HA    = WP_EKIN_EV / HA_TO_EV;
static const     double WP_K0         = std::sqrt(2.0 * WP_EKIN_HA);

static const double WP_CX = L_BOHR / 2.0;
static const double WP_CY = L_BOHR / 2.0;
static const double WP_CZ = 5.0 * WP_SIGMA_BOHR;

static constexpr int    N_STEPS     = 417;
static constexpr double DT_AU       = 0.02;
static constexpr int    WRITE_EVERY = 2;

static void add_field_inplace(inqkit::fields::RealField3D & a,
                               inqkit::fields::RealField3D const& b) {
    for (std::size_t i = 0; i < a.values.size(); i++) a.values[i] += b.values[i];
}

int main() {
    std::cout << "\n=== jellium run_01_base: N=40, 200 eV, sigma=0.53A, +z ===\n";

    auto cell = systems::cell::cubic(L_BOHR*1.0_b).periodic();
    auto ions = systems::ions(cell);

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(0.50*1.0_b)
            .extra_electrons(38)
            .extra_states(3)
            .temperature(0.00862*1.0_eV),
        input::kpoints::gamma());

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(1e-4_Ha)
            .max_steps(300)
            .broyden_mixing()
            .mixing_ndim(8)
            .mixing(0.1));

   // TODO: Save ground state 
   // 1. Total density
   // 2. Orbital density

    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    auto wp = inqkit::WavePacket{}
        .center(WP_CX, WP_CY, WP_CZ)
        .sigma(WP_SIGMA_BOHR)
        .k0(0.0, 0.0, WP_K0)
        .orthogonalise_against_occupied(electrons);

    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    std::cout << "  state_index = " << report.state_index
              << "  norm_after = "  << report.norm_after
              << "  max_overlap = " << report.max_overlap << "\n";


    inqkit::io::RealField3DWriter density_writer("results/density_rt",
        {.field_name = "density", .include_meta = true}, {.overwrite = true});
    // t=0 density: jellium GS (40 e-) + WP orbital
    auto rho_t0 = inqkit::fields::density::total(electrons);
    auto rho_wp = inqkit::fields::density::orbital(electrons, report.state_index);
    


    // TODO: Calculate each oritals overlap with initial ground state orbitals
    // and save them in the results folder, under a subfolder with an appropriate name
    // This should be written as a class in inqview library (I will review the code myself)
    // to ensure I've understood what is being done clearly. For every KS orbital at every timestep, 
    // the overlap of the said orbital with all the KS ground state orbitals is calculated. This is
    // repeated for the entire duration of the simulation run. Then, a GIF is made for each orbital using the overlap
    //  matrices. The range on the x and y axis must clearly identify similar/degenerate (or near degenerate) orbitals
    
    add_field_inplace(rho_t0, rho_wp);
    density_writer.write(rho_t0, 0.0, 0);

    // TODO: Write the total electronic density, the target system electornic density (given by fields::density::total, 
    // as currently, it only adds up the KS orbitals that were occupied before the wavefunction was added) and
    // wave packet density. 
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x  = sel.dipole_y  = sel.dipole_z  = true;
    inqkit::io::ObservablesWriter obs_writer("results/observables.csv", sel);
    obs_writer.write_header();

    // Add 20 screens (almost evenly spaced, with a small random offset, so as to be secure from falling 
    // into some resonance effects due to the size of the box.
    // Screens at z=10, 20, 30 bohr (L/4, L/2, 3L/4)
    inqkit::screens::LeedPatternAccumulator sc1(inqkit::screens::PlaneScreen{10.0, "screen_z10"});
    inqkit::screens::LeedPatternAccumulator sc2(inqkit::screens::PlaneScreen{20.0, "screen_z20"});
    inqkit::screens::LeedPatternAccumulator sc3(inqkit::screens::PlaneScreen{30.0, "screen_z30"});

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
    sc1.save("results/screens/screen_z10.dat");
    sc2.save("results/screens/screen_z20.dat");
    sc3.save("results/screens/screen_z30.dat");

    std::cout << "Done. Output in results/\n";
    return 0;
}
