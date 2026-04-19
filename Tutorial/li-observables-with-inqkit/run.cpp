#include <inq/inq.hpp>

#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/io/observables_writer.hpp>

using namespace inq;
using namespace inq::magnitude;

// Ionic kick velocity in a.u. (bohr / a.u. time).
// Reference: QBall Li series, smallest kick — linear-response regime.
// Santervás-Arranz, Stengel, Artacho, Phys. Rev. Research 7, 033292 (2025).
static constexpr double KICK_VEL_AU = 0.0123;

int main() {
    // ── System ────────────────────────────────────────────────────────────────
    // BCC Li, a = 3.51 Å, 2×2×2 supercell (L = 7.02 Å), 16 atoms.
    // Gamma-point only; ionic velocity kick along +x after GS.
    auto cell = systems::cell::cubic(7.02_angstrom).periodic();
    auto ions = systems::ions::parse("li_bcc_2x2x2.xyz", cell);

    // Metal: needs Fermi smearing + extra empty states
    systems::electrons electrons(
        ions,
        options::electrons{}
            .cutoff(30.0_Ry)
            .extra_states(4)
            .temperature(0.001_Ha));

    // ── Ground state ──────────────────────────────────────────────────────────
    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}.energy_tolerance(1e-6_Ha));

    auto gs_energy = gs.energy.total();
    std::cout << "GS total energy = " << gs_energy << " Ha\n";
    std::cout << "Kicking all ions: v_x = " << KICK_VEL_AU << " a.u.\n";

    // ── Ionic kick ────────────────────────────────────────────────────────────
    for (int i = 0; i < ions.size(); ++i)
        ions.velocities()[i] = inq::vector3<double>{KICK_VEL_AU, 0.0, 0.0};

    // ── Real-time setup ───────────────────────────────────────────────────────
    const int write_every = 10;
    inqkit::RealTimeSession rt(ions, electrons, write_every);

    // All observables enabled (dipole_x included for spectral analysis)
    inqkit::io::ObservableSelection sel;
    sel.dipole_x = true;
    sel.dipole_y = true;
    sel.dipole_z = true;
    sel.energy_hartree = true;
    sel.energy_xc      = true;

    inqkit::io::ObservablesWriter obs_writer("results/observables.csv", sel);
    obs_writer.write_header();

    rt.add([&](inqkit::StepContext const& ctx) {
        obs_writer.append(ctx);
    });
    rt.add([&](inqkit::StepContext const& ctx) {
        if (ctx.step % 200 == 0) {
            double drift = (ctx.energy_total - gs_energy) / std::abs(gs_energy) * 100.0;
            std::cout << "  step=" << ctx.step
                      << "  t=" << ctx.time_au << " au"
                      << "  E=" << ctx.energy_total << " Ha"
                      << "  drift=" << drift << "%\n";
        }
    });

    // ── Propagation ───────────────────────────────────────────────────────────
    // .impulsive(): ions move at constant velocity; no force evaluation
    // .observables_current() + .observables_dipole(): enable observables in viewables
    real_time::propagate(
        ions, electrons,
        [&](auto const& data) { rt.step(data); },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(2000)
            .dt(0.05_atomictime)
            .impulsive()
            .observables_current()
            .observables_dipole());

    std::cout << "Done. Results written to results/observables.csv\n";
    return 0;
}
