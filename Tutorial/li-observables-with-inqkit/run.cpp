#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>

using namespace inq;
using namespace inq::magnitude;

// Ionic kick velocity: 6.56 Ang/fs = 0.3000 a.u. (bohr/a.u.-time).
// 1 a.u. velocity = a0/t0 = 0.529177 Ang / 0.024189 fs = 21.877 Ang/fs.
static constexpr double KICK_VEL_AU = 0.3000;

int main() {
    // ── System ────────────────────────────────────────────────────────────────
    // BCC Li, a = 3.51 Å, 3×3×3 supercell (L = 10.53 Å), 54 atoms.
    // Gamma-point only; ionic velocity kick along +x after GS.
    auto cell = systems::cell::cubic(10.53_angstrom).periodic();
    auto ions = systems::ions::parse("li_bcc_3x3x3.xyz", cell);

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

    // ── Ground-state density ──────────────────────────────────────────────────
    auto rho_gs = inqkit::fields::density::total(electrons);
    inqkit::io::RealField3DWriter gs_density_writer(
        "results/gs_density",
        {.field_name = "gs_total_density", .include_meta = true},
        {.overwrite = true});
    gs_density_writer.write(rho_gs, "gs_density");
    std::cout << "GS density written to results/gs_density/\n";

    std::cout << "Kicking all ions: v_x = " << KICK_VEL_AU << " a.u. (6.56 Ang/fs)\n";

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
