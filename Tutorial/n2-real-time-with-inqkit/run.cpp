#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>

using namespace inq;
using namespace inq::magnitude;

int main() {
  // ── System ─────────────────────────────────────────────────────────────────
  // N2 slightly stretched (1.15 Å vs equilibrium 1.098 Å), centred in a
  // 20-bohr finite cubic cell.  Atoms at (L/2, L/2, L/2 ± half_bond).
  // INQ uses (0,0,0) as the cell corner, so L/2 = 10 bohr is the true
  // centre.

  auto L = 20.0_bohr;
  auto cell = systems::cell::cubic(L).finite();
  auto cx = 5.29177_angstrom;

  auto bond_length = 1.11946_angstrom;
  auto half_bond = (bond_length) / 2; // 1.15 Å / 2
  auto zero = 0.0_bohr;

  systems::ions ions(cell);
  ions.insert("N", {cx, cx, cx-half_bond});
  ions.insert("N", {cx, cx, cx+half_bond});

  // E_cut = 20 Ha (40 Ry): coarse but stable at dt = 0.1 au.
  systems::electrons electrons(ions, options::electrons{}.cutoff(30.0_Ha));

  // ── Ground state ───────────────────────────────────────────────────────────
  ground_state::initial_guess(ions, electrons);
  auto gs = ground_state::calculate(
      ions, electrons, options::theory{}.pbe(),
      options::ground_state{}.energy_tolerance(1e-8_Ha).mixing(0.01).max_steps(
          1000));

  std::cout << "GS total energy = " << gs.energy.total() << " Ha\n";

  // ── Real-time setup ────────────────────────────────────────────────────────
  // Write total density every 100 steps → 61 frames (steps 0, 100, …, 6000).
  /* The idea of the session and context manager is simple yet
   * powerful. The RealTimeSession is something to which we give
   * access to the ground state ions, electrons and write_every
   * argument. The instance of RealTimeSession has function named
   * .step(). The step command, re-assigns the ions and electrons
   * attributes to the context class. To the real time session,
   * we add a task, which is essentially a function that handles
   * what we want to do with the data, like writing it. RealTimeSessi
   * -on takes in an array of tasks that are completed one after the
   * other.
   * */

  const int write_every = 100;
  inqkit::RealTimeSession rt(ions, electrons, write_every);

  inqkit::io::RealField3DWriter density_writer(
      "results/real_time/density",
      {.field_name = "total_density", .include_meta = true});

  rt.add([&density_writer](inqkit::StepContext const &ctx) {
    auto rho = inqkit::fields::density::total(*ctx.electrons);
    density_writer.write(rho, ctx.time_au, ctx.step);
    std::cout << "  wrote density t=" << ctx.time_au << " au  step=" << ctx.step
              << "\n";
  });

  // ── Propagation ────────────────────────────────────────────────────────────
  // Electronic delta-kick in z (along molecular axis) with amplitude 0.01
  // bohr^-1. Free propagation afterwards reveals electronic density
  // oscillations.
  auto kick = perturbations::kick{cell, {0.0, 0.0, 0.01}};

  real_time::propagate(
      ions, electrons, [&](auto const &data) { rt.step(data); },
      options::theory{}.lda(),
      options::real_time{}.num_steps(6000).dt(0.05_atomictime), kick);

  std::cout << "Done. Wrote 61 density frames (t=0..300 au) to "
               "results/real_time/density/\n";

  return 0;
}
