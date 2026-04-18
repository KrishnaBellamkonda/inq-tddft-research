// Tutorial: Li BCC — Ground State + TDDFT Ionic Velocity Kick
// =============================================================
// Reproduces the type of calculation in:
//   QuantumKickExtension/qball-codebase/Li/td_kicks/
// using INQ instead of QBall.
//
// SYSTEM
//   Li body-centred cubic (BCC), conventional 2-atom cubic unit cell.
//   The QBall production runs use a 54-atom 3×3×3 supercell with
//   Gamma-only k-sampling. This tutorial is equivalent: the 2-atom
//   primitive cell with a 4×4×4 k-grid gives the same physics
//   at a fraction of the cost.
//
// PROTOCOL
//   1. Converge the ground state (SCF) with Fermi smearing.
//      Li is a metal — smearing is essential for SCF convergence.
//   2. Apply an ionic velocity kick: all atoms get velocity v in +x.
//   3. Propagate in real time with impulsive ionic dynamics:
//      ions move at constant velocity while electrons respond via TDDFT.
//   4. Record energy and electronic current at each time step.
//
// QBall → INQ translation:
//   atoms_dyn IMPULSIVE          →  ionic::propagator::impulsive
//   wf_dyn FORKTD                →  real_time::propagate (ETRS, default)
//   set_velocity Li_i v 0 0      →  ions.velocities()[i] = {v, 0, 0}
//   smearing FERMI 1000 kelvin   →  .temperature(0.086_eV) (kB × 1000 K)
//   nempty 20                    →  .extra_states(4)       (small cell)
//   ecut 74 rydberg              →  .cutoff(40.0_Ry)       (tutorial)
//
// Reference:
//   Santervás-Arranz, Stengel, Artacho, Phys. Rev. Research 7, 033292 (2025)

#include <inq/inq.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>

int main(){
    using namespace inq;
    using namespace inq::magnitude;

    // ── 1. Crystal geometry ──────────────────────────────────────────────────
    // BCC Li, lattice constant a = 3.51 Å.
    // Conventional cubic cell: 2 atoms at fractional coords (0,0,0) and (½,½,½).

    auto a = 3.51_angstrom;

    systems::ions ions(systems::cell::cubic(a).periodic());
    ions.insert_fractional("Li", {0.0, 0.0, 0.0});
    ions.insert_fractional("Li", {0.5, 0.5, 0.5});

    // ── 2. Electronic structure ──────────────────────────────────────────────
    // Li is metallic → Fermi smearing is required for SCF convergence.
    // kT = kB × 1000 K ≈ 0.086 eV  (matches QBall's smearing_width).
    //
    // Cutoff: 40 Ry (tutorial). QBall production: 74 Ry.
    // k-grid: 4×4×4 Monkhorst-Pack (replaces QBall's Gamma + 54-atom supercell).

    auto kT_1000K = 0.086_eV;

    systems::electrons electrons(ions,
        options::electrons{}
            .cutoff(40.0_Ry)
            .extra_states(4)
            .temperature(kT_1000K),
        input::kpoints::grid({4, 4, 4}, true));

    // ── 3. Ground state SCF ──────────────────────────────────────────────────
    ground_state::initial_guess(ions, electrons);

    auto gs = ground_state::calculate(ions, electrons,
        options::theory{}.pbe(),
        options::ground_state{}.energy_tolerance(1e-6_Ha));

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n=== Ground State (BCC Li, 2 atoms) ===\n";
    std::cout << "  Total energy  : " << gs.energy.total()   << " Ha\n";
    std::cout << "  Kinetic       : " << gs.energy.kinetic() << " Ha\n";
    std::cout << "  XC (PBE)      : " << gs.energy.xc()     << " Ha\n";
    std::cout << "  Ion-ion       : " << gs.energy.ion()     << " Ha\n";
    std::cout << "  SCF iters     : " << gs.total_iter       << "\n\n";

    // Save the ground state so the TD run can be cheaply restarted.
    electrons.save("results/li_bcc_gs");

    // ── 4. Ionic velocity kick ────────────────────────────────────────────────
    // All atoms receive an impulsive velocity in the +x direction.
    // Units: bohr / a.u. of time  (same as QBall's 'atomicvelocity').
    // v = 0.04 a.u. is a moderate kick, typical of QBall production runs.

    double v_kick = 0.04;
    for(int i = 0; i < ions.size(); i++)
        ions.velocities()[i] = vector3<double>{v_kick, 0.0, 0.0};

    std::cout << "=== TDDFT Propagation ===\n";
    std::cout << "  Kick velocity  : " << v_kick << " a.u. (+x)\n";
    std::cout << "  Time step      : 0.04 a.u. ≈ 0.97 as\n";
    std::cout << "  Steps          : 2000  →  ~1.94 fs total\n\n";

    // ── 5. Real-time TDDFT propagation ───────────────────────────────────────
    // ionic::propagator::impulsive: ions move at constant velocity,
    //   no forces applied to them. Matches QBall's 'atoms_dyn IMPULSIVE'.
    //
    // The electrons respond self-consistently to the moving ions via TDDFT
    //   using the ETRS propagator (default, efficient for metals).
    //
    // Columns written to results/tddft.dat:
    //   step | time[atu] | E_total[Ha] | E_kin_elec[Ha] | J_x | J_y | J_z | x_atom0[bohr]

    std::ofstream out("results/tddft.dat");
    out << "# step   time[atu]   E_total[Ha]   E_kin[Ha]   J_x   J_y   J_z   x_atom0[bohr]\n";

    auto output = [&](auto data){
        if(data.every(10)){
            auto cur = data.current();
            auto pos = data.positions();
            out << std::setw(5) << data.iter()           << "  "
                << std::fixed << std::setprecision(5)
                << data.time()                           << "  "
                << data.energy().total()                 << "  "
                << data.energy().kinetic()               << "  "
                << cur[0]                                << "  "
                << cur[1]                                << "  "
                << cur[2]                                << "  "
                << pos[0][0]                             << "\n";
            out.flush();
        }
    };

    real_time::propagate(
        ions, electrons, output,
        options::theory{}.pbe(),
        options::real_time{}
            .num_steps(2000)
            .dt(0.04_atomictime)
            .impulsive()
            .observables_current()
    );

    std::cout << "\nTDDFT complete. Results written to results/tddft.dat\n";

    return 0;
}
