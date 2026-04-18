// Include <inq/inq.hpp> for adding the inq library
#include <inq/inq.hpp>
#include <iostream>
#include <iomanip>

int main() {
    using namespace inq;
    using namespace inq::magnitude;

    auto half_bond = 0.37_angstrom;   // experimental H-H bond = 0.741 Å

    // ── 1. Cell & geometry ──────────────────────────────────────────────────
    // A 20 bohr cube with open (non-periodic) boundary conditions.
    // The cell must be large enough so the molecule doesn't interact with
    // its periodic images (~5 Å vacuum on each side is usually sufficient).
    systems::ions ions(systems::cell::cubic(20.0_bohr).finite());
    ions.insert("H", {0.0_angstrom, 0.0_angstrom, -half_bond});
    ions.insert("H", {0.0_angstrom, 0.0_angstrom,  half_bond});

    // ── 2. Electronic structure ─────────────────────────────────────────────
    // 80 Ry plane-wave cutoff — well converged for H ONCV pseudopotentials.
    systems::electrons electrons(ions, options::electrons{}.cutoff(80.0_Ry));

    // ── 3. SCF ──────────────────────────────────────────────────────────────
    ground_state::initial_guess(ions, electrons);
    auto result = ground_state::calculate(ions, electrons, options::theory{}.pbe(),
                                         options::ground_state{}.calculate_forces());

    // ── 4. Results ──────────────────────────────────────────────────────────
    // The << operator prints a compact summary; we also pull out individual
    // terms to understand the energy decomposition.

    std::cout << std::fixed << std::setprecision(6);

    std::cout << "\n=== Energy breakdown (Ha) ===\n";
    std::cout << "  Total        : " << result.energy.total()       << "\n";
    std::cout << "  Kinetic      : " << result.energy.kinetic()     << "\n";
    std::cout << "  Hartree      : " << result.energy.hartree()     << "\n";
    std::cout << "  XC (PBE)     : " << result.energy.xc()         << "\n";
    std::cout << "  External     : " << result.energy.external()    << "\n";
    std::cout << "  Non-local PP : " << result.energy.non_local()   << "\n";
    std::cout << "  Ion-ion      : " << result.energy.ion()         << "\n";
    std::cout << "  SCF iters    : " << result.total_iter           << "\n";

    // Convert total energy to eV for reference
    double total_eV = result.energy.total() / in_atomic_units(1.0_eV);
    std::cout << "\n  Total        : " << total_eV << " eV\n";

    // ── 5. Forces ───────────────────────────────────────────────────────────
    // Forces are in Ha/bohr. For H2 along z, F_z on atom 0 should be
    // positive (pushed away) and equal/opposite on atom 1 at equilibrium.
    std::cout << "\n=== Forces (Ha/bohr) ===\n";
    for(int i = 0; i < ions.size(); i++){
        std::cout << "  H[" << i << "]  Fx=" << result.forces[i][0]
                            << "  Fy=" << result.forces[i][1]
                            << "  Fz=" << result.forces[i][2] << "\n";
    }

    // ── 6. Dipole ───────────────────────────────────────────────────────────
    // H2 is homonuclear so the dipole should be ~0. Good sanity check.
    std::cout << "\n=== Dipole moment (a.u.) ===\n";
    std::cout << "  d = (" << result.dipole[0] << ", "
                           << result.dipole[1] << ", "
                           << result.dipole[2] << ")\n";

    return 0;
}
