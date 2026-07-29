// ============================================================================
// dryrun.cpp — pre-GS dry-run for the v2 cubic 50x50x50 / N=162 / dx=0.248
// setup, BEFORE committing to the full SCF.
//
// Same idea as v1 dryrun: construct cell + ions + electrons, print grid
// dimensions, exit. Skips ground_state::calculate().
//
// What it checks:
//   1. INQ accepts cell::cubic(50).periodic() with spacing 0.248.
//   2. systems::electrons() can be constructed at this size without OOM.
//      Expected dimensions ~200^3 = 8.0M points (INQ may round to FFT-
//      friendly sizes).
//   3. State count = N_e/2 + EXTRA_STATES = 81 + 20 = 101.
// ============================================================================

#include <inq/inq.hpp>
#include "../../../shared/configs/electron_proj_E1500_L50_cubic.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::Common_E1500_L50_cubic;

int main() {
    std::cout << std::setprecision(8);
    std::cout << "=== dryrun (v2): 50^3 cubic / N=162 / dx=0.248 ===\n\n";

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    std::cout << cell;

    auto ions = systems::ions(cell);
    std::cout << "ions.size() = " << ions.size() << " (expect 0)\n";

    std::cout << "Constructing electrons (this allocates the full grid)...\n";
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());

    std::cout << "  num_states     = " << electrons.states().num_states()
              << "  (expect " << (Cfg::N_ELECTRONS/2 + Cfg::EXTRA_STATES) << ")\n";
    std::cout << "  num_electrons  = " << electrons.states().num_electrons()
              << "  (expect " << Cfg::N_ELECTRONS << ")\n";

    int n_est = int(std::round(Cfg::L_BOHR / Cfg::SPACING_BOHR));
    std::cout << "  grid_estimate  = " << n_est << "^3 = "
              << (n_est*n_est*n_est) << " points (INQ may round up)\n";

    std::cout << "\n[dryrun v2] PASS — cell + electrons constructed without errors.\n";
    return 0;
}
