// ============================================================================
// dryrun.cpp — pre-GS dry-run for the orthorhombic 40 x 40 x 150 / N=162 /
// dx=0.30 setup, BEFORE committing to the full SCF (which is hours of GPU).
//
// What it checks:
//   1. INQ accepts cell::orthorhombic(40, 40, 150).periodic().
//   2. systems::electrons() can be constructed at this size without OOM
//      or grid-shape errors. Dimensions reported should be ~134 x 134 x 500.
//   3. The state count (n_states = N_e/2 + EXTRA_STATES = 81 + 20 = 101)
//      is what we expect.
//
// What it does NOT do: ground_state::calculate(). Skipping the SCF saves
// hours of compute. If this dryrun passes, we know the setup is sane and
// can launch the real save_gs/gs_L40x40x150_orth_N162_dx0p30/run.cpp.
// ============================================================================

#include <inq/inq.hpp>
#include "../../../shared/configs/electron_proj_E1000_L40x40x150.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::Common_E1000_L40x40x150;

int main() {
    // Note: do NOT call `input::environment{}` here — INQ's
    // systems::electrons() constructor initialises MPI itself, and a manual
    // `input::environment{}` call would cause MPI_Init to be invoked twice,
    // which OpenMPI rejects as a fatal error.
    std::cout << std::setprecision(8);
    std::cout << "=== dryrun: 40x40x150 / N=162 / dx=0.30 setup check ===\n\n";

    auto cell = systems::cell::orthorhombic(
        Cfg::LX_BOHR * 1.0_b,
        Cfg::LY_BOHR * 1.0_b,
        Cfg::LZ_BOHR * 1.0_b).periodic();
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

    // Estimated grid points (independent of INQ — just a calculation).
    int nx_est = int(std::round(Cfg::LX_BOHR / Cfg::SPACING_BOHR));
    int ny_est = int(std::round(Cfg::LY_BOHR / Cfg::SPACING_BOHR));
    int nz_est = int(std::round(Cfg::LZ_BOHR / Cfg::SPACING_BOHR));
    std::cout << "  grid_estimate  = " << nx_est << " x " << ny_est
              << " x " << nz_est << " = " << (nx_est*ny_est*nz_est)
              << " points\n";

    std::cout << "\n[dryrun] PASS — cell + electrons constructed without errors.\n";
    return 0;
}
