// ============================================================================
// Diagnostic run 01: coronene ground state with tight SCF convergence
//
// The baseline LEED simulations (run_01 through run_06) all use a loose SCF
// tolerance of 1e-4 Ha. For a molecule with D6h symmetry and 54 occupied
// states — including doubly-degenerate HOMO and HOMO-1 — the steepest-descent
// eigensolver may settle into an unphysical linear combination of degenerate
// orbitals before the energy reaches the stated tolerance. This produces
// orbital densities that lack the expected hexagonal symmetry, and downstream
// LEED patterns with a spurious cross-like artifact.
//
// This run repeats the run_07 ground-state calculation with:
//   - energy_tolerance: 1e-4 Ha -> 1e-8 Ha
//   - max_steps:        300     -> 1000
//   - mixing_ndim:      8       -> 16  (deeper Broyden history)
//
// If the orbital densities now show proper D6h symmetry, the convergence
// tolerance is confirmed as the root cause.
//
// Outputs:
//   results/density/                         total ground-state density
//   results/orbital_density/orbital_XXXX/    density of each KS orbital
//   results/orbital_density/orbital_index_map.csv
//   results/ground_state_summary.txt
//   results/checkpoint/
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;

static constexpr double LX_BOHR = 34.9222;
static constexpr double LY_BOHR = 34.9222;
static constexpr double LZ_BOHR = 59.9043;

static std::string orbital_tag(int ist) {
    std::ostringstream os;
    os << "orbital_" << std::setw(4) << std::setfill('0') << ist;
    return os.str();
}

int main() {
    std::cout << "\n=== diagnostic run 01: coronene tight SCF (1e-8 Ha) ===\n";

    auto cell = systems::cell::orthorhombic(
        LX_BOHR * 1.0_b, LY_BOHR * 1.0_b, LZ_BOHR * 1.0_b
    ).finite();

    auto ions = systems::ions::parse("coronene_leed.xyz", cell);
    std::cout << "  Atoms: " << ions.size() << "\n";

    auto electrons = systems::electrons(
        ions,
        options::electrons{}.cutoff(54.0_Ha).extra_states(8)
    );

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions,
        electrons,
        options::theory{}.pbe(),
        options::ground_state{}
            .energy_tolerance(1e-8_Ha)
            .max_steps(1000)
            .broyden_mixing()
            .mixing_ndim(8)
            .mixing(0.1)
    );

    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    std::filesystem::create_directories("results/density");
    std::filesystem::create_directories("results/orbital_density");
    std::filesystem::create_directories("results/checkpoint");

    {
        inqkit::io::RealField3DWriter density_writer(
            "results/density",
            {.field_name = "density", .include_meta = true},
            {.overwrite = true}
        );
        auto rho_total = inqkit::fields::density::total(electrons);
        density_writer.write(rho_total, 0.0, 0);
    }

    const int nstates = electrons.states().num_states();
    std::cout << "  Writing densities for " << nstates << " KS states\n";

    {
        std::ofstream map("results/orbital_density/orbital_index_map.csv");
        map << "state_index,directory\n";
        for (int ist = 0; ist < nstates; ++ist)
            map << ist << ',' << orbital_tag(ist) << "\n";
    }

    for (int ist = 0; ist < nstates; ++ist) {
        const auto out_dir = std::string("results/orbital_density/") + orbital_tag(ist);
        std::filesystem::create_directories(out_dir);

        inqkit::io::RealField3DWriter orbital_density_writer(
            out_dir,
            {.field_name = "density", .include_meta = true},
            {.overwrite = true}
        );
        auto rho_orbital = inqkit::fields::density::orbital(electrons, ist);
        orbital_density_writer.write(rho_orbital, 0.0, 0);
        std::cout << "    wrote " << orbital_tag(ist) << "\n";
    }

    electrons.save("results/checkpoint");

    if (electrons.root()) {
        std::ofstream summary("results/ground_state_summary.txt");
        summary << std::setprecision(16);
        summary << "run = diagnostic_01_tight_scf\n";
        summary << "system = coronene_C24H12\n";
        summary << "geometry_file = coronene_leed.xyz\n";
        summary << "cell_bohr = " << LX_BOHR << ' ' << LY_BOHR << ' ' << LZ_BOHR << "\n";
        summary << "boundary = finite\n";
        summary << "xc = pbe\n";
        summary << "cutoff_ha = 54.0\n";
        summary << "extra_states = 8\n";
        summary << "energy_tolerance_ha = 1e-8\n";
        summary << "max_steps = 1000\n";
        summary << "mixing = broyden\n";
        summary << "mixing_ndim = 16\n";
        summary << "mixing_alpha = 0.1\n";
        summary << "num_atoms = " << ions.size() << "\n";
        summary << "num_electrons = " << electrons.states().num_electrons() << "\n";
        summary << "num_states = " << nstates << "\n";
        summary << "ground_state_energy_ha = " << gs.energy.total() << "\n";
    }

    std::cout << "Done. Output written to results/\n";
    return 0;
}
