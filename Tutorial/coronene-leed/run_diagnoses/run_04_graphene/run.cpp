// ============================================================================
// Diagnostic run 04: small graphene nanoflake — extended honeycomb reference
//
// Coronene is the smallest fully aromatic PAH with 7 hexagonal rings. To test
// whether the cross artifact is specific to this molecule or appears for any
// in-plane carbon system with the graphene honeycomb lattice, this run uses a
// small graphene nanoflake: all C sites within r < 5.0 Ang of the patch centre,
// taken from the graphene primitive lattice (a = 2.46 Ang, no hydrogen
// passivation). The patch contains 31 carbon atoms.
//
// Bare graphene edges carry dangling bonds and metallic-like edge states.
// extra_states is increased to 12 to accommodate these near-Fermi unoccupied
// levels. This run is a diagnostic, not a model of physical graphene.
//
// Geometry: 31 C atoms, centred at (9.2408, 9.2408, 15.8496) Ang.
// Lattice vectors used to generate the patch:
//   a1 = (2.46, 0.000) Ang  (armchair direction)
//   a2 = (1.23, 2.131) Ang
//   Basis: A at (0,0), B at (1.23, 0.711) per unit cell
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
    std::cout << "\n=== diagnostic run 04: graphene nanoflake (31 C, no H) ===\n";

    auto cell = systems::cell::orthorhombic(
        LX_BOHR * 1.0_b, LY_BOHR * 1.0_b, LZ_BOHR * 1.0_b
    ).finite();

    auto ions = systems::ions::parse("graphene_nanoflake.xyz", cell);
    std::cout << "  Atoms: " << ions.size() << "\n";

    // extra_states=12: graphene edge states produce metallic-like near-Fermi levels
    auto electrons = systems::electrons(
        ions,
        options::electrons{}.cutoff(54.0_Ha).extra_states(12)
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
        summary << "run = diagnostic_04_graphene_nanoflake\n";
        summary << "system = graphene_nanoflake_31C\n";
        summary << "geometry_file = graphene_nanoflake.xyz\n";
        summary << "cell_bohr = " << LX_BOHR << ' ' << LY_BOHR << ' ' << LZ_BOHR << "\n";
        summary << "boundary = finite\n";
        summary << "xc = pbe\n";
        summary << "cutoff_ha = 54.0\n";
        summary << "extra_states = 12\n";
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
