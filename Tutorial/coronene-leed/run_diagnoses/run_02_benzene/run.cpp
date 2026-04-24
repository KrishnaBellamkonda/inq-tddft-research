// ============================================================================
// Diagnostic run 02: benzene ground state — cross-run reference system
//
// If the orbital density cross artifact in the coronene LEED simulations is
// caused by near-degenerate frontier orbitals specific to coronene's 54-state
// D6h system, the artifact should be absent (or weaker) for benzene.
//
// Benzene shares the same D6h point group as coronene but has only 30 electrons
// and 15 occupied states. Its HOMO-LUMO gap is ~6 eV under LDA — far larger
// than coronene's — so degenerate-orbital mixing at loose convergence is much
// less likely. Using the same tight SCF settings as run_01_tight_scf provides
// a clean like-for-like comparison.
//
// If benzene's orbital densities show the same cross artifact -> grid artefact
// If benzene's orbital densities look physically correct -> coronene-specific issue
//
// Geometry: C6H6 D6h, C-C = 1.397 Ang, C-H = 1.086 Ang,
//           centred at (LX/2, LY/2, LZ/2) = (9.2408, 9.2408, 15.8496) Ang.
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
    std::cout << "\n=== diagnostic run 02: benzene C6H6 tight SCF (1e-8 Ha) ===\n";

    auto cell = systems::cell::orthorhombic(
        LX_BOHR * 1.0_b, LY_BOHR * 1.0_b, LZ_BOHR * 1.0_b
    ).finite();

    auto ions = systems::ions::parse("benzene.xyz", cell);
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
        summary << "run = diagnostic_02_benzene\n";
        summary << "system = benzene_C6H6\n";
        summary << "geometry_file = benzene.xyz\n";
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
