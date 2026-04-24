// ============================================================================
// Diagnostic run 03: coronene ground state with hardcoded atom coordinates
//
// The previous runs load coronene geometry from an .xyz file via
// systems::ions::parse(). Although the INQ parser is verified correct
// (parse/xyz.hpp reads Angstroms and converts to Bohr via in_atomic_units()),
// this run eliminates the parser entirely by inserting all 36 coronene atoms
// directly in C++ using the _angstrom magnitude suffix.
//
// If the orbital densities from this run match run_01_tight_scf exactly, the
// xyz parser is confirmed innocent. If they differ, the parser is introducing
// an offset or misassignment.
//
// Coordinates taken verbatim from coronene_leed.xyz (C24H12 D6h, centred at
// (9.2408, 9.2408, 15.8496) Ang in the 18.48 x 18.48 x 31.7 Ang cell).
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
    std::cout << "\n=== diagnostic run 03: coronene hardcoded coordinates ===\n";

    auto cell = systems::cell::orthorhombic(
        LX_BOHR * 1.0_b, LY_BOHR * 1.0_b, LZ_BOHR * 1.0_b
    ).finite();

    // Coronene C24H12 — all 36 atoms from coronene_leed.xyz hardcoded.
    // Coordinates in Angstroms; INQ converts to Bohr internally.
    auto ions = systems::ions(cell);

    // Carbons (C1-C24)
    ions.insert("C", {10.661840_angstrom,  9.240840_angstrom, 15.849600_angstrom});
    ions.insert("C", {12.082840_angstrom,  9.240840_angstrom, 15.849600_angstrom});
    ions.insert("C", { 9.951340_angstrom, 10.471462_angstrom, 15.849600_angstrom});
    ions.insert("C", {10.661840_angstrom, 11.702084_angstrom, 15.849600_angstrom});
    ions.insert("C", { 8.530340_angstrom, 10.471462_angstrom, 15.849600_angstrom});
    ions.insert("C", { 7.819840_angstrom, 11.702084_angstrom, 15.849600_angstrom});
    ions.insert("C", { 7.819840_angstrom,  9.240840_angstrom, 15.849600_angstrom});
    ions.insert("C", { 6.398840_angstrom,  9.240840_angstrom, 15.849600_angstrom});
    ions.insert("C", { 8.530340_angstrom,  8.010218_angstrom, 15.849600_angstrom});
    ions.insert("C", { 7.819840_angstrom,  6.779596_angstrom, 15.849600_angstrom});
    ions.insert("C", { 9.951340_angstrom,  8.010218_angstrom, 15.849600_angstrom});
    ions.insert("C", {10.661840_angstrom,  6.779596_angstrom, 15.849600_angstrom});
    ions.insert("C", {12.793340_angstrom, 10.471462_angstrom, 15.849600_angstrom});
    ions.insert("C", {12.082840_angstrom, 11.702084_angstrom, 15.849600_angstrom});
    ions.insert("C", { 9.951340_angstrom, 12.932706_angstrom, 15.849600_angstrom});
    ions.insert("C", { 8.530340_angstrom, 12.932706_angstrom, 15.849600_angstrom});
    ions.insert("C", { 6.398840_angstrom, 11.702084_angstrom, 15.849600_angstrom});
    ions.insert("C", { 5.688340_angstrom, 10.471462_angstrom, 15.849600_angstrom});
    ions.insert("C", { 5.688340_angstrom,  8.010218_angstrom, 15.849600_angstrom});
    ions.insert("C", { 6.398840_angstrom,  6.779596_angstrom, 15.849600_angstrom});
    ions.insert("C", { 8.530340_angstrom,  5.548974_angstrom, 15.849600_angstrom});
    ions.insert("C", { 9.951340_angstrom,  5.548974_angstrom, 15.849600_angstrom});
    ions.insert("C", {12.082840_angstrom,  6.779596_angstrom, 15.849600_angstrom});
    ions.insert("C", {12.793340_angstrom,  8.010218_angstrom, 15.849600_angstrom});

    // Hydrogens (H1-H12)
    ions.insert("H", {13.819514_angstrom, 10.826939_angstrom, 15.849600_angstrom});
    ions.insert("H", {12.903779_angstrom, 12.413038_angstrom, 15.849600_angstrom});
    ions.insert("H", {10.156575_angstrom, 13.999137_angstrom, 15.849600_angstrom});
    ions.insert("H", { 8.325105_angstrom, 13.999137_angstrom, 15.849600_angstrom});
    ions.insert("H", { 5.577901_angstrom, 12.413038_angstrom, 15.849600_angstrom});
    ions.insert("H", { 4.662166_angstrom, 10.826939_angstrom, 15.849600_angstrom});
    ions.insert("H", { 4.662166_angstrom,  7.654741_angstrom, 15.849600_angstrom});
    ions.insert("H", { 5.577901_angstrom,  6.068642_angstrom, 15.849600_angstrom});
    ions.insert("H", { 8.325105_angstrom,  4.482543_angstrom, 15.849600_angstrom});
    ions.insert("H", {10.156575_angstrom,  4.482543_angstrom, 15.849600_angstrom});
    ions.insert("H", {12.903779_angstrom,  6.068642_angstrom, 15.849600_angstrom});
    ions.insert("H", {13.819514_angstrom,  7.654741_angstrom, 15.849600_angstrom});

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
        summary << "run = diagnostic_03_coronene_hardcoded\n";
        summary << "system = coronene_C24H12_hardcoded\n";
        summary << "geometry_source = hardcoded_in_cpp\n";
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
