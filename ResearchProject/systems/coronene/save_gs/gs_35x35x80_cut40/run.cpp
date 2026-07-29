// ============================================================================
// save_gs/gs_35x35x80_cut40/run.cpp
//
// Coronene GS at cell = 35 x 35 x 80 Bohr, cutoff 40 Ha.
// Loaded by run_b18_35x35x80 and run_b6_35x35x80.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include "../../shared/configs/cell_35x35x80.hpp"
#include "../../shared/cpp/eigenvalues_writer.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = coronene::config::cell_35x35x80;

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

int main() {
    const std::string GEOMETRY_XYZ =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/shared/geometry/coronene.xyz";
    const std::string CHECKPOINT_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/checkpoints/gs_35x35x80_cut40";

    std::cout << "\n=== save_gs/gs_35x35x80_cut40 ===\n"
              << "  cell = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR
              << " x " << Cfg::LZ_BOHR << " Bohr\n"
              << "  cutoff = " << Cfg::CUTOFF_HA << " Ha\n"
              << "  geometry   = " << GEOMETRY_XYZ << "\n"
              << "  checkpoint = " << CHECKPOINT_DIR << "\n";

    auto cell = systems::cell::orthorhombic(
        Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b
    ).finite();
    auto ions = systems::ions::parse(GEOMETRY_XYZ, cell);
    std::cout << "  Atoms: " << ions.size() << "\n";

    {
        const double hx = 0.5 * Cfg::LX_BOHR, hy = 0.5 * Cfg::LY_BOHR,
                     hz = 0.5 * Cfg::LZ_BOHR;
        for (int i = 0; i < static_cast<int>(ions.size()); ++i) {
            auto const &p = ions.positions()[i];
            if (std::fabs(p[0]) > hx || std::fabs(p[1]) > hy ||
                std::fabs(p[2]) > hz) {
                std::cerr << "FATAL: atom " << i << " out of cell.\n";
                return 2;
            }
        }
    }

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .cutoff(Cfg::CUTOFF_HA * 1.0_Ha)
            .extra_states(Cfg::EXTRA_STATES));

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(Cfg::SCF_TOL_HA * 1.0_Ha)
            .max_steps(Cfg::SCF_MAX_STEPS)
            .broyden_mixing()
            .mixing_ndim(Cfg::SCF_MIX_NDIM)
            .mixing(Cfg::SCF_MIX_ALPHA));
    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    const int n_occupied  = n_electrons / 2;
    std::cout << "  num_states = " << n_states
              << "  num_electrons = " << n_electrons
              << "  n_occupied = " << n_occupied << "\n";

    std::filesystem::create_directories(CHECKPOINT_DIR);
    electrons.save(CHECKPOINT_DIR);
    std::cout << "  Wrote checkpoint to " << CHECKPOINT_DIR << "\n";

    coronene::eigenvalues::dump(electrons, CHECKPOINT_DIR);
    coronene::eigenvalues::dump(electrons,
                                "results/raw/observables/eigenvalues");
    std::cout << "  Wrote eigenvalues + occupations CSVs\n";

    std::filesystem::create_directories("results/density_gs_system");
    {
        inqkit::io::RealField3DWriter gs_wr("results/density_gs_system",
            { .field_name = "density",
              .include_meta = false,
              .emit_raw = false,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
            { .overwrite = true });
        gs_wr.write(inqkit::fields::density::total(electrons),
                    "density_gs_system");
    }
    std::filesystem::create_directories("results/density_gs_orbitals");
    {
        inqkit::io::RealField3DWriter orb_wr("results/density_gs_orbitals",
            { .field_name = "density",
              .include_meta = false,
              .emit_raw = false,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
            { .overwrite = true });
        for (int i = 0; i < n_states; ++i) {
            orb_wr.write(inqkit::fields::density::orbital(electrons, i),
                         "orbital_" + zero_pad(i, 4));
        }
    }

    if (electrons.root()) {
        std::ofstream summary("results/run_summary.txt");
        summary << std::setprecision(16);
        summary << "run = save_gs/gs_35x35x80_cut40\n"
                << "system = coronene_C24H12\n"
                << "geometry_file = " << GEOMETRY_XYZ << "\n"
                << "checkpoint_dir = " << CHECKPOINT_DIR << "\n"
                << "cell_bohr = " << Cfg::LX_BOHR << ' ' << Cfg::LY_BOHR << ' '
                                  << Cfg::LZ_BOHR << "\n"
                << "boundary = finite\n"
                << "xc = LDA\n"
                << "cutoff_ha = " << Cfg::CUTOFF_HA << "\n"
                << "extra_states = " << Cfg::EXTRA_STATES << "\n"
                << "scf_tol_ha = " << Cfg::SCF_TOL_HA << "\n"
                << "ground_state_energy_ha = " << gs.energy.total() << "\n"
                << "num_states = " << n_states << "\n"
                << "num_electrons = " << n_electrons << "\n"
                << "n_occupied = " << n_occupied << "\n"
                << "vti_format = binary\n"
                << "raw_emitted = no\n";
    }

    std::cout << "Done.\n";
    return 0;
}
