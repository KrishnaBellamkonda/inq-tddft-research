// ============================================================================
// run_save_gs_paper_replica
//
// Computes the Tsubonoya 2014 coronene ground state and saves it to a shared
// checkpoint directory via electrons.save(). Intended to be run *once* — the
// resulting checkpoint is then loaded by every downstream RT run that uses
// the same paper-replica configuration.
//
// Outputs:
//   <GS_CHECKPOINT_DIR>/                  - electrons.save() directory
//   results/density_gs/density_t000000.vti          (binary VTI)
//   results/density_gs_orbitals/orbital_NNNN/...    (binary VTI)
//   results/run_summary.txt                          (text)
//
// No .raw / .meta.txt sidecars: emit_raw=false, emit_vti=true.
// Binary VTI is the format ParaView 6.1 (and any vtkXMLImageDataReader)
// reads natively — confirmed end-to-end against run_09 outputs.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/config/tsubonoya_2014_coronene.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include "/local/data/public/skcb2/tddft/ResearchProject/systems/coronene/configurations/tsubonoya_2014_paper_replica/paths.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace cfg = inqkit::config::tsubonoya_2014;
namespace paths = coronene::paper_replica;

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

int main() {
    std::cout << "\n=== run_save_gs_paper_replica ===\n";
    std::cout << "  cell = " << cfg::LX_BOHR << " x " << cfg::LY_BOHR
              << " x " << cfg::LZ_BOHR << " Bohr\n";
    std::cout << "  geometry = " << paths::GEOMETRY_XYZ << "\n";
    std::cout << "  checkpoint = " << paths::GS_CHECKPOINT_DIR << "\n";

    // ----- Cell + atoms ----------------------------------------------------
    auto cell = systems::cell::orthorhombic(
        cfg::LX_BOHR * 1.0_b, cfg::LY_BOHR * 1.0_b, cfg::LZ_BOHR * 1.0_b
    ).finite();
    auto ions = systems::ions::parse(paths::GEOMETRY_XYZ, cell);
    std::cout << "  Atoms: " << ions.size() << "\n";

    // Defensive cell-bounds check.
    {
        const double half_lx = 0.5 * cfg::LX_BOHR;
        const double half_ly = 0.5 * cfg::LY_BOHR;
        const double half_lz = 0.5 * cfg::LZ_BOHR;
        for (int iatom = 0; iatom < static_cast<int>(ions.size()); ++iatom) {
            auto const & p = ions.positions()[iatom];
            if (std::fabs(p[0]) > half_lx ||
                std::fabs(p[1]) > half_ly ||
                std::fabs(p[2]) > half_lz) {
                std::cerr << "FATAL: atom " << iatom << " outside [-L/2, +L/2]: "
                          << "(" << p[0] << ", " << p[1] << ", " << p[2] << ")\n";
                return 2;
            }
        }
    }

    // ----- Electrons + GS SCF ---------------------------------------------
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .cutoff(cfg::CUTOFF_HA * 1.0_Ha)
            .extra_states(cfg::EXTRA_STATES)
    );

    ground_state::initial_guess(ions, electrons);
    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(cfg::SCF_TOL_HA * 1.0_Ha)
            .max_steps(cfg::SCF_MAX_STEPS)
            .broyden_mixing()
            .mixing_ndim(cfg::SCF_MIX_NDIM)
            .mixing(cfg::SCF_MIX_ALPHA)
    );
    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    const int n_states    = electrons.states().num_states();
    const int n_electrons = electrons.states().num_electrons();
    const int n_occupied  = n_electrons / 2;
    std::cout << "  num_states = " << n_states
              << "  num_electrons = " << n_electrons
              << "  n_occupied = " << n_occupied << "\n";

    // ----- Save the electronic system ------------------------------------
    std::filesystem::create_directories(paths::GS_CHECKPOINT_DIR);
    electrons.save(paths::GS_CHECKPOINT_DIR);
    std::cout << "  Wrote checkpoint to " << paths::GS_CHECKPOINT_DIR << "\n";

    // ----- VTI-only density outputs ---------------------------------------
    std::filesystem::create_directories("results");

    // GS total density (binary VTI; readable by ParaView's bundled VTK)
    {
        inqkit::io::RealField3DWriter gs_wr("results/density_gs",
            { .field_name = "density",
              .include_meta = false,
              .emit_raw = false,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
            { .overwrite = true });
        gs_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
    }

    // Per-orbital GS densities (binary VTI)
    std::filesystem::create_directories("results/density_gs_orbitals");
    for (int i = 0; i < n_states; ++i) {
        const auto out_dir = std::string("results/density_gs_orbitals/orbital_") + zero_pad(i, 4);
        std::filesystem::create_directories(out_dir);
        inqkit::io::RealField3DWriter orb_wr(out_dir,
            { .field_name = "density",
              .include_meta = false,
              .emit_raw = false,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
            { .overwrite = true });
        orb_wr.write(inqkit::fields::density::orbital(electrons, i), 0.0, 0);
    }
    std::cout << "  Wrote GS density and " << n_states << " orbital VTIs\n";

    // ----- Run summary ----------------------------------------------------
    if (electrons.root()) {
        std::ofstream summary("results/run_summary.txt");
        summary << std::setprecision(16);
        summary << "run = run_save_gs_paper_replica\n";
        summary << "system = coronene_C24H12\n";
        summary << "geometry_file = " << paths::GEOMETRY_XYZ << "\n";
        summary << "checkpoint_dir = " << paths::GS_CHECKPOINT_DIR << "\n";
        summary << "cell_bohr = " << cfg::LX_BOHR << ' ' << cfg::LY_BOHR << ' '
                << cfg::LZ_BOHR << "\n";
        summary << "boundary = finite\n";
        summary << "xc = ALDA\n";
        summary << "cutoff_ha = " << cfg::CUTOFF_HA << "\n";
        summary << "extra_states = " << cfg::EXTRA_STATES << "\n";
        summary << "scf_tol_ha = " << cfg::SCF_TOL_HA << "\n";
        summary << "ground_state_energy_ha = " << gs.energy.total() << "\n";
        summary << "num_states = " << n_states << "\n";
        summary << "num_electrons = " << n_electrons << "\n";
        summary << "n_occupied = " << n_occupied << "\n";
        summary << "vti_format = binary\n";
        summary << "raw_emitted = no\n";
    }

    std::cout << "Done. Checkpoint + VTI density outputs written.\n";
    return 0;
}
