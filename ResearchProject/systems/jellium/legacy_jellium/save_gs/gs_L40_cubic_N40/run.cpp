// ============================================================================
// save_gs/gs_L40_cubic_N40/run.cpp
//
// Jellium GS at cell = 40^3 Bohr cubic-periodic, N = 40 open-shell.
// Loaded by run_E200_s0p53_N40 only.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include "../../shared/configs/E200_s0p53_N40.hpp"
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
using Cfg = jellium::config::E200_s0p53_N40;

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

int main() {
    const std::string CHECKPOINT_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L40_cubic_N40";

    std::cout << "\n=== save_gs/gs_L40_cubic_N40 ===\n"
              << "  cell = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
              << "  N_electrons = " << Cfg::N_ELECTRONS << " (open-shell)\n"
              << "  spacing = " << Cfg::SPACING_BOHR << " bohr\n"
              << "  temperature = " << Cfg::TEMPERATURE_EV << " eV\n"
              << "  checkpoint = " << CHECKPOINT_DIR << "\n";

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    std::cout << "  Atoms: " << ions.size() << " (jellium — no nuclei)\n";

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());

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

    jellium::eigenvalues::dump(electrons, CHECKPOINT_DIR);
    jellium::eigenvalues::dump(electrons,
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
        summary << "run = save_gs/gs_L40_cubic_N40\n"
                << "system = jellium_N40_open_shell\n"
                << "geometry_file = (none, jellium)\n"
                << "checkpoint_dir = " << CHECKPOINT_DIR << "\n"
                << "cell_bohr = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
                << "boundary = periodic\n"
                << "xc = LDA\n"
                << "spacing_bohr = " << Cfg::SPACING_BOHR << "\n"
                << "temperature_ev = " << Cfg::TEMPERATURE_EV << "\n"
                << "extra_electrons = " << Cfg::N_ELECTRONS << "\n"
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
