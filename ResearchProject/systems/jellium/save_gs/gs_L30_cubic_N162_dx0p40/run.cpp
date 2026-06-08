// ============================================================================
// save_gs/gs_L30_cubic_N162_dx0p40/run.cpp — high-density GS for Run-6.
//
// Cubic 30^3 Bohr periodic jellium at N=162 closed shell. Density
// 6.000e-3 e/Bohr^3, r_s=3.41 Bohr (between Li and Al). Used by both
// run_wp_n162_L30_E100_highdens/ and run_classical_n162_L30_E100_highdens/.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include "../../shared/configs/highdens_n162_L30_E100.hpp"
#include "../../shared/cpp/eigenvalues_writer.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::HighDens_N162_L30_E100_Classical_dx0p40;

int main() {
    const std::string CHECKPOINT_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L30_cubic_N162_dx0p40";

    std::cout << "\n=== save_gs/gs_L30_cubic_N162_dx0p40 (Run-6 GS) ===\n"
              << "  cell = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
              << "  volume = " << (Cfg::L_BOHR * Cfg::L_BOHR * Cfg::L_BOHR)
              << " Bohr^3\n"
              << "  N_electrons = " << Cfg::N_ELECTRONS << "\n"
              << "  density   = "
              << (Cfg::N_ELECTRONS /
                  (Cfg::L_BOHR * Cfg::L_BOHR * Cfg::L_BOHR))
              << " e/Bohr^3\n"
              << "  r_s       = "
              << std::cbrt(3.0 / (4.0 * M_PI * Cfg::N_ELECTRONS /
                                  (Cfg::L_BOHR * Cfg::L_BOHR * Cfg::L_BOHR)))
              << " Bohr\n"
              << "  spacing = " << Cfg::SPACING_BOHR << " Bohr\n"
              << "  cutoff_Ha = " << Cfg::CUTOFF_HA << "\n"
              << "  checkpoint = " << CHECKPOINT_DIR << "\n\n";

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
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

    std::filesystem::create_directories(CHECKPOINT_DIR);
    electrons.save(CHECKPOINT_DIR);
    jellium::eigenvalues::dump(electrons, CHECKPOINT_DIR);
    std::filesystem::create_directories("results/raw/observables/eigenvalues");
    jellium::eigenvalues::dump(electrons,
                               "results/raw/observables/eigenvalues");

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

    if (electrons.root()) {
        std::ofstream s("results/run_summary.txt");
        s << std::setprecision(16);
        s << "run                = save_gs/gs_L30_cubic_N162_dx0p40\n"
          << "checkpoint_dir     = " << CHECKPOINT_DIR << "\n"
          << "cell_bohr          = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
          << "spacing_bohr       = " << Cfg::SPACING_BOHR << "\n"
          << "xc_functional      = LDA\n"
          << "n_electrons        = " << Cfg::N_ELECTRONS << "\n"
          << "ground_state_energy_ha = " << gs.energy.total() << "\n";
    }

    std::cout << "Done.\n";
    return 0;
}
