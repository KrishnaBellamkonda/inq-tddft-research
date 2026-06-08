// ============================================================================
// save_gs/gs_L50_cubic_N138_dx1p0/run.cpp
//
// Jellium GS at cell = 50^3 Bohr cubic-periodic, N = 138 (closed-shell magic
// |G|^2 <= 6 — same shell as the L=60 and L=30 runs, intermediate box).
// Density 1.104e-3 e/bohr^3, r_s ~= 5.97 bohr (lithium-like).
// Spacing 1.0 bohr (k_Nyquist = pi/1.0 = 3.14 Bohr^-1; resolves WP up to
// k_0 + 3 sigma_k = 0.332 + 0.600 = 0.932 Bohr^-1 with ~3.4x margin).
// EXTRA_STATES = 20 (project standard).
// SCF tolerance 1e-6 Ha (project standard).
// Cfg = jellium::config::Base_N138_L50_E1p5.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include "../../shared/configs/base_n138_L50_E1p5.hpp"
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
using Cfg = jellium::config::Base_N138_L50_E1p5;

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

int main() {
    const std::string CHECKPOINT_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/checkpoints/gs_L50_cubic_N138_dx1p0";

    std::cout << "\n=== save_gs/gs_L50_cubic_N138_dx1p0 ===\n"
              << "  cell = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
              << "  N_electrons = " << Cfg::N_ELECTRONS
              << " (closed shell |G|^2 <= 6)\n"
              << "  spacing = " << Cfg::SPACING_BOHR << " bohr\n"
              << "  k_Nyquist = " << M_PI / Cfg::SPACING_BOHR << " Bohr^-1\n"
              << "  WP_kinetic = " << Cfg::WP_EKIN_EV << " eV (k0="
              << Cfg::WP_K0 << " Bohr^-1)\n"
              << "  WP_sigma = " << Cfg::WP_SIGMA_BOHR << " bohr\n"
              << "  scf_tol_ha = " << Cfg::SCF_TOL_HA << "\n"
              << "  extra_states = " << Cfg::EXTRA_STATES << "\n"
              << "  checkpoint = " << CHECKPOINT_DIR << "\n";

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    std::cout << "  Atoms: " << ions.size() << " (jellium - no nuclei)\n";

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

    jellium::eigenvalues::dump(electrons, CHECKPOINT_DIR);
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
        summary << "run = save_gs/gs_L50_cubic_N138_dx1p0\n"
                << "system = jellium_N138_L50_closed_shell_E1p5_wide_WP\n"
                << "checkpoint_dir = " << CHECKPOINT_DIR << "\n"
                << "cell_bohr = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
                << "boundary = periodic\n"
                << "xc = LDA\n"
                << "spacing_bohr = " << Cfg::SPACING_BOHR << "\n"
                << "k_nyquist_bohr_inv = " << M_PI / Cfg::SPACING_BOHR << "\n"
                << "wp_kinetic_ev = " << Cfg::WP_EKIN_EV << "\n"
                << "wp_k0_bohr_inv = " << Cfg::WP_K0 << "\n"
                << "wp_sigma_bohr = " << Cfg::WP_SIGMA_BOHR << "\n"
                << "temperature_ev = " << Cfg::TEMPERATURE_EV << "\n"
                << "extra_electrons = " << Cfg::N_ELECTRONS << "\n"
                << "extra_states = " << Cfg::EXTRA_STATES << "\n"
                << "scf_tol_ha = " << Cfg::SCF_TOL_HA << "\n"
                << "ground_state_energy_ha = " << gs.energy.total() << "\n"
                << "num_states = " << n_states << "\n"
                << "num_electrons = " << n_electrons << "\n"
                << "n_occupied = " << n_occupied << "\n";
    }

    std::cout << "Done.\n";
    return 0;
}
