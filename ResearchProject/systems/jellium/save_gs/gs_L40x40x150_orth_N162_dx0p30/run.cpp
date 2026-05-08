// ============================================================================
// save_gs/gs_L40x40x150_orth_N162_dx0p30/run.cpp
//
// Jellium GS for the elongated 40 x 40 x 150 Bohr orthorhombic-periodic box,
// N=162 electrons, spacing 0.30 Bohr (Nyquist-safe at k_0 up to ~10.5 Bohr^-1
// — see electron_proj_E1000_L40x40x150.hpp for the derivation). EXTRA_STATES
// = 20 (project standard). LDA xc, T = 100 K Fermi-Dirac smearing.
//
// Density n = 162 / 240000 = 6.75e-4 e/Bohr^3, r_s = 7.07 Bohr — much more
// dilute than the L=50 N=162 base run (r_s = 5.69), but holds the same
// closed-shell magic number for ease of comparison. (See
// docs/plans/the-objective-in-this-dapper-moon.md "Open issues" §1.)
//
// Reused by both run dirs:
//   * run_wp_e1000_L40x40x150/         (Gaussian wave-packet projectile)
//   * run_classical_e1000_L40x40x150/  (classical-electron projectile)
//
// Cost: grid is 134 x 134 x 500 ≈ 9.0M points. With 91 occupied + 20 extra
// = 111 spatial states this GS is *substantial* — expect order-of-hours
// wall on a single A30 GPU. SCF tolerance 1e-6 Ha (project standard).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include "../../shared/configs/electron_proj_E1000_L40x40x150.hpp"
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

// Use the Common base struct — neither the WP nor the classical Cfg
// specialises any field this GS depends on; either would work.
using Cfg = jellium::config::Common_E1000_L40x40x150;

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

int main() {
    // Do NOT instantiate input::environment{} here — INQ's
    // systems::electrons() constructor initialises MPI itself, and a manual
    // input::environment{} would cause OpenMPI to abort with double-init
    // (verified in scripts/classical_electron_smoke/C_pre_gs_dryrun).
    const std::string CHECKPOINT_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L40x40x150_orth_N162_dx0p30";

    std::cout << "\n=== save_gs/gs_L40x40x150_orth_N162_dx0p30 ===\n"
              << "  cell = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR
              << " x " << Cfg::LZ_BOHR << " Bohr (orthorhombic, periodic)\n"
              << "  volume = " << (Cfg::LX_BOHR * Cfg::LY_BOHR * Cfg::LZ_BOHR)
              << " Bohr^3\n"
              << "  N_electrons = " << Cfg::N_ELECTRONS
              << " (closed shell from L=50 cubic — re-derive in orthorhombic)\n"
              << "  spacing = " << Cfg::SPACING_BOHR << " Bohr\n"
              << "  k_Nyquist = " << M_PI / Cfg::SPACING_BOHR << " Bohr^-1\n"
              << "  WP_kinetic = " << Cfg::WP_EKIN_EV
              << " eV (k0 = " << Cfg::WP_K0 << " Bohr^-1)\n"
              << "  WP_sigma = " << Cfg::WP_SIGMA_BOHR << " Bohr\n"
              << "  scf_tol_ha = " << Cfg::SCF_TOL_HA << "\n"
              << "  extra_states = " << Cfg::EXTRA_STATES << "\n"
              << "  checkpoint = " << CHECKPOINT_DIR << "\n\n";

    auto cell = systems::cell::orthorhombic(
        Cfg::LX_BOHR * 1.0_b,
        Cfg::LY_BOHR * 1.0_b,
        Cfg::LZ_BOHR * 1.0_b).periodic();

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

    // Single t=0 density VTI (the run_*/ scripts will re-emit at every step).
    std::filesystem::create_directories("results/density_gs_system");
    {
        inqkit::io::RealField3DWriter gs_wr("results/density_gs_system",
            { .field_name = "density",
              .include_meta = false,
              .emit_raw     = false,
              .emit_vti     = true,
              .vti_format   = inqkit::io::VTIWriteOptions::Format::binary },
            { .overwrite = true });
        gs_wr.write(inqkit::fields::density::total(electrons),
                    "density_gs_system");
    }

    if (electrons.root()) {
        std::ofstream summary("results/run_summary.txt");
        summary << std::setprecision(16);
        summary << "run = save_gs/gs_L40x40x150_orth_N162_dx0p30\n"
                << "system = jellium_N162_L40x40x150_orth_E1000_runs\n"
                << "checkpoint_dir = " << CHECKPOINT_DIR << "\n"
                << "cell_bohr_x = " << Cfg::LX_BOHR << "\n"
                << "cell_bohr_y = " << Cfg::LY_BOHR << "\n"
                << "cell_bohr_z = " << Cfg::LZ_BOHR << "\n"
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
