// ============================================================================
// save_gs/gs_L50_cubic_N162_dx0p30/run.cpp  (v2 of the comparison plan)
//
// Jellium GS at cell = 50^3 Bohr cubic-periodic, N = 162 (closed-shell magic
// |G|^2 <= 6 — same lineage as the L=50 / N=162 base run, the plasmon
// variants, and the positive-ion companion). Density n = 162/125000 =
// 1.296e-3 e/Bohr^3, r_s = 5.69 Bohr, Lithium-like.
//
// Spacing 0.30 Bohr (relaxed from user-spec 0.248 to fit the eigensolver
// workspaces in 24 GB GPU memory). Nyquist k = pi/0.30 = 10.47 Bohr^-1
// sits just below WP k_0 = 10.50; the WP centre is at the edge of the
// resolved k-space, and its 3 sigma_k = 0.6 high-k tail aliases by ~6 %.
// For Bethe-regime stopping at 1500 eV (k_0 dominates) this is acceptable;
// a publication-quality run would tighten back to 0.248 on a workstation
// with >24 GB GPU memory.
//
// Reused by both v2 run dirs:
//   * run_wp_e1500_L50_cubic/         (Gaussian wave-packet projectile)
//   * run_classical_e1500_L50_cubic/  (classical-electron projectile)
//
// Cost estimate: grid lands at ~167^3 = 4.66M points (half of dx=0.248).
// With 91 occupied + 20 extra = 111 spatial states this GS should fit
// comfortably on one A30 (eigensolver buffers ~ 6 x 7.5 GB = 45 GB
// nominal but only 2-3 are simultaneously resident, fitting in 24 GB).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include "../../shared/configs/electron_proj_E1500_L50_cubic.hpp"
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
using Cfg = jellium::config::Common_E1500_L50_cubic;

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

int main() {
    // Do NOT call input::environment{} here — INQ's systems::electrons()
    // initialises MPI itself; calling environment{} would cause a fatal
    // double MPI_Init (verified in v1 pre-GS dryrun).

    const std::string CHECKPOINT_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p30";

    std::cout << "\n=== save_gs/gs_L50_cubic_N162_dx0p30 ===\n"
              << "  cell = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
              << "  volume = " << (Cfg::L_BOHR * Cfg::L_BOHR * Cfg::L_BOHR)
              << " Bohr^3\n"
              << "  N_electrons = " << Cfg::N_ELECTRONS
              << " (closed shell |G|^2 <= 6)\n"
              << "  spacing = " << Cfg::SPACING_BOHR << " Bohr\n"
              << "  cutoff_Ha = " << Cfg::CUTOFF_HA << " (= " << (2*Cfg::CUTOFF_HA)
              << " Ry)\n"
              << "  k_Nyquist = " << M_PI / Cfg::SPACING_BOHR << " Bohr^-1\n"
              << "  WP_kinetic = " << Cfg::WP_EKIN_EV
              << " eV (k0 = " << Cfg::WP_K0 << " Bohr^-1)\n"
              << "  WP_sigma = " << Cfg::WP_SIGMA_BOHR << " Bohr\n"
              << "  scf_tol_ha = " << Cfg::SCF_TOL_HA << "\n"
              << "  extra_states = " << Cfg::EXTRA_STATES << "\n"
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

    if (electrons.root()) {
        std::ofstream summary("results/run_summary.txt");
        summary << std::setprecision(16);
        summary << "run = save_gs/gs_L50_cubic_N162_dx0p30\n"
                << "system = jellium_N162_L50_cubic_E1500_runs\n"
                << "checkpoint_dir = " << CHECKPOINT_DIR << "\n"
                << "cell_bohr = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
                << "boundary = periodic\n"
                << "xc = LDA\n"
                << "spacing_bohr = " << Cfg::SPACING_BOHR << "\n"
                << "cutoff_Ha = " << Cfg::CUTOFF_HA << "\n"
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
