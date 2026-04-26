// ============================================================================
// run_09_gs_vti_writer
//
// Same GS-only configuration as run_08_gs_only_wp_check, but every
// RealField3DWriter is opted in to native VTI emission via the new
// emit_vti / vti_format options on RealField3DLayout.
//
// What we expect to see compared to run_08:
//   - Each output directory contains <basename>.raw + <basename>.meta.txt
//     (back-compat) AND a new <basename>.vti written directly by the C++
//     side (no Python conversion needed).
//   - The .vti files are byte-equivalent (numerically) to the
//     inqview-converted .vti files in run_08_gs_only_wp_check_results/, up
//     to the writer's 17-digit ASCII precision / exact binary copy.
//
// Storage policy:
//   - Density of interest (GS total, WP pre/post normalisation): ASCII VTI
//     for diff-friendly verification.
//   - Per-orbital GS densities (62 fields): binary VTI to keep disk usage
//     manageable (~31 MB each instead of ~70 MB).
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/config/tsubonoya_2014_coronene.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

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

static std::string zero_pad(int n, int width) {
    std::ostringstream ss;
    ss << std::setfill('0') << std::setw(width) << n;
    return ss.str();
}

int main() {
    std::cout << "\n=== run_09: GS-only with native VTI writer ===\n";
    std::cout << "  cell = " << cfg::LX_BOHR << " x " << cfg::LY_BOHR
              << " x " << cfg::LZ_BOHR << " Bohr\n";
    std::cout << "  WP sigma = " << cfg::WP_SIGMA_BOHR
              << " Bohr, |k| = " << cfg::WP_K0 << " Bohr^-1\n";

    // ----- Cell + atoms ----------------------------------------------------
    auto cell = systems::cell::orthorhombic(
        cfg::LX_BOHR * 1.0_b, cfg::LY_BOHR * 1.0_b, cfg::LZ_BOHR * 1.0_b
    ).finite();
    auto ions = systems::ions::parse("coronene_centred.xyz", cell);
    std::cout << "  Atoms: " << ions.size() << "\n";

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

    std::filesystem::create_directories("results");

    // ----- GS total density (raw + meta + ASCII VTI) ----------------------
    {
        inqkit::io::RealField3DWriter gs_wr("results/density_gs",
            { .field_name = "density",
              .include_meta = true,
              .emit_raw = true,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::ascii },
            { .overwrite = true });
        gs_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
    }

    // ----- Per-orbital GS densities (raw + meta + binary VTI) -------------
    std::filesystem::create_directories("results/density_gs_orbitals");
    for (int i = 0; i < n_states; ++i) {
        const auto out_dir = std::string("results/density_gs_orbitals/orbital_") + zero_pad(i, 4);
        std::filesystem::create_directories(out_dir);
        inqkit::io::RealField3DWriter orb_wr(out_dir,
            { .field_name = "density",
              .include_meta = true,
              .emit_raw = true,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
            { .overwrite = true });
        orb_wr.write(inqkit::fields::density::orbital(electrons, i), 0.0, 0);
    }
    std::cout << "  Wrote GS density and " << n_states << " orbital densities\n";

    // ----- WP pre-normalisation (raw analytical Gaussian) -----------------
    auto wp_raw = inqkit::WavePacket{}
        .center(cfg::WP_CX_BOHR, cfg::WP_CY_BOHR, cfg::WP_CZ_BOHR)
        .sigma(cfg::WP_SIGMA_BOHR)
        .k0(cfg::WP_KX, cfg::WP_KY, cfg::WP_KZ);
    auto report_raw = wp_raw.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report_raw.state_index;
    std::cout << "  WP pass 1 (pre-normalisation, no ortho):"
              << " state_index = " << wp_idx
              << "  grid_norm = " << report_raw.norm_after << "\n";
    {
        inqkit::io::RealField3DWriter w("results/density_wp_pre_normalisation",
            { .field_name = "density",
              .include_meta = true,
              .emit_raw = true,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::ascii },
            { .overwrite = true });
        w.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    }

    // ----- WP post-normalisation (MGS + renorm) ---------------------------
    auto wp_norm = inqkit::WavePacket{}
        .center(cfg::WP_CX_BOHR, cfg::WP_CY_BOHR, cfg::WP_CZ_BOHR)
        .sigma(cfg::WP_SIGMA_BOHR)
        .k0(cfg::WP_KX, cfg::WP_KY, cfg::WP_KZ)
        .orthogonalise_against_occupied(electrons);
    auto report = wp_norm.inject_into_last_extra_state(electrons, 1.0);
    std::cout << "  WP pass 2 (post-normalisation, ortho+renorm):"
              << " norm_before = " << report.norm_before
              << "  norm_after = "  << report.norm_after
              << "  max_overlap = " << report.max_overlap << "\n";
    {
        inqkit::io::RealField3DWriter w("results/density_wp_post_normalisation",
            { .field_name = "density",
              .include_meta = true,
              .emit_raw = true,
              .emit_vti = true,
              .vti_format = inqkit::io::VTIWriteOptions::Format::ascii },
            { .overwrite = true });
        w.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    }

    // ----- Run summary ----------------------------------------------------
    if (electrons.root()) {
        std::ofstream summary("results/run_summary.txt");
        summary << std::setprecision(16);
        summary << "run = run_09_gs_vti_writer\n";
        summary << "system = coronene_C24H12\n";
        summary << "geometry_file = coronene_centred.xyz\n";
        summary << "cell_bohr = " << cfg::LX_BOHR << ' ' << cfg::LY_BOHR << ' ' << cfg::LZ_BOHR << "\n";
        summary << "boundary = finite\n";
        summary << "xc = ALDA\n";
        summary << "cutoff_ha = " << cfg::CUTOFF_HA << "\n";
        summary << "extra_states = " << cfg::EXTRA_STATES << "\n";
        summary << "scf_tol_ha = " << cfg::SCF_TOL_HA << "\n";
        summary << "ground_state_energy_ha = " << gs.energy.total() << "\n";
        summary << "num_states = " << n_states << "\n";
        summary << "num_electrons = " << n_electrons << "\n";
        summary << "n_occupied = " << n_occupied << "\n";
        summary << "wp_state_index = " << wp_idx << "\n";
        summary << "wp_pre_norm_grid_norm = " << report_raw.norm_after << "\n";
        summary << "wp_post_norm_norm_before_renorm = " << report.norm_before << "\n";
        summary << "wp_post_norm_norm_after_renorm = " << report.norm_after << "\n";
        summary << "wp_post_norm_max_overlap = " << report.max_overlap << "\n";
        summary << "wp_sigma_bohr = " << cfg::WP_SIGMA_BOHR << "\n";
        summary << "wp_k0 = " << cfg::WP_K0 << "\n";
        summary << "wp_offset_bohr = " << cfg::WP_OFFSET_BOHR << "\n";
        summary << "wp_ekin_ev = " << cfg::WP_EKIN_EV << "\n";
        summary << "rt_propagation = DISABLED\n";
        summary << "vti_native_writer = ENABLED\n";
    }

    std::cout << "Done. Output in results/\n";
    return 0;
}
