// ============================================================================
// run_02_wide_wp: coronene + wide Gaussian wavepacket
//
// System: same cell and geometry as run_01.
//
// Wavepacket:
//   sigma = 2.0 Å = 3.779 bohr    (4× wider than paper value)
//   E_kin = 200 eV → k0 = 3.834 bohr⁻¹
//   center: cell centre in x,y; D = 6.35 Å above flake in z
//   direction: -z
//
// Purpose: compare WP orbital density shape for broad vs narrow packet.
// ============================================================================

#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;

static constexpr double ANG_TO_BOHR = 1.8897259886;
static constexpr double HA_TO_EV    = 27.21138625;

static constexpr double LX_BOHR = 34.771;
static constexpr double LY_BOHR = 34.771;
static constexpr double LZ_BOHR = 89.856;

static constexpr double WP_SIGMA_ANG  = 2.0;
static constexpr double WP_SIGMA_BOHR = WP_SIGMA_ANG * ANG_TO_BOHR;   // 3.779 bohr
static constexpr double WP_EKIN_EV    = 200.0;
static constexpr double WP_EKIN_HA    = WP_EKIN_EV / HA_TO_EV;
static const     double WP_K0_BOHR    = std::sqrt(2.0 * WP_EKIN_HA);  // 3.834 bohr⁻¹

static constexpr double D_IMPACT_ANG  = 6.35;
static constexpr double D_IMPACT_BOHR = D_IMPACT_ANG * ANG_TO_BOHR;
static const     double WP_CX_BOHR    = LX_BOHR / 2.0;
static const     double WP_CY_BOHR    = LY_BOHR / 2.0;
static const     double WP_CZ_BOHR    = LZ_BOHR / 2.0 + D_IMPACT_BOHR;

int main() {
    auto cell = systems::cell::orthorhombic(
        LX_BOHR*1.0_b, LY_BOHR*1.0_b, LZ_BOHR*1.0_b).finite();

    auto ions = systems::ions::parse("coronene_centered.xyz", cell);

    std::cout << "\n=== run_02_wide_wp ===\n";
    std::cout << "  Atoms: " << ions.size() << "  (expect 36)\n";

    auto electrons = systems::electrons(
        ions, options::electrons{}.cutoff(40.0_Ha).extra_states(3));

    ground_state::initial_guess(ions, electrons);

    auto gs = ground_state::calculate(
        ions, electrons,
        options::theory{}.lda(),
        options::ground_state{}
            .energy_tolerance(1e-4_Ha)
            .max_steps(300)
            .broyden_mixing()
            .mixing_ndim(8)
            .mixing(0.1));

    std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

    std::cout << "Writing GS total density...\n";
    auto rho_gs = inqkit::fields::density::total(electrons);
    inqkit::io::RealField3DWriter(
        "results/density",
        {.field_name = "total_density", .include_meta = true},
        {.overwrite = true})
        .write(rho_gs, "density_total");

    std::cout << "Injecting wavepacket...\n";
    auto wp = inqkit::WavePacket{}
        .center(WP_CX_BOHR, WP_CY_BOHR, WP_CZ_BOHR)
        .sigma(WP_SIGMA_BOHR)
        .k0(0.0, 0.0, -WP_K0_BOHR)
        .orthogonalise_against_occupied(electrons);

    auto report = wp.inject_into_last_extra_state(electrons, 1.0);

    std::cout << "  state_index   = " << report.state_index   << "\n";
    std::cout << "  norm_before   = " << report.norm_before   << "\n";
    std::cout << "  max_overlap   = " << report.max_overlap   << "\n";
    std::cout << "  norm_after    = " << report.norm_after    << "  (expect ≈ 1.0)\n";
    std::cout << "  orthogonalised = " << (report.orthogonalised ? "yes" : "no") << "\n";

    std::cout << "Writing WP orbital density...\n";
    auto rho_wp = inqkit::fields::density::orbital(electrons, report.state_index);
    inqkit::io::RealField3DWriter(
        "results/wp_density",
        {.field_name = "wp_density", .include_meta = true},
        {.overwrite = true})
        .write(rho_wp, "wp_density");

    std::filesystem::create_directories("results");
    std::ofstream pf("results/wp_params.txt");
    pf << std::fixed << std::setprecision(6);
    pf << "run              = run_02_wide_wp\n";
    pf << "description      = Wide wavepacket: sigma = 2.0 Ang\n";
    pf << "sigma_ang        = " << WP_SIGMA_ANG   << "\n";
    pf << "sigma_bohr       = " << WP_SIGMA_BOHR  << "\n";
    pf << "ekin_ev          = " << WP_EKIN_EV     << "\n";
    pf << "ekin_ha          = " << WP_EKIN_HA     << "\n";
    pf << "k0_bohr_inv      = " << WP_K0_BOHR     << "\n";
    pf << "kx_bohr_inv      = " << 0.0            << "\n";
    pf << "ky_bohr_inv      = " << 0.0            << "\n";
    pf << "kz_bohr_inv      = " << -WP_K0_BOHR   << "\n";
    pf << "center_x_bohr    = " << WP_CX_BOHR     << "\n";
    pf << "center_y_bohr    = " << WP_CY_BOHR     << "\n";
    pf << "center_z_bohr    = " << WP_CZ_BOHR     << "\n";
    pf << "d_impact_ang     = " << D_IMPACT_ANG   << "\n";
    pf << "d_impact_bohr    = " << D_IMPACT_BOHR  << "\n";
    pf << "occupation       = " << 1.0            << "\n";
    pf << "state_index      = " << report.state_index  << "\n";
    pf << "norm_before      = " << report.norm_before  << "\n";
    pf << "norm_after       = " << report.norm_after   << "\n";
    pf << "max_overlap      = " << report.max_overlap  << "\n";
    pf << "orthogonalised   = " << (report.orthogonalised ? 1 : 0) << "\n";
    pf << "passed_tolerance = " << (report.passed_tolerance ? 1 : 0) << "\n";

    std::cout << "\nDone. Output in results/\n";
    return 0;
}
