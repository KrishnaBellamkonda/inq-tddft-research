// ============================================================================
// localised_jellium / scripts/campaign_autorun / gs_ghost / run.cpp
//
// B2 (energy book-keeping campaign): SCF WITH the projectile present.
//  - LJ_GHOST=1: stationary classical Gaussian ghost (electron_gaussian_wpsigma0p5,
//    z_valence 0) at z=LJ_LAUNCH_Z, SCF converges the slab AROUND it -> the
//    self-consistently screened classical projectile.
//  - LJ_GHOST=0 + LJ_N=83: the "83-electron SCF" illustration (extra electron
//    relaxes unconstrained; background still from LJ_N_BG=82).
// Derived from the validated gs/run.cpp (same SCF settings) + the ion-insertion
// lines of classical/run.cpp. Saves density VTI + full energy decomposition.
//
// Env: as gs/run.cpp, plus LJ_GHOST(0) LJ_LAUNCH_Z(-16.5) LJ_PROJ_UPF
//      LJ_N_BG(defaults to LJ_N) — background n0 = LJ_N_BG/(LX*LY*2*HALF).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/analytics.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }
static const char* PROJ_PSEUDO =
  "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main() {
	constexpr double M_PROJ = 1.0/1822.8885;
	const double LX = env_d("LJ_LX", 50.0), LY = env_d("LJ_LY", 50.0), LZ = env_d("LJ_LZ", 120.0);
	const double HALF = env_d("LJ_HALF", 12.5);
	const int    N = env_i("LJ_N", 82);
	const int    N_BG = env_i("LJ_N_BG", N);
	const double EDGE_W = env_d("LJ_EDGE_W", 0.0);
	const int    PERIODICITY = env_i("LJ_PERIODICITY", 2);
	const double SPACING = env_d("LJ_SPACING", 0.5);
	const int    GHOST = env_i("LJ_GHOST", 0);
	const double LAUNCH_Z = env_d("LJ_LAUNCH_Z", -16.5);
	const std::string GS_DIR = env_s("LJ_GS_DIR", "");
	const std::string TAG = env_s("LJ_TAG", "gs_ghost");
	if (GS_DIR.empty()) { std::cerr << "FATAL: LJ_GS_DIR required\n"; return 2; }

	const double N0 = double(N_BG) / (LX * LY * (2.0 * HALF));

	std::cout << "\n=== B2 gs_ghost [" << TAG << "] ghost=" << GHOST << " z=" << LAUNCH_Z
	          << " N=" << N << " N_bg=" << N_BG << " per=" << PERIODICITY << " ===\n";

	auto cell0 = systems::cell::orthorhombic(LX * 1.0_b, LY * 1.0_b, LZ * 1.0_b);
	auto cell = (PERIODICITY == 2) ? cell0.periodicity(2) : cell0.periodic();
	auto ions = systems::ions(cell);
	if (GHOST) {
		auto sp = ionic::species("H").pseudo_file(env_s("LJ_PROJ_UPF", PROJ_PSEUDO)).mass(M_PROJ);
		ions.insert(sp, {0.0 * 1.0_b, 0.0 * 1.0_b, LAUNCH_Z * 1.0_b});
	}

	auto electrons = systems::electrons(
		ions,
		options::electrons{}
			.spacing(SPACING * 1.0_b)
			.extra_electrons(N)
			.extra_states(20)
			.temperature(0.00862 * 1.0_eV),
		input::kpoints::gamma());

	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = N0; bg.half_width = HALF; bg.slab_axis = 2;
	bg.center = {0.0, 0.0, 0.0}; bg.edge_width = EDGE_W;
	inqkit::jellium::localised_background_perturbation pert(bg);

	ground_state::initial_guess(ions, electrons);
	auto gs = ground_state::calculate(
		ions, electrons,
		options::theory{}.lda(),
		options::ground_state{}
			.energy_tolerance(1.0e-4 * 1.0_Ha)
			.max_steps(300)
			.broyden_mixing().mixing_ndim(8).mixing(0.1),
		pert);

	const double e_tot = gs.energy.total();
	std::cout << "  GS energy = " << std::setprecision(12) << e_tot << " Ha\n";

	std::filesystem::create_directories(GS_DIR);
	electrons.save(GS_DIR);

	std::filesystem::create_directories("results/density_gs_system");
	{ inqkit::io::RealField3DWriter wr("results/density_gs_system",
		{.field_name="density", .include_meta=false, .emit_raw=false, .emit_vti=true,
		 .vti_format=inqkit::io::VTIWriteOptions::Format::binary}, {.overwrite=true});
	  wr.write(inqkit::fields::density::total(electrons), "density_gs_system"); }

	if (electrons.root()) {
		std::filesystem::create_directories("results");
		std::ofstream s("results/run_summary.txt");
		s << std::setprecision(16)
		  << "run = campaign_autorun/gs_ghost/" << TAG << "\nengine = inq-study\n"
		  << "cell_bohr = " << LX << " x " << LY << " x " << LZ << "\n"
		  << "periodicity = " << PERIODICITY << "\nxc = LDA\n"
		  << "slab_half_width = " << HALF << "\nedge_width = " << EDGE_W << "\n"
		  << "n_electrons_requested = " << N << "\nn_bg_electrons = " << N_BG << "\nn0_a0m3 = " << N0 << "\n"
		  << "ghost = " << GHOST << "\nghost_z = " << LAUNCH_Z << "\n"
		  << "proj_upf = " << (GHOST ? env_s("LJ_PROJ_UPF", PROJ_PSEUDO) : std::string("none")) << "\n"
		  << "spacing_bohr = " << SPACING << "\n"
		  << "ground_state_energy_ha = " << e_tot << "\n"
		  << "energy_kinetic_ha = " << gs.energy.kinetic() << "\n"
		  << "energy_hartree_ha = " << gs.energy.hartree() << "\n"
		  << "energy_xc_ha = " << gs.energy.xc() << "\n"
		  << "energy_external_ha = " << gs.energy.external() << "\n"
		  << "num_electrons = " << electrons.states().num_electrons() << "\n"
		  << "checkpoint_dir = " << GS_DIR << "\nrun_completed = true\n";
	}
	std::cout << "Done.\n";
	return 0;
}
