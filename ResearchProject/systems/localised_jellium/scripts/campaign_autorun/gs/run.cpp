// ============================================================================
// localised_jellium / scripts/campaign_autorun / gs / run.cpp
//
// Fully ENV-PARAMETRISED localised-jellium slab GS, for the autonomous H1/H2/H3
// ladder (edge-w, Lz, thickness sweeps) and the periodicity-2 (open-z) GS used by
// H2/H4/H5. Background injected as a static perturbation during SCF. Saves the
// checkpoint to LJ_GS_DIR. Prints GS energy + writes a run_summary.
//
// Env (all optional, defaults = baseline slab_n82):
//   LJ_LX(50) LJ_LY(50) LJ_LZ(90) LJ_HALF(12.5) LJ_N(82) LJ_EDGE_W(0)
//   LJ_PERIODICITY(3|2) LJ_SPACING(0.5) LJ_GS_DIR(REQUIRED) LJ_TAG(gs)
//   n0 is computed = LJ_N / (LX*LY*2*HALF).
// Build against INQ_SOURCE=inq-study; runtime shares from inq/install/share.
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

int main() {
	const double LX = env_d("LJ_LX", 50.0), LY = env_d("LJ_LY", 50.0), LZ = env_d("LJ_LZ", 90.0);
	const double HALF = env_d("LJ_HALF", 12.5);
	const int    N = env_i("LJ_N", 82);
	const double EDGE_W = env_d("LJ_EDGE_W", 0.0);
	const int    PERIODICITY = env_i("LJ_PERIODICITY", 3);
	const double SPACING = env_d("LJ_SPACING", 0.5);
	const std::string GS_DIR = env_s("LJ_GS_DIR", "");
	const std::string TAG = env_s("LJ_TAG", "gs");
	if (GS_DIR.empty()) { std::cerr << "FATAL: LJ_GS_DIR required\n"; return 2; }

	const double N0 = double(N) / (LX * LY * (2.0 * HALF));

	std::cout << "\n=== campaign_autorun GS [" << TAG << "] ===\n"
	          << "  cell = " << LX << " x " << LY << " x " << LZ << "  periodicity=" << PERIODICITY << "\n"
	          << "  slab half-width = " << HALF << "  N = " << N << "  n0 = " << N0
	          << "  (r_s=" << inqkit::jellium::rs_from_n0(N0) << ")  edge_w=" << EDGE_W
	          << "  spacing=" << SPACING << "\n  checkpoint = " << GS_DIR << "\n\n";

	auto cell0 = systems::cell::orthorhombic(LX * 1.0_b, LY * 1.0_b, LZ * 1.0_b);
	auto cell = (PERIODICITY == 2) ? cell0.periodicity(2) : cell0.periodic();
	auto ions = systems::ions(cell);

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
	const int n_states = electrons.states().num_states();
	const int n_el = electrons.states().num_electrons();
	std::cout << "  GS energy = " << std::setprecision(12) << e_tot << " Ha   N=" << n_el << "\n";

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
		  << "run = campaign_autorun/gs/" << TAG << "\nengine = inq-study\n"
		  << "cell_bohr = " << LX << " x " << LY << " x " << LZ << "\n"
		  << "periodicity = " << PERIODICITY << "\nxc = LDA\n"
		  << "slab_half_width = " << HALF << "\nedge_width = " << EDGE_W << "\n"
		  << "n_electrons_requested = " << N << "\nn0_a0m3 = " << N0 << "\n"
		  << "r_s = " << inqkit::jellium::rs_from_n0(N0) << "\nspacing_bohr = " << SPACING << "\n"
		  << "ground_state_energy_ha = " << e_tot << "\n"
		  << "num_states = " << n_states << "\nnum_electrons = " << n_el << "\n"
		  << "checkpoint_dir = " << GS_DIR << "\nrun_completed = true\n";
	}
	std::cout << "Done.\n";
	return 0;
}
