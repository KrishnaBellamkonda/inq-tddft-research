// ============================================================================
// localised_jellium / pbc_open_z / gs / run.cpp
//
// Ground state of the slab_n52 witness system at SWITCHABLE z-periodicity —
// the task-0 GS for the PBC-vs-open-z (Arm B) energy-oscillation campaign
// (docs/campaigns/localised_jellium/pbc-open-z-oscillation.md).
//
//   EM_PERIODICITY=2 (default) : periodic x,y + OPEN z (slab-truncated Poisson;
//                                no electrostatic images along z)
//   EM_PERIODICITY=3           : fully periodic (the convention of ALL prior
//                                slab_n52 oscillation runs)
//
// Clone of muon_mass_fork/effmass_sigma1/gs/run.cpp (the producer of the p3
// witness GS shared_gs/slab_n52_L40x40x80_dx0p333) with two additions only:
// the periodicity knob and a GS-density VTI dump for the n(z) sanity gate.
// Checkpoint default: shared_gs/slab_n52_L40x40x80_dx0p333_per2.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/analytics.hpp>
#include "../../../shared/configs/slab_n52_L40x40x80.hpp"
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN52_L40x40x80;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
	const double SPACING     = env_d("EM_SPACING", 0.33333);
	const int    PERIODICITY = env_i("EM_PERIODICITY", 2);
	const std::string CHECKPOINT_DIR = env_s("EM_GS_DIR",
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n52_L40x40x80_dx0p333_per2");

	if (PERIODICITY != 2 && PERIODICITY != 3) {
		std::cerr << "FATAL: EM_PERIODICITY must be 2 or 3, got " << PERIODICITY << "\n";
		return 2;
	}

	std::cout << "\n=== pbc_open_z GS: slab_n52 40x40x80 dx=" << SPACING
	          << " periodicity=" << PERIODICITY << " ===\n"
	          << "  N = " << Cfg::N_ELECTRONS << "  r_s = " << inqkit::jellium::rs_from_n0(Cfg::N0)
	          << "  checkpoint = " << CHECKPOINT_DIR << "\n\n";

	auto cell0 = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b,
	                                          Cfg::LZ_BOHR * 1.0_b);
	auto cell = (PERIODICITY == 2) ? cell0.periodicity(2) : cell0.periodic();
	auto ions = systems::ions(cell);
	auto electrons = systems::electrons(
		ions,
		options::electrons{}
			.spacing(SPACING * 1.0_b)
			.extra_electrons(Cfg::N_ELECTRONS)
			.extra_states(Cfg::EXTRA_STATES)
			.temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
		input::kpoints::gamma());

	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation pert(bg);

	ground_state::initial_guess(ions, electrons);
	auto gs = ground_state::calculate(
		ions, electrons, options::theory{}.lda(),
		options::ground_state{}
			.energy_tolerance(Cfg::SCF_TOL_HA * 1.0_Ha)
			.max_steps(Cfg::SCF_MAX_STEPS)
			.broyden_mixing().mixing_ndim(Cfg::SCF_MIX_NDIM).mixing(Cfg::SCF_MIX_ALPHA),
		pert);

	const int n_states = electrons.states().num_states();
	std::cout << "  GS energy = " << std::setprecision(12) << gs.energy.total()
	          << " Ha   num_states = " << n_states << "\n";

	std::filesystem::create_directories(CHECKPOINT_DIR);
	electrons.save(CHECKPOINT_DIR);

	// GS density VTI for the n(z) sanity gate (loaded via inqview.load_vti)
	{
		inqkit::io::RealField3DWriter gs_wr("results/density_gs",
			{ .field_name = "density", .include_meta = false, .emit_raw = false,
			  .emit_vti = true,
			  .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
			{ .overwrite = true });
		gs_wr.write(inqkit::fields::density::total(electrons), "density_gs");
	}

	if (electrons.root()) {
		std::filesystem::create_directories("results");
		std::ofstream s("results/run_summary.txt");
		s << std::setprecision(14)
		  << "run = localised_jellium/pbc_open_z/gs\n"
		  << "engine = inq-study\nxc = LDA\n"
		  << "cell_bohr = 40x40x80 (orthorhombic)\n"
		  << "periodicity = " << PERIODICITY << "\n"
		  << "n0_a0m3 = " << Cfg::N0 << "\nr_s = " << inqkit::jellium::rs_from_n0(Cfg::N0) << "\n"
		  << "spacing_bohr = " << SPACING << "\n"
		  << "extra_electrons = " << Cfg::N_ELECTRONS << "\nextra_states = " << Cfg::EXTRA_STATES << "\n"
		  << "ground_state_energy_ha = " << gs.energy.total() << "\n"
		  << "num_states = " << n_states << "\n"
		  << "checkpoint_dir = " << CHECKPOINT_DIR << "\nrun_completed = true\n";
	}
	std::cout << "Done.\n";
	return 0;
}
