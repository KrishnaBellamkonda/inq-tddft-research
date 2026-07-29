// ============================================================================
// localised_jellium / 01_slab_validation / static_2au / run.cpp
//
// Phase-2 STATIC run + T3.4 gate: load the slab GS and propagate 2 a.u. with the
// background well ON but NO projectile and NO CAP. Two purposes:
//   (1) T3.4 — the background perturbation is static/Hermitian, so total energy
//       must be conserved and the density must stay stationary.
//   (2) the user-requested "density as a function of time" (xz gif) baseline.
//
// dt=0.02, 100 steps = 2.0 a.u.; density frame every 2 steps (50 frames).
// Energy logged every step. Engine: stock inq (no CAP).
// ============================================================================
#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include "../../../shared/configs/slab_n234_L50.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN234L50;

static std::string pad4(int n) {
	std::ostringstream s; s << std::setfill('0') << std::setw(4) << n; return s.str();
}

int main() {
	const std::string GS_DIR =
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n234_L50";

	constexpr double DT_AU       = 0.02;
	constexpr int    N_STEPS     = 100;     // 2.0 a.u.
	constexpr int    WRITE_EVERY = 2;       // 50 frames

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
	electrons.load(GS_DIR);
	std::cout << "Loaded GS from " << GS_DIR << "\n";

	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation pert(bg);

	std::filesystem::create_directories("results/density_frames");
	std::filesystem::create_directories("results");
	inqkit::io::RealField3DWriter frame_wr("results/density_frames",
		{ .field_name = "density", .include_meta = false, .emit_raw = false,
		  .emit_vti = true,
		  .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
		{ .overwrite = true });

	std::ofstream elog("results/energy_vs_time.csv");
	elog << std::setprecision(12) << "time_au,total_ha\n";
	int frame = 0;

	real_time::propagate(
		ions, electrons,
		[&](auto const & data) {
			if (data.root())
				elog << data.time() << "," << data.energy().total() << "\n";
			if (data.iter() % WRITE_EVERY == 0) {
				frame_wr.write(inqkit::fields::density::total(electrons),
				               "density_" + pad4(frame));
				++frame;
			}
		},
		options::theory{}.lda(),
		options::real_time{}
			.num_steps(N_STEPS)
			.dt(DT_AU * 1.0_atomictime)
			.observables_dipole(),
		pert);

	std::cout << "Static run done: " << frame << " density frames.\n";
	return 0;
}
