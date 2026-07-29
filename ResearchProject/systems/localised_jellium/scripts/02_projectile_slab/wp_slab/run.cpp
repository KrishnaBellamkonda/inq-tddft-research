// ============================================================================
// localised_jellium / 02_projectile_slab / wp_slab / run.cpp
//
// Phase 3: a Gaussian wave-packet electron projectile (σ=0.5 Bohr, E=100 eV,
// +z) fired through the localised jellium SLAB. No CAP. Loads the slab GS,
// injects the WP into the last extra state outside the slab (z=−23, 4σ from the
// −z wall), and propagates it through the slab (faces ±12.5) to near the +z
// wall (boundary-rule stop ~1σ short).
//
// Observables: total-density frames (xz gif), 20 evenly-spaced plane screens
// (normal z) for LEED + transmitted/reflected flux, dipole + current, total
// energy. Engine: stock inq.
// ============================================================================
#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/screens/leed_pattern_accumulator.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include "../../../shared/configs/slab_n234_L50.hpp"

#include <array>
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
	constexpr int    N_STEPS     = 880;     // 17.6 a.u.: run-up + traverse + exit
	constexpr int    WRITE_EVERY = 10;      // 88 density frames; screen cadence
	constexpr int    N_SCREENS   = 20;

	// ----- System + GS -----------------------------------------------------
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

	// ----- WP projectile injection ----------------------------------------
	auto wp = inqkit::WavePacket{}
	              .center(0.0, 0.0, Cfg::WP_CZ_BOHR)
	              .sigma(Cfg::WP_SIGMA_BOHR)
	              .k0(0.0, 0.0, Cfg::WP_KZ)
	              .orthogonalise_against_occupied(electrons);
	auto report = wp.inject_into_last_extra_state(electrons, 1.0);
	std::cout << "WP injected: norm_after=" << report.norm_after
	          << " orthogonalised=" << (report.orthogonalised ? "yes" : "no") << "\n";

	// ----- Localised background well --------------------------------------
	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation pert(bg);

	// ----- 20 evenly-spaced plane screens (centred box, normal z) ---------
	std::array<double, N_SCREENS> screen_z;
	for (int k = 0; k < N_SCREENS; ++k)
		screen_z[k] = -24.0 + k * (48.0 / (N_SCREENS - 1));   // −24 … +24

	using inqkit::screens::LeedPatternAccumulator;
	using inqkit::screens::PlaneScreen;
	std::array<LeedPatternAccumulator, N_SCREENS> acc;
	for (int k = 0; k < N_SCREENS; ++k)
		acc[k] = LeedPatternAccumulator(PlaneScreen{screen_z[k], "screen_" + pad4(k)});

	// ----- Output writers --------------------------------------------------
	std::filesystem::create_directories("results/density_frames");
	std::filesystem::create_directories("results/screens");
	std::filesystem::create_directories("results");
	inqkit::io::RealField3DWriter frame_wr("results/density_frames",
		{ .field_name = "density", .include_meta = false, .emit_raw = false,
		  .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
		{ .overwrite = true });
	std::ofstream elog("results/energy_dipole_vs_time.csv");
	elog << std::setprecision(12) << "time_au,total_ha,dipole_z\n";

	const double SCREEN_DT = WRITE_EVERY * DT_AU;
	int frame = 0;

	real_time::propagate(
		ions, electrons,
		[&](auto const & data) {
			if (data.root())
				elog << data.time() << "," << data.energy().total() << ","
				     << data.dipole()[2] << "\n";
			if (data.iter() % WRITE_EVERY == 0) {
				frame_wr.write(inqkit::fields::density::total(electrons),
				               "density_" + pad4(frame));
				for (int k = 0; k < N_SCREENS; ++k) acc[k].accumulate(electrons, SCREEN_DT);
				++frame;
			}
		},
		options::theory{}.lda(),
		options::real_time{}
			.num_steps(N_STEPS)
			.dt(DT_AU * 1.0_atomictime)
			.observables_current()
			.observables_dipole(),
		pert);

	for (int k = 0; k < N_SCREENS; ++k)
		acc[k].save("results/screens/screen_" + pad4(k) + ".dat");

	if (electrons.root()) {
		std::ofstream s("results/run_summary.txt");
		s << std::setprecision(12)
		  << "run = localised_jellium/02_projectile_slab/wp_slab\n"
		  << "engine = inq (stock)\n"
		  << "projectile = WP sigma " << Cfg::WP_SIGMA_BOHR << " E "
		  << Cfg::WP_EKIN_EV << " eV k0 " << Cfg::WP_K0 << " launch_z "
		  << Cfg::WP_CZ_BOHR << "\n"
		  << "wp_norm_after = " << report.norm_after << "\n"
		  << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS
		  << "  total_time_au = " << (DT_AU * N_STEPS) << "\n"
		  << "n_screens = " << N_SCREENS << "\n"
		  << "density_frames = " << frame << "\n";
	}
	std::cout << "Phase-3 WP run done: " << frame << " frames.\n";
	return 0;
}
