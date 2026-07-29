// ============================================================================
// localised_jellium / 03_cap_stopping / wp_cap / run.cpp
//
// Phase 5 (WP): WP projectile through the slab WITH a two-sided sin² CAP, on
// inq-study (the CAP's imaginary potential needs the complexified scalar
// potential). Perturbation = sum(background, sum(CAP_-z, CAP_+z)).
//
//   CAP: eta=-0.5 Ha, 7.5 Bohr each side. −z window z∈[−25,−17.5] (mid −0.425,
//        width 0.15 fractional); +z window z∈[17.5,25] (mid +0.425).
//   WP : σ=0.5, E=100 eV, +z, launched at z=−15.5 (4σ from the −z CAP inner
//        edge at −17.5; 3 Bohr to the slab face at −12.5).
//
// Measured: total norm ∫n(t) = num_electrons (cumulative absorbed = N0−∫n), and
// total energy(t). After the WP is fully absorbed and CAP drain dies, the bath
// energy increase gives the deposited (stopping) energy. Engine: inq-study.
// ============================================================================
#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
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
	constexpr int    N_STEPS     = 900;     // 18 a.u.: traverse + absorb + settle
	constexpr int    WRITE_EVERY = 10;
	constexpr double WP_LAUNCH_Z = -15.5;   // 4σ from −z CAP edge (−17.5)
	constexpr double CAP_ETA_HA  = -0.5;
	constexpr double CAP_MID     = 0.425;   // |mid| fractional (z=±21.25 Bohr)
	constexpr double CAP_WIDTH   = 0.15;    // 7.5 Bohr / 50

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
	std::cout << "Loaded GS (on inq-study) from " << GS_DIR << "\n";

	auto wp = inqkit::WavePacket{}
	              .center(0.0, 0.0, WP_LAUNCH_Z)
	              .sigma(Cfg::WP_SIGMA_BOHR)
	              .k0(0.0, 0.0, Cfg::WP_KZ)
	              .orthogonalise_against_occupied(electrons);
	auto report = wp.inject_into_last_extra_state(electrons, 1.0);
	std::cout << "WP injected: norm_after=" << report.norm_after << "\n";

	// Background well + two-sided sin² CAP.
	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation bg_pert(bg);

	perturbations::absorbing cap_lo(CAP_ETA_HA * 1.0_Ha, -CAP_MID, CAP_WIDTH);
	perturbations::absorbing cap_hi(CAP_ETA_HA * 1.0_Ha,  CAP_MID, CAP_WIDTH);
	auto pert = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

	std::filesystem::create_directories("results/density_frames");
	std::filesystem::create_directories("results");
	inqkit::io::RealField3DWriter frame_wr("results/density_frames",
		{ .field_name = "density", .include_meta = false, .emit_raw = false,
		  .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
		{ .overwrite = true });
	std::ofstream log("results/cap_trace.csv");
	log << std::setprecision(12) << "time_au,total_ha,num_electrons,dipole_z\n";
	int frame = 0;

	real_time::propagate(
		ions, electrons,
		[&](auto const & data) {
			if (data.root())
				log << data.time() << "," << data.energy().total() << ","
				    << data.num_electrons() << "," << data.dipole()[2] << "\n";
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
			.observables_current()
			.observables_dipole(),
		pert);

	if (electrons.root()) {
		std::ofstream s("results/run_summary.txt");
		s << std::setprecision(12)
		  << "run = localised_jellium/03_cap_stopping/wp_cap\n"
		  << "engine = inq-study\n"
		  << "cap = sin2 eta " << CAP_ETA_HA << " width_frac " << CAP_WIDTH
		  << " mid_frac +/-" << CAP_MID << " (7.5 Bohr/side)\n"
		  << "wp_launch_z = " << WP_LAUNCH_Z << "  wp_norm_after = " << report.norm_after << "\n"
		  << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "\n"
		  << "density_frames = " << frame << "\n";
	}
	std::cout << "Phase-5 WP+CAP run done: " << frame << " frames.\n";
	return 0;
}
