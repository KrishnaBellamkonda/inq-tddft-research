// ============================================================================
// localised_jellium / qsp_phase1 / sie / run.cpp   (campaign: quantum-stopping-power, P1.2)
//
// SIE diagnostic short RT. Loads the slab_n82_L50x50x70 GS, injects the σ=0.5 /
// 100 eV wavepacket FAR from the slab (z=−32, CAP OFF), and propagates a few
// steps. We need only the t=0 quantities:
//   * E_total(0)  from observables.csv (energy_total)
//   * KE_WP       from wp_momentum_stats.csv (e_kin_ha = ⟨p²⟩/2)
// SIE_a = E_total(0) − (E_GS_slab + 100 eV)        [user ref; = SIE + zero-point]
// SIE_b = E_total(0) − E_GS_slab − KE_WP           [= SIE]
// cross-check: SIE_a − SIE_b ≈ zero-point KE 3/(4σ²) = 81.6 eV.
//
// CAP off ⇒ stock inq is sufficient (no absorbing perturbation). A few steps also
// confirm E_total(0)/KE_WP are stable (guard rail) and the WP norm = 1 at launch.
// ============================================================================
#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>

#include "../../../shared/configs/slab_n82_L50x50x70.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN82_L50x50x70;

int main() {
	const std::string GS_DIR =
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n82_L50x50x70";
	const std::string OUT = "results/sie_wp_far";
	const double DT_AU   = 0.02;
	const int    N_STEPS = 10;          // a few steps; only t=0 is needed
	const double LAUNCH_Z = Cfg::WP_CZ_BOHR;   // −32 (far, CAP off)

	std::cout << "\n=== qsp_phase1 SIE short-RT (WP far, CAP off) ===\n"
	          << "  GS=" << GS_DIR << "  launch_z=" << LAUNCH_Z
	          << "  N_STEPS=" << N_STEPS << " dt=" << DT_AU << "\n";

	auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b,
	                                         Cfg::LZ_BOHR * 1.0_b).periodic();
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
	std::cout << "  Loaded GS from " << GS_DIR << "\n";

	std::filesystem::create_directories(OUT + "/raw/observables");

	// ----- WP injection (far from slab) ---------------------------------
	auto wp = inqkit::WavePacket{}
	              .center(0.0, 0.0, LAUNCH_Z).sigma(Cfg::WP_SIGMA_BOHR)
	              .k0(0.0, 0.0, Cfg::WP_KZ).orthogonalise_against_occupied(electrons);
	auto report = wp.inject_into_last_extra_state(electrons, 1.0);
	const int wp_idx = report.state_index;
	std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << report.norm_after
	          << " max_overlap=" << report.max_overlap << "\n";

	// ----- background well only (NO CAP) --------------------------------
	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation bg_pert(bg);

	// ----- observables (energies) + WP momentum/real-space stats --------
	inqkit::io::ObservableSelection sel;
	sel.step = sel.time_au = true;
	sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
	inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
	obs_writer.write_header();
	inqkit::observables::WPMomentumStats wp_momentum_stats(
		OUT + "/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every=1});
	inqkit::observables::WPRealSpaceStats wp_real_space_stats(
		OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every=1});

	inqkit::RealTimeSession rt_obs(ions, electrons, 1);
	rt_obs.add([&](inqkit::StepContext const& ctx) { obs_writer.append(ctx); });

	auto step_fn = [&](auto const& data) {
		rt_obs.step(data);
		wp_momentum_stats.maybe_accumulate(data);
		wp_real_space_stats.maybe_accumulate(data);
	};

	auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
	real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, bg_pert);

	if (electrons.root()) {
		std::ofstream s(OUT + "/raw/run_summary.txt");
		s << std::setprecision(16);
		s << "run = localised_jellium/qsp_phase1/sie (quantum-stopping-power P1.2)\n"
		  << "gs_dir = " << GS_DIR << "\n"
		  << "cell_bohr = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x " << Cfg::LZ_BOHR << "\n"
		  << "wp_sigma_bohr = " << Cfg::WP_SIGMA_BOHR << "\n"
		  << "wp_ekin_ev = " << Cfg::WP_EKIN_EV << "  wp_k0 = " << Cfg::WP_K0 << "\n"
		  << "launch_z_bohr = " << LAUNCH_Z << "\n"
		  << "wp_norm_after = " << report.norm_after << "\n"
		  << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "\n"
		  << "cap = OFF\n";
	}
	std::cout << "Done.\n";
	return 0;
}
