// ============================================================================
// localised_jellium / 01_slab_validation / gs_slab / run.cpp
//
// Ground state of the localised jellium SLAB: N=234 electrons confined by a
// static positive background well (inqkit::jellium::localised_background_
// perturbation), 50 Bohr cubic periodic box, LDA. The background is injected
// into the KS potential via the perturbation hook and is present throughout the
// SCF, so the electrons self-consistently localise inside the slab.
//
// Engine: stock inq (no CAP here). The checkpoint is verified to load on
// inq-study before the CAP phase (Phase 5).
//
// Validation this run feeds: T1 (SCF converges; density peaks inside slab;
// interior density flat to a few %) — see docs/plans/localised-jellium.md.
// ============================================================================
#include <inq/inq.hpp>

#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/analytics.hpp>

#include "../../../shared/configs/slab_n234_L50.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN234L50;

int main() {
	const std::string CHECKPOINT_DIR =
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n234_L50";

	std::cout << "\n=== localised_jellium GS: slab_n234_L50 ===\n"
	          << "  cell      = " << Cfg::L_BOHR << "^3 Bohr (cubic, periodic)\n"
	          << "  slab      = half-width " << Cfg::SLAB_HALF_WIDTH
	          << " Bohr along axis " << Cfg::SLAB_AXIS << " (25 Bohr thick)\n"
	          << "  N         = " << Cfg::N_ELECTRONS
	          << "  n0 = " << Cfg::N0 << " a0^-3  (r_s = "
	          << inqkit::jellium::rs_from_n0(Cfg::N0) << ")\n"
	          << "  spacing   = " << Cfg::SPACING_BOHR << " Bohr\n"
	          << "  checkpoint= " << CHECKPOINT_DIR << "\n\n";

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

	// ---- Localised positive background (the confining well) ----------------
	inqkit::jellium::localised_background_params bg;
	bg.shape       = inqkit::jellium::background_shape::slab;
	bg.n0          = Cfg::N0;
	bg.half_width  = Cfg::SLAB_HALF_WIDTH;
	bg.slab_axis   = Cfg::SLAB_AXIS;
	bg.center      = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR};
	bg.edge_width  = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation pert(bg);

	ground_state::initial_guess(ions, electrons);
	auto gs = ground_state::calculate(
		ions, electrons,
		options::theory{}.lda(),
		options::ground_state{}
			.energy_tolerance(Cfg::SCF_TOL_HA * 1.0_Ha)
			.max_steps(Cfg::SCF_MAX_STEPS)
			.broyden_mixing()
			.mixing_ndim(Cfg::SCF_MIX_NDIM)
			.mixing(Cfg::SCF_MIX_ALPHA),
		pert);

	std::cout << "  GS energy = " << gs.energy.total() << " Ha\n";

	const int n_states    = electrons.states().num_states();
	const int n_electrons = electrons.states().num_electrons();
	std::cout << "  num_states = " << n_states
	          << "  num_electrons = " << n_electrons << "\n";

	std::filesystem::create_directories(CHECKPOINT_DIR);
	electrons.save(CHECKPOINT_DIR);

	// GS electron density for T1 (interior flat? density peaks in slab?).
	std::filesystem::create_directories("results/density_gs_system");
	{
		inqkit::io::RealField3DWriter gs_wr("results/density_gs_system",
			{ .field_name = "density", .include_meta = false, .emit_raw = false,
			  .emit_vti = true,
			  .vti_format = inqkit::io::VTIWriteOptions::Format::binary },
			{ .overwrite = true });
		gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");
	}

	if (electrons.root()) {
		std::filesystem::create_directories("results");
		std::ofstream s("results/run_summary.txt");
		s << std::setprecision(16);
		s << "run = localised_jellium/01_slab_validation/gs_slab\n"
		  << "engine = inq (stock)\n"
		  << "cell_bohr = " << Cfg::L_BOHR << "^3 (cubic, periodic)\n"
		  << "xc = LDA\n"
		  << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH
		  << " axis " << Cfg::SLAB_AXIS << " edge_width " << Cfg::EDGE_WIDTH_BOHR << "\n"
		  << "n0_a0m3 = " << Cfg::N0 << "\n"
		  << "r_s = " << inqkit::jellium::rs_from_n0(Cfg::N0) << "\n"
		  << "spacing_bohr = " << Cfg::SPACING_BOHR << "\n"
		  << "extra_electrons = " << Cfg::N_ELECTRONS << "\n"
		  << "extra_states = " << Cfg::EXTRA_STATES << "\n"
		  << "temperature_ev = " << Cfg::TEMPERATURE_EV << "\n"
		  << "scf_tol_ha = " << Cfg::SCF_TOL_HA << "\n"
		  << "ground_state_energy_ha = " << gs.energy.total() << "\n"
		  << "num_states = " << n_states << "\n"
		  << "num_electrons = " << n_electrons << "\n"
		  << "checkpoint_dir = " << CHECKPOINT_DIR << "\n";
	}

	std::cout << "Done.\n";
	return 0;
}
