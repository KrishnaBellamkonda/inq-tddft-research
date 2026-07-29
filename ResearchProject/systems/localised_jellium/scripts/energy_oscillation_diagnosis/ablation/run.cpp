// ============================================================================
// localised_jellium / energy_oscillation_diagnosis / ablation / run.cpp
//
// ABLATION-LADDER probe binary for the ΔE_total energy-oscillation diagnosis
// campaign (docs/campaigns/localised_jellium/energy-oscillation-diagnosis.md).
//
// It is a diagnostic clone of muon_mass_fork/effmass_sigma1/wp/run.cpp with TWO
// changes only (no new numerics):
//   1. Turns ON the FULL energy decomposition on ObservableSelection
//      (energy_external / nonlocal / ion / ion_kinetic / eigenvalues / nvxc /
//      exact_exchange) and populates every StepContext energy field from
//      data.energy() — so we can see WHICH functional term carries the rise.
//   2. Ablation knobs (env, no code edit needed per probe):
//        EM_CAP=0/1     two-sided CAP off/on   (hypothesis a control)
//        EM_WP=0/1      inject the WP electron off/on (pure-GS / +v_bg floor)
//        EM_WRITE_EVERY (de-alias the oscillation; default 5 here)
//        EM_N_STEPS     tiny probe length (default 200)
//   Everything else (cell, GS, mass fork, ETRS, LDA, CAP geometry) is identical
//   to the reference run so the probe reproduces the SAME physics/numerics.
//
// Engine: inq-study (mass fork) — set INQ_SOURCE=.../inq-study before inq-run.
// Reuses the effmass_sigma1 GS checkpoint (shared_gs/slab_n52_L40x40x80_dx0p333).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/analytics.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include "../../../shared/configs/slab_n52_L40x40x80.hpp"
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN52_L40x40x80;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
	auto t0 = std::chrono::steady_clock::now();
	const std::string OUT     = "results/" + env_s("EM_OUT", "ablation");
	const double SPACING      = env_d("EM_SPACING", 0.33333);
	const double SIGMA_WP     = env_d("EM_SIGMA_WP", 1.0);
	const double K0           = env_d("EM_K0", 5.693);
	const double INV_MASS     = env_d("EM_INV_MASS", 0.476190);   // 1/2.10 (m_eff)
	const double LAUNCH_Z     = env_d("EM_LAUNCH_Z", -16.5);
	const double DT_AU        = env_d("EM_DT", 0.04);
	const int    N_STEPS      = env_i("EM_N_STEPS", 200);
	const int    WRITE_EVERY  = env_i("EM_WRITE_EVERY", 5);
	const bool   USE_CAP      = env_i("EM_CAP", 1) != 0;
	const bool   USE_WP       = env_i("EM_WP", 1) != 0;   // 0 => pure-GS / +v_bg floor
	const std::string GS_DIR  = env_s("EM_GS_DIR",
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n52_L40x40x80_dx0p333");
	const double CAP_ETA      = env_d("EM_CAP_ETA", -1.0);
	const double CAP_CENTER_B = env_d("EM_CAP_CENTER_BOHR", 32.5);
	const double CAP_WIDTH_B  = env_d("EM_CAP_WIDTH_BOHR", 15.0);
	const bool   USE_BG       = env_i("EM_BG", 1) != 0;   // 0 => drop the v_bg term (pure numerics floor)
	const double CAP_MID   = CAP_CENTER_B / Cfg::LZ_BOHR;
	const double CAP_WIDTH = CAP_WIDTH_B  / Cfg::LZ_BOHR;

	if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

	auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b).periodic();
	auto ions = systems::ions(cell);
	auto electrons = systems::electrons(
		ions,
		options::electrons{}
			.spacing(SPACING * 1.0_b)
			.extra_electrons(Cfg::N_ELECTRONS)
			.extra_states(Cfg::EXTRA_STATES)
			.temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
		input::kpoints::gamma());

	electrons.load(GS_DIR);

	int wp_idx = -1;
	if (USE_WP) {
		const double FOCUS_DIST = env_d("EM_FOCUS_DIST", 4.0);
		auto wp = inqkit::WavePacket{}
			          .center(0.0, 0.0, LAUNCH_Z)
			          .sigma(SIGMA_WP)
			          .k0(0.0, 0.0, K0)
			          .focus_z(FOCUS_DIST, 1.0/INV_MASS)
			          .orthogonalise_against_occupied(electrons);
		auto report = wp.inject_into_last_extra_state(electrons, 1.0);
		wp_idx = report.state_index;
		electrons.inverse_mass()[0][wp_idx] = INV_MASS;   // the mass fork
	}
	const int n_states = electrons.states().num_states();

	if (electrons.root())
		std::cout << std::setprecision(8)
			<< "\n=== energy_oscillation_diagnosis ABLATION out=" << OUT << " ===\n"
			<< "  spacing=" << SPACING << " dt=" << DT_AU << " n_steps=" << N_STEPS
			<< " write_every=" << WRITE_EVERY << "\n"
			<< "  WP=" << (USE_WP?"on":"off") << " CAP=" << (USE_CAP?"on":"off")
			<< " BG=" << (USE_BG?"on":"off") << " cap_eta=" << CAP_ETA << "\n"
			<< "  wp_idx=" << wp_idx << " n_states=" << n_states << "\n";

	// ----- background well (+ CAP), with ablation term-drops -----------------
	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation bg_pert(bg);
	const double eta_bg = USE_BG ? 1.0 : 0.0;   // scale the bg perturbation off if requested
	const double eta    = USE_CAP ? CAP_ETA : 0.0;
	perturbations::absorbing cap_lo(eta * 1.0_Ha, -CAP_MID, CAP_WIDTH);
	perturbations::absorbing cap_hi(eta * 1.0_Ha,  CAP_MID, CAP_WIDTH);
	// bg_pert is scaled by amplitude via its own perturbation; when EM_BG=0 we
	// still include a zero-strength CAP-only sum below. (bg_pert has no scalar
	// knob, so EM_BG=0 simply drops it from the sum.)
	auto cap_sum = perturbations::sum(cap_lo, cap_hi);
	// Compose the perturbation depending on the ablation:
	//   BG on : sum(bg_pert, cap_sum)   (cap_sum is zero-strength if CAP off)
	//   BG off: cap_sum only
	// We build both and select at propagate time via a lambda-free branch.
	(void)eta_bg;

	// ----- observables: FULL energy decomposition ON ------------------------
	std::filesystem::create_directories(OUT + "/raw/observables");

	inqkit::io::ObservableSelection sel;
	sel.step = sel.time_au = true;
	sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
	sel.energy_external = sel.energy_nonlocal = sel.energy_ion = true;
	sel.energy_ion_kinetic = sel.energy_exact_exchange = true;
	sel.energy_nvxc = sel.energy_eigenvalues = true;
	sel.density_l2 = true;
	inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
	obs_writer.write_header();

	// snapshot writer for the WP density difference (l2) — cheap, keeps parity
	inqkit::observables::DensityDelta density_delta(
		OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
		{.emit_raw_vti=false, .emit_coarse_vti=false, .compute_l2=true, .coarse_bin_bohr=3.0});
	{ auto s0 = inqkit::fields::density::total(electrons); density_delta.snapshot(s0, 0.0, 0); }

	auto record = [&](auto const& data){
		const int    step = data.iter();
		const double t    = data.time();
		if (step % WRITE_EVERY != 0) return;
		auto sys_f = inqkit::fields::density::total(electrons);
		const double l2 = density_delta.snapshot(sys_f, t, step);
		auto const& en = data.energy();
		inqkit::StepContext c; c.step = step; c.time_au = t;
		c.energy_total    = en.total();
		c.energy_kinetic  = en.kinetic();
		c.energy_hartree  = en.hartree();
		c.energy_xc       = en.xc();
		c.energy_external = en.external();
		c.energy_nonlocal = en.non_local();
		c.energy_ion      = en.ion();
		c.energy_ion_kinetic    = en.ion_kinetic();
		c.energy_exact_exchange = en.exact_exchange();
		c.energy_nvxc           = en.nvxc();
		c.energy_eigenvalues    = en.eigenvalues();
		c.density_l2 = l2;
		obs_writer.append(c);
	};

	auto opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
	if (USE_BG) {
		auto pert = perturbations::sum(bg_pert, cap_sum);
		real_time::propagate(ions, electrons, record, options::theory{}.lda(), opts, pert, 0);
	} else {
		real_time::propagate(ions, electrons, record, options::theory{}.lda(), opts, cap_sum, 0);
	}

	const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	if (electrons.root()) {
		std::ofstream s(OUT + "/run_summary.txt");
		s << std::setprecision(12)
		  << "run = localised_jellium/energy_oscillation_diagnosis/ablation\n"
		  << "engine = inq-study (mass fork)\nxc = LDA  propagator = ETRS  gamma-only\n"
		  << "cell_bohr = 40x40x80  spacing = " << SPACING << "\n"
		  << "WP = " << (USE_WP?"on":"off") << "  CAP = " << (USE_CAP?"on":"off")
		  << "  BG = " << (USE_BG?"on":"off") << "\n"
		  << "cap_eta = " << CAP_ETA << "  cap_center_bohr = " << CAP_CENTER_B
		  << "  cap_width_bohr = " << CAP_WIDTH_B << "\n"
		  << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
		  << "wp_idx = " << wp_idx << "  n_states = " << n_states << "\n"
		  << "gs_dir = " << GS_DIR << "\nwall_time_s = " << wall << "\n"
		  << "full_energy_decomposition = true\n"
		  << "run_completed = true\n";
		std::cout << "  ablation run done (" << wall << " s), " << N_STEPS << " steps.\n";
	}
	return 0;
}
