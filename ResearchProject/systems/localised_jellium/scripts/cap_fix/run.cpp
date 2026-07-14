// ============================================================================
// localised_jellium / cap_fix / run.cpp
//
// EXPERIMENT binary for the CAP energy-artifact FIX campaign
// (docs/campaigns/localised_jellium/cap-fix-experimentation.md).
//
// Clone of energy_oscillation_diagnosis/ablation/run.cpp (which is itself a
// diagnostic clone of muon_mass_fork/effmass_sigma1/wp/run.cpp) with THREE
// additions and no other numerics changes:
//   1. EM_CAP_MODE=two|wrap — CAP topology:
//        two  : the standard two sin² bumps at ±EM_CAP_CENTER_BOHR (inq
//               perturbations::absorbing). NOTE these fall to W=0 exactly at
//               the periodic boundary z=±L/2.
//        wrap : ONE smooth cos² bump of full width EM_WRAP_WIDTH_BOHR centred
//               ON the periodic boundary plane (inqkit absorbing_wrap). With
//               EM_WRAP_WIDTH_BOHR = 2×15 = 30 this has the SAME footprint
//               (|z|>25) and SAME ∫W dz as the default two-sided CAP —
//               topology is the only difference.
//   2. charge.csv — ∫n dV (total electron count of electrons.density()) per
//      write step. Closes the diagnosis Part-IV instrumentation gap; feeds
//      the absorbed_frac metric. (Caveat: whether the WP state is included in
//      electrons.density() is config-dependent — read the t=0 value.)
//   3. Full energy decomposition ON (as in the ablation binary).
//
// Ablation knobs kept: EM_CAP, EM_WP, EM_BG, EM_CAP_ETA, EM_CAP_CENTER_BOHR,
// EM_CAP_WIDTH_BOHR, EM_N_STEPS, EM_WRITE_EVERY, EM_DT, EM_OUT, EM_GS_DIR.
//
// Engine: inq-study (mass fork) — set INQ_SOURCE=.../inq-study before inq-run.
// Reuses the effmass_sigma1 GS checkpoint (shared_gs/slab_n52_L40x40x80_dx0p333).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/perturbations/absorbing_wrap.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/analytics.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include "../../shared/configs/slab_n52_L40x40x80.hpp"
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN52_L40x40x80;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
	auto t0 = std::chrono::steady_clock::now();
	const std::string OUT     = "results/" + env_s("EM_OUT", "cap_fix");
	const double SPACING      = env_d("EM_SPACING", 0.33333);
	const double SIGMA_WP     = env_d("EM_SIGMA_WP", 1.0);
	const double K0           = env_d("EM_K0", 5.693);
	const double INV_MASS     = env_d("EM_INV_MASS", 0.476190);   // 1/2.10 (m_eff)
	const double LAUNCH_Z     = env_d("EM_LAUNCH_Z", -16.5);
	const double DT_AU        = env_d("EM_DT", 0.04);
	const int    N_STEPS      = env_i("EM_N_STEPS", 700);
	const int    WRITE_EVERY  = env_i("EM_WRITE_EVERY", 5);
	const bool   USE_CAP      = env_i("EM_CAP", 1) != 0;
	const bool   USE_WP       = env_i("EM_WP", 1) != 0;
	const std::string GS_DIR  = env_s("EM_GS_DIR",
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n52_L40x40x80_dx0p333");
	const double CAP_ETA      = env_d("EM_CAP_ETA", -1.0);
	const double CAP_CENTER_B = env_d("EM_CAP_CENTER_BOHR", 32.5);
	const double CAP_WIDTH_B  = env_d("EM_CAP_WIDTH_BOHR", 15.0);
	const std::string CAP_MODE = env_s("EM_CAP_MODE", "two");     // two | wrap
	const double WRAP_WIDTH_B = env_d("EM_WRAP_WIDTH_BOHR", 30.0);
	const bool   USE_BG       = env_i("EM_BG", 1) != 0;
	// z-periodicity of the cell (Arm-B campaign, 2026-07-14): 3 = fully
	// periodic (ALL prior oscillation runs), 2 = periodic x,y + open z
	// (slab-truncated Poisson — no electrostatic images along z). The GS
	// loaded via EM_GS_DIR MUST have been converged at the SAME periodicity.
	const int    PERIODICITY  = env_i("EM_PERIODICITY", 3);
	const double CAP_MID   = CAP_CENTER_B / Cfg::LZ_BOHR;
	const double CAP_WIDTH = CAP_WIDTH_B  / Cfg::LZ_BOHR;
	const double WRAP_WIDTH = WRAP_WIDTH_B / Cfg::LZ_BOHR;

	if (CAP_MODE != "two" && CAP_MODE != "wrap") {
		std::cerr << "FATAL: EM_CAP_MODE must be 'two' or 'wrap', got: " << CAP_MODE << "\n";
		return 2;
	}
	if (PERIODICITY != 2 && PERIODICITY != 3) {
		std::cerr << "FATAL: EM_PERIODICITY must be 2 or 3, got " << PERIODICITY << "\n";
		return 2;
	}
	if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

	auto cell0 = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b);
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
			<< "\n=== cap_fix EXPERIMENT out=" << OUT << " ===\n"
			<< "  spacing=" << SPACING << " dt=" << DT_AU << " n_steps=" << N_STEPS
			<< " write_every=" << WRITE_EVERY << "\n"
			<< "  WP=" << (USE_WP?"on":"off") << " CAP=" << (USE_CAP?"on":"off")
			<< " BG=" << (USE_BG?"on":"off") << " cap_eta=" << CAP_ETA
			<< " cap_mode=" << CAP_MODE << "\n"
			<< "  wp_idx=" << wp_idx << " n_states=" << n_states << "\n";

	// ----- background well + CAP (topology-selectable) -----------------------
	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation bg_pert(bg);

	// Exactly one topology is active: the other gets amplitude 0 (its kernel
	// then adds nothing). EM_CAP=0 zeroes both — the CAP-off control.
	const double eta      = USE_CAP ? CAP_ETA : 0.0;
	const double eta_two  = (CAP_MODE == "two")  ? eta : 0.0;
	const double eta_wrap = (CAP_MODE == "wrap") ? eta : 0.0;
	perturbations::absorbing cap_lo(eta_two * 1.0_Ha, -CAP_MID, CAP_WIDTH);
	perturbations::absorbing cap_hi(eta_two * 1.0_Ha,  CAP_MID, CAP_WIDTH);
	inqkit::perturbations::absorbing_wrap cap_wrap(eta_wrap * 1.0_Ha, WRAP_WIDTH);
	auto cap_sum = perturbations::sum(perturbations::sum(cap_lo, cap_hi), cap_wrap);

	// ----- observables: FULL energy decomposition + charge -------------------
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

	// charge ledger: total electron count ∫n dV per write step
	std::ofstream charge_csv(OUT + "/raw/observables/charge.csv");
	charge_csv << std::setprecision(15) << "step,time_au,n_total\n";
	auto integrate = [](inqkit::fields::RealField3D const& f){
		const double dV = f.dx_bohr * f.dy_bohr * f.dz_bohr;
		long double s = 0.0L;
		for (auto v : f.values) s += v;
		return static_cast<double>(s * dV);
	};

	inqkit::observables::DensityDelta density_delta(
		OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
		{.emit_raw_vti=false, .emit_coarse_vti=false, .compute_l2=true, .coarse_bin_bohr=3.0});
	{
		auto s0 = inqkit::fields::density::total(electrons);
		density_delta.snapshot(s0, 0.0, 0);
		charge_csv << 0 << "," << 0.0 << "," << integrate(s0) << "\n";
	}

	auto record = [&](auto const& data){
		const int    step = data.iter();
		const double t    = data.time();
		if (step % WRITE_EVERY != 0) return;
		auto sys_f = inqkit::fields::density::total(electrons);
		const double l2 = density_delta.snapshot(sys_f, t, step);
		charge_csv << step << "," << t << "," << integrate(sys_f) << "\n";
		charge_csv.flush();
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
		  << "run = localised_jellium/cap_fix\n"
		  << "engine = inq-study (mass fork)\nxc = LDA  propagator = ETRS  gamma-only\n"
		  << "cell_bohr = 40x40x80  spacing = " << SPACING << "\n"
		  << "WP = " << (USE_WP?"on":"off") << "  CAP = " << (USE_CAP?"on":"off")
		  << "  BG = " << (USE_BG?"on":"off") << "\n"
		  << "cap_mode = " << CAP_MODE << "\n"
		  << "periodicity = " << PERIODICITY << "\n"
		  << "cap_eta = " << CAP_ETA << "  cap_center_bohr = " << CAP_CENTER_B
		  << "  cap_width_bohr = " << CAP_WIDTH_B << "\n"
		  << "wrap_width_bohr = " << WRAP_WIDTH_B << "\n"
		  << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
		  << "wp_idx = " << wp_idx << "  n_states = " << n_states << "\n"
		  << "gs_dir = " << GS_DIR << "\nwall_time_s = " << wall << "\n"
		  << "full_energy_decomposition = true\ncharge_csv = true\n"
		  << "run_completed = true\n";
		std::cout << "  cap_fix run done (" << wall << " s), " << N_STEPS << " steps.\n";
	}
	return 0;
}
