// ============================================================================
// localised_jellium / muon_mass_fork / effmass_sigma2 / quantum / run.cpp
//
// QUANTUM effective-mass WAVEPACKET projectile through the r_s=5.69 slab,
// RE-PLANNED for <=12 h (N3/A, 2026-07-08). sigma_WP=2, m_eff=2.506 m_e,
// v=2.711 a.u. (= 100 eV electron velocity -> same S(v)), E=251 eV, k0=6.7933,
// dx=0.40, dt=0.05, N=846 steps on the shrunk 36x36x80 (N=42) slab.
// Gaussian WP injected into the last extra state with a tuned inverse_mass
// (the inq-study mass fork). Stopping read from the n(k,t) coherent peak.
// Loads shared_gs/slab_n82_L50x50x90_dx0p40. Adapted from effmass_pair/quantum.
//
// CHECKPOINT / RESUME (pause & continue):
//   * Every EM_CKPT_EVERY steps, electrons.save(EM_RT_CKPT_DIR) + rt_state.txt.
//   * EM_RESUME=1 loads the RT checkpoint, re-applies the inverse_mass fork
//     (save/load does NOT persist it), reads start_step from rt_state.txt, and
//     hands it to real_time::propagate's native start_step. Static jellium ions
//     + static background/CAP make the ETRS restart bit-faithful. On resume the
//     CSVs are written to `.from<START>` segment files (concatenate in analysis);
//     VTIs continue in the shared per-step dirs.
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
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include "../../../../shared/configs/slab_n82_L50x50x90.hpp"
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN82_L50x50x90;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

// Read "last_step" from rt_state.txt (returns -1 if absent/unreadable).
static int read_last_step(const std::string& path) {
	std::ifstream f(path); if (!f) return -1;
	std::string line; int last = -1;
	while (std::getline(f, line)) {
		auto p = line.find("last_step=");
		if (p == 0) last = std::atoi(line.c_str() + 10);
	}
	return last;
}
static int read_kv_int(const std::string& path, const std::string& key, int dflt) {
	std::ifstream f(path); if (!f) return dflt;
	std::string line;
	while (std::getline(f, line)) {
		auto p = line.find(key + "=");
		if (p == 0) return std::atoi(line.c_str() + key.size() + 1);
	}
	return dflt;
}

int main() {
	auto t0 = std::chrono::steady_clock::now();
	const std::string OUT     = "results/" + env_s("EM_OUT", "quantum");
	const double SPACING      = env_d("EM_SPACING", 0.40);
	const double SIGMA_WP     = env_d("EM_SIGMA_WP", 2.0);
	const double K0           = env_d("EM_K0", 6.793);
	const double INV_MASS     = env_d("EM_INV_MASS", 0.3991);   // 1/2.506
	const double LAUNCH_Z     = env_d("EM_LAUNCH_Z", -19.21);
	const double DT_AU        = env_d("EM_DT", 0.05);
	const int    N_STEPS      = env_i("EM_N_STEPS", 819);
	const int    WRITE_EVERY  = env_i("EM_WRITE_EVERY", 20);
	const int    CKPT_EVERY   = env_i("EM_CKPT_EVERY", 234);
	const bool   RESUME       = env_i("EM_RESUME", 0) != 0;
	const bool   USE_CAP      = env_i("EM_CAP", 1) != 0;
	const std::string GS_DIR  = env_s("EM_GS_DIR",
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n82_L50x50x90_dx0p40");
	const std::string RT_CKPT_DIR  = env_s("EM_RT_CKPT_DIR", OUT + "/rt_ckpt");
	const std::string RT_STATE_TXT = RT_CKPT_DIR + "/rt_state.txt";
	// CAP retuned for Lz=80: centre |z|=32 (0.40*Lz), 8 Bohr/side (0.10*Lz), inner face ~28.
	const double CAP_ETA = -1.0, CAP_MID = 40.0/90.0, CAP_WIDTH = 10.0/90.0;

	if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

	// ----- resume bookkeeping -----------------------------------------------
	int START = 0;
	if (RESUME) {
		if (!std::filesystem::exists(RT_CKPT_DIR)) { std::cerr << "FATAL: resume but no RT ckpt: " << RT_CKPT_DIR << "\n"; return 2; }
		START = read_last_step(RT_STATE_TXT);
		if (START < 0) { std::cerr << "FATAL: resume but unreadable " << RT_STATE_TXT << "\n"; return 2; }
		if (START >= N_STEPS) { std::cout << "Already at/after target (" << START << ">=" << N_STEPS << "); nothing to do.\n"; return 0; }
	}
	const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

	auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b).periodic();
	auto ions = systems::ions(cell);                     // jellium: no nuclei (static ions -> restart OK)
	auto electrons = systems::electrons(
		ions,
		options::electrons{}
			.spacing(SPACING * 1.0_b)
			.extra_electrons(Cfg::N_ELECTRONS)
			.extra_states(Cfg::EXTRA_STATES)
			.temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
		input::kpoints::gamma());

	int wp_idx = -1;
	if (RESUME) {
		electrons.load(RT_CKPT_DIR);                     // already holds the propagated WP
		wp_idx = read_kv_int(RT_STATE_TXT, "wp_idx", electrons.states().num_states() - 1);
	} else {
		electrons.load(GS_DIR);
		auto wp = inqkit::WavePacket{}
		              .center(0.0, 0.0, LAUNCH_Z)
		              .sigma(SIGMA_WP)
		              .k0(0.0, 0.0, K0)
		              .orthogonalise_against_occupied(electrons);
		auto report = wp.inject_into_last_extra_state(electrons, 1.0);
		wp_idx = report.state_index;
	}
	const int n_states = electrons.states().num_states();
	electrons.inverse_mass()[0][wp_idx] = INV_MASS;      // the mass fork (re-applied on resume)

	if (electrons.root())
		std::cout << std::setprecision(8)
			<< "\n=== effmass_sigma2 QUANTUM out=" << OUT << (RESUME?"  [RESUME]":"  [FRESH]") << " ===\n"
			<< "  spacing=" << SPACING << " sigma_WP=" << SIGMA_WP << " k0=" << K0
			<< " inv_mass=" << INV_MASS << " (m=" << 1.0/INV_MASS << ")\n"
			<< "  launch_z=" << LAUNCH_Z << " dt=" << DT_AU << " start_step=" << START
			<< " n_steps=" << N_STEPS << " ckpt_every=" << CKPT_EVERY << " cap=" << (USE_CAP?"on":"off") << "\n"
			<< "  WP idx=" << wp_idx << " n_states=" << n_states << "\n";

	// ----- background well (+ CAP) ------------------------------------------
	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation bg_pert(bg);
	const double eta = USE_CAP ? CAP_ETA : 0.0;
	perturbations::absorbing cap_lo(eta * 1.0_Ha, -CAP_MID, CAP_WIDTH);
	perturbations::absorbing cap_hi(eta * 1.0_Ha,  CAP_MID, CAP_WIDTH);
	auto pert = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

	// ----- observables (segment-suffixed on resume) -------------------------
	std::filesystem::create_directories(OUT + "/raw/observables");
	for (auto sub : {"density_total","density_wp","density_gs_system","density_delta","density_delta_coarse"})
		std::filesystem::create_directories(OUT + "/raw/vti/" + sub);

	inqkit::io::ObservableSelection sel;
	sel.step = sel.time_au = true;
	sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
	sel.density_l2 = true;
	inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables" + SEG + ".csv", sel);
	obs_writer.write_header();
	inqkit::observables::WPRealSpaceStats wp_rs(
		OUT + "/raw/observables/wp_real_space_stats" + SEG + ".csv", wp_idx, {.write_every = WRITE_EVERY});
	inqkit::observables::WPMomentumStats wp_ms(
		OUT + "/raw/observables/wp_momentum_stats" + SEG + ".csv", wp_idx, {.write_every = WRITE_EVERY});

	inqkit::io::RealField3DLayout vti_layout{
		.field_name = "density", .include_meta = false, .emit_raw = false,
		.emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
	if (!RESUME) { inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
	               gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system"); }
	inqkit::io::RealField3DWriter total_wr(OUT + "/raw/vti/density_total", vti_layout, {.overwrite = !RESUME});
	inqkit::io::RealField3DWriter wp_wr   (OUT + "/raw/vti/density_wp",    vti_layout, {.overwrite = !RESUME});
	inqkit::observables::DensityDelta density_delta(
		OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
		{.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});
	if (!RESUME) { auto s0 = inqkit::fields::density::total(electrons);
	  total_wr.write(s0,0.0,0); density_delta.snapshot(s0,0.0,0);
	  wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0); }

	std::filesystem::create_directories(RT_CKPT_DIR);

	// ----- propagate (LDA, ETRS default; native start_step restart) ---------
	auto func = [&](auto const& data){
		const int    step = data.iter();     // authoritative global step (resume-safe)
		const double t    = data.time();
		wp_rs.maybe_accumulate(data);
		wp_ms.maybe_accumulate(data);
		if (step % WRITE_EVERY == 0) {
			auto sys_f = inqkit::fields::density::total(electrons);
			total_wr.write(sys_f, t, step);
			wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), t, step);
			const double l2 = density_delta.snapshot(sys_f, t, step);
			inqkit::StepContext c; c.step = step; c.time_au = t;
			c.energy_total = data.energy().total(); c.energy_kinetic = data.energy().kinetic();
			c.energy_hartree = data.energy().hartree(); c.energy_xc = data.energy().xc();
			c.density_l2 = l2;
			obs_writer.append(c);
		}
		if (step > START && step % CKPT_EVERY == 0 && step < N_STEPS) {
			electrons.save(RT_CKPT_DIR);      // collective; overwrites the checkpoint
			if (electrons.root()) {
				std::ofstream st(RT_STATE_TXT, std::ios::trunc);
				st << "last_step=" << step << "\ntime_au=" << t << "\nwp_idx=" << wp_idx
				   << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS
				   << "\ndt_au=" << DT_AU << "\n";
				std::cout << "  [ckpt] saved at step " << step << " (t=" << t << ")\n";
			}
		}
	};
	auto opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
	real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts, pert, START);

	// final checkpoint at the end (so a follow-on resume is a clean no-op)
	electrons.save(RT_CKPT_DIR);
	if (electrons.root()) {
		std::ofstream st(RT_STATE_TXT, std::ios::trunc);
		st << "last_step=" << N_STEPS << "\ntime_au=" << N_STEPS*DT_AU << "\nwp_idx=" << wp_idx
		   << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS << "\ndt_au=" << DT_AU << "\n";
	}

	const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	if (electrons.root()) {
		std::ofstream s(OUT + "/run_summary" + SEG + ".txt");
		s << std::setprecision(12)
		  << "run = localised_jellium/muon_mass_fork/effmass_sigma2/quantum\n"
		  << "engine = inq-study (mass fork)\nxc = LDA  propagator = ETRS  gamma-only\n"
		  << "resume = " << (RESUME?"true":"false") << "  start_step = " << START << "\n"
		  << "cell_bohr = 50x50x90  spacing = " << SPACING << "\n"
		  << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS << "\n"
		  << "n0_a0m3 = " << Cfg::N0 << "  r_s = " << inqkit::jellium::rs_from_n0(Cfg::N0) << "\n"
		  << "projectile = quantum WP  sigma_WP = " << SIGMA_WP << "  k0 = " << K0
		  << "  inverse_mass = " << INV_MASS << "  m_eff = " << 1.0/INV_MASS << "\n"
		  << "velocity_au = " << K0*INV_MASS << "  E_eV = " << 0.5*(1.0/INV_MASS)*std::pow(K0*INV_MASS,2)*27.211386 << "\n"
		  << "launch_z = " << LAUNCH_Z << "  cap = " << (USE_CAP?"on":"off") << "\n"
		  << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY
		  << "  ckpt_every = " << CKPT_EVERY << "\n"
		  << "wp_idx = " << wp_idx << "  n_states = " << n_states << "\n"
		  << "gs_dir = " << GS_DIR << "\nrt_ckpt_dir = " << RT_CKPT_DIR << "\nwall_time_s = " << wall << "\n"
		  << "run_completed = true\n";
		std::cout << "  quantum run done (" << wall << " s), final step " << N_STEPS << ".\n";
	}
	return 0;
}
