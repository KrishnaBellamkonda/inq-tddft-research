// ============================================================================
// localised_jellium / mass_pair_n162 / wp / run.cpp
//
// Matched quantum-wavepacket projectile through the GENUINE 162-electron
// localised jellium slab. ONE binary, dispatched twice (only the projectile
// mass changes): m1 (EM_INV_MASS=1.0, EM_K0=2.711) and m2 (EM_INV_MASS=0.5,
// EM_K0=3.834). Bath electrons stay mass 1 in BOTH (the per-orbital inq-study
// mass fork is applied to the WP state ONLY). E=100 eV, sigma_WP=1, plain
// (unchirped) launch at z=-16.5, two-sided CAP 10 Bohr/side eta=-1.0, 100 a.u.
// Plan: docs/plans/mass-pair-n162-sigma1-cap.md (user spec 2026-07-19).
//
// Merges the fullsuite_wp EXTENSIVE observable suite with the muon_mass_fork
// mass fork + checkpoint/resume. Adds a PER-STEP full energy decomposition
// (energy_decomp.csv: all 11 INQ energy terms every timestep, user spec).
//
// CHECKPOINT / RESUME (user: checkpoint every 500 steps):
//   * Every EM_CKPT_EVERY(500) steps + at the end: electrons.save(EM_RT_CKPT_DIR)
//     + rt_state.txt (last_step/wp_idx/dt). EM_RESUME=1 loads the RT checkpoint,
//     re-applies the inverse_mass fork (save/load does NOT persist it), reads
//     start_step, hands it to real_time::propagate's native start_step, and writes
//     CSVs to `.from<START>` segment files (concatenate in analysis).
//
// Env: EM_OUT EM_SPACING(0.40) EM_SIGMA_WP(1.0) EM_K0(2.711) EM_INV_MASS(1.0)
//      EM_LAUNCH_Z(-16.5) EM_DT(0.04) EM_N_STEPS(2500) EM_WRITE_EVERY(8)
//      EM_WF_EVERY(40) EM_CKPT_EVERY(500) EM_RESUME(0) EM_CAP(1)
//      EM_CAP_ETA(-1.0) EM_CAP_CENTER_BOHR(26.2) EM_CAP_WIDTH_BOHR(10.0) EM_GS_DIR.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/analytics.hpp>

#include "../../../shared/configs/slab_n120_L60x60x62.hpp"
#include "../../../../jellium/shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN120_L60x60x62;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static int read_last_step(const std::string& path) {
	std::ifstream f(path); if (!f) return -1;
	std::string line; int last = -1;
	while (std::getline(f, line)) { auto p = line.find("last_step="); if (p == 0) last = std::atoi(line.c_str() + 10); }
	return last;
}
static int read_kv_int(const std::string& path, const std::string& key, int dflt) {
	std::ifstream f(path); if (!f) return dflt;
	std::string line;
	while (std::getline(f, line)) { auto p = line.find(key + "="); if (p == 0) return std::atoi(line.c_str() + key.size() + 1); }
	return dflt;
}

int main() {
	auto t0 = std::chrono::steady_clock::now();
	const std::string OUT     = "results/" + env_s("EM_OUT", "m1");
	const double SPACING      = env_d("EM_SPACING", 0.40);
	const double SIGMA_WP     = env_d("EM_SIGMA_WP", 1.0);
	const double K0           = env_d("EM_K0", 2.711);
	const double INV_MASS     = env_d("EM_INV_MASS", 1.0);      // 1.0 => m=1; 0.5 => m=2
	const double LAUNCH_Z     = env_d("EM_LAUNCH_Z", -16.5);
	const double DT_AU        = env_d("EM_DT", 0.04);
	const int    N_STEPS      = env_i("EM_N_STEPS", 2500);      // 100 a.u.
	const int    WRITE_EVERY  = env_i("EM_WRITE_EVERY", 8);     // ~313 frames
	const int    WF_EVERY     = env_i("EM_WF_EVERY", 40);       // complex WF VTI (heavy)
	const int    CKPT_EVERY   = env_i("EM_CKPT_EVERY", 500);    // user: every 500 steps
	const bool   RESUME       = env_i("EM_RESUME", 0) != 0;
	const bool   USE_CAP      = env_i("EM_CAP", 1) != 0;
	// Observable toggles (default on) — for bisecting the step-0 GPU crash on the
	// large orthorhombic grid. The fullsuite-only extras (momentum/WF/overlap/state,
	// current/dipole) were only ever validated on cubic boxes.
	const bool OBS_MOM    = env_i("EM_OBS_MOM", 1) != 0;   // momentum_distribution (FFT of all states)
	const bool OBS_WF     = env_i("EM_OBS_WF", 1) != 0;    // complex wavefunction VTI
	const bool OBS_OVL    = env_i("EM_OBS_OVL", 1) != 0;   // orbital-overlap matrices
	const bool OBS_STATE  = env_i("EM_OBS_STATE", 1) != 0; // state_energies + occupations
	const bool OBS_DIPCUR = env_i("EM_OBS_DIPCUR", 1) != 0;// current + dipole observables
	const std::string GS_DIR  = env_s("EM_GS_DIR",
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n120_L60x60x62_dx0p40");
	const std::string RT_CKPT_DIR  = env_s("EM_RT_CKPT_DIR", OUT + "/rt_ckpt");
	const std::string RT_STATE_TXT = RT_CKPT_DIR + "/rt_state.txt";

	// CAP: 10 Bohr/side at the cell boundaries. For Lz=62.4 the outer face is at
	// |z|=31.2 (=Lz/2), inner face |z|=21.2, centre |z|=26.2, full width 10.0.
	const double CAP_ETA      = env_d("EM_CAP_ETA", -1.0);
	const double CAP_CENTER_B = env_d("EM_CAP_CENTER_BOHR", 26.2);
	const double CAP_WIDTH_B  = env_d("EM_CAP_WIDTH_BOHR", 10.0);
	const double CAP_MID   = CAP_CENTER_B / Cfg::LZ_BOHR;
	const double CAP_WIDTH = CAP_WIDTH_B  / Cfg::LZ_BOHR;

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
	auto ions = systems::ions(cell);
	auto electrons = systems::electrons(
		ions,
		options::electrons{}
			.spacing(SPACING * 1.0_b)
			.extra_electrons(Cfg::N_ELECTRONS)
			.extra_states(Cfg::EXTRA_STATES)
			.temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
		input::kpoints::gamma());

	int wp_idx = -1;
	inqkit::InjectionReport report{};
	if (RESUME) {
		electrons.load(RT_CKPT_DIR);
		wp_idx = read_kv_int(RT_STATE_TXT, "wp_idx", electrons.states().num_states() - 1);
	} else {
		electrons.load(GS_DIR);
		jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OUT + "/raw/observables/eigenvalues");
		auto wp = inqkit::WavePacket{}
		              .center(0.0, 0.0, LAUNCH_Z)
		              .sigma(SIGMA_WP)
		              .k0(0.0, 0.0, K0)              // plain drifting Gaussian (no chirp; dx=0.40)
		              .orthogonalise_against_occupied(electrons);
		report = wp.inject_into_last_extra_state(electrons, 1.0);
		wp_idx = report.state_index;
	}
	const int n_states = electrons.states().num_states();
	electrons.inverse_mass()[0][wp_idx] = INV_MASS;      // the mass fork (re-applied on resume)

	if (electrons.root())
		std::cout << std::setprecision(8)
			<< "\n=== mass_pair_n162 WP out=" << OUT << (RESUME?"  [RESUME]":"  [FRESH]") << " ===\n"
			<< "  spacing=" << SPACING << " sigma_WP=" << SIGMA_WP << " k0=" << K0
			<< " inv_mass=" << INV_MASS << " (m=" << 1.0/INV_MASS << ")  v=" << K0*INV_MASS << "\n"
			<< "  launch_z=" << LAUNCH_Z << " dt=" << DT_AU << " start_step=" << START
			<< " n_steps=" << N_STEPS << " ckpt_every=" << CKPT_EVERY << " cap=" << (USE_CAP?"on":"off") << "\n"
			<< "  WP idx=" << wp_idx << " n_states=" << n_states
			<< (RESUME? "" : ("  norm_after=" + std::to_string(report.norm_after)
			   + "  max_overlap=" + std::to_string(report.max_overlap))) << "\n";

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

	// ----- output skeleton --------------------------------------------------
	std::filesystem::create_directories(OUT + "/raw/observables/overlap");
	std::filesystem::create_directories(OUT + "/raw/observables/overlap_full");
	for (auto sub : {"density_total","density_system","density_gs_system","density_wp",
	                 "wavefunction_wp","density_delta","density_delta_coarse"})
		std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
	std::filesystem::create_directories(RT_CKPT_DIR);

	inqkit::io::RealField3DLayout vti_layout{
		.field_name = "density", .include_meta = false, .emit_raw = false,
		.emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};

	if (!RESUME) { inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
	               gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system"); }

	inqkit::io::RealField3DWriter total_wr (OUT + "/raw/vti/density_total",  vti_layout, {.overwrite=!RESUME});
	inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite=!RESUME});
	inqkit::io::RealField3DWriter wp_density_wr(OUT + "/raw/vti/density_wp", vti_layout, {.overwrite=!RESUME});
	inqkit::io::ComplexField3DWriter wp_wf_wr(
		OUT + "/raw/vti/wavefunction_wp",
		{.field_name="wavefunction", .include_meta=false, .emit_raw=false,
		 .emit_vti=true, .vti_format=inqkit::io::VTIWriteOptions::Format::binary},
		{.overwrite=!RESUME});
	inqkit::observables::DensityDelta density_delta(
		OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
		{.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});
	if (!RESUME) { auto s0 = inqkit::fields::density::total(electrons);
		total_wr.write(s0,0.0,0); system_wr.write(s0,0.0,0); density_delta.snapshot(s0,0.0,0);
		wp_density_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0); }

	// ----- scalar observables -----------------------------------------------
	// (a) extensive obs.csv at WRITE_EVERY: energies + current + dipole + L2
	inqkit::io::ObservableSelection sel;
	sel.step = sel.time_au = true;
	sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
	sel.current_x = sel.current_y = sel.current_z = OBS_DIPCUR;
	sel.dipole_x = sel.dipole_y = sel.dipole_z = OBS_DIPCUR;
	sel.density_l2 = true;
	inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables" + SEG + ".csv", sel);
	obs_writer.write_header();

	// (b) PER-STEP full energy decomposition -> energy_decomp.csv (user spec:
	//     all energy components recorded EACH timestep). The 8 non-diagnostic
	//     terms (kinetic+hartree+xc+external+nonlocal+ion+ion_kinetic+exact_exch)
	//     sum to total; nvxc and eigenvalues are diagnostics.
	inqkit::io::ObservableSelection esel;
	esel.step = esel.time_au = true;
	esel.energy_total = esel.energy_kinetic = esel.energy_hartree = esel.energy_xc = true;
	esel.energy_external = esel.energy_nonlocal = esel.energy_ion = esel.energy_ion_kinetic = true;
	esel.energy_exact_exchange = esel.energy_nvxc = esel.energy_eigenvalues = true;
	esel.current_x = esel.current_y = esel.current_z = false;
	inqkit::io::ObservablesWriter energy_wr(OUT + "/raw/observables/energy_decomp" + SEG + ".csv", esel);
	energy_wr.write_header();

	inqkit::observables::StateEnergyWriter state_energy_wr(OUT + "/raw/observables/state_energies" + SEG + ".csv", true);
	inqkit::observables::OccupationsWriter occupations_wr(OUT + "/raw/observables/occupations_vs_time" + SEG + ".csv");

	// ----- WP-specific observables (overlap matrices guarded — they store a full
	//        GS-orbital reference on GPU; suspected step-0 OOM on the large grid) --
	std::unique_ptr<inqkit::observables::OrbitalOverlapMatrix> overlap_obs, overlap_full_obs;
	if (OBS_OVL) {
		overlap_obs = std::make_unique<inqkit::observables::OrbitalOverlapMatrix>(electrons, wp_idx, OUT + "/raw/observables/overlap");
		overlap_full_obs = std::make_unique<inqkit::observables::OrbitalOverlapMatrix>(electrons, n_states - 1, OUT + "/raw/observables/overlap_full");
		if (!RESUME) { overlap_full_obs->snapshot(electrons, 0.0, 0); overlap_obs->snapshot_wp_only(electrons, 0.0, 0); }
	}
	inqkit::observables::MomentumDistribution momentum_dist(
		OUT + "/raw/observables/momentum_distribution" + SEG + ".csv", wp_idx, Cfg::LZ_BOHR,
		{.n_bins=64, .k_max_bohr_inv=0.0, .write_every=WRITE_EVERY});
	inqkit::observables::WPMomentumStats wp_momentum_stats(
		OUT + "/raw/observables/wp_momentum_stats" + SEG + ".csv", wp_idx, {.write_every=WRITE_EVERY});
	inqkit::observables::WPRealSpaceStats wp_real_space_stats(
		OUT + "/raw/observables/wp_real_space_stats" + SEG + ".csv", wp_idx, {.write_every=WRITE_EVERY});

	std::ofstream nlog(OUT + "/raw/observables/electron_number" + SEG + ".csv");
	nlog << std::setprecision(12) << "step,time_au,N_total\n";

	// ----- density/observable session (WRITE_EVERY cadence) -----------------
	inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
	rt_obs.add([&](inqkit::StepContext const& ctx) {
		auto sys_f = inqkit::fields::density::total(*ctx.electrons);
		system_wr.write(sys_f, ctx.time_au, ctx.step);
		total_wr.write (sys_f, ctx.time_au, ctx.step);
		const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
		inqkit::StepContext c = ctx; c.density_l2 = l2; obs_writer.append(c);
		wp_density_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, ctx.step);
		if (OBS_WF && ctx.step % WF_EVERY == 0) {
			char nm[64]; std::snprintf(nm, sizeof(nm), "wavefunction_t%06d", ctx.step);
			wp_wf_wr.write(inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx), std::string(nm));
		}
	});

	// ----- propagate (LDA, ETRS default; native start_step restart) ---------
	auto step_fn = [&](auto const& data) {
		const int    it = data.iter();
		const double t  = data.time();
		// per-step FULL energy decomposition
		inqkit::StepContext ec; ec.step = it; ec.time_au = t;
		ec.energy_total = data.energy().total(); ec.energy_kinetic = data.energy().kinetic();
		ec.energy_hartree = data.energy().hartree(); ec.energy_xc = data.energy().xc();
		ec.energy_external = data.energy().external(); ec.energy_nonlocal = data.energy().non_local();
		ec.energy_ion = data.energy().ion(); ec.energy_ion_kinetic = data.energy().ion_kinetic();
		ec.energy_exact_exchange = data.energy().exact_exchange(); ec.energy_nvxc = data.energy().nvxc();
		ec.energy_eigenvalues = data.energy().eigenvalues();
		energy_wr.append(ec);
		// WRITE_EVERY-cadence heavy observables
		rt_obs.step(data);
		if (OBS_MOM) momentum_dist.maybe_accumulate(data);
		wp_momentum_stats.maybe_accumulate(data);
		wp_real_space_stats.maybe_accumulate(data);
		if (it % (5 * WRITE_EVERY) == 0) {
			if (OBS_STATE) { state_energy_wr.snapshot(data); occupations_wr.snapshot(data); }
			if (OBS_OVL) overlap_obs->snapshot_wp_only(electrons, t, it);
		}
		if (data.root()) nlog << it << "," << t << "," << data.num_electrons() << "\n";
		if (it > START && it % CKPT_EVERY == 0 && it < N_STEPS) {
			electrons.save(RT_CKPT_DIR);
			if (electrons.root()) {
				std::ofstream st(RT_STATE_TXT, std::ios::trunc);
				st << "last_step=" << it << "\ntime_au=" << t << "\nwp_idx=" << wp_idx
				   << "\nn_states=" << n_states << "\nn_steps_target=" << N_STEPS
				   << "\ndt_au=" << DT_AU << "\n";
				std::cout << "  [ckpt] saved at step " << it << " (t=" << t << ")\n";
			}
		}
	};
	auto rt_opts = OBS_DIPCUR
		? options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime).observables_current().observables_dipole()
		: options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
	real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, pert, START);

	if (OBS_OVL) overlap_full_obs->snapshot(electrons, DT_AU * N_STEPS, N_STEPS);

	// final checkpoint (so a follow-on resume is a clean no-op)
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
		  << "run = localised_jellium/mass_pair_n162/wp/" << env_s("EM_OUT","m1") << "\n"
		  << "engine = inq-study (mass fork)\nxc = LDA  propagator = ETRS  gamma-only\n"
		  << "resume = " << (RESUME?"true":"false") << "  start_step = " << START << "\n"
		  << "cell_bohr = 60.8x60.8x62.4  spacing = " << SPACING << "\n"
		  << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " edge_width " << Cfg::EDGE_WIDTH_BOHR << " axis " << Cfg::SLAB_AXIS << "\n"
		  << "n0_a0m3 = " << Cfg::N0 << "  r_s = " << inqkit::jellium::rs_from_n0(Cfg::N0) << "\n"
		  << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "  wp_state_index = " << wp_idx << "\n"
		  << "wp_norm_after = " << report.norm_after << "  wp_max_overlap = " << report.max_overlap << "\n"
		  << "projectile = quantum WP  sigma_WP = " << SIGMA_WP << "  k0 = " << K0
		  << "  inverse_mass = " << INV_MASS << "  m_eff = " << 1.0/INV_MASS << "\n"
		  << "velocity_au = " << K0*INV_MASS << "  E_eV = " << 0.5*(1.0/INV_MASS)*std::pow(K0*INV_MASS,2)*27.211386 << "\n"
		  << "launch_z = " << LAUNCH_Z << "  cap = " << (USE_CAP?"on":"off") << "\n"
		  << "cap_eta = " << CAP_ETA << "  cap_center_bohr = " << CAP_CENTER_B << "  cap_width_bohr = " << CAP_WIDTH_B
		  << "  cap_region = [" << (CAP_CENTER_B-CAP_WIDTH_B/2) << "," << (CAP_CENTER_B+CAP_WIDTH_B/2) << "]"
		  << "  slab_cap_gap = " << (CAP_CENTER_B-CAP_WIDTH_B/2-Cfg::SLAB_HALF_WIDTH) << "\n"
		  << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY
		  << "  wf_every = " << WF_EVERY << "  ckpt_every = " << CKPT_EVERY << "\n"
		  << "gs_dir = " << GS_DIR << "\nrt_ckpt_dir = " << RT_CKPT_DIR << "\nwall_time_s = " << wall << "\n"
		  << "run_completed = true\n";
		std::cout << "  WP run done (" << wall << " s), final step " << N_STEPS << ".\n";
	}
	return 0;
}
