// ============================================================================
// localised_jellium / muon_mass_fork / effmass_pair / quantum / run.cpp
//
// QUANTUM effective-mass WAVEPACKET projectile through the r_s=5.665 jellium
// slab. Contender D: sigma_WP=2, m_eff=3.085 m_e, v=2.711 a.u. (= 100 eV
// electron velocity -> same S(v)), E=309 eV, k0=8.364, dx=0.333, dt=0.04.
// The projectile is a Gaussian WP injected into the last extra state with a
// tuned inverse_mass (the inq-study mass fork). Bath = mass-1 electrons.
// Stopping read from the n(k,t) coherent peak (robust to WP spreading).
// Pairs 1:1 with ../classical (same slab, momentum, mass, launch, dt, CAP).
// Loads the shared GS from shared_gs/slab_n82_L50x50x90_dx0p333.
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

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN82_L50x50x90;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
	auto t0 = std::chrono::steady_clock::now();
	const std::string OUT     = "results/" + env_s("EM_OUT", "quantum");
	const double SPACING      = env_d("EM_SPACING", 0.33333);
	const double SIGMA_WP     = env_d("EM_SIGMA_WP", 2.0);
	const double K0           = env_d("EM_K0", 8.3641);
	const double INV_MASS     = env_d("EM_INV_MASS", 0.324127);   // 1/3.0852
	const double LAUNCH_Z     = env_d("EM_LAUNCH_Z", -16.743);
	const double DT_AU        = env_d("EM_DT", 0.04);
	const int    N_STEPS      = env_i("EM_N_STEPS", 2000);
	const int    WRITE_EVERY  = env_i("EM_WRITE_EVERY", 20);
	const bool   USE_CAP      = env_i("EM_CAP", 1) != 0;
	const std::string GS_DIR  = env_s("EM_GS_DIR",
		"/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
		"shared_gs/slab_n82_L50x50x90_dx0p333");
	const double CAP_ETA = -0.7, CAP_MID = 40.0/90.0, CAP_WIDTH = 10.0/90.0;

	if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

	auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b).periodic();
	auto ions = systems::ions(cell);                     // jellium: no nuclei
	auto electrons = systems::electrons(
		ions,
		options::electrons{}
			.spacing(SPACING * 1.0_b)
			.extra_electrons(Cfg::N_ELECTRONS)
			.extra_states(Cfg::EXTRA_STATES)
			.temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
		input::kpoints::gamma());
	electrons.load(GS_DIR);
	const int n_states = electrons.states().num_states();

	// ----- inject the effective-mass WP into the last extra state -----------
	auto wp = inqkit::WavePacket{}
	              .center(0.0, 0.0, LAUNCH_Z)
	              .sigma(SIGMA_WP)
	              .k0(0.0, 0.0, K0)
	              .orthogonalise_against_occupied(electrons);
	auto report = wp.inject_into_last_extra_state(electrons, 1.0);
	const int wp_idx = report.state_index;
	electrons.inverse_mass()[0][wp_idx] = INV_MASS;      // the mass fork

	if (electrons.root())
		std::cout << std::setprecision(8)
			<< "\n=== effmass_pair QUANTUM out=" << OUT << " ===\n"
			<< "  spacing=" << SPACING << " sigma_WP=" << SIGMA_WP << " k0=" << K0
			<< " inv_mass=" << INV_MASS << " (m=" << 1.0/INV_MASS << ")\n"
			<< "  launch_z=" << LAUNCH_Z << " dt=" << DT_AU << " n_steps=" << N_STEPS
			<< " cap=" << (USE_CAP?"on":"off") << "\n"
			<< "  WP idx=" << wp_idx << " norm_after=" << report.norm_after
			<< " n_states=" << n_states << "\n";

	// ----- background well (+ CAP) ------------------------------------------
	inqkit::jellium::localised_background_params bg;
	bg.shape = inqkit::jellium::background_shape::slab;
	bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
	bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
	inqkit::jellium::localised_background_perturbation bg_pert(bg);
	const double eta = USE_CAP ? CAP_ETA : 0.0;   // eta=0 -> CAP inert (keeps one pert type)
	perturbations::absorbing cap_lo(eta * 1.0_Ha, -CAP_MID, CAP_WIDTH);
	perturbations::absorbing cap_hi(eta * 1.0_Ha,  CAP_MID, CAP_WIDTH);
	auto pert = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

	// ----- observables -------------------------------------------------------
	std::filesystem::create_directories(OUT + "/raw/observables");
	for (auto sub : {"density_total","density_wp","density_gs_system","density_delta","density_delta_coarse"})
		std::filesystem::create_directories(OUT + "/raw/vti/" + sub);

	inqkit::io::ObservableSelection sel;
	sel.step = sel.time_au = true;
	sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
	sel.density_l2 = true;
	inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
	obs_writer.write_header();
	inqkit::observables::WPRealSpaceStats wp_rs(
		OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every = WRITE_EVERY});
	inqkit::observables::WPMomentumStats wp_ms(
		OUT + "/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every = WRITE_EVERY});

	inqkit::io::RealField3DLayout vti_layout{
		.field_name = "density", .include_meta = false, .emit_raw = false,
		.emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
	{ inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
	  gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system"); }
	inqkit::io::RealField3DWriter total_wr(OUT + "/raw/vti/density_total", vti_layout, {.overwrite=true});
	inqkit::io::RealField3DWriter wp_wr   (OUT + "/raw/vti/density_wp",    vti_layout, {.overwrite=true});
	inqkit::observables::DensityDelta density_delta(
		OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
		{.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});
	{ auto s0 = inqkit::fields::density::total(electrons);
	  total_wr.write(s0,0.0,0); density_delta.snapshot(s0,0.0,0);
	  wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0); }

	// ----- propagate (LDA, ETRS default; background+CAP perturbation) --------
	int g = 0;
	auto func = [&](auto const& data){
		const double t = g * DT_AU;
		wp_rs.maybe_accumulate(data);
		wp_ms.maybe_accumulate(data);
		if (g % WRITE_EVERY == 0) {
			auto sys_f = inqkit::fields::density::total(electrons);
			total_wr.write(sys_f, t, g);
			wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), t, g);
			const double l2 = density_delta.snapshot(sys_f, t, g);
			inqkit::StepContext c; c.step = g; c.time_au = t;
			c.energy_total = data.energy().total(); c.energy_kinetic = data.energy().kinetic();
			c.energy_hartree = data.energy().hartree(); c.energy_xc = data.energy().xc();
			c.density_l2 = l2;
			obs_writer.append(c);
		}
		++g;
	};
	auto opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
	real_time::propagate(ions, electrons, func, options::theory{}.lda(), opts, pert);

	const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
	if (electrons.root()) {
		std::ofstream s(OUT + "/run_summary.txt");
		s << std::setprecision(12)
		  << "run = localised_jellium/muon_mass_fork/effmass_pair/quantum\n"
		  << "engine = inq-study (mass fork)\nxc = LDA  propagator = ETRS  gamma-only\n"
		  << "cell_bohr = 50x50x90  spacing = " << SPACING << "\n"
		  << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS << "\n"
		  << "n0_a0m3 = " << Cfg::N0 << "  r_s = " << inqkit::jellium::rs_from_n0(Cfg::N0) << "\n"
		  << "projectile = quantum WP  sigma_WP = " << SIGMA_WP << "  k0 = " << K0
		  << "  inverse_mass = " << INV_MASS << "  m_eff = " << 1.0/INV_MASS << "\n"
		  << "velocity_au = " << K0*INV_MASS << "  E_eV = " << 0.5*(1.0/INV_MASS)*std::pow(K0*INV_MASS,2)*27.211386 << "\n"
		  << "launch_z = " << LAUNCH_Z << "  cap = " << (USE_CAP?"on":"off") << "\n"
		  << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
		  << "wp_idx = " << wp_idx << "  n_states = " << n_states << "\n"
		  << "gs_dir = " << GS_DIR << "\nwall_time_s = " << wall << "\n"
		  << "run_completed = true\n";
		std::cout << "  quantum run done (" << wall << " s).\n";
	}
	return 0;
}
