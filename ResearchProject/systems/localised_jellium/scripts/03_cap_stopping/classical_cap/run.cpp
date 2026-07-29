// ============================================================================
// localised_jellium / 03_cap_stopping / classical_cap / run.cpp
//
// Phase 5 (CLASSICAL twin): a classical Gaussian-electron projectile fired
// through the slab WITH the two-sided sin² CAP, on inq-study. The projectile is
// an Ehrenfest ion (Gaussian electron pseudopotential, mass = 1 a.u.); the CAP
// absorbs the bath wake. The ion's kinetic-energy loss ΔKE_ion gives the CLEAN
// stopping power S = ΔKE_ion / x (x = 25 Bohr) — the cross-check to the WP run's
// bath-energy estimate. σ-matching: σ_pot=0.35 = σ_WP/√2 (density std matched).
//
// Geometry matches the WP+CAP run: launch z=−15.5, CAP eta=−0.5, 7.5 Bohr/side.
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
	const std::string PROJ_UPF =
		"/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/shared/"
		"pseudopotentials/electron_gaussian_sigma0p35.upf";

	constexpr double DT_AU       = 0.02;
	constexpr int    N_STEPS     = 900;
	constexpr int    WRITE_EVERY = 10;
	constexpr double LAUNCH_Z    = -15.5;
	constexpr double PROJ_MASS_AMU = 1.0 / 1822.8885;   // electron mass (1 a.u.)
	const     double VEL_Z       = Cfg::WP_K0;          // v=sqrt(2E/m)=2.71
	constexpr double CAP_ETA_HA  = -0.5;
	constexpr double CAP_MID     = 0.425;
	constexpr double CAP_WIDTH   = 0.15;

	auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();
	auto ions = systems::ions(cell);
	auto sp_e = ionic::species("H").pseudo_file(PROJ_UPF).mass(PROJ_MASS_AMU);
	ions.insert(sp_e, {0.0_b, 0.0_b, LAUNCH_Z * 1.0_b});
	std::cout << "Projectile ion mass_au = " << ions.species(0).mass() << "\n";

	auto electrons = systems::electrons(
		ions,
		options::electrons{}
			.spacing(Cfg::SPACING_BOHR * 1.0_b)
			.extra_electrons(Cfg::N_ELECTRONS)
			.extra_states(Cfg::EXTRA_STATES)
			.temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
		input::kpoints::gamma());
	electrons.load(GS_DIR);

	ions.velocities()[0] = vector3<double>{0.0, 0.0, VEL_Z};
	std::cout << "KE_ion(0) = " << ions.kinetic_energy() << " Ha (expect ~3.675)\n";

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
	std::ofstream log("results/classical_trace.csv");
	log << std::setprecision(12)
	    << "time_au,total_ha,num_electrons,ion_z,ion_vz,ke_ion_ha\n";
	int frame = 0;

	real_time::propagate(
		ions, electrons,
		[&](auto const & data) {
			if (data.root()) {
				auto z  = data.positions()[0][2];
				auto vz = data.velocities()[0][2];
				log << data.time() << "," << data.energy().total() << ","
				    << data.num_electrons() << "," << z << "," << vz << ","
				    << 0.5 * (1.0) * vz * vz << "\n";   // m=1 a.u.
			}
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
			.ehrenfest()
			.observables_current()
			.observables_dipole(),
		pert);

	if (electrons.root()) {
		std::ofstream s("results/run_summary.txt");
		s << std::setprecision(12)
		  << "run = localised_jellium/03_cap_stopping/classical_cap\n"
		  << "engine = inq-study\n"
		  << "projectile = classical Gaussian-e ion (sigma_pot=0.35), mass_au "
		  << ions.species(0).mass() << "\n"
		  << "cap = sin2 eta " << CAP_ETA_HA << " width " << CAP_WIDTH
		  << " mid +/-" << CAP_MID << "\n"
		  << "launch_z = " << LAUNCH_Z << " vel_z = " << VEL_Z << "\n"
		  << "ke_ion_initial_ha = " << (0.5 * VEL_Z * VEL_Z) << "\n"
		  << "n_steps = " << N_STEPS << " density_frames = " << frame << "\n";
	}
	std::cout << "Phase-5 classical+CAP run done: " << frame << " frames.\n";
	return 0;
}
