// ============================================================================
// localised_jellium / scripts/fullsuite_classical / run.cpp
//
// FULL-SUITE classical-projectile re-run (Phase 5 headline, D2 2026-06-22).
// Mirrors the proven cap_baselines b2 path (Ehrenfest Gaussian-electron ion) +
// the localised background well + two-sided sin² CAP. Build ONCE against
// inq-study; reuse shared_gs/slab_n82_L50x50x101_h0p40 with the background re-applied.
//
// The projectile is a classical Gaussian electron (electron-Gaussian UPF,
// sigma_pot = sigma_WP/sqrt2 = 0.35; mass = m_e), launched +z. Its kinetic-energy
// loss across the slab faces gives the stopping power S = ΔKE_ion / x.
//
// Emits: density VTIs total/system/gs_system + delta(+coarse); observables
// (energies/current/dipole/L2), state_energies, occupations, eigenvalues,
// overlap_full (bath excitations); electron_track.csv (z,v,KE every step);
// electron_number.csv (bath norm -> over-drain check).
//
// Env: LJ_OUT(p5_classical) LJ_CAP(1) LJ_N_STEPS(900) LJ_WRITE_EVERY(10)
//      LJ_LAUNCH_Z(-15.5) LJ_DT(0.02).
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/orbital_overlap.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include "../../../shared/configs/slab_n82_L50x50x111.hpp"
#include "../../../../jellium/shared/cpp/eigenvalues_writer.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN82_L50x50x111;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static const char* PROJ_PSEUDO =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
    "shared/pseudopotentials/electron_gaussian_wpsigma3p5.upf";

int main() {
    auto t0 = std::chrono::steady_clock::now();
    constexpr double M_PROJ = 1.0 / 1822.8885;   // electron mass in amu

    const bool   USE_CAP     = env_i("LJ_CAP", 1) != 0;
    const std::string OUT    = "results/" + env_s("LJ_OUT", "p5_classical");
    const double DT_AU       = env_d("LJ_DT", 0.02);
    const int    N_STEPS     = env_i("LJ_N_STEPS", 900);
    const int    WRITE_EVERY = env_i("LJ_WRITE_EVERY", 10);
    const double LAUNCH_Z    = env_d("LJ_LAUNCH_Z", -26.5);
    const double V0          = env_d("LJ_K0", Cfg::WP_K0);
    // Matched to the WP run: two-sided sin^2 CAP, eta -1.0 Ha, 14 Bohr/side, region
    // [+/-41.5, +/-55.5]. mid = 48.5/111, width = 14.0/111.
    const double CAP_ETA = -1.0, CAP_MID = 48.5/111.0, CAP_WIDTH = 14.0/111.0;

    const std::string GS_DIR = env_s("LJ_GS_DIR",
        "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "shared_gs/slab_n82_L50x50x111_h0p40_per2");
    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

    std::cout << "\n=== fullsuite_classical (cap=" << (USE_CAP?"on":"off") << ", out=" << OUT << ") ===\n"
              << "  N_STEPS=" << N_STEPS << " dt=" << DT_AU << " launch_z=" << LAUNCH_Z << " v0=" << V0 << "\n";

    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b).periodicity(2);  // open-z: kills the point-charge z self-image drag
    auto ions = systems::ions(cell);
    auto sp = ionic::species("H").pseudo_file(PROJ_PSEUDO).mass(M_PROJ);
    ions.insert(sp, {0.0 * 1.0_b, 0.0 * 1.0_b, LAUNCH_Z * 1.0_b});

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);
    std::cout << "  Loaded GS (inq-study) from " << GS_DIR << "\n";
    jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OUT + "/raw/observables/eigenvalues");
    const int n_states = electrons.states().num_states();

    ions.velocities()[0] = vector3<double>{0.0, 0.0, V0};

    // ----- background well (+ optional CAP) ------------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::slab;
    bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
    bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);
    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID, CAP_WIDTH);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID, CAP_WIDTH);
    auto pert_bg  = bg_pert;
    auto pert_cap = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

    // ----- output skeleton + writers ------------------------------------
    for (auto sub : {"density_total","density_system","density_gs_system",
                     "density_delta","density_delta_coarse"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
    std::filesystem::create_directories(OUT + "/raw/observables/overlap_full");

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    { inqkit::io::RealField3DWriter gs_wr(OUT + "/raw/vti/density_gs_system", vti_layout, {.overwrite=true});
      gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system"); }
    inqkit::io::RealField3DWriter total_wr (OUT + "/raw/vti/density_total",  vti_layout, {.overwrite=true});
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite=true});
    { auto s0 = inqkit::fields::density::total(electrons); total_wr.write(s0,0.0,0); system_wr.write(s0,0.0,0); }

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
    obs_writer.write_header();
    inqkit::observables::StateEnergyWriter state_energy_wr(OUT + "/raw/observables/state_energies.csv", true);
    inqkit::observables::OccupationsWriter occupations_wr(OUT + "/raw/observables/occupations_vs_time.csv");
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});
    inqkit::observables::OrbitalOverlapMatrix overlap_full_obs(electrons, n_states - 1, OUT + "/raw/observables/overlap_full");
    overlap_full_obs.snapshot(electrons, 0.0, 0);

    // KE in atomic units uses the electron mass = 1 a.u. (NOT the amu value
    // M_PROJ passed to INQ's .mass(), which is 1/1822.8885 amu). Logging with
    // M_PROJ would under-scale KE by ~1822x.
    auto ke_ion = [&](){ auto const& v = ions.velocities()[0];
        return 0.5*1.0*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); };
    std::ofstream trk(OUT + "/raw/observables/electron_track.csv");
    trk << std::setprecision(16) << "step,time_au,x,y,z,vx,vy,vz,ke_ion_ha\n";
    { auto const& p = ions.positions()[0]; auto const& v = ions.velocities()[0];
      trk << 0 << ",0," << p[0] << "," << p[1] << "," << p[2] << "," << v[0] << "," << v[1] << "," << v[2] << "," << ke_ion() << "\n"; }
    std::ofstream nlog(OUT + "/raw/observables/electron_number.csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    // density/observable session (every WRITE_EVERY) + per-step ion track.
    inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        total_wr.write (sys_f, ctx.time_au, ctx.step);
        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext c = ctx; c.density_l2 = l2; obs_writer.append(c);
    });

    auto step_fn = [&](auto const& data) {
        rt_obs.step(data);
        const int it = data.iter();
        if (it % (5 * WRITE_EVERY) == 0) { state_energy_wr.snapshot(data); occupations_wr.snapshot(data); }
        if (data.root()) {
            nlog << it << "," << (it*DT_AU) << "," << data.num_electrons() << "\n";
            auto const& p = ions.positions()[0]; auto const& v = ions.velocities()[0];
            trk << it << "," << (it*DT_AU) << "," << p[0] << "," << p[1] << "," << p[2] << ","
                << v[0] << "," << v[1] << "," << v[2] << "," << ke_ion() << "\n";
        }
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime)
                       .observables_current().observables_dipole().ehrenfest();
    if (USE_CAP) real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, pert_cap);
    else         real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, pert_bg);

    overlap_full_obs.snapshot(electrons, DT_AU * N_STEPS, N_STEPS);

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/wide_wp_classical/" << env_s("LJ_OUT","p5_classical") << "\n"
          << "engine = inq-study\n"
          << "projectile = classical Gaussian-e ion (sigma_pot 2.475 = sigma_WP 3.5, mass m_e)\n"
          << "boundary = periodicity 2 (open-z; kills point-charge z self-image)\n"
          << "cap = " << (USE_CAP?"on (sin2 eta -1.0 Ha, 14 Bohr/side, region +/-41.5..+/-55.5)":"off") << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "  spacing = " << Cfg::SPACING_BOHR << "\n"
          << "background = slab half_width " << Cfg::SLAB_HALF_WIDTH << " axis " << Cfg::SLAB_AXIS << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n_states = " << n_states << "\n"
          << "launch_z = " << LAUNCH_Z << "  v0 = " << V0 << "  ke_ion_initial_ha = " << ke_ion() << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
          << "final_z = " << ions.positions()[0][2] << "  final_vz = " << ions.velocities()[0][2] << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
    }
    std::cout << "  done. final_z=" << ions.positions()[0][2] << " wall=" << wall << "s\n";
    return 0;
}
