// ============================================================================
// cylindrical_jellium / scripts/annular_sv / classical / run.cpp
//
// Gliding classical-electron projectile down the bore of a PERIODIC annular
// jellium tube — the production S(v) run. Loads a bare-tube GS (CJ_GS_DIR),
// inserts a classical Gaussian electron (electron_gaussian_wpsigma0p5.upf,
// sigma_pot=0.354, mass m_e) on-axis near the −z face, launches it +z at PROJ_V0,
// and propagates with EHRENFEST ion dynamics + the annulus background re-applied.
//
// S is extracted two ways downstream: (a) ΔE_system(t) regression vs path
// (PRIMARY), (b) ΔKE_ion/Δz from the track. Emits: electron_track.csv (every
// step: z,v,KE_ion), observables.csv (energies + current_z = wall current, every
// WRITE_EVERY), density_system VTI + induced density_delta(raw+coarse) at the
// 300-frame cadence, electron_number.csv (bath-norm conservation).
//
// Env (orchestrator supplies all):
//   geometry: CJ_LXY(40) CJ_LZ(48) CJ_RIN(5) CJ_ROUT(13) CJ_N(24) CJ_EDGE_W(1.0)
//             CJ_SPACING(0.5)
//   run:      CJ_GS_DIR(REQUIRED) PROJ_V0(0.30) SV_N_STEPS(2000)
//             SV_WRITE_EVERY(=N_STEPS/300) SV_OUT_SUBDIR(v_run)
//             CJ_LAUNCH_Z(=-(L_z/2)+1)
// Build against INQ_SOURCE=inq-study; runtime shares from inq/install/share.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static const char* PROJ_PSEUDO =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/cylindrical_jellium/"
    "shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main() {
    auto t0 = std::chrono::steady_clock::now();
    constexpr double M_PROJ = 1.0 / 1822.8885;   // electron mass in amu (INQ wants amu)

    const double LXY = env_d("CJ_LXY", 40.0), LZ = env_d("CJ_LZ", 48.0);
    const double RIN = env_d("CJ_RIN", 5.0), ROUT = env_d("CJ_ROUT", 13.0);
    const int    N = env_i("CJ_N", 24);
    const double EDGE_W = env_d("CJ_EDGE_W", 1.0);
    const double SPACING = env_d("CJ_SPACING", 0.5);
    const double V0 = env_d("PROJ_V0", 0.30);
    const int    N_STEPS = env_i("SV_N_STEPS", 2000);
    const int    WRITE_EVERY = env_i("SV_WRITE_EVERY", std::max(1, (int)std::lround(N_STEPS/300.0)));
    const double DT_AU = 0.020;
    const double LAUNCH_Z = env_d("CJ_LAUNCH_Z", -(LZ/2.0) + 1.0);
    const std::string OUT = "results/" + env_s("SV_OUT_SUBDIR", "v_run");
    const std::string GS_DIR = env_s("CJ_GS_DIR", "");
    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2;
    }

    const double V_ann = M_PI * (ROUT*ROUT - RIN*RIN) * LZ;
    const double N0 = double(N) / V_ann;

    std::cout << "\n=== annular_sv classical OUT=" << OUT << " ===\n"
              << "  tube R_in=" << RIN << " R_out=" << ROUT << " L_z=" << LZ << " (periodic)\n"
              << "  N=" << N << " n0=" << N0 << " v0=" << V0 << " launch_z=" << LAUNCH_Z << "\n"
              << "  N_STEPS=" << N_STEPS << " dt=" << DT_AU << " write_every=" << WRITE_EVERY << "\n";

    auto cell = systems::cell::orthorhombic(LXY*1.0_b, LXY*1.0_b, LZ*1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto sp = ionic::species("H").pseudo_file(PROJ_PSEUDO).mass(M_PROJ);
    ions.insert(sp, {0.0*1.0_b, 0.0*1.0_b, LAUNCH_Z*1.0_b});

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING*1.0_b)
            .extra_electrons(N)
            .extra_states(20)
            .temperature(0.00862*1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);
    const int n_states = electrons.states().num_states();

    ions.velocities()[0] = vector3<double>{0.0, 0.0, V0};

    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::annulus;
    bg.n0 = N0; bg.half_width = ROUT; bg.inner_radius = RIN; bg.slab_axis = 2;
    bg.center = {0.0, 0.0, 0.0}; bg.edge_width = EDGE_W;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    // ----- output skeleton + writers ------------------------------------
    for (auto sub : {"density_system","density_delta","density_delta_coarse"})
        std::filesystem::create_directories(OUT + "/raw/vti/" + sub);
    std::filesystem::create_directories(OUT + "/raw/observables");

    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite=true});
    { auto s0 = inqkit::fields::density::total(electrons); system_wr.write(s0, 0.0, 0); }

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;   // current_z = wall current
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
    obs_writer.write_header();
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti=true, .emit_coarse_vti=true, .compute_l2=true, .coarse_bin_bohr=3.0});

    // KE of the ion in atomic units uses electron mass = 1.0 a.u. (NOT M_PROJ amu,
    // which would under-scale KE by ~1822x — see fullsuite_classical note).
    auto ke_ion = [&](){ auto const& v = ions.velocities()[0];
        return 0.5*1.0*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); };
    std::ofstream trk(OUT + "/raw/observables/electron_track.csv");
    trk << std::setprecision(16) << "step,time_au,x,y,z,vx,vy,vz,ke_ion_ha\n";
    { auto const& p = ions.positions()[0]; auto const& v = ions.velocities()[0];
      trk << 0 << ",0," << p[0] << "," << p[1] << "," << p[2] << "," << v[0] << "," << v[1] << "," << v[2] << "," << ke_ion() << "\n"; }
    std::ofstream nlog(OUT + "/raw/observables/electron_number.csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    // density/observable session (every WRITE_EVERY): density VTI + induced delta + energies.
    inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext c = ctx; c.density_l2 = l2; obs_writer.append(c);
    });

    auto step_fn = [&](auto const& data) {
        rt_obs.step(data);
        const int it = data.iter();
        if (data.root()) {
            nlog << it << "," << (it*DT_AU) << "," << data.num_electrons() << "\n";
            auto const& p = ions.positions()[0]; auto const& v = ions.velocities()[0];
            trk << it << "," << (it*DT_AU) << "," << p[0] << "," << p[1] << "," << p[2] << ","
                << v[0] << "," << v[1] << "," << v[2] << "," << ke_ion() << "\n";
        }
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU*1.0_atomictime)
                       .observables_current().observables_dipole().ehrenfest();
    real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, bg_pert);

    trk.flush();
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = annular_sv/classical/" << env_s("SV_OUT_SUBDIR","v_run") << "\nengine = inq-study\n"
          << "projectile = classical Gaussian-e ion (sigma_pot 0.354, mass m_e, ehrenfest)\n"
          << "geometry = annular_tube  R_in=" << RIN << " R_out=" << ROUT << " L_z=" << LZ << " (periodic)\n"
          << "cell_bohr = " << LXY << " x " << LXY << " x " << LZ << "  spacing = " << SPACING << "\n"
          << "n_electrons = " << N << "  n0 = " << N0 << "  r_s_eff = " << std::cbrt(3.0/(4.0*M_PI*N0)) << "\n"
          << "n_states = " << n_states << "\n"
          << "launch_z = " << LAUNCH_Z << "  v0 = " << V0 << "  ke_ion_initial_ha = " << (0.5*V0*V0) << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
          << "gs_dir = " << GS_DIR << "\n"
          << "final_z = " << ions.positions()[0][2] << "  final_vz = " << ions.velocities()[0][2] << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
    }
    std::cout << "  done. final_z=" << ions.positions()[0][2]
              << " final_vz=" << ions.velocities()[0][2] << " wall=" << wall << "s\n";
    return 0;
}
