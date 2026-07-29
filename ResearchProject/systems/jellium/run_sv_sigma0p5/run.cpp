// run_sv_sigma0p5/run.cpp — minimal free-Ehrenfest S(v) run.
//
// One classical sigma=0.5 erf-smoothed electron (m_e), launched +z on-axis,
// decelerating through the r_s=5.69 jellium bath. Minimal observables:
// electron_track.csv (pos/vel every step -> KE_proj and instantaneous v) and
// observables.csv (energy components -> dE_electrons cross-check). The heavy
// density VTI / momentum / overlap writers are intentionally omitted: at
// ~14 s/step the propagation dominates, and S(v) needs only the track + energy.
//
// Runtime env (one build serves the whole ladder):
//   PROJ_V0         initial +z velocity in a.u.            (default 1.0)
//   SV_N_STEPS      number of propagation steps            (default 400)
//   SV_WRITE_EVERY  energy-observable cadence              (default 50)
//   SV_OUT_SUBDIR   results/<subdir>/ output location      (default v_run)
//   SV_PSEUDO       override psp path (sigma=0.4 sibling)  (default sigma0p5)
// ----------------------------------------------------------------------------
#include <inq/inq.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>

#include "../shared/configs/sv_ladder_L50_sigma0p5.hpp"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = jellium::config::SV_Ladder_L50_sigma0p5;

static double env_d(const char* k, double d) {
    const char* v = std::getenv(k); return v ? std::atof(v) : d;
}
static int env_i(const char* k, int d) {
    const char* v = std::getenv(k); return v ? std::atoi(v) : d;
}
static std::string env_s(const char* k, const std::string& d) {
    const char* v = std::getenv(k); return v ? std::string(v) : d;
}

int main() {
    const double V0          = env_d("PROJ_V0", Cfg::PROJ_VEL_Z_DEFAULT);
    const int    N_STEPS     = env_i("SV_N_STEPS", Cfg::N_STEPS_DEFAULT);
    const int    WRITE_EVERY = env_i("SV_WRITE_EVERY", Cfg::WRITE_EVERY_DEFAULT);
    const std::string OUT    = "results/" + env_s("SV_OUT_SUBDIR", "v_run");
    const std::string PSEUDO = env_s("SV_PSEUDO", Cfg::PROJ_PSEUDO_PATH);

    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p40";

    std::cout << "\n=== run_sv_sigma0p5 ===\n"
              << "  v0=" << V0 << " a.u.  N_STEPS=" << N_STEPS
              << "  dt=" << Cfg::DT_AU << "  out=" << OUT << "\n"
              << "  psp=" << PSEUDO << "\n"
              << "  launch_z=" << Cfg::PROJ_LAUNCH_Z << " Bohr\n";

    if (!std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: GS checkpoint missing\n"; return 2;
    }
    std::filesystem::create_directories(OUT);

    auto cell = systems::cell::cubic(Cfg::L_BOHR * 1.0_b).periodic();

    auto ions = systems::ions(cell);
    auto sp_e = ionic::species(Cfg::PROJ_SPECIES_SYMBOL)
                    .pseudo_file(PSEUDO)
                    .mass(Cfg::PROJ_MASS_AMU);              // m_e
    ions.insert(sp_e, {Cfg::PROJ_LAUNCH_X * 1.0_b,
                       Cfg::PROJ_LAUNCH_Y * 1.0_b,
                       Cfg::PROJ_LAUNCH_Z * 1.0_b});

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);

    ions.velocities()[0] = vector3<double>{0.0, 0.0, V0};
    std::cout << "  mass_au=" << ions.species(0).mass()
              << " (expect 1.0)  KE0=" << ions.kinetic_energy() << " Ha\n";

    // stub run_summary
    {
        std::ofstream s(OUT + "/run_summary.txt");
        s << "run_name       = run_sv_sigma0p5\n"
          << "projectile     = classical electron, erf sigma=0.5, m_e, ehrenfest\n"
          << "rs             = 5.69 (N=162, L=50, dx=0.40)\n"
          << "v0_au          = " << V0 << "\n"
          << "dt_au          = " << Cfg::DT_AU << "\n"
          << "n_steps        = " << N_STEPS << "\n"
          << "psp            = " << PSEUDO << "\n"
          << "run_completed  = false\n";
    }

    // electron_track.csv (every step): step,time,pos,vel  -> KE_proj, v(t)
    std::ofstream trk(OUT + "/electron_track.csv");
    trk << std::setprecision(16)
        << "step,time_au,x,y,z,vx,vy,vz\n";
    {
        auto const& p = ions.positions()[0];
        auto const& v = ions.velocities()[0];
        trk << 0 << ",0," << p[0] << "," << p[1] << "," << p[2] << ","
            << v[0] << "," << v[1] << "," << v[2] << "\n";
    }

    // observables.csv (energy / current / dipole) every WRITE_EVERY
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    inqkit::io::ObservablesWriter obs(OUT + "/observables.csv", sel);
    obs.write_header();

    inqkit::RealTimeSession rt_track(ions, electrons, 1);
    rt_track.add([&](inqkit::StepContext const& ctx) {
        auto const& p = ctx.ions->positions()[0];
        auto const& v = ctx.ions->velocities()[0];
        trk << ctx.step << "," << ctx.time_au << ","
            << p[0] << "," << p[1] << "," << p[2] << ","
            << v[0] << "," << v[1] << "," << v[2] << "\n";
    });

    inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const& ctx) { obs.append(ctx); });

    real_time::propagate(
        ions, electrons,
        [&](auto const& data) { rt_track.step(data); rt_obs.step(data); },
        options::theory{}.lda(),
        options::real_time{}
            .num_steps(N_STEPS)
            .dt(Cfg::DT_AU * 1.0_atomictime)
            .ehrenfest()
            .observables_current()
            .observables_dipole());

    trk.flush();
    {
        std::ofstream s(OUT + "/run_summary.txt", std::ios::app);
        s << "run_completed  = true\n"
          << "final_z        = " << ions.positions()[0][2] << "\n"
          << "final_vz       = " << ions.velocities()[0][2] << "\n";
    }
    std::cout << "  done. final_z=" << ions.positions()[0][2]
              << " final_vz=" << ions.velocities()[0][2] << "\n";
    return 0;
}
