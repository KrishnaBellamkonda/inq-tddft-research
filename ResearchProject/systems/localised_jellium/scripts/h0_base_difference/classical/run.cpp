// ============================================================================
// localised_jellium / scripts/h0_base_difference / classical / run.cpp
//
// H0 (base WP-vs-classical E_total(0) gap), classical half. Loads the 50x50x120
// slab GS, inserts a STATIONARY matched classical Gaussian-charge ghost
// (electron_gaussian_wpsigma0p5.upf: z_valence 0, charge std sigma_WP/sqrt2 =
// 0.354, the SAME cloud as the sigma_WP = 0.5 WP) at z = LJ_LAUNCH_Z with
// velocity 0, re-applies the background well, and propagates a few short steps to
// read E_total(t=0) from the standard observables writer (records step 0).
//
// NOTE: the chargeless ghost contributes the BARE, unscreened int v_ghost*n_GS to
// E_total; the compensating int v_ghost*n_+ is OMITTED here (z_valence 0) and must
// be re-added in the H0 notebook for a fair WP-vs-classical comparison.
//
// Env: LJ_OUT(h0_cl_r4) LJ_LAUNCH_Z(-16.5) LJ_N_STEPS(2) LJ_DT(0.01).
// Build against INQ_SOURCE=inq-study; runtime shares from inq/install/share.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>

#include "../../../shared/configs/slab_n82_L50x50x120.hpp"

#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN82_L50x50x120;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

static const char* PROJ_PSEUDO =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
    "shared/pseudopotentials/electron_gaussian_wpsigma0p5.upf";

int main() {
    auto t0 = std::chrono::steady_clock::now();
    constexpr double M_PROJ = 1.0 / 1822.8885;   // electron mass in amu (irrelevant at v=0)

    const std::string OUT    = "results/" + env_s("LJ_OUT", "h0_cl_r4");
    const double DT_AU       = env_d("LJ_DT", 0.01);
    const int    N_STEPS     = env_i("LJ_N_STEPS", 2);
    const double LAUNCH_Z    = env_d("LJ_LAUNCH_Z", -16.5);

    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/localised_jellium/"
        "shared_gs/slab_n82_L50x50x120";
    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing: " << GS_DIR << "\n"; return 2; }

    std::cout << "\n=== H0 classical (out=" << OUT << ") launch_z=" << LAUNCH_Z
              << " ghost=wpsigma0p5 N_STEPS=" << N_STEPS << " dt=" << DT_AU << " ===\n";

    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b, Cfg::LY_BOHR * 1.0_b, Cfg::LZ_BOHR * 1.0_b).periodic();
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
    std::cout << "  Loaded GS from " << GS_DIR << "\n";

    ions.velocities()[0] = vector3<double>{0.0, 0.0, 0.0};   // stationary

    // ----- background well (re-applied during propagation) ---------------
    inqkit::jellium::localised_background_params bg;
    bg.shape = inqkit::jellium::background_shape::slab;
    bg.n0 = Cfg::N0; bg.half_width = Cfg::SLAB_HALF_WIDTH; bg.slab_axis = Cfg::SLAB_AXIS;
    bg.center = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR}; bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    // ----- observables (records step 0 => E_total(0)) --------------------
    std::filesystem::create_directories(OUT + "/raw/observables");
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.density_l2 = false;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
    obs_writer.write_header();

    std::ofstream nlog(OUT + "/raw/observables/electron_number.csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    inqkit::RealTimeSession rt_obs(ions, electrons, 1);
    rt_obs.add([&](inqkit::StepContext const& ctx) { obs_writer.append(ctx); });

    auto step_fn = [&](auto const& data) {
        rt_obs.step(data);
        const int it = data.iter();
        if (data.root()) nlog << it << "," << (it*DT_AU) << "," << data.num_electrons() << "\n";
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
    real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, bg_pert);

    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/h0_base_difference/classical/" << env_s("LJ_OUT","h0_cl_r4") << "\n"
          << "engine = inq-study\n"
          << "projectile = classical Gaussian ghost wpsigma0p5 (z_valence 0, charge std 0.354), stationary\n"
          << "launch_z = " << LAUNCH_Z << "  (r_from_face = " << (-LAUNCH_Z - Cfg::SLAB_HALF_WIDTH) << " Bohr)\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "  spacing = " << Cfg::SPACING_BOHR << "\n"
          << "ghost_background_term_omitted = true  (re-add int v_ghost*n_+ in analysis)\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
    }
    std::cout << "  done. wall=" << wall << "s\n";
    return 0;
}
