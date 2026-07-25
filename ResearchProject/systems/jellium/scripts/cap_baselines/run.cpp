// scripts/cap_baselines/run.cpp — CAP-in-jellium baselines (b1/b2/b3, full suite).
//
// First use of the inq-study built-in CAP (complexified scalar potential) in an
// INTERACTING (LDA) jellium bath. Build ONCE against inq-study.
//
//   CAP_MODE = b1 : CAP on, NO projectile (bath-drainage reference)
//              b2 : CAP + classical sigma=0.5 e- (electron-Gaussian UPF, m_e), 100 eV
//              b3 : CAP + sigma=0.5 Gaussian WP, 100 eV
//
// Geometry (grill 2026-06-17): 50 Bohr cubic, two-sided sin^2 CAP 10 Bohr/side,
// eta=-0.5 -> slabs |z| in [15,25]; free [-15,+15]. Fractional absorbing(eta,+-0.4,0.2)
// via sum. ETRS. Pure-jellium GS reused. Projectile/WP launched z0=-13 (4 sigma
// inside the -z CAP edge), +z, exits through the far CAP.
//
// Env: CAP_MODE(b1) CAP_N_STEPS(7000) CAP_WRITE_EVERY(23) CAP_OUT_SUBDIR
//      CAP_ETA(-0.5) CAP_WIDTH_BOHR(10) CAP_LAUNCH_Z(-13) CAP_V0(2.7111).
// PROVISIONAL until Task #7.
// ----------------------------------------------------------------------------
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/occupations_writer.hpp>
#include <inqkit/observables/state_energy_writer.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include "../../shared/configs/electron_proj_E100_L50_cubic.hpp"
#include "../../shared/cpp/eigenvalues_writer.hpp"

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
using Cfg = jellium::config::Common_E100_L50_cubic;

static double env_d(const char* k, double d) { const char* v = std::getenv(k); return v ? std::atof(v) : d; }
static int env_i(const char* k, int d) { const char* v = std::getenv(k); return v ? std::atoi(v) : d; }
static std::string env_s(const char* k, const std::string& d) { const char* v = std::getenv(k); return v ? std::string(v) : d; }

static const char* DEFAULT_PROJ_PSEUDO =
    "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
    "shared/pseudopotentials/electron_gaussian_sigma0p5.upf";

int main() {
    auto t0 = std::chrono::steady_clock::now();

    const std::string MODE        = env_s("CAP_MODE", "b1");
    const int         N_STEPS     = env_i("CAP_N_STEPS", 7000);
    const int         WRITE_EVERY = env_i("CAP_WRITE_EVERY", 23);
    const std::string OUT         = "results/" + env_s("CAP_OUT_SUBDIR", MODE);
    const double      ETA         = env_d("CAP_ETA", -0.5);
    const double      CAP_W_BOHR  = env_d("CAP_WIDTH_BOHR", 10.0);
    const double      LAUNCH_Z    = env_d("CAP_LAUNCH_Z", -13.0);
    const double      V0          = env_d("CAP_V0", std::sqrt(2.0 * 100.0 / 27.211386));

    const bool is_b2 = (MODE == "b2"), is_b3 = (MODE == "b3"), is_b1 = (MODE == "b1");
    if (!is_b1 && !is_b2 && !is_b3) { std::cerr << "FATAL: CAP_MODE must be b1|b2|b3\n"; return 3; }

    const double L          = Cfg::L_BOHR;
    const double width_frac = CAP_W_BOHR / L;                 // 0.2
    const double mid_frac   = (L / 2 - CAP_W_BOHR / 2) / L;   // 0.4
    const double WP_SIGMA   = env_d("CAP_WP_SIGMA", 0.5), WP_K0 = V0;
    const std::string PROJ_PSEUDO = env_s("CAP_PROJ_PSEUDO", DEFAULT_PROJ_PSEUDO);

    const std::string GS_DIR =
        "/local/data/public/skcb2/tddft/ResearchProject/systems/jellium/"
        "checkpoints/gs_L50_cubic_N162_dx0p40";

    std::cout << "\n=== cap_baselines (" << MODE << ") ===\n"
              << "  CAP eta=" << ETA << " width/side=" << CAP_W_BOHR
              << " (mid_frac=" << mid_frac << " width_frac=" << width_frac << ")\n"
              << "  launch_z=" << LAUNCH_Z << " v0=" << V0
              << "  N_STEPS=" << N_STEPS << " WRITE_EVERY=" << WRITE_EVERY << "\n";

    if (!std::filesystem::exists(GS_DIR)) { std::cerr << "FATAL: GS missing\n"; return 2; }
    std::filesystem::create_directories(OUT + "/raw/observables");
    std::filesystem::create_directories(OUT + "/raw/vti/density_system");
    std::filesystem::create_directories(OUT + "/raw/vti/density_delta");
    std::filesystem::create_directories(OUT + "/raw/vti/density_delta_coarse");

    auto cell = systems::cell::cubic(L * 1.0_b).periodic();
    auto ions = systems::ions(cell);

    // b2: insert the classical sigma=0.5 electron projectile (electron-Gaussian
    // UPF, mass overridden to m_e). Same proven pattern as run_classical/run_sv.
    if (is_b2) {
        auto sp = ionic::species("H").pseudo_file(PROJ_PSEUDO).mass(1.0 / 1822.8885);
        ions.insert(sp, {0.0 * 1.0_b, 0.0 * 1.0_b, LAUNCH_Z * 1.0_b});
    }

    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(Cfg::SPACING_BOHR * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    electrons.load(GS_DIR);
    jellium::eigenvalues::copy_from_checkpoint(GS_DIR, OUT + "/raw/observables/eigenvalues");

    if (is_b2) ions.velocities()[0] = vector3<double>{0.0, 0.0, V0};

    // b3: inject the sigma=0.5 Gaussian WP into the last extra state.
    int wp_idx = -1;
    if (is_b3) {
        auto wp = inqkit::WavePacket{}
                      .center(0.0, 0.0, LAUNCH_Z).sigma(WP_SIGMA).k0(0.0, 0.0, WP_K0)
                      .orthogonalise_against_occupied(electrons);
        auto rep = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = rep.state_index;
        std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << rep.norm_after
                  << " max_overlap=" << rep.max_overlap << "\n";
    }

    auto electron_number = [&electrons]() -> double {
        auto rho = observables::density::calculate(electrons);
        return operations::integral_partial_sum(rho, std::min(2, (int)rho.set_size()));
    };
    const double N0 = electron_number();
    std::cout << "  N(t=0) = " << std::setprecision(10) << N0 << "\n";

    auto cap = perturbations::absorbing(ETA * 1.0_Ha,  mid_frac, width_frac)
             + perturbations::absorbing(ETA * 1.0_Ha, -mid_frac, width_frac);

    // ----- writers (shared suite) ---------------------------------------
    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter system_wr(OUT + "/raw/vti/density_system", vti_layout, {.overwrite = true});
    { auto sys0 = inqkit::fields::density::total(electrons); system_wr.write(sys0, 0.0, 0); }

    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.current_x = sel.current_y = sel.current_z = true;
    sel.dipole_x = sel.dipole_y = sel.dipole_z = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
    obs_writer.write_header();

    inqkit::observables::StateEnergyWriter state_energy_wr(
        OUT + "/raw/observables/state_energies.csv", true);
    inqkit::observables::OccupationsWriter occupations_wr(
        OUT + "/raw/observables/occupations_vs_time.csv");
    inqkit::observables::DensityDelta density_delta(
        OUT + "/raw/vti/density_delta", OUT + "/raw/vti/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true, .compute_l2 = true, .coarse_bin_bohr = 3.0});

    // b3: total-WP momentum distribution (proven writer).
    std::unique_ptr<inqkit::observables::MomentumDistribution> momentum_dist;
    if (is_b3) momentum_dist = std::make_unique<inqkit::observables::MomentumDistribution>(
        OUT + "/raw/observables/momentum_distribution.csv", wp_idx, Cfg::L_BOHR,
        inqkit::observables::MomentumDistributionConfig{.n_bins = 64, .k_max_bohr_inv = 0.0,
                                                        .write_every = WRITE_EVERY});

    // b2: classical projectile track.
    std::ofstream trk;
    if (is_b2) {
        trk.open(OUT + "/raw/observables/electron_track.csv");
        trk << std::setprecision(16) << "step,time_au,x,y,z,vx,vy,vz\n";
        auto const& p = ions.positions()[0]; auto const& v = ions.velocities()[0];
        trk << 0 << ",0," << p[0] << "," << p[1] << "," << p[2] << ","
            << v[0] << "," << v[1] << "," << v[2] << "\n";
    }

    std::ofstream nlog(OUT + "/raw/observables/electron_number.csv");
    nlog << std::setprecision(12) << "step,time_au,N_total,absorbed_frac\n";
    nlog << 0 << ",0," << N0 << ",0\n";

    bool nan_seen = false;

    inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const& ctx) {
        auto sys_f = inqkit::fields::density::total(*ctx.electrons);
        system_wr.write(sys_f, ctx.time_au, ctx.step);
        const double l2 = density_delta.snapshot(sys_f, ctx.time_au, ctx.step);
        inqkit::StepContext c = ctx; c.density_l2 = l2;
        obs_writer.append(c);
    });

    // b2: per-step projectile track (separate session, every step).
    inqkit::RealTimeSession rt_trk(ions, electrons, 1);
    if (is_b2) rt_trk.add([&](inqkit::StepContext const& ctx) {
        auto const& p = ctx.ions->positions()[0]; auto const& v = ctx.ions->velocities()[0];
        trk << ctx.step << "," << ctx.time_au << "," << p[0] << "," << p[1] << "," << p[2]
            << "," << v[0] << "," << v[1] << "," << v[2] << "\n";
    });

    auto step_fn = [&](auto const& data) {
        rt_obs.step(data);
        if (is_b2) rt_trk.step(data);
        const int it = data.iter();
        if (it % WRITE_EVERY == 0) { state_energy_wr.snapshot(data); occupations_wr.snapshot(data); }
        if (is_b3 && momentum_dist) momentum_dist->maybe_accumulate(data);
        const double N = electron_number();
        if (!std::isfinite(N)) nan_seen = true;
        nlog << it << "," << (it * Cfg::DT_AU) << "," << N << ","
             << (N0 > 0 ? (N0 - N) / N0 : 0.0) << "\n";
    };

    auto rt_opts = options::real_time{}
                       .num_steps(N_STEPS)
                       .dt(Cfg::DT_AU * 1.0_atomictime)
                       .observables_current()
                       .observables_dipole();
    if (is_b2) {
        real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(),
                             rt_opts.ehrenfest(), cap);
    } else {
        real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(), rt_opts, cap);
    }

    nlog.flush(); if (is_b2) trk.flush();
    const double Nf = electron_number();
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    std::ofstream s(OUT + "/run_summary.txt");
    s << std::setprecision(12)
      << "run            = cap_baselines/" << MODE << "\n"
      << "rs             = 5.69 (N=162, L=50, dx=0.40)\n"
      << "cap_eta_ha     = " << ETA << "\ncap_width_bohr = " << CAP_W_BOHR << " (per side)\n"
      << "launch_z       = " << LAUNCH_Z << "\nv0_au          = " << V0 << "\n"
      << "dt_au          = " << Cfg::DT_AU << "\nn_steps        = " << N_STEPS << "\n"
      << "write_every    = " << WRITE_EVERY << "\npropagator     = etrs\n"
      << "wp_state_index = " << wp_idx << "\n"
      << "N0             = " << N0 << "\nN_final        = " << Nf << "\n"
      << "absorbed_frac  = " << (N0 > 0 ? (N0 - Nf) / N0 : 0.0) << "\n"
      << "nan_seen       = " << (nan_seen ? "true" : "false") << "\n"
      << "wall_time_s    = " << wall << "\nrun_completed  = true\n";
    if (is_b2) s << "final_z        = " << ions.positions()[0][2]
                 << "\nfinal_vz       = " << ions.velocities()[0][2] << "\n";

    std::cout << "  done. N0=" << N0 << " N_final=" << Nf
              << " absorbed=" << (N0 > 0 ? (N0 - Nf) / N0 : 0.0)
              << " nan=" << nan_seen << " wall=" << wall << "s\n";
    return nan_seen ? 4 : 0;
}
