// ============================================================================
// localised_jellium / muon_mass_fork / sigma1_massonly / vacuum / run.cpp
//
// VACUUM CALIBRATION for the σ=1 mass-only slab run (2026-07-09). Same box,
// grid, CAP, launch and mass as sigma1_massonly/wp — but NO jellium slab
// (empty box, 1 spectator + WP, NON_INTERACTING). It isolates and measures:
//   (1) CAP reflection — WP norm drainage + whether a reflected −z lobe appears
//       in the per-step WP density VTI (the direct test of the σ=2 reflection
//       complaint the user raised).
//   (2) True free spreading — σ_ρ(t) vs the oracle √(σ_ρ0² + (t/2mσ_ρ0)²).
// If reflection is non-trivial, retune EM_CAP_ETA / band before trusting the
// slab stopping. Cheap: 2 orbitals vs 82 — minutes. Adapted from vacuum_wp.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

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

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
    auto t0 = std::chrono::steady_clock::now();

    const std::string OUT     = "results/" + env_s("WP_OUT", "vac_sigma1");
    const double LX_BOHR      = env_d("WP_LX", 50.0);
    const double LY_BOHR      = env_d("WP_LY", 50.0);
    const double LZ_BOHR      = env_d("WP_LZ", 64.0);
    const double SPACING      = env_d("WP_SPACING", 0.40);
    const double DT_AU        = env_d("WP_DT", 0.05);
    const int    N_STEPS      = env_i("WP_TSTEPS", 1192);
    const double SIGMA_WP     = env_d("WP_SIGMA", 1.0);          // sigma_WP (wavefunction Gaussian width)
    const double K0           = env_d("WP_K0", 5.0356);         // launch momentum along z (fork k0)
    const double INV_MASS     = env_d("WP_INV_MASS", 0.289855); // 1/3.45
    const double LAUNCH_Z     = env_d("WP_LAUNCH_Z", -16.5);    // match the slab run
    const int    WRITE_EVERY  = env_i("WP_WRITE_EVERY", 4);
    const bool   EMIT_VTI     = env_i("WP_EMIT_VTI", 1) != 0;   // per-step WP density (reflection movie)
    const bool   USE_CAP      = env_i("WP_CAP", 1) != 0;
    // Same CAP as the slab run: gentler η=−0.6, sin² band peak |z|=26.25 (0.410·Lz),
    // width 11.5 Bohr (0.180·Lz) → band [±20.5,±32].
    const double CAP_ETA = env_d("WP_CAP_ETA", -0.6), CAP_MID = 26.25/64.0, CAP_WIDTH = 11.5/64.0;

    const double MASS         = 1.0 / INV_MASS;
    const double SIGMA_RHO0   = SIGMA_WP / std::sqrt(2.0);       // density std at t=0
    const double L_HALF       = 0.5 * LZ_BOHR;

    // ----- near-empty box: 1 spectator electron + 1 extra_state for the WP ---
    auto cell = systems::cell::orthorhombic(LX_BOHR * 1.0_b, LY_BOHR * 1.0_b, LZ_BOHR * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(1.0)
            .extra_states(1),
        input::kpoints::gamma());
    ground_state::initial_guess(ions, electrons);   // populate the spectator orbital

    if (electrons.root()) {
        std::cout << std::setprecision(8)
            << "\n=== vacuum calibration (sigma1_massonly) out=" << OUT << " ===\n"
            << "  box=" << LX_BOHR << "x" << LY_BOHR << "x" << LZ_BOHR << " spacing=" << SPACING
            << " dt=" << DT_AU << " n_steps=" << N_STEPS << " (t_max=" << N_STEPS*DT_AU << " a.u.)\n"
            << "  sigma_WP=" << SIGMA_WP << " sigma_rho0=" << SIGMA_RHO0
            << " k0=" << K0 << " inv_mass=" << INV_MASS << " mass=" << MASS << "\n"
            << "  launch_z=" << LAUNCH_Z << " cap=" << (USE_CAP?"on":"off") << " cap_eta=" << CAP_ETA << "\n";
    }

    // ----- inject the Gaussian WP at the slab-run launch point ---------------
    auto wp = inqkit::WavePacket{}
                  .center(0.0, 0.0, LAUNCH_Z)
                  .sigma(SIGMA_WP)
                  .k0(0.0, 0.0, K0);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;
    electrons.inverse_mass()[0][wp_idx] = INV_MASS;      // the mass fork

    if (electrons.root())
        std::cout << "  WP injected: idx=" << wp_idx << " norm_after=" << report.norm_after
                  << "  inverse_mass[" << wp_idx << "]=" << INV_MASS << "\n";

    // ----- observables -------------------------------------------------------
    std::filesystem::create_directories(OUT + "/raw/observables");
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = true;
    inqkit::io::ObservablesWriter obs_writer(OUT + "/raw/observables/observables.csv", sel);
    obs_writer.write_header();
    inqkit::observables::WPRealSpaceStats wp_rs(
        OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every = WRITE_EVERY});
    inqkit::observables::WPMomentumStats wp_ms(
        OUT + "/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every = WRITE_EVERY});

    // ----- per-step WP density VTI (reflection movie) ------------------------
    std::unique_ptr<inqkit::io::RealField3DWriter> vti_wr;
    if (EMIT_VTI) {
        std::filesystem::create_directories(OUT + "/raw/vti/density_wp");
        inqkit::io::RealField3DLayout lay{};
        lay.field_name = "density"; lay.emit_raw = false; lay.emit_vti = true;
        lay.vti_format = inqkit::io::VTIWriteOptions::Format::binary;
        vti_wr = std::make_unique<inqkit::io::RealField3DWriter>(
            OUT + "/raw/vti/density_wp", lay, inqkit::io::RealField3DWriteOptions{.overwrite = true});
        vti_wr->write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);   // t=0
    }

    inqkit::RealTimeSession rt_obs(ions, electrons, WRITE_EVERY);
    rt_obs.add([&](inqkit::StepContext const& ctx){ obs_writer.append(ctx); });

    double sigma_rho_pred_final = 0.0;
    auto step_fn = [&](auto const& data){
        rt_obs.step(data);
        wp_rs.maybe_accumulate(data);
        wp_ms.maybe_accumulate(data);
        const int    step = data.iter();
        const double t    = step * DT_AU;
        sigma_rho_pred_final = std::sqrt(SIGMA_RHO0*SIGMA_RHO0 + std::pow(t/(2.0*MASS*SIGMA_RHO0), 2.0));
        if (EMIT_VTI && vti_wr && step % WRITE_EVERY == 0 && step > 0)
            vti_wr->write(inqkit::fields::density::orbital(electrons, wp_idx), t, step);
    };

    // ----- propagation: non_interacting free WP + CAP ------------------------
    const double eta = USE_CAP ? CAP_ETA : 0.0;
    perturbations::absorbing cap_lo(eta * 1.0_Ha, -CAP_MID, CAP_WIDTH);
    perturbations::absorbing cap_hi(eta * 1.0_Ha,  CAP_MID, CAP_WIDTH);
    auto pert = perturbations::sum(cap_lo, cap_hi);
    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
    real_time::propagate(ions, electrons, step_fn, options::theory{}.non_interacting(), rt_opts, pert);

    if (EMIT_VTI && vti_wr)
        vti_wr->write(inqkit::fields::density::orbital(electrons, wp_idx), N_STEPS*DT_AU, N_STEPS);

    // ----- provenance --------------------------------------------------------
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/muon_mass_fork/sigma1_massonly/vacuum/" << env_s("WP_OUT","vac_sigma1") << "\n"
          << "engine = inq-study (per-state mass fork)\n"
          << "theory = non_interacting (free particle + CAP; no slab)\n"
          << "cell_bohr = " << LX_BOHR << "x" << LY_BOHR << "x" << LZ_BOHR << "  spacing = " << SPACING << "\n"
          << "sigma_WP = " << SIGMA_WP << "  sigma_rho0 = " << SIGMA_RHO0 << "\n"
          << "k0 = " << K0 << "  inverse_mass = " << INV_MASS << "  mass = " << MASS << "\n"
          << "velocity_au = " << K0*INV_MASS << "  launch_z = " << LAUNCH_Z << "\n"
          << "cap = " << (USE_CAP?"on":"off") << "  cap_eta = " << CAP_ETA
          << "  cap_mid = " << CAP_MID << "  cap_width = " << CAP_WIDTH << "\n"
          << "wp_state_index = " << wp_idx << "  wp_norm_after = " << report.norm_after << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
          << "predicted_sigma_rho_final_free = " << sigma_rho_pred_final << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
        std::cout << "  done. wall=" << wall << "s\n";
    }
    return 0;
}
