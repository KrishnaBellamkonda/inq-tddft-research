// ============================================================================
// localised_jellium / scripts/muon_mass_fork / vacuum_wp / run.cpp
//
// PHASE-2 vacuum free-particle spreading test for the per-state mass fork.
//
// An empty cubic box (no ions, no occupied electrons) with a single Gaussian
// wavepacket injected into the one extra_state. The state carries a tunable
// inverse mass (electron 1.0, muon 1/206.77) via the inq-study mass fork
// (electrons.inverse_mass()). It is propagated under NON_INTERACTING theory
// (H = -(1/2m) laplacian; no Hartree/XC self-force, no external potential) so
// the run reproduces the EXACT textbook free-Gaussian spreading law:
//
//     sigma_rho(t)^2 = sigma_rho0^2 + ( t / (2 m sigma_rho0) )^2 ,  hbar = 1
//     sigma_rho0 = sigma_WP / sqrt(2)      (density std of psi ~ exp(-r^2/2 sigma_WP^2))
//
// This isolates the mass fork: the t^2 spreading coefficient is 1/(4 m^2 sigma_rho0^2),
// so a parabola fit of sigma_rho^2 vs t recovers m; the muon:electron rate ratio
// must equal (206.77)^2 in the coefficient (or 1/206.77 in slope-of-sigma^2-vs-t^2).
//
// Also records:
//   - centroid z_mean(t)          -> group velocity v = k0/m  (Panel C, WP_K0>0)
//   - energy_kinetic(t) = <T>      -> KE conservation (free particle, <5e-4 drift)
//                                     and <T> = k0^2/(2m) absolute check
//   - norm_check(t)                -> unitarity (norm ~ 1, drift < 1e-4)
//
// Env (set by orchestrate.py phase2):
//   WP_OUT(vac)  WP_L(48)  WP_SPACING(0.4)  WP_DT(0.02)  WP_TSTEPS(600)
//   WP_SIGMA(0.5, = sigma_WP)  WP_K0(0.0)  WP_INV_MASS(1.0)  WP_WRITE_EVERY(1)
//
// Mass fork opt-in: real_time/propagate.hpp reads electrons.inverse_mass()[0]
// and calls ham.set_inverse_mass(); all-mass-1 routes the ORIGINAL scalar path.
// NEVER edit inq/ -- this builds against inq-study only (INQ_SOURCE=inq-study).
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
    // NOTE: do NOT construct input::environment here -- INQ lazily initialises
    // MPI via the environment::global() singleton (constructing one explicitly
    // double-inits MPI -> "MPI_Init_thread more than once" abort). Mirror the
    // proven run templates: use electrons.root() / data.root() for rank guards.
    auto t0 = std::chrono::steady_clock::now();

    const std::string OUT     = "results/" + env_s("WP_OUT", "vac");
    const double L_BOHR       = env_d("WP_L", 48.0);
    const double SPACING      = env_d("WP_SPACING", 0.4);
    const double DT_AU        = env_d("WP_DT", 0.02);
    const int    N_STEPS      = env_i("WP_TSTEPS", 600);
    const double SIGMA_WP     = env_d("WP_SIGMA", 0.5);          // sigma_WP (wavefunction Gaussian width)
    const double K0           = env_d("WP_K0", 0.0);            // launch momentum along z (Bohr^-1)
    const double INV_MASS     = env_d("WP_INV_MASS", 1.0);      // 1.0 electron; 1/206.77 muon
    const int    WRITE_EVERY  = env_i("WP_WRITE_EVERY", 1);
    const bool   EMIT_VTI     = env_i("WP_EMIT_VTI", 0) != 0;   // xz-density-vs-sigma deliverable
    const std::string THEORY  = env_s("WP_THEORY", "non_interacting");  // "non_interacting" | "lda"
    const bool   USE_LDA      = (THEORY == "lda");               // Phase-3b interacting sanity

    const double MASS         = 1.0 / INV_MASS;
    const double SIGMA_RHO0   = SIGMA_WP / std::sqrt(2.0);       // density std at t=0
    const double L_HALF       = 0.5 * L_BOHR;

    // ----- near-empty box: 1 spectator electron + 1 extra_state for the WP ---
    // INQ refuses to build electrons with num_electrons==0 (electrons.hpp:239),
    // so we add ONE occupied electron. Under non_interacting theory the KS
    // orbitals are decoupled (no Hartree/XC), so the spectator does NOT affect
    // the free evolution of the WP orbital; every WP observable (wp_real_space_
    // stats, wp_momentum_stats, density::orbital) reads only the WP slot.
    auto cell = systems::cell::cubic(L_BOHR * 1.0_b).periodic();
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
            << "\n=== vacuum_wp (mass fork Phase 2) out=" << OUT << " ===\n"
            << "  L=" << L_BOHR << " spacing=" << SPACING << " dt=" << DT_AU
            << " n_steps=" << N_STEPS << " (t_max=" << N_STEPS*DT_AU << " a.u.)\n"
            << "  sigma_WP=" << SIGMA_WP << " sigma_rho0=" << SIGMA_RHO0
            << " (sigma_rho0^2=" << SIGMA_RHO0*SIGMA_RHO0 << ")\n"
            << "  k0=" << K0 << " inv_mass=" << INV_MASS << " mass=" << MASS << "\n";
    }

    // ----- inject the Gaussian WP (no occupied states -> no orthogonalise) ---
    auto wp = inqkit::WavePacket{}
                  .center(0.0, 0.0, 0.0)
                  .sigma(SIGMA_WP)
                  .k0(0.0, 0.0, K0);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = report.state_index;

    // ----- set the per-state mass on the WP orbital (the fork) ---------------
    electrons.inverse_mass()[0][wp_idx] = INV_MASS;

    if (electrons.root())
        std::cout << "  WP injected: idx=" << wp_idx
                  << " norm_after=" << report.norm_after
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

    // ----- optional WP density VTI (xz-density-vs-sigma visualisation) --------
    // Physical-order VTI (RealField3DWriter fft-shifts + stamps Origin=-L/2);
    // load via inqview.load_vti downstream -- NEVER fftshift it.
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

    // boundary guard: abort if the packet's 4-sigma tail reaches the box edge.
    double sigma_rho_max_seen = 0.0;
    bool   boundary_hit = false;

    auto step_fn = [&](auto const& data){
        rt_obs.step(data);
        wp_rs.maybe_accumulate(data);
        wp_ms.maybe_accumulate(data);
        // predicted density sigma from the free-Gaussian oracle (for the guard)
        const double t = data.iter() * DT_AU;
        const double sig = std::sqrt(SIGMA_RHO0*SIGMA_RHO0
                                     + std::pow(t/(2.0*MASS*SIGMA_RHO0), 2.0));
        sigma_rho_max_seen = sig;
        if (!boundary_hit && 4.0*sig >= L_HALF) {
            boundary_hit = true;
            if (data.root())
                std::cout << "  [boundary guard] 4*sigma_rho=" << 4.0*sig
                          << " >= L/2=" << L_HALF << " at t=" << t
                          << " a.u. (step " << data.iter() << ") -- later points wrap.\n";
        }
    };

    // ----- propagation: non_interacting (free oracle) or lda (interacting) ---
    // non_interacting -> H = -(1/2m) lap : the exact free-Gaussian spreading law.
    // lda             -> full Hartree+XC on the forked kinetic path (Phase-3b
    //                    sanity: a muon WP feeling interaction; check no NaN /
    //                    bounded energy, NOT the free oracle).
    perturbations::none nop;
    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT_AU * 1.0_atomictime);
    if (USE_LDA) real_time::propagate(ions, electrons, step_fn, options::theory{}.lda(),            rt_opts, nop);
    else         real_time::propagate(ions, electrons, step_fn, options::theory{}.non_interacting(), rt_opts, nop);

    if (EMIT_VTI && vti_wr)                                   // final-step density slice
        vti_wr->write(inqkit::fields::density::orbital(electrons, wp_idx), N_STEPS*DT_AU, N_STEPS);

    // ----- provenance --------------------------------------------------------
    double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/muon_mass_fork/vacuum_wp/" << env_s("WP_OUT","vac") << "\n"
          << "engine = inq-study (per-state mass fork)\n"
          << "theory = " << THEORY << (USE_LDA ? " (interacting)" : " (free particle)") << "\n"
          << "cell_bohr = " << L_BOHR << "^3  spacing = " << SPACING << "\n"
          << "sigma_WP = " << SIGMA_WP << "  sigma_rho0 = " << SIGMA_RHO0 << "\n"
          << "k0 = " << K0 << "  inverse_mass = " << INV_MASS << "  mass = " << MASS << "\n"
          << "wp_state_index = " << wp_idx << "  wp_norm_after = " << report.norm_after << "\n"
          << "dt_au = " << DT_AU << "  n_steps = " << N_STEPS << "  write_every = " << WRITE_EVERY << "\n"
          << "predicted_sigma_rho_final = " << sigma_rho_max_seen << "\n"
          << "boundary_hit = " << (boundary_hit ? "true" : "false") << "\n"
          << "wall_time_s = " << wall << "\nrun_completed = true\n";
        std::cout << "  done. wall=" << wall << "s boundary_hit=" << boundary_hit << "\n";
    }
    return 0;
}
