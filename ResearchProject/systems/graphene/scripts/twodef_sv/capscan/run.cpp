// ============================================================================
// systems/graphene/scripts/twodef_sv/capscan/run.cpp
//
// PHASE 1 — CAP reflection/transmission scan for the two-definitions campaign.
// Plan: docs/plans/real-material-stopping-comparison.md (Stage D / Phase 1).
// Clone of localised_jellium/scripts/sigma56_sv/vac/run.cpp (the CAP-definition
// replica), adapted to the graphene production cell + cutoff-matched grid.
//
// A FREE electron WP in the EMPTY graphene box (no ions, non_interacting, one
// electron) with the production two-sided sin^2 CAP. Scans (eta, W, v): the
// residue (final norm) and reflection diagnostics pick the production CAP.
// Analytic transmission through one hump: exp(-2|eta|(W/2)/v) — reflection at
// low v is what the numerics must measure (no closed form).
//
// GRID NOTE: production graphene runs use .cutoff(50 Ha) (not .spacing), so
// this scan uses the SAME option — the grid must be the production grid or the
// CAP verdict does not transfer.
//
// Aliasing check (recorded, not gated — cutoff_guard rule): at sigma_WP=0.5,
// sigma_p = 1/(sqrt2*0.5) = 1.414; k_Nyq = sqrt(2*50) = 10 a.u.; worst case
// k0 = 4.70 (300 eV) -> (k_Nyq-k0)/sigma_p = 3.75 -> ~0.009% aliased: PASS.
//
// Env: CS_K0(2.71) CS_ETA(-1.0) CS_CAP_L(16) CS_LZ(80) CS_CUTOFF_HA(50)
//      CS_SIGMA(0.5) CS_LAUNCH_Z(-12) CS_DT(0.025) CS_NSTEPS(auto) CS_OUT(auto)
//      CS_SAVE_EVERY(0 = no VTI)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include "../../../shared/configs/twodef_gs.hpp"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace fs = std::filesystem;
namespace obs_ = inqkit::observables;
namespace Cfg = graphene_twodef;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
    const double HA_TO_EV = 27.211386245988;

    const double LX = Cfg::LX_BOHR, LY = Cfg::LY_BOHR;
    const double LZ     = env_d("CS_LZ", 80.0);
    const double CUTOFF = env_d("CS_CUTOFF_HA", Cfg::CUTOFF_HA);

    const double SIGMA    = env_d("CS_SIGMA", 0.5);       // sigma_WP (sigma-wp-convention)
    const double K0       = env_d("CS_K0", 2.71);         // = v (m = 1)
    const double LAUNCH_Z = env_d("CS_LAUNCH_Z", -12.0);

    const double CAP_L    = env_d("CS_CAP_L", 16.0);      // Bohr per face
    const double ETA      = env_d("CS_ETA", -1.0);        // Ha; <0 absorbs; 0 = control
    const double CAP_WIDTH_FRAC = CAP_L / LZ;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;
    const bool   CAP_ON   = (ETA != 0.0);

    const double DT = env_d("CS_DT", 0.025);
    // auto steps: launch -> +z face -> through CAP; buffer covers the slow
    // momentum tail (k0 - 2 sigma_p components). CS_BUFFER=2.0 for sigma=2 runs.
    const double BUFFER = env_d("CS_BUFFER", 1.4);
    const int N_STEPS_AUTO = int(std::ceil(BUFFER * (std::abs(LAUNCH_Z) + LZ/2.0) / K0 / DT));
    const int N_STEPS    = env_i("CS_NSTEPS", N_STEPS_AUTO);
    const int SAVE_EVERY = env_i("CS_SAVE_EVERY", 0);

    std::ostringstream outdef;
    outdef << "eta" << std::abs(ETA) << "_W" << int(CAP_L) << "_Lz" << int(LZ)
           << "_k" << K0 << "_s" << SIGMA;
    const std::string TAG = env_s("CS_OUT", outdef.str());
    const std::string OUT = "results/" + TAG;

    const double E_eV     = 0.5 * K0 * K0 * HA_TO_EV;
    const double z_cap_in = LZ/2.0 - CAP_L;
    const double survival = std::exp(-2.0 * std::abs(ETA) * (CAP_L/2.0) / K0);

    fs::create_directories(OUT + "/raw/observables");

    std::cout << std::setprecision(8)
              << "\n=== twodef_sv capscan [" << TAG << "] ===\n"
              << "  cell    = " << LX << " x " << LY << " x " << LZ
              << " Bohr, periodicity(2), cutoff=" << CUTOFF << " Ha (production grid)\n"
              << "  WP      : sigma_WP=" << SIGMA << "  k0=" << K0 << "  E=" << E_eV
              << " eV  launch_z=" << LAUNCH_Z << "\n"
              << "  CAP     : " << (CAP_ON ? "ON" : "OFF (CONTROL)") << "  eta=" << ETA
              << " Ha  W=" << CAP_L << " Bohr/face  bands +/-[" << z_cap_in
              << "," << LZ/2 << "]\n"
              << "  analytic one-hump transmission = " << survival << "\n"
              << "  dt=" << DT << "  N_STEPS=" << N_STEPS
              << "  t_total=" << DT*N_STEPS << " a.u.\n"
              << "  spread rate = " << 1.0/(std::sqrt(2.0)*SIGMA) << " Bohr/a.u.\n\n";

    // TRUE VACUUM: WP is the only electron (sigma56 vac construction)
    auto cell = systems::cell::orthorhombic(LX*1.0_b, LY*1.0_b, LZ*1.0_b).periodicity(2);
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions, options::electrons{}.cutoff(CUTOFF*1.0_Ha).extra_states(0).extra_electrons(1.0));

    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons, options::theory{}.non_interacting(),
                            options::ground_state{}.energy_tolerance(1.0e-8_Ha).max_steps(200));

    auto rep = inqkit::WavePacket{}.center(0.0, 0.0, LAUNCH_Z).sigma(SIGMA).k0(0.0, 0.0, K0)
                   .inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = rep.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx
              << "  norm_after=" << rep.norm_after << "\n";

    perturbations::absorbing cap_lo(ETA * 1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(ETA * 1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
    auto cap_both = perturbations::sum(cap_lo, cap_hi);

    obs_::WPMomentumStats  wp_mom(OUT + "/raw/observables/wp_momentum_stats.csv",  wp_idx, {.write_every = 1});
    obs_::WPRealSpaceStats wp_pos(OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every = 1});

    std::ofstream en(OUT + "/raw/observables/energies.csv");
    en << std::setprecision(12) << "step,time_au,total,kinetic\n";

    // t = 0 analytic gates (sigma56 vac pattern; relative tolerances — sigma=0.5
    // on this grid has ~2 points per sigma, discretisation error is real)
    const double sigma_p2 = 1.0 / (2.0 * SIGMA * SIGMA);
    {
        auto m0 = wp_mom.compute(electrons);
        auto r0 = wp_pos.compute(electrons);
        int fails = 0;
        auto gate_rel = [&](char const* nm, double got, double want, double relpc){
            const double rel = (want != 0.0) ? 100.0*(got-want)/std::abs(want) : 0.0;
            const bool ok = std::abs(rel) <= relpc;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << " (expect " << want << ", dev " << rel << " %)\n";
            if (!ok) ++fails;
        };
        auto gate_abs = [&](char const* nm, double got, double want, double tol){
            const bool ok = std::abs(got-want) <= tol;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << " (expect " << want << " +/- " << tol << ")\n";
            if (!ok) ++fails;
        };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        gate_abs("norm",                    r0.N,  1.0, 0.02);
        gate_rel("<p_z> = k0",              m0.pz, K0,  2.0);
        gate_rel("sigma_pz^2 = 1/(2 s^2)",  m0.sz2, sigma_p2, 10.0);
        gate_abs("centroid z",              r0.zc, LAUNCH_Z, 0.05);
        {
            const double sp   = 1.0/(std::sqrt(2.0)*SIGMA);
            const double knyq = std::sqrt(2.0*CUTOFF);
            const double tail = 0.5*std::erfc((knyq - K0)/sp/std::sqrt(2.0));
            std::cout << "  [info] k_Nyq=" << knyq << " sigma_p=" << sp
                      << " aliased tail=" << 100.0*tail << " %\n";
        }
        if (fails > 0) {
            std::cerr << "\nFATAL: " << fails << " t=0 gate(s) failed. Aborting.\n";
            return 4;
        }
        std::cout << "  all t=0 gates PASSED\n\n";
    }

    auto step_fn = [&](auto const& data) {
        const int step = data.iter();
        if (data.root()) {
            auto e = data.energy();
            en << step << ',' << data.time() << ',' << e.total() << ',' << e.kinetic() << '\n';
        }
        wp_mom.maybe_accumulate(data);
        wp_pos.maybe_accumulate(data);
    };

    auto theory = options::theory{}.non_interacting();
    auto opts   = options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    if (CAP_ON) real_time::propagate(ions, electrons, step_fn, theory, opts, cap_both);
    else        real_time::propagate(ions, electrons, step_fn, theory, opts);
    en.close();

    // final norm = the residue this scan exists to measure
    double norm_final = -1.0, pz_final = 0.0;
    {
        auto rF = wp_pos.compute(electrons);
        auto mF = wp_mom.compute(electrons);
        norm_final = rF.N; pz_final = mF.pz;
    }

    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = graphene/twodef_sv/capscan/" << TAG << "\n"
          << "purpose = Phase 1 CAP scan (free WP, production graphene box+grid)\n"
          << "plan = docs/plans/real-material-stopping-comparison.md\n"
          << "engine = inq-study\ntheory = non_interacting\n"
          << "cell_bohr = " << LX << " x " << LY << " x " << LZ
          << "  periodicity = 2  cutoff_ha = " << CUTOFF << "\n"
          << "wp_sigma_bohr = " << SIGMA << "  wp_k0 = " << K0
          << "  wp_energy_ev = " << E_eV << "  launch_z = " << LAUNCH_Z << "\n"
          << "cap_eta_ha = " << ETA << "  cap_width_bohr = " << CAP_L
          << "  cap_mid_frac = " << CAP_MID_FRAC
          << "  cap_width_frac = " << CAP_WIDTH_FRAC << "\n"
          << "analytic_one_hump_transmission = " << survival << "\n"
          << "dt_au = " << DT << "  n_steps = " << N_STEPS
          << "  total_time_au = " << DT*N_STEPS << "\n"
          << "norm_final = " << norm_final << "\n"
          << "pz_final = " << pz_final << "\n"
          << "run_completed = true\n";
    }
    std::cout << "  done [" << TAG << "]  norm_final=" << norm_final
              << "  pz_final=" << pz_final << "\n";
    return 0;
}
