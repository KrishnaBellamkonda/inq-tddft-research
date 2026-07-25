// ============================================================================
// cap_engine_probe / run.cpp
//
// Minimal EMPIRICAL test: does the stock inq `perturbations::absorbing` CAP
// actually absorb a wavepacket, or does only the inq-study fork?
//
// Free Gaussian WP launched +z in a vacuum box, with a one-sided sin^2 CAP near
// the +z edge. We record the total electron number N(t) = |WP|^2 each step.
//   * If the CAP works  -> N(t) decays from 1.0 toward 0 as the WP enters the CAP.
//   * If the CAP is dropped (imaginary term lost on a REAL potential field) ->
//     N(t) stays ~1.0 (WP just crosses / wraps, norm conserved).
//
// Build against EITHER engine by setting INQ_SOURCE before inq-run:
//   export INQ_SOURCE=/local/data/public/skcb2/tddft/inq        (stock)
//   export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study  (fork)
// Output: results/<PROBE_TAG>/norm_vs_time.csv
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

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

int main() {
    const std::string TAG   = env_s("PROBE_TAG", "stock");
    const double L          = env_d("PROBE_L", 36.0);
    const double SPACING    = env_d("PROBE_SPACING", 0.6);
    const double KZ         = env_d("PROBE_KZ", 3.0);
    const double SIGMA      = env_d("PROBE_SIGMA", 1.0);
    const double LAUNCH_Z   = env_d("PROBE_LAUNCH_Z", -10.0);
    const double CAP_ETA    = env_d("PROBE_CAP_ETA", -0.5);   // Ha
    const double CAP_MID    = env_d("PROBE_CAP_MID", 0.333);  // fractional (z=+12)
    const double CAP_WIDTH  = env_d("PROBE_CAP_WIDTH", 0.222);// fractional (8 Bohr)
    const int    N_STEPS    = env_i("PROBE_N_STEPS", 600);
    const double DT         = env_d("PROBE_DT", 0.02);
    const std::string OUT   = "results/" + TAG;
    std::filesystem::create_directories(OUT);

    std::cout << "\n=== cap_engine_probe (tag=" << TAG << ") ===\n"
              << "  L=" << L << " h=" << SPACING << " kz=" << KZ << " sigma=" << SIGMA
              << " cap_eta=" << CAP_ETA << " N_STEPS=" << N_STEPS << "\n";

    auto cell = systems::cell::cubic(L * 1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(2)
            .extra_states(3),
        input::kpoints::gamma());

    // Inject a free Gaussian WP into the last extra state, occupation 1.
    auto wp = inqkit::WavePacket{}
                  .center(0.0, 0.0, LAUNCH_Z).sigma(SIGMA)
                  .k0(0.0, 0.0, KZ);
    auto report = wp.inject_into_last_extra_state(electrons, 1.0);
    std::cout << "  WP injected: idx=" << report.state_index
              << " norm_after=" << report.norm_after << "\n";

    // One-sided sin^2 CAP near the +z edge (stock inq class; inq-study fork
    // makes its imaginary part actually reach the orbitals).
    perturbations::absorbing cap(CAP_ETA * 1.0_Ha, CAP_MID, CAP_WIDTH);

    std::ofstream nlog(OUT + "/norm_vs_time.csv");
    nlog << std::setprecision(12) << "step,time_au,N_total\n";

    auto step_fn = [&](auto const& data) {
        const int it = data.iter();
        if (data.root()) nlog << it << "," << (it*DT) << "," << data.num_electrons() << "\n";
    };

    auto rt_opts = options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime);
    real_time::propagate(ions, electrons, step_fn, options::theory{}.non_interacting(), rt_opts, cap);

    std::cout << "  done. wrote " << OUT << "/norm_vs_time.csv\n";
    return 0;
}
