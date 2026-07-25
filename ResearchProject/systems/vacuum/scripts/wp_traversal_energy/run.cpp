// ============================================================================
// vacuum / scripts/wp_traversal_energy / run.cpp
//
// QUICK warm-up experiment (precedes the localised-jellium CAP campaign):
// a SINGLE Gaussian wavepacket traversing a VACUUM box, NON-INTERACTING, run
// WITH and WITHOUT a CAP. Purpose: watch total energy(t) and its DECOMPOSITION
//   * no CAP  (WP_ETA=0)      -> periodic box, WP wraps, energy_total CONSERVED.
//   * CAP     (WP_ETA=-0.7)   -> WP absorbed at the +z edge, energy_total DROPS
//                               to ~0 as the packet leaves the box.
//
// Non-interacting: E_total = E_kinetic; hartree=external=non_local=xc=0 (a
// bookkeeping check — all recorded each step in energies.csv).
//
// CAP = INQ's in-built perturbations::absorbing (region-restricted sin^2 i-pot),
// ONE-SIDED at the +z end only (a SINGLE absorbing band; two-sided needs two
// summed bands, as in the jellium runs). Functional ONLY on inq-study
// (scalar-potential complexification). Build:
//   export INQ_SOURCE=/local/data/public/skcb2/tddft/inq-study
//   inq-run --reconfig
//
// GEOMETRY (corrected 2026-07-24): the WP MUST launch >=5 sigma from every CAP
// boundary. The +z CAP occupies physical z in [LZ/2-CAP_L, LZ/2]; because the box
// is periodic the -z wall (z=-LZ/2) coincides with the CAP's OUTER edge (+LZ/2),
// so the launch must also be >=5 sigma from -LZ/2. Defaults: LZ=80 box [-40,40],
// CAP z in [30,40], WP launch z=-30 (10 sigma from the -40 wall, 60 Bohr traversal
// to the +30 CAP inner edge). The old LZ=60 / launch=-26 put the WP only 4 sigma
// from the wrapped CAP edge -- superseded.
//
// Env: WP_K0(2.711) WP_SIGMA(1.0) WP_ETA(0=no CAP; -0.7=CAP) WP_OUT
//      WP_LZ(80) WP_LPERP(12) WP_H(0.5) WP_DT(0.02) WP_NSTEPS(1600)
//      WP_LAUNCH_Z(-30) WP_CAP_L(10) WP_WF_EVERY(20) WP_MOM_EVERY(1)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
namespace fs = std::filesystem;
namespace obs_ = inqkit::observables;

static double env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int    env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

int main() {
    const double HA_TO_EV = 27.211386245988;
    const double K0       = env_d("WP_K0", 2.711);          // E=100 eV, m=1
    const double SIGMA    = env_d("WP_SIGMA", 1.0);
    const double ETA      = env_d("WP_ETA", 0.0);           // 0 -> no CAP
    const double LZ       = env_d("WP_LZ", 80.0);
    const double LPERP    = env_d("WP_LPERP", 12.0);
    const double H        = env_d("WP_H", 0.5);
    const double DT       = env_d("WP_DT", 0.02);
    const int    N_STEPS  = env_i("WP_NSTEPS", 1600);
    const double LAUNCH_Z = env_d("WP_LAUNCH_Z", -30.0);
    const double CAP_L    = env_d("WP_CAP_L", 10.0);        // one-sided, +z end
    const int    WF_EVERY = env_i("WP_WF_EVERY", 20);
    const int    MOM_EVERY= env_i("WP_MOM_EVERY", 1);
    const std::string OUT = "results/" + env_s("WP_OUT", (ETA==0.0? "nocap":"cap"));

    const double E_eV      = 0.5 * K0 * K0 * HA_TO_EV;
    const double width_frac= CAP_L / LZ;                    // +z last CAP_L Bohr
    const double mid_frac  = 0.5 - width_frac / 2.0;        // CAP z in [LZ/2-CAP_L, LZ/2]
    const double z_cap0    = LZ/2.0 - CAP_L;
    const bool   CAP_ON    = (ETA != 0.0);

    fs::create_directories(OUT + "/raw/observables");
    fs::create_directories(OUT + "/raw/vti/density_wp");
    fs::create_directories(OUT + "/raw/vti/density_total");

    std::cout << std::setprecision(6)
              << "\n=== wp_traversal_energy (out=" << OUT << ") ===\n"
              << "  k0=" << K0 << " E=" << E_eV << " eV sigma=" << SIGMA
              << " eta=" << ETA << " (CAP " << (CAP_ON?"ON":"OFF") << ")\n"
              << "  Lz=" << LZ << " Lperp=" << LPERP << " h=" << H
              << " launch_z=" << LAUNCH_Z << " z_cap0=" << z_cap0
              << " N_STEPS=" << N_STEPS << " dt=" << DT << "\n";

    auto cell = systems::cell::orthorhombic(LPERP*1.0_b, LPERP*1.0_b, LZ*1.0_b).periodic();
    auto ions = systems::ions(cell);
    auto electrons = systems::electrons(
        ions, options::electrons{}.spacing(H*1.0_b).extra_states(1).extra_electrons(2.0));
    // Relax the 2 "bath" electrons to their k=0 eigenstate (kinetic ~ 0) so the
    // reported total energy is the WP's alone (initial_guess would leave the bath
    // with tens of Ha of spurious kinetic that swamps the ~4 Ha WP signal).
    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons, options::theory{}.non_interacting(),
                            options::ground_state{}.energy_tolerance(1.0e-8_Ha).max_steps(200));

    auto rep = inqkit::WavePacket{}.center(0.0,0.0,LAUNCH_Z).sigma(SIGMA).k0(0.0,0.0,K0)
                   .inject_into_last_extra_state(electrons, 1.0);
    const long wp_idx = rep.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx << " norm_after=" << rep.norm_after << "\n";

    // one-sided +z sin^2 CAP (inq-study makes its imaginary part reach orbitals)
    perturbations::absorbing cap(ETA*1.0_Ha, mid_frac, width_frac);

    // ---- energy ledger: ALL components each step ----------------------------
    std::ofstream en(OUT + "/raw/observables/energies.csv");
    en << std::setprecision(12)
       << "step,time_au,total,kinetic,hartree,external,non_local,xc,exact_exchange,ion,ion_kinetic,wp_norm\n";

    // ---- momentum distribution (each step) + WP stats -----------------------
    obs_::MomentumDistribution mom_dist(OUT + "/raw/observables/momentum_distribution.csv", wp_idx, LZ,
                                        {.n_bins=64, .k_max_bohr_inv=0.0, .write_every=MOM_EVERY});
    obs_::WPMomentumStats  wp_mom(OUT + "/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every=MOM_EVERY});
    obs_::WPRealSpaceStats wp_rs (OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every=WF_EVERY});

    // ---- density writers for the density GIF --------------------------------
    using RLay = inqkit::io::RealField3DLayout;
    const auto vti = inqkit::io::VTIWriteOptions::Format::binary;
    RLay lay_wp{.field_name="density_wp", .include_meta=false, .emit_raw=false, .emit_vti=true, .vti_format=vti};
    RLay lay_tot{.field_name="density_total", .include_meta=false, .emit_raw=false, .emit_vti=true, .vti_format=vti};
    inqkit::io::RealField3DWriter wp_wr (OUT + "/raw/vti/density_wp",   lay_wp,  {.overwrite=true});
    inqkit::io::RealField3DWriter tot_wr(OUT + "/raw/vti/density_total",lay_tot, {.overwrite=true});
    wp_wr.write (inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    tot_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);

    real_time::propagate(
        ions, electrons,
        [&](auto const& data) {
            const int step = data.iter();
            auto e = data.energy();
            // WP norm = integral of |psi_wp|^2 (drops when CAP absorbs it)
            double wpn = std::numeric_limits<double>::quiet_NaN();
            if (data.root()) {
                en << step << ',' << data.time() << ',' << e.total() << ',' << e.kinetic()
                   << ',' << e.hartree() << ',' << e.external() << ',' << e.non_local()
                   << ',' << e.xc() << ',' << e.exact_exchange() << ',' << e.ion()
                   << ',' << e.ion_kinetic() << ',' << wpn << '\n';
            }
            mom_dist.maybe_accumulate(data);
            wp_mom.maybe_accumulate(data);
            wp_rs.maybe_accumulate(data);
            if (step % WF_EVERY == 0) {
                wp_wr.write (inqkit::fields::density::orbital(electrons, wp_idx), data.time(), step);
                tot_wr.write(inqkit::fields::density::total(electrons), data.time(), step);
            }
        },
        options::theory{}.non_interacting(),
        options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime),   // ETRS
        cap);
    en.close();

    std::ofstream sum(OUT + "/run_summary.txt");
    sum << std::setprecision(12)
        << "run = vacuum/wp_traversal_energy/" << env_s("WP_OUT",(ETA==0.0?"nocap":"cap")) << "\n"
        << "engine = inq-study\ntheory = non_interacting\n"
        << "cap = " << (CAP_ON? "on":"off") << "  eta_Ha = " << ETA
        << "  cap_L = " << CAP_L << "  z_cap0 = " << z_cap0
        << "  cap_z = [" << z_cap0 << "," << LZ/2.0 << "] (one-sided +z)\n"
        << "wp = gaussian sigma " << SIGMA << " k0 " << K0 << " E " << E_eV << " eV mass 1\n"
        << "launch_z = " << LAUNCH_Z << "  wp_state_index = " << wp_idx
        << "  launch_clearance_sigma = " << (LZ/2.0 + LAUNCH_Z)/SIGMA
        << " (from -z wall / wrapped CAP edge)\n"
        << "cell_bohr = " << LPERP << "x" << LPERP << "x" << LZ << "  spacing = " << H << "\n"
        << "dt_au = " << DT << "  n_steps = " << N_STEPS << "  wf_every = " << WF_EVERY
        << "  mom_every = " << MOM_EVERY << "\nrun_completed = true\n";
    sum.close();
    std::cout << "  done -> " << OUT << "/\n";
    return 0;
}
