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
// TRUE VACUUM (corrected 2026-07-27): extra_electrons(0) -> the WP is the ONLY
// electron (no uniform background gas). density_total == density_wp. The earlier
// extra_electrons(2.0) added a k=0 background that, though kinetically inert in
// non_interacting theory, polluted density_total plots -- removed.
//
// LOW-SPREADING DESIGN (CORRECTED 2026-07-27): the free Gaussian DENSITY width is
//   sigma_dens(t) = sqrt(sigma0^2/2 + t^2/(2 sigma0^2))   [atomic units, m_e=1]
// i.e. expansion factor R(t) = sqrt(1 + (t/sigma0^2)^2), spreading time tau =
// sigma0^2 (NOT 2 sigma0^2 -- an earlier note had a 2x error, verified vs data).
// Over a travel L = 5 sigma0 (box clearance) at v = k0 = sqrt(2E), the transit
// expansion is R = sqrt(1 + (5/(k0 sigma0))^2): dispersion is controlled ONLY by the
// product k0*sigma0. R<=1.05 needs k0*sigma0 >= 16. sigma0=3/100 eV gave k0*sigma0=8.1
// -> ~17% (too much). PRODUCTION (2026-07-27): sigma0=3, E=400 eV (k0=5.421,
// k0*sigma0=16.3) -> transit expansion ~5%. Grid h=0.4 so k_max=pi/h=7.85 > k0+4dk
// (dk=1/(sigma0 sqrt2)=0.236) -- cutoff-aliasing guard PASS. dt=0.01 (finer for the
// higher energy). Box 30x30x45 [-22.5,22.5], one-sided +z CAP z in [7.5,22.5] (WIDE
// W=15, eta=-1.0: full absorption exp(-2|eta|W/v)~4e-3, adiabatic => low reflection),
// launch z=-7.5 (5 sigma0 from CAP inner edge AND the wrapped back wall). The CAP run
// (NSTEPS=800) is fully absorbed -> NO wrapping remnant; the no-CAP control is run
// SHORT (NSTEPS~350) so the WP stops in the CAP region BEFORE it can wrap the box.
// Supersedes sigma0=3/100 eV (17% spread) and sigma0=1 (dispersing) runs.
//
// CAP strength (validated 2026-07-27): eta=-3.5 Ha gives survival
// exp(-|eta|W/v)=exp(-3.5*15/5.421)~6e-5 -> the wrapped remnant is BELOW the log-GIF
// floor (invisible); reflection stays 0.000 (adiabatic W=15). The no-CAP CONTROL is
// run SHORT (NSTEPS=350) so its WP stops in the CAP region before it can wrap.
//
// Env: WP_K0(5.421) WP_SIGMA(3.0) WP_ETA(0=no CAP; -3.5=CAP) WP_OUT
//      WP_LZ(45) WP_LPERP(30) WP_H(0.4) WP_DT(0.01) WP_NSTEPS(800 cap / 350 nocap)
//      WP_LAUNCH_Z(-7.5) WP_CAP_L(15) WP_WF_EVERY(20) WP_MOM_EVERY(1)
//
// EXTENSIVE-KINETIC TEST ADDITIONS (2026-07-29):
//   WP_CAP2=1    -> DOUBLE-SIDED CAP: a second absorbing band at the -z end
//                   (frac mid=-0.5+w/2, same eta/width), composed with the +z
//                   band via perturbations::sum. Geometry for this mode:
//                   LZ=60 box [-30,30], CAP_L=15 both ends, launch z=0
//                   (5 sigma0 clearance to BOTH CAP inner edges).
//   WP_EXTKIN=1  -> OrbitalKineticStats: per-orbital BARE kinetic + norm each
//                   step (orbital_kinetic_stats.csv) — the extensive-kinetic
//                   fix for the norm-divided energy_kinetic (energy.hpp:55).
//                   WP_EXTKIN_EVERY(1) sets its cadence. Its wall_ms column +
//                   the propagate wall-time in run_summary give the overhead.
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/momentum_distribution.hpp>
#include <inqkit/observables/orbital_kinetic_stats.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/absorbers/mask_absorber.hpp>

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <optional>
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
    const double K0       = env_d("WP_K0", 5.421);          // E=400 eV, m=1 (k0*sigma0=16.3)
    const double SIGMA    = env_d("WP_SIGMA", 3.0);        // low-spreading (~5% transit)
    const double ETA      = env_d("WP_ETA", 0.0);           // 0 -> no CAP
    const double LZ       = env_d("WP_LZ", 45.0);
    const double LPERP    = env_d("WP_LPERP", 30.0);
    const double H        = env_d("WP_H", 0.4);            // k_max=pi/h=7.85 > k0+4dk (guard)
    const double DT       = env_d("WP_DT", 0.01);
    const int    N_STEPS  = env_i("WP_NSTEPS", 800);
    const double LAUNCH_Z = env_d("WP_LAUNCH_Z", -7.5);
    const double CAP_L    = env_d("WP_CAP_L", 15.0);        // one-sided, +z end (wide, adiabatic)
    const int    WF_EVERY = env_i("WP_WF_EVERY", 20);
    const int    MOM_EVERY= env_i("WP_MOM_EVERY", 1);
    const int    CAP2     = env_i("WP_CAP2", 0);            // 1 -> double-sided CAP
    const int    EXTKIN   = env_i("WP_EXTKIN", 0);          // 1 -> OrbitalKineticStats
    const int    EK_EVERY = env_i("WP_EXTKIN_EVERY", 1);
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
        ions, options::electrons{}.spacing(H*1.0_b).extra_states(0).extra_electrons(1.0));
    // TRUE VACUUM: exactly ONE electron and ONE state -- and the WP REPLACES it
    // (rather than being added on top of a background). INQ needs >=1 electron for
    // the ground state and validates the count in propagate, so extra_electrons(0)
    // throws "no electrons". With extra_states(0).extra_electrons(1.0) the single
    // state holds the one electron; ground_state relaxes it to a uniform k=0 plane
    // wave (kinetic ~0), and inject_into_last_extra_state OVERWRITES that state (the
    // last = only state) with the Gaussian, occupation 1.0. The WP thus BECOMES the
    // sole electron: density_total == density_wp, no uniform background gas.
    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons, options::theory{}.non_interacting(),
                            options::ground_state{}.energy_tolerance(1.0e-8_Ha).max_steps(200));

    auto rep = inqkit::WavePacket{}.center(0.0,0.0,LAUNCH_Z).sigma(SIGMA).k0(0.0,0.0,K0)
                   .inject_into_last_extra_state(electrons, 1.0);
    const long wp_idx = rep.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx << " norm_after=" << rep.norm_after << "\n";

    // one-sided +z sin^2 CAP (inq-study makes its imaginary part reach orbitals);
    // WP_CAP2 adds the mirror -z band (frac mid = -0.5+w/2 -> physical
    // [-LZ/2, -LZ/2+CAP_L]); absorbing works in CONTRAVARIANT (cell-fraction)
    // coordinates, mid_pos in [-0.5, 0.5), matching rvector.
    perturbations::absorbing cap  (ETA*1.0_Ha,  mid_frac, width_frac);
    perturbations::absorbing cap_m(ETA*1.0_Ha, -mid_frac, width_frac);
    perturbations::sum       dcap (cap, cap_m);

    // --- absorber / propagator selection (energy-normalization investigation) ---
    // WP_ABS = cap (non-Hermitian CAP, default) | mask (spatial sin^2 mask on the WP
    //          orbital, SAME +z band [z_cap0, LZ/2]).
    // WP_PROP = etrs (default) | cn (Crank-Nicolson). Decisive test: mask+ETRS loses
    //          norm (reproduces the artifact); mask+CN renormalises each step -> norm
    //          held ~1 -> should show NO energy rise if the hypothesis holds.
    const std::string ABS  = env_s("WP_ABS",  "cap");
    const std::string PROP = env_s("WP_PROP", "etrs");
    const bool USE_MASK = (ABS == "mask");
    const bool USE_CN   = (PROP == "cn");
    inqkit::absorbers::MaskAbsorber mask(2, z_cap0, CAP_L, wp_idx);  // +z band, WP orbital
    std::cout << "  absorber=" << ABS << " propagator=" << PROP
              << (USE_MASK ? " (mask band [" : " (cap band [") << z_cap0 << "," << LZ/2.0 << "])\n";

    // ---- energy ledger: ALL components each step ----------------------------
    std::ofstream en(OUT + "/raw/observables/energies.csv");
    en << std::setprecision(12)
       << "step,time_au,total,kinetic,hartree,external,non_local,xc,exact_exchange,ion,ion_kinetic,wp_norm\n";

    // ---- momentum distribution (each step) + WP stats -----------------------
    obs_::MomentumDistribution mom_dist(OUT + "/raw/observables/momentum_distribution.csv", wp_idx, LZ,
                                        {.n_bins=64, .k_max_bohr_inv=0.0, .write_every=MOM_EVERY});
    obs_::WPMomentumStats  wp_mom(OUT + "/raw/observables/wp_momentum_stats.csv", wp_idx, {.write_every=MOM_EVERY});
    obs_::WPRealSpaceStats wp_rs (OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every=WF_EVERY});

    // extensive per-orbital kinetic + norm (the norm-division fix, measured in-run)
    std::optional<obs_::OrbitalKineticStats> extkin;
    if (EXTKIN)
        extkin.emplace(OUT + "/raw/observables/orbital_kinetic_stats.csv",
                       obs_::OrbitalKineticStatsConfig{.write_every=EK_EVERY});

    // ---- density writers for the density GIF --------------------------------
    using RLay = inqkit::io::RealField3DLayout;
    const auto vti = inqkit::io::VTIWriteOptions::Format::binary;
    RLay lay_wp{.field_name="density_wp", .include_meta=false, .emit_raw=false, .emit_vti=true, .vti_format=vti};
    RLay lay_tot{.field_name="density_total", .include_meta=false, .emit_raw=false, .emit_vti=true, .vti_format=vti};
    inqkit::io::RealField3DWriter wp_wr (OUT + "/raw/vti/density_wp",   lay_wp,  {.overwrite=true});
    inqkit::io::RealField3DWriter tot_wr(OUT + "/raw/vti/density_total",lay_tot, {.overwrite=true});
    wp_wr.write (inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    tot_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);

    auto step_fn = [&](auto const& data) {
        const int step = data.iter();
        auto e = data.energy();
        // WP norm = integral of |psi_wp|^2 (drops when CAP/mask absorbs it)
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
        if (extkin) extkin->maybe_accumulate(data);
        if (step % WF_EVERY == 0) {
            wp_wr.write (inqkit::fields::density::orbital(electrons, wp_idx), data.time(), step);
            tot_wr.write(inqkit::fields::density::total(electrons), data.time(), step);
        }
        // Spatial mask applied AFTER recording (each recorded step is the propagated
        // state; the mask acts between steps). ETRS keeps the removal -> norm decays;
        // CN renormalises the density next step -> norm held ~1.
        if (USE_MASK && step > 0) mask.apply(electrons);
    };

    auto theory = options::theory{}.non_interacting();
    const auto wall0 = std::chrono::steady_clock::now();
    if (USE_MASK) {
        if (USE_CN)
            real_time::propagate(ions, electrons, step_fn, theory,
                options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime).crank_nicolson());
        else
            real_time::propagate(ions, electrons, step_fn, theory,
                options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime));   // ETRS
    } else if (CAP2) {
        real_time::propagate(ions, electrons, step_fn, theory,
            options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime), dcap); // 2-sided CAP (ETRS)
    } else {
        real_time::propagate(ions, electrons, step_fn, theory,
            options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime), cap);  // CAP (ETRS)
    }
    const double wall_s = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - wall0).count();
    en.close();

    // final checkpoint (rule: every run checkpoints its last timestep)
    electrons.save(OUT + "/checkpoint");
    {
        std::ofstream rt(OUT + "/rt_state.txt");
        rt << std::setprecision(17)
           << "last_step=" << N_STEPS << "\ntime_au=" << N_STEPS*DT
           << "\ndt=" << DT << "\nwp_idx=" << wp_idx << "\n";
    }

    std::ofstream sum(OUT + "/run_summary.txt");
    sum << std::setprecision(12)
        << "run = vacuum/wp_traversal_energy/" << env_s("WP_OUT",(ETA==0.0?"nocap":"cap")) << "\n"
        << "engine = inq-study\ntheory = non_interacting\n"
        << "cap = " << (CAP_ON? "on":"off") << "  eta_Ha = " << ETA
        << "  cap_L = " << CAP_L << "  z_cap0 = " << z_cap0
        << "  cap_z = [" << z_cap0 << "," << LZ/2.0 << "]"
        << (CAP2 ? " + [-" : " (one-sided +z)")
        << (CAP2 ? std::to_string(LZ/2.0) + ",-" + std::to_string(z_cap0) + "] (two-sided)" : "")
        << "\nextkin = " << (EXTKIN? "on":"off") << "  extkin_every = " << EK_EVERY
        << "\npropagate_wall_s = " << wall_s
        << "  per_step_ms = " << 1000.0*wall_s/N_STEPS << "\n"
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
