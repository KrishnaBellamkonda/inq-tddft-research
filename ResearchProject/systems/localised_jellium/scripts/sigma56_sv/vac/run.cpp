// ============================================================================
// localised_jellium / scripts/sigma56_sv / vac / run.cpp
//
// CAP DEFINITION VALIDATION REPLICA (user request 2026-07-30: "make a smaller
// replica of a free wavepacket and a CAP to ensure you know how to define it").
//
// A FREE Gaussian electron wavepacket in an EMPTY box with the EXACT production
// geometry and the EXACT production CAP, so the CAP definition is proven before
// six 1.5x-length production runs are launched.
//
// Plan: docs/plans/sigma56-sv-twin.md   (sweep sigma56_sv)
// Supplies BOTH the CAP-definition check and the per-(sigma, v) CAP-only baselines.
//
// WHAT IS IDENTICAL TO PRODUCTION
//   cell 35 x 35 x 105 Bohr, orthorhombic, periodicity(2), dx = 0.40
//   sigma_WP = 5 or 6 Bohr, launch z = -27.5, dt = 0.04
//   CAP: TWO sin^2 absorbing bands, 12.5 Bohr wide per z face, |eta| = 1 Ha
//
// WHAT IS DELIBERATELY REMOVED (this is the "smaller replica")
//   no jellium slab, no background well, theory = non_interacting, ONE electron
//   (the WP itself). So the ONLY thing acting on the packet is free dispersion
//   plus the CAP -> any norm loss is unambiguously the CAP.
//
// ---------------------------------------------------------------------------
// THE UNIT TRAP THIS RUN EXISTS TO CATCH.  perturbations::absorbing takes
// (amplitude, mid_pos, width) where mid_pos and width are FRACTIONAL CELL
// COORDINATES, not Bohr: it compares point_op.rvector()[2], which uses the
// point_operator's CONTRAVARIANT spacing (inq/src/basis/real_space.hpp:105,129)
// and therefore lies in [-0.5, 0.5). The constructor's
// assert(mid_pos_ >= -0.5 and mid_pos_ < 0.5) is the tell. Passing Bohr would
// put a 12.5-Bohr CAP at z ~ 0.4 Bohr -- straight through the slab centre.
//
//   CAP_WIDTH_FRAC = 12.5 / 105 = 0.119047619048
//   CAP_MID_FRAC   = 0.5 - WIDTH/2 = 0.440476190476  ( = 46.25 Bohr )
//   +z band: z in [ +40.0, +52.5 ] Bohr      (peak |eta| at +46.25)
//   -z band: z in [ -52.5, -40.0 ] Bohr      (peak |eta| at -46.25)
//
// SIGN: eta < 0 absorbs. exp(-iVt) with V = i*eta*sin^2 gives exp(eta*sin^2*t),
// which decays only for eta < 0. The user specified "eta is 1", i.e. STRENGTH 1
// Ha; it is applied as -1.0 Ha because +1.0 would be an exponentially GROWING
// (gain) potential. Prior localised-jellium runs used -0.5 Ha.
//
// ENGINE: inq-study, NOT stock inq. Stock inq keeps the scalar potential real
// (field<real_space,double>) so the CAP's `vk[...] += complex(0.0, ...)` does
// not even compile; inq-study complexifies it (self_consistency.hpp:176,
// ks_hamiltonian.hpp). Set INQ_SOURCE=<repo>/inq-study before inq-run.
//
// ---------------------------------------------------------------------------
// WHAT TO LOOK FOR (the pass criteria, checked in the notebook)
//   1. GEOMETRY. norm stays ~1 while the packet is inside |z| < 30 and only
//      starts falling when its density reaches z = 30. A Bohr/fractional mix-up
//      would absorb it instantly at launch (or never).
//   2. ABSORPTION. With the CAP on, norm -> ~0 after the packet crosses +30.
//      Analytic survival through ONE sin^2 hump is exp(-2|eta|(W/2)/v)
//      = exp(-12.5/v): 1.9e-3 at v=2.0 and 6.2e-2 at v=4.5 (times a second
//      factor if it wraps through the -z band).
//   3. NO WRAP RE-ENTRY. With the CAP on, no density reappears at z ~ -42.5.
//      The WPC_ETA=0 CONTROL is the contrast: it MUST show the packet wrapping
//      the +z face and re-entering at -z, because periodicity(2) switches only
//      the electrostatics -- the orbital basis stays a 3-D FFT periodic in z
//      (verified in inq/src/basis/fourier_space.hpp, hamiltonian/ks_hamiltonian.hpp:200).
//   4. LOW REFLECTION. <p_z> must not go negative before absorption; a too-abrupt
//      CAP reflects. |eta|=1 over 12.5 Bohr is the adiabatic regime that the
//      vacuum study (vacuum/scripts/wp_traversal_energy) ran at eta=-1.0/W=15.
//
// WIDTH NOTE. At sigma_WP = 5/6 the packet spreads at 1/(sqrt2 sigma) = 0.141 /
// 0.118 Bohr per a.u. — an order of magnitude slower than the sigma = 0.5 runs
// this replica was first written for. The -z CAP inner edge is 12.5 Bohr behind
// the launch point (2.95 sigma_d at sigma = 6), so the t = 0 weight inside the
// absorber is 0.16 %, and the backward tail does not add materially to it. That
// one-off loss is exactly what this control measures, so the slab runs can quote
// it rather than absorb it silently into S.
//
// Env: WPC_K0(2.0) WPC_ETA(-1.0) WPC_NSTEPS(4360) WPC_DT(0.04) WPC_H(0.40)
//      WPC_SIGMA(6.0) WPC_LAUNCH_Z(-27.5) WPC_LX/LY/LZ(35/35/105)
//      WPC_SAVE_EVERY(14) WPC_CAP_L(12.5) WPC_OUT
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
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

    // ---- geometry: EVERYTHING env-driven, defaults = SlabN100_L35x35x105 -----
    // The ancestor hard-coded LZ = 85. It must be a knob here, because a vacuum
    // control is only subtractable from its slab twin if it shares the box: the
    // CAP fractions, and therefore WHEN the packet meets the absorber, are set by
    // L_z. Passing the slab run's L_z is the whole point of this binary.
    const double LX = env_d("WPC_LX", 35.0), LY = env_d("WPC_LY", 35.0);
    const double LZ = env_d("WPC_LZ", 105.0);
    const double H  = env_d("WPC_H", 0.40);       // dx — production grid
    const int    PER = 2;                         // periodicity(2)

    // ---- wavepacket -------------------------------------------------------
    const double SIGMA    = env_d("WPC_SIGMA", 6.0);      // sigma_WP (sigma-wp-convention)
    const double K0       = env_d("WPC_K0", 2.0);         // = v (m = 1)
    const double LAUNCH_Z = env_d("WPC_LAUNCH_Z", -27.5); // must match the slab twin

    // ---- CAP --------------------------------------------------------------
    const double CAP_L    = env_d("WPC_CAP_L", 12.5);     // Bohr, per face
    const double ETA      = env_d("WPC_ETA", -1.0);       // Ha; <0 absorbs; 0 = control
    const double CAP_WIDTH_FRAC = CAP_L / LZ;                        // 0.119047619048 at LZ=105
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;        // 0.440476190476
    const bool   CAP_ON   = (ETA != 0.0);

    // ---- propagation ------------------------------------------------------
    const double DT       = env_d("WPC_DT", 0.04);
    // Default matches the v = 2.0 production point so the baseline is subtractable
    // step for step; the dispatcher passes the per-velocity value.
    const int    N_STEPS  = env_i("WPC_NSTEPS", 4360);
    const int    SAVE_EVERY = env_i("WPC_SAVE_EVERY", 14);
    const std::string OUT = "results/" + env_s("WPC_OUT", (CAP_ON ? "cap" : "nocap"));

    const double E_eV     = 0.5 * K0 * K0 * HA_TO_EV;
    const double z_cap_in = LZ/2.0 - CAP_L;               // +40.0 Bohr at LZ = 105
    // survival of |psi|^2 through ONE sin^2 hump: integral W dz = |eta| * W/2
    const double survival = std::exp(-2.0 * std::abs(ETA) * (CAP_L/2.0) / K0);

    fs::create_directories(OUT + "/raw/observables");
    fs::create_directories(OUT + "/raw/vti/density_wp");

    std::cout << std::setprecision(8)
              << "\n=== sigma56_sv vac / CAP control (out=" << OUT << ") ===\n"
              << "  cell      = " << LX << " x " << LY << " x " << LZ
              << " Bohr, periodicity(" << PER << "), dx=" << H << "\n"
              << "  WP        : sigma_WP=" << SIGMA << " (density std "
              << SIGMA/std::sqrt(2.0) << ")  k0=" << K0
              << "  E=" << E_eV << " eV  launch_z=" << LAUNCH_Z << "\n"
              << "  CAP       : " << (CAP_ON ? "ON" : "OFF (CONTROL)")
              << "  eta=" << ETA << " Ha  width=" << CAP_L << " Bohr/face\n"
              << "  CAP frac  : mid=" << CAP_MID_FRAC << "  width=" << CAP_WIDTH_FRAC
              << "   -> +z band [" << z_cap_in << ", " << LZ/2.0 << "] Bohr,"
              << " -z band [" << -LZ/2.0 << ", " << -z_cap_in << "] Bohr\n"
              << "  predicted survival through one hump = " << survival << "\n"
              << "  dt=" << DT << "  N_STEPS=" << N_STEPS
              << "  t_total=" << (DT*N_STEPS) << " a.u.\n"
              << "  spreading : sigma_d(t) = sqrt(sigma^2/2 + t^2/(2 sigma^2)),"
              << " rate " << 1.0/(std::sqrt(2.0)*SIGMA) << " Bohr/a.u.\n\n";

    // ---- system: TRUE VACUUM, the WP is the only electron ------------------
    // extra_states(0).extra_electrons(1.0): INQ needs >= 1 electron, and
    // inject_into_last_extra_state then OVERWRITES that single state with the
    // Gaussian, so density_total == density_wp (no background gas polluting the
    // picture). Same construction as vacuum/scripts/wp_traversal_energy/run.cpp.
    auto cell0 = systems::cell::orthorhombic(LX*1.0_b, LY*1.0_b, LZ*1.0_b);
    auto cell  = (PER == 2) ? cell0.periodicity(2) : cell0.periodic();
    auto ions  = systems::ions(cell);
    auto electrons = systems::electrons(
        ions, options::electrons{}.spacing(H*1.0_b).extra_states(0).extra_electrons(1.0));

    ground_state::initial_guess(ions, electrons);
    ground_state::calculate(ions, electrons, options::theory{}.non_interacting(),
                            options::ground_state{}.energy_tolerance(1.0e-8_Ha).max_steps(200));

    auto rep = inqkit::WavePacket{}.center(0.0, 0.0, LAUNCH_Z).sigma(SIGMA).k0(0.0, 0.0, K0)
                   .inject_into_last_extra_state(electrons, 1.0);
    const int wp_idx = rep.state_index;
    std::cout << "  WP injected: state_index=" << wp_idx
              << "  norm_after=" << rep.norm_after << "\n";

    // ---- the two CAPs (production definition) ------------------------------
    perturbations::absorbing cap_lo(ETA * 1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(ETA * 1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
    auto cap_both = perturbations::sum(cap_lo, cap_hi);

    // ---- CAP geometry record (fractional -> Bohr, for the notebook) --------
    {
        std::ofstream cp(OUT + "/raw/observables/cap_profile.csv");
        cp << std::setprecision(12) << "z_bohr,W_ha\n";
        const int nz = int(std::lround(LZ / H));
        for (int i = 0; i < nz; ++i) {
            const int ii = (i < nz/2) ? i : i - nz;          // to_symmetric_range
            const double zf = double(ii) / double(nz);       // contravariant coord
            double W = 0.0;
            for (double mid : {-CAP_MID_FRAC, CAP_MID_FRAC}) {
                const double lo = mid - CAP_WIDTH_FRAC/2.0, hi = mid + CAP_WIDTH_FRAC/2.0;
                if (zf > lo && zf < hi)
                    W += std::abs(ETA) * std::pow(std::sin((zf - lo)*M_PI/CAP_WIDTH_FRAC), 2);
            }
            cp << zf*LZ << "," << W << "\n";                 // z in Bohr
        }
    }

    // ---- observables -------------------------------------------------------
    obs_::WPMomentumStats  wp_mom(OUT + "/raw/observables/wp_momentum_stats.csv",  wp_idx, {.write_every = 1});
    obs_::WPRealSpaceStats wp_pos(OUT + "/raw/observables/wp_real_space_stats.csv", wp_idx, {.write_every = 1});

    std::ofstream en(OUT + "/raw/observables/energies.csv");
    en << std::setprecision(12)
       << "step,time_au,total,kinetic,hartree,external,non_local,xc\n";

    inqkit::io::RealField3DLayout lay{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter wp_wr(OUT + "/raw/vti/density_wp", lay, {.overwrite = true});
    wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);

    // ---- t = 0 analytic gates ---------------------------------------------
    // sigma_p^2 = 1/(2 sigma^2) for psi ~ exp(-r^2/2 sigma^2).
    const double sigma_p2 = 1.0 / (2.0 * SIGMA * SIGMA);
    {
        auto m0 = wp_mom.compute(electrons);
        auto r0 = wp_pos.compute(electrons);
        const double T1 = m0.ekin;
        const double T2 = 0.5*(m0.px*m0.px + m0.py*m0.py + m0.pz*m0.pz);
        int fails = 0;
        // RELATIVE tolerances. sigma_WP = 0.5 on a dx = 0.5 grid is only ONE grid
        // point per sigma, so the momentum moments carry a real O(1%)
        // discretisation error (measured: sigma_pz^2 +1.6%, <p_z> -0.5% at k0=2).
        // The gates must accept that while still catching a genuine factor-2
        // blunder — hence percent-level relative bounds, not the 0.02 ABSOLUTE
        // bound that aborted the first cap_check (job 32416846).
        auto gate_rel = [&](char const* nm, double got, double want, double relpc){
            const double rel = (want != 0.0) ? 100.0*(got-want)/std::abs(want) : 0.0;
            const bool ok = std::abs(rel) <= relpc;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << ", dev " << rel << " %, tol +/-"
                      << relpc << " %)\n";
            if (!ok) ++fails;
        };
        auto gate_abs = [&](char const* nm, double got, double want, double tol){
            const bool ok = std::abs(got-want) <= tol;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << " +/- " << tol << ")\n";
            if (!ok) ++fails;
        };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        gate_abs("norm (real space)",        r0.N,  1.0, 0.02);
        gate_rel("<p_z> = k0",               m0.pz, K0,  2.0);
        gate_rel("sigma_pz^2 = 1/(2 s^2)",   m0.sz2, sigma_p2, 10.0);
        gate_rel("T1 = (k0^2+3 sp2)/2 (Ha)", T1, 0.5*(K0*K0 + 3.0*sigma_p2), 3.0);
        gate_rel("T1-T2 = 3/(4 s^2) (Ha)",   T1-T2, 3.0/(4.0*SIGMA*SIGMA), 5.0);
        gate_abs("centroid z (circular)",    r0.zc, LAUNCH_Z, 0.05);
        gate_rel("density std = s/sqrt2",    std::sqrt(r0.sz2), SIGMA/std::sqrt(2.0), 5.0);

        // Momentum-space aliasing diagnostic (.claude memory
        // reference_cutoff_aliasing_guard). The WP's k-distribution is centred at
        // k0 with std sigma_p = 1/(sqrt2 sigma); whatever lies beyond
        // k_Nyq = pi/dx folds back and corrupts <p^2> (i.e. T1) worst of all.
        // A classical Gaussian CHARGE has no such tail — this is WP-specific.
        {
            const double sp   = 1.0/(std::sqrt(2.0)*SIGMA);
            const double knyq = M_PI / H;
            const double zsc  = (knyq - K0)/sp;
            const double tail = 0.5*std::erfc(zsc/std::sqrt(2.0));
            std::cout << "  [info] k_Nyq=" << knyq << "  sigma_p=" << sp
                      << "  (k_Nyq-k0)/sigma_p=" << zsc
                      << "  -> " << 100.0*tail << " % of the z-momentum weight is "
                      << "BEYOND Nyquist and aliases\n";
        }
        if (fails > 0) {
            std::cerr << "\nFATAL: " << fails << " t=0 gate(s) failed — the injected "
                         "packet is not the one this run claims. Aborting.\n";
            return 4;
        }
        std::cout << "  all t=0 gates PASSED\n\n";
    }

    auto step_fn = [&](auto const& data) {
        const int step = data.iter();
        auto e = data.energy();
        if (data.root())
            en << step << ',' << data.time() << ',' << e.total() << ',' << e.kinetic()
               << ',' << e.hartree() << ',' << e.external() << ',' << e.non_local()
               << ',' << e.xc() << '\n';
        wp_mom.maybe_accumulate(data);
        wp_pos.maybe_accumulate(data);
        if (SAVE_EVERY > 0 && step % SAVE_EVERY == 0)
            wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), data.time(), step);
    };

    auto theory = options::theory{}.non_interacting();
    auto opts   = options::real_time{}.num_steps(N_STEPS).dt(DT*1.0_atomictime);
    if (CAP_ON) real_time::propagate(ions, electrons, step_fn, theory, opts, cap_both);
    else        real_time::propagate(ions, electrons, step_fn, theory, opts);
    en.close();

    if (electrons.root()) {
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(12)
          << "run = localised_jellium/wp_highdensity_sv/cap_check/" << env_s("WPC_OUT", CAP_ON?"cap":"nocap") << "\n"
          << "purpose = CAP definition validation replica (free WP, no slab)\n"
          << "plan = docs/plans/wavepacket-highdensity-sv-twin.md\n"
          << "engine = inq-study (scalar-potential complexification required for CAP)\n"
          << "theory = non_interacting\n"
          << "cell_bohr = " << LX << " x " << LY << " x " << LZ
          << "  periodicity = " << PER << "  spacing = " << H << "\n"
          << "wp_sigma_bohr = " << SIGMA << "  wp_sigma_density = " << SIGMA/std::sqrt(2.0) << "\n"
          << "wp_k0 = " << K0 << "  wp_energy_ev = " << E_eV
          << "  launch_z = " << LAUNCH_Z << "  wp_state_index = " << wp_idx << "\n"
          << "cap = " << (CAP_ON ? "on" : "off (CONTROL)")
          << "  eta_ha = " << ETA << "  cap_width_bohr = " << CAP_L << " per face\n"
          << "cap_mid_frac = " << CAP_MID_FRAC << "  cap_width_frac = " << CAP_WIDTH_FRAC << "\n"
          << "cap_band_hi_bohr = [" << z_cap_in << "," << LZ/2.0 << "]\n"
          << "cap_band_lo_bohr = [" << -LZ/2.0 << "," << -z_cap_in << "]\n"
          << "predicted_survival_one_hump = " << survival << "\n"
          << "dt_au = " << DT << "  n_steps = " << N_STEPS
          << "  total_time_au = " << (DT*N_STEPS) << "  save_every = " << SAVE_EVERY << "\n"
          << "spread_rate_bohr_per_au = " << 1.0/(std::sqrt(2.0)*SIGMA) << "\n"
          << "run_completed = true\n";
    }
    std::cout << "  done -> " << OUT << "/\n";
    return 0;
}
