// ============================================================================
// localised_jellium / scripts/wp_highdensity_sv / wp / run.cpp
//
// WAVEPACKET twin of the high-density classical S(v) benchmark.
// Campaign: docs/campaigns/localised_jellium/classical-highdensity-sv-benchmark.md
// Plan:     docs/plans/wavepacket-highdensity-sv-twin.md
//
// A Gaussian electron wavepacket replaces the classical Gaussian-charge
// projectile, with EVERY other physical parameter held at the classical
// campaign's values, so the classical curve (S = 1.087 -> 0.283 eV/Bohr over
// v = 2.0 -> 4.5) becomes the like-for-like reference for a quantum S(v).
//
// ---------------------------------------------------------------------------
// THE WIDTH MAPPING (.claude/rules/sigma-wp-convention.md)
// sigma always means sigma_WP, the psi-width. The classical campaign ran a
// Gaussian CHARGE of std sigma_pot = 0.35355 = 0.5/sqrt(2), i.e. the sigma_pot
// of sigma_WP = 0.5. This run therefore uses sigma_WP = 0.5, whose DENSITY std
// sigma_WP/sqrt(2) = 0.35355 equals the classical charge std exactly. Both
// halves are labelled sigma = 0.5. LJ_SIGMA is sigma_WP in BOTH binaries —
// the sqrt(2) lives inside each one, never at the call site.
//
// ---------------------------------------------------------------------------
// WHY THERE IS A CAP HERE WHEN THE CLASSICAL CAMPAIGN IS CAP-FREE
// The classical campaign is CAP-free by design: periodicity(2) makes z open for
// the ELECTROSTATICS, so the moving Gaussian CHARGE leaves the box and energy is
// conserved, giving an exact post-exit plateau. That mechanism does NOT transfer
// to a wavepacket. periodicity(2) is consulted only by the Poisson solver
// (inq/src/solvers/poisson.hpp:189,206), ionic replicas and the kick gauge; the
// wavefunction basis and kinetic operator are a plain 3-D FFT periodic in ALL
// three directions (inq/src/basis/fourier_space.hpp:60-151,
// inq/src/hamiltonian/ks_hamiltonian.hpp:200-204). A KS ORBITAL travelling +z
// therefore WRAPS and re-enters at -z, arriving behind the slab and
// contaminating the bath. (Confirmed in-repo: docs/handovers/pbc-open-z-
// oscillation.md:20 "wavefunction always wraps on the FFT grid (p2 switches
// electrostatics only)".)
//
// So two absorbing bands are added at the z faces (user decision 2026-07-30):
//   width 12.5 Bohr per face, |eta| = 1 Ha
// perturbations::absorbing takes FRACTIONAL cell coordinates (it compares
// point_op.rvector()[2], which uses the contravariant spacing and lies in
// [-0.5,0.5) — see real_space.hpp:105,129; the constructor's
// assert(mid_pos in [-0.5,0.5)) is the tell):
//   CAP_WIDTH_FRAC = 12.5/85 = 0.147058823529
//   CAP_MID_FRAC   = 0.5 - W/2 = 0.426470588235   ( = 36.25 Bohr )
//   +z band z in [+30.0,+42.5];  -z band z in [-42.5,-30.0]
// eta < 0 ABSORBS (exp(-iVt) = exp(eta sin^2 t)); +1.0 would be a gain medium.
// The CAP geometry and strength are validated independently, before production,
// by scripts/wp_highdensity_sv/cap_check/ (free WP, no slab).
//
// CONSEQUENCES OF THE CAP, which the notebooks must state:
//   * energy_total is NO LONGER CONSERVED (the CAP is non-Hermitian), so the
//     classical correctness gate is replaced by norm/absorption monitoring.
//   * the WP norm decays; WPMomentumStats divides every moment by the CURRENT
//     norm, so T1/T2 stay valid expectation values OF THE SURVIVING PACKET —
//     but the CAP preferentially removes the slow, spread tail, which biases
//     <p_z> upward. This is the qsp5 lesson (docs/handovers/qsp5-momentum-
//     stopping.md) and is why norm_check is recorded every step.
//
// ---------------------------------------------------------------------------
// WHY E_absorbed HAS NO DIRECT WP ANALOGUE
// In the classical run the projectile is an EXTERNAL moving perturbation doing
// work, so energy_total climbs and plateaus and E_absorbed = plateau - E_GS.
// Here the projectile IS part of the system, so (CAP aside) the Hamiltonian is
// time-independent and energy_total is conserved: there is no plateau to read.
// The deposited energy must come from the pairwise Coulomb ledger instead —
// hence interactions.csv is written EVERY step via compute_coulomb_wp, whose
// documented WP closure is
//     E_hartree = E_SS + E_PS + E_PP ,  E_external = E_SB + E_PB
// and which emits e_hartree_check / e_external_check as closure gates against
// INQ's own scalars.
//
// ---------------------------------------------------------------------------
// STOPPING POWER (definitions locked in docs/plans/bulk-jellium-ks-stopping.md §4)
//   T1 = <p^2>/2m   -> wp_momentum_stats.csv, e_kin_ha
//   T2 = <p>^2/2m   -> 0.5*(px_mean^2+py_mean^2+pz_mean^2), same file
//   s3 = WP density centroid -> wp_real_space_stats.csv, z_mean_circ (CIRCULAR:
//        the naive centroid is discontinuous across the periodic z face)
//   s4 = integral of <p_z> dt -> cumulative trapezoid of pz_mean
// S_ij = -dT_i/ds_j. Both stats files are written EVERY step: they ARE the
// measurement, and one extra FFT per step is negligible against ~74 orbitals.
//
// ANALYSIS WINDOW CAVEAT (sigma_WP = 0.5 disperses fast). sigma_d(t) =
// sqrt(sigma^2/2 + t^2/(2 sigma^2)) spreads at 1/(sqrt2 sigma) = 1.414 Bohr per
// a.u. The transverse periodic images overlap (6 sigma_d = L_xy = 35) at
// t = 4.12 a.u. = step 103, independent of velocity; the leading 3-sigma tail
// reaches the +z CAP by t ~ 7.6-10.7 a.u. The packet is projectile-like only for
// the first ~100-260 of 1611-3623 steps. Later steps are still recorded (they
// carry the density GIFs and the absorption physics) but must be excluded from
// the slope fit. T1 - T2 = 3/(4 sigma^2) = 3.0 Ha = 81.6 eV of localisation
// energy is likewise inherent to the matched width.
//
// ---------------------------------------------------------------------------
// CHECKPOINTS (.claude/rules/final-timestep-checkpoint.md + checkpoint-dont-block)
// The user asked for AT LEAST 4 resumable checkpoints per run. This writes:
//   * results/checkpoint          — rolling latest; what LJ_RESUME=1 loads
//   * results/ckpt_step<N>        — RETAINED numbered snapshots, one per
//                                   LJ_CKPT_EVERY (default N_STEPS/5 -> 5 of them)
//   * results/rt_state.txt        — last_step, time_au, dt, wp_idx
// so a killed job loses at most one interval and any quarter can be resumed from.
//
// ENGINE: inq-study, NOT stock inq (the CAP's complex scalar potential does not
// compile, let alone propagate, against stock inq).
//
// Env: LJ_K0(2.0) LJ_N_STEPS(3623) LJ_DT(0.04) LJ_SIGMA(0.5) LJ_LAUNCH_Z(-24)
//      LJ_CAP_ETA(-1.0) LJ_CAP_L(12.5) LJ_SAVE_EVERY(12) LJ_WF_EVERY(36)
//      LJ_CKPT_EVERY(0=auto N/5) LJ_RESUME(0) LJ_OUT(REQUIRED) LJ_GS_DIR(REQUIRED)
// ============================================================================
#include <inq/inq.hpp>
#include <inqkit/fields/density.hpp>
#include <inqkit/fields/orbital.hpp>
#include <inqkit/io/complex_field_3d_writer.hpp>
#include <inqkit/io/observables_writer.hpp>
#include <inqkit/io/real_field_3d_writer.hpp>
#include <inqkit/observables/density_delta.hpp>
#include <inqkit/observables/wp_momentum_stats.hpp>
#include <inqkit/observables/wp_real_space_stats.hpp>
#include <inqkit/real_time/real_time_session.hpp>
#include <inqkit/real_time/step_context.hpp>
#include <inqkit/wavepacket/injection_report.hpp>
#include <inqkit/wavepacket/wavepacket.hpp>
#include <inqkit/jellium/localised_background.hpp>
#include <inqkit/jellium/background_perturbation.hpp>
#include <inqkit/jellium/interaction_energies.hpp>
#include <inqkit/jellium/analytics.hpp>

#include "../../../shared/configs/slab_n100_L35x35x85.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

using namespace inq;
using namespace inq::magnitude;
using Cfg = localised_jellium::config::SlabN100_L35x35x85;

static double      env_d(const char* k, double d){ const char* v=std::getenv(k); return v?std::atof(v):d; }
static int         env_i(const char* k, int d){ const char* v=std::getenv(k); return v?std::atoi(v):d; }
static std::string env_s(const char* k, const std::string& d){ const char* v=std::getenv(k); return v?std::string(v):d; }

// rt_state.txt is written as `key=value` lines; tolerate several per line.
static double read_state_d(const std::string& path, const char* key, double def){
    std::ifstream f(path); std::string line; const std::string k = std::string(key) + "=";
    while (std::getline(f, line)) {
        auto p = line.find(k);
        if (p != std::string::npos) return std::atof(line.substr(p + k.size()).c_str());
    }
    return def;
}
static std::string tag6(int n){ std::ostringstream o; o << std::setw(6) << std::setfill('0') << n; return o.str(); }
static std::string iso_now(){
    auto t = std::time(nullptr); auto tm = *std::localtime(&t);
    char b[64]; std::strftime(b, sizeof(b), "%Y-%m-%dT%H:%M:%S", &tm); return std::string(b);
}

int main() {
    auto t_wall0 = std::chrono::steady_clock::now();
    const double HA = 27.211386245988;

    // ---- parameters (geometry pinned to the classical campaign's Cfg) -------
    const double SIGMA_WP = env_d("LJ_SIGMA", Cfg::WP_SIGMA_BOHR);   // 0.5
    // PRODUCTION GRID = 0.40, NOT the classical campaign's 0.50 (user decision
    // 2026-07-30). dx is a numerical, not physical, parameter: at dx = 0.5 the
    // sigma_WP = 0.5 packet is ONE grid point per sigma and aliases 10.4 % of its
    // z-momentum weight past k_Nyq = pi/dx at v = 4.5, biasing T1 = <p^2>/2 — the
    // primary stopping observable. At 0.40 that falls to <= 0.9 % everywhere, for
    // 2x the grid. The campaign's own <resolved_decisions> called dx = 0.50
    // PROVISIONAL and pre-authorised 0.35-0.40. The GS_DIR must match this dx.
    const double SPACING  = env_d("LJ_SPACING", 0.40);
    const double K0       = env_d("LJ_K0", 2.0);                     // = v (m = 1)
    const double LAUNCH_Z = env_d("LJ_LAUNCH_Z", -24.0);
    const double DT       = env_d("LJ_DT", 0.04);
    const int    N_STEPS  = env_i("LJ_N_STEPS", 3623);
    const double CAP_ETA  = env_d("LJ_CAP_ETA", -1.0);               // Ha, <0 absorbs
    const double CAP_L    = env_d("LJ_CAP_L", 12.5);                 // Bohr per face
    const int    SAVE_EVERY = env_i("LJ_SAVE_EVERY", 12);
    const int    WF_EVERY   = env_i("LJ_WF_EVERY", 36);
    const bool   RESUME     = env_i("LJ_RESUME", 0) != 0;
    const std::string OUT    = "results/" + env_s("LJ_OUT", "wp");
    const std::string GS_DIR = env_s("LJ_GS_DIR", "");

    int CKPT_EVERY = env_i("LJ_CKPT_EVERY", 0);
    if (CKPT_EVERY <= 0) CKPT_EVERY = std::max(1, N_STEPS / 5);      // >= 4 retained

    if (GS_DIR.empty() || !std::filesystem::exists(GS_DIR)) {
        std::cerr << "FATAL: LJ_GS_DIR missing or unset: '" << GS_DIR << "'\n";
        return 2;
    }

    const double CAP_WIDTH_FRAC = CAP_L / Cfg::LZ_BOHR;
    const double CAP_MID_FRAC   = 0.5 - CAP_WIDTH_FRAC / 2.0;
    const double z_cap_in       = Cfg::LZ_BOHR/2.0 - CAP_L;          // +30.0
    const bool   CAP_ON         = (CAP_ETA != 0.0);
    const double sigma_p2       = 1.0 / (2.0 * SIGMA_WP * SIGMA_WP); // = 2.0
    const double rs             = inqkit::jellium::rs_from_n0(Cfg::N0);

    // ---- resume state -------------------------------------------------------
    const std::string CKPT = OUT + "/checkpoint", RT_STATE = OUT + "/rt_state.txt";
    int START = 0, wp_idx = -1;
    if (RESUME) {
        START  = (int)read_state_d(RT_STATE, "last_step", -1);
        wp_idx = (int)read_state_d(RT_STATE, "wp_idx", -1);
        if (START < 0 || wp_idx < 0) {
            std::cerr << "FATAL: LJ_RESUME=1 but no readable " << RT_STATE << "\n"; return 2;
        }
        if (START >= N_STEPS) {
            std::cout << "Already at/after target (" << START << " >= " << N_STEPS
                      << "); nothing to do.\n"; return 0;
        }
    }
    const std::string SEG = (START > 0) ? (".from" + std::to_string(START)) : std::string("");

    std::cout << std::setprecision(10)
              << "\n=== wp_highdensity_sv  OUT=" << OUT << " ===\n"
              << "  cell      = " << Cfg::LX_BOHR << " x " << Cfg::LY_BOHR << " x "
              << Cfg::LZ_BOHR << " Bohr, periodicity(2), dx=" << SPACING << "\n"
              << "  slab      = 25 Bohr (half " << Cfg::SLAB_HALF_WIDTH << ", edge "
              << Cfg::EDGE_WIDTH_BOHR << "), N=" << Cfg::N_ELECTRONS
              << ", n0=" << Cfg::N0 << ", r_s=" << rs << "\n"
              << "  WP        = sigma_WP " << SIGMA_WP << " (density std "
              << SIGMA_WP/std::sqrt(2.0) << " == classical sigma_pot)"
              << "  k0=" << K0 << " (v=" << K0 << ")  E_drift="
              << 0.5*K0*K0*HA << " eV\n"
              << "  launch_z  = " << LAUNCH_Z << "\n"
              << "  CAP       = " << (CAP_ON ? "ON" : "OFF") << "  eta=" << CAP_ETA
              << " Ha  width=" << CAP_L << " Bohr/face  bands +/-["
              << z_cap_in << "," << Cfg::LZ_BOHR/2.0 << "]\n"
              << "  dt=" << DT << "  START=" << START << " -> N_STEPS=" << N_STEPS
              << "  t_total=" << (DT*N_STEPS) << " a.u." << (RESUME ? "  [RESUME]" : "") << "\n"
              << "  cadence   : density/" << SAVE_EVERY << "  wavefn/" << WF_EVERY
              << "  stats/1  ckpt/" << CKPT_EVERY << "\n"
              << "  GS        = " << GS_DIR << "\n"
              << "  spreading : rate " << 1.0/(std::sqrt(2.0)*SIGMA_WP)
              << " Bohr/a.u.;  T1-T2 = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA << " eV\n\n";

    // ---- system -------------------------------------------------------------
    auto cell = systems::cell::orthorhombic(Cfg::LX_BOHR * 1.0_b,
                                            Cfg::LY_BOHR * 1.0_b,
                                            Cfg::LZ_BOHR * 1.0_b).periodicity(2);
    auto ions = systems::ions(cell);                    // jellium: no nuclei
    auto electrons = systems::electrons(
        ions,
        options::electrons{}
            .spacing(SPACING * 1.0_b)
            .extra_electrons(Cfg::N_ELECTRONS)
            .extra_states(Cfg::EXTRA_STATES)
            .temperature(Cfg::TEMPERATURE_EV * 1.0_eV),
        input::kpoints::gamma());
    const int n_states = electrons.states().num_states();

    // ---- output skeleton ----------------------------------------------------
    const std::string OBS = OUT + "/raw/observables";
    const std::string VTI = OUT + "/raw/vti";
    for (auto const& d : {OBS, VTI + "/density_total", VTI + "/density_wp",
                          VTI + "/density_delta", VTI + "/density_delta_coarse",
                          VTI + "/density_gs_system", VTI + "/wavefunction_wp"})
        std::filesystem::create_directories(d);

    inqkit::InjectionReport report{};

    if (RESUME) {
        electrons.load(CKPT);
        std::cout << "  RESUMED from step " << START << " (wp_idx=" << wp_idx << ")\n";
    } else {
        electrons.load(GS_DIR);
        std::cout << "  Loaded GS from " << GS_DIR << "\n";

        // t=0 BATH density BEFORE the WP goes in — the baseline every induced
        // density in post-processing is measured against.
        inqkit::io::RealField3DLayout lay{
            .field_name = "density", .include_meta = false, .emit_raw = false,
            .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
        inqkit::io::RealField3DWriter gs_wr(VTI + "/density_gs_system", lay, {.overwrite = true});
        gs_wr.write(inqkit::fields::density::total(electrons), "density_gs_system");

        auto wp = inqkit::WavePacket{}
                      .center(0.0, 0.0, LAUNCH_Z)
                      .sigma(SIGMA_WP)
                      .k0(0.0, 0.0, K0)
                      .orthogonalise_against_occupied(electrons);
        report = wp.inject_into_last_extra_state(electrons, 1.0);
        wp_idx = report.state_index;
        std::cout << "  WP injected: state_index=" << wp_idx
                  << "  norm_after=" << report.norm_after
                  << "  max_overlap=" << report.max_overlap << "\n";

        std::ofstream f(OBS + "/wp_config.txt");
        f << std::setprecision(16)
          << "wp_center_bohr = 0 0 " << LAUNCH_Z << "\n"
          << "wp_sigma_bohr  = " << SIGMA_WP << "\n"
          << "wp_sigma_density = " << SIGMA_WP/std::sqrt(2.0) << "\n"
          << "wp_k0_bohr_inv = " << K0 << "\n"
          << "wp_state_index = " << wp_idx << "\n"
          << "norm_after     = " << report.norm_after << "\n"
          << "max_overlap    = " << report.max_overlap << "\n";
    }

    // ---- background well + the two CAPs ------------------------------------
    inqkit::jellium::localised_background_params bg;
    bg.shape      = inqkit::jellium::background_shape::slab;
    bg.n0         = Cfg::N0;
    bg.half_width = Cfg::SLAB_HALF_WIDTH;
    bg.slab_axis  = Cfg::SLAB_AXIS;
    bg.center     = {0.0, 0.0, Cfg::SLAB_CENTER_BOHR};
    bg.edge_width = Cfg::EDGE_WIDTH_BOHR;
    inqkit::jellium::localised_background_perturbation bg_pert(bg);

    perturbations::absorbing cap_lo(CAP_ETA * 1.0_Ha, -CAP_MID_FRAC, CAP_WIDTH_FRAC);
    perturbations::absorbing cap_hi(CAP_ETA * 1.0_Ha,  CAP_MID_FRAC, CAP_WIDTH_FRAC);
    auto pert_with_cap = perturbations::sum(bg_pert, perturbations::sum(cap_lo, cap_hi));

    // constant background fields for the pairwise ledger
    auto basis   = electrons.density().basis();
    auto nplus   = bg_pert.background_density(basis);
    auto phiplus = solvers::poisson::solve(nplus);
    const double E_BB = inqkit::jellium::background_self_energy(nplus, phiplus);

    // ---- writers (segment-suffixed on resume) -------------------------------
    inqkit::io::RealField3DLayout vti_layout{
        .field_name = "density", .include_meta = false, .emit_raw = false,
        .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary};
    inqkit::io::RealField3DWriter total_wr(VTI + "/density_total", vti_layout, {.overwrite = (START == 0)});
    inqkit::io::RealField3DWriter wp_wr   (VTI + "/density_wp",    vti_layout, {.overwrite = (START == 0)});
    inqkit::io::ComplexField3DWriter wf_wr(
        VTI + "/wavefunction_wp",
        {.field_name = "wavefunction", .include_meta = false, .emit_raw = false,
         .emit_vti = true, .vti_format = inqkit::io::VTIWriteOptions::Format::binary},
        {.overwrite = (START == 0)});

    // FULL energy decomposition — the user asked for every component.
    inqkit::io::ObservableSelection sel;
    sel.step = sel.time_au = true;
    sel.energy_total = sel.energy_kinetic = sel.energy_hartree = sel.energy_xc = true;
    sel.energy_external = sel.energy_nonlocal = sel.energy_ion = true;
    sel.energy_ion_kinetic = sel.energy_exact_exchange = true;
    sel.energy_nvxc = sel.energy_eigenvalues = true;
    sel.density_l2 = true;
    inqkit::io::ObservablesWriter obs(OBS + "/observables" + SEG + ".csv", sel);
    obs.write_header();

    inqkit::observables::WPMomentumStats  wp_mom(OBS + "/wp_momentum_stats"  + SEG + ".csv", wp_idx, {.write_every = 1});
    inqkit::observables::WPRealSpaceStats wp_pos(OBS + "/wp_real_space_stats" + SEG + ".csv", wp_idx, {.write_every = 1});

    // snapshot() is called EVERY step (the density_delta L2 is an every-step
    // scalar in observables.csv), but the 18 MB delta FIELD is only wanted at the
    // density cadence. Without .emit_every this wrote 3624 frames per run instead
    // of 302 -- 66 GB each, which filled /rds on 2026-07-31 and killed three
    // sigma=3 runs mid-flight. Aligning on SAVE_EVERY also matches delta frames
    // one-to-one with density_total / density_wp, which is what the GIF battery
    // pairs them against.
    inqkit::observables::DensityDelta density_delta(
        VTI + "/density_delta", VTI + "/density_delta_coarse",
        {.emit_raw_vti = true, .emit_coarse_vti = true,
         .compute_l2 = true, .coarse_bin_bohr = 3.0,
         .emit_every = (SAVE_EVERY > 0 ? SAVE_EVERY : 1)});

    std::ofstream ix;   // pairwise Coulomb ledger, Hartree units
    if (electrons.root()) {
        ix.open(OBS + "/interactions" + SEG + ".csv");
        ix << std::setprecision(12)
           << "step,time_au,e_ss,e_pp,e_ps,e_sb,e_pb,e_bb,"
              "e_hartree_check,e_external_check,norm_wp,norm_total\n";
    }

    // ---- t = 0 analytic gates (abort before burning GPU time) ---------------
    if (START == 0) {
        auto m0 = wp_mom.compute(electrons);
        auto r0 = wp_pos.compute(electrons);
        const double T1 = m0.ekin;
        const double T2 = 0.5*(m0.px*m0.px + m0.py*m0.py + m0.pz*m0.pz);
        int fails = 0;
        // RELATIVE tolerances: sigma_WP = 0.5 on a dx = 0.5 grid is ONE grid point
        // per sigma, so the momentum moments carry a real O(1%) discretisation
        // error (measured in cap_check: sigma_pz^2 +1.6%, <p_z> -0.5% at k0 = 2).
        // Percent-level bounds accept that while still catching a factor-2 blunder.
        auto gate_rel = [&](char const* nm, double got, double want, double relpc){
            const double rel = (want != 0.0) ? 100.0*(got-want)/std::abs(want) : 0.0;
            const bool ok = std::abs(rel) <= relpc;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << ", dev " << rel << " %, tol +/-"
                      << relpc << " %)\n";
            if (!ok) ++fails;
        };
        auto gate_abs = [&](char const* nm, double got, double want, double tol){
            const bool ok = std::abs(got - want) <= tol;
            std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << nm << ": " << got
                      << "  (expect " << want << " +/- " << tol << ")\n";
            if (!ok) ++fails;
        };
        std::cout << "\n  --- t=0 analytic gates ---\n";
        gate_abs("norm (real space)",         r0.N,  1.0, 0.02);
        gate_rel("<p_z> = k0",                m0.pz, K0,  2.0);
        gate_rel("sigma_pz^2 = 1/(2 s^2)",    m0.sz2, sigma_p2, 10.0);
        gate_rel("T1 = (k0^2+3 sp2)/2 (Ha)",  T1, 0.5*(K0*K0 + 3.0*sigma_p2), 3.0);
        gate_rel("T1-T2 = 3/(4 s^2) (Ha)",    T1-T2, 3.0/(4.0*SIGMA_WP*SIGMA_WP), 5.0);
        gate_abs("centroid z (circular)",     r0.zc, LAUNCH_Z, 0.05);
        gate_rel("density std = s/sqrt2",     std::sqrt(r0.sz2), SIGMA_WP/std::sqrt(2.0), 5.0);
        std::cout << "  [info] max_overlap with occupied manifold = "
                  << report.max_overlap << " (want < 1e-3)\n";

        // Momentum-space aliasing diagnostic. The WP k-distribution is centred at
        // k0 with std sigma_p = 1/(sqrt2 sigma) = 1.414; whatever lies beyond
        // k_Nyq = pi/dx folds back and corrupts <p^2> (T1) worst of all. At
        // dx = 0.5 this is 0.12 % at v = 2.0 but 10.4 % at v = 4.5 — a
        // WP-SPECIFIC problem the classical Gaussian CHARGE never had, and the
        // reason the campaign flagged dx = 0.5 as PROVISIONAL. Reported, not
        // gated, so the number lands in every run log and notebook.
        {
            const double sp   = 1.0/(std::sqrt(2.0)*SIGMA_WP);
            const double knyq = M_PI / SPACING;
            const double zsc  = (knyq - K0)/sp;
            const double tail = 0.5*std::erfc(zsc/std::sqrt(2.0));
            std::cout << "  [info] ALIASING: k_Nyq=" << knyq << " sigma_p=" << sp
                      << " (k_Nyq-k0)/sigma_p=" << zsc << " -> " << 100.0*tail
                      << " % of the z-momentum weight is beyond Nyquist\n";
            if (tail > 0.02)
                std::cout << "  [WARN] > 2 % aliased: T1 = <p^2>/2 is biased at this "
                             "velocity on dx=" << SPACING << ".\n";
        }
        if (fails > 0) {
            std::cerr << "\nFATAL: " << fails << " t=0 gate(s) failed — the injected "
                         "packet is not the one this run claims. Aborting.\n";
            return 4;
        }
        std::cout << "  all t=0 gates PASSED\n\n";

        total_wr.write(inqkit::fields::density::total(electrons), 0.0, 0);
        wp_wr.write(inqkit::fields::density::orbital(electrons, wp_idx), 0.0, 0);
    }

    // ---- real-time -----------------------------------------------------------
    inqkit::RealTimeSession rt(ions, electrons, 1);     // callback EVERY step
    rt.add([&](inqkit::StepContext const& ctx) {
        const int step = ctx.step;

        auto n_tot = inqkit::fields::density::total(*ctx.electrons);
        const double l2 = density_delta.snapshot(n_tot, ctx.time_au, step);
        inqkit::StepContext out = ctx; out.density_l2 = l2;
        obs.append(out);

        // pairwise Coulomb ledger with the WP as charge group P (Poisson-linearity
        // form, so the terms close EXACTLY against INQ's own E_hartree/E_external)
        auto n_wp_f = inqkit::jellium::orbital_density_field(*ctx.electrons, wp_idx);
        auto ct = inqkit::jellium::compute_coulomb_wp(ctx.electrons->density(), n_wp_f, phiplus);
        if (ctx.electrons->root())
            ix << step << ',' << ctx.time_au << ',' << ct.e_ss << ',' << ct.e_pp << ','
               << ct.e_ps << ',' << ct.e_sb << ',' << ct.e_pb << ',' << E_BB << ','
               << ct.e_hartree_check << ',' << ct.e_external_check << ','
               << ct.norm_wp << ',' << ct.norm_total << '\n';

        if (SAVE_EVERY > 0 && step % SAVE_EVERY == 0) {
            total_wr.write(n_tot, ctx.time_au, step);
            wp_wr.write(inqkit::fields::density::orbital(*ctx.electrons, wp_idx), ctx.time_au, step);
        }
        if (WF_EVERY > 0 && step % WF_EVERY == 0) {
            auto wf = inqkit::fields::orbital::wavefunction(*ctx.electrons, wp_idx);
            wf_wr.write(wf, "wavefunction_t" + tag6(step));
        }
    });

    auto write_rt_state = [&](int last){
        if (!electrons.root()) return;
        std::ofstream st(RT_STATE);
        st << std::setprecision(12)
           << "last_step=" << last << "\ntime_au=" << (last*DT) << "\ndt=" << DT
           << "\nwp_idx=" << wp_idx << "\nk0=" << K0 << "\nsigma_wp=" << SIGMA_WP << "\n";
    };

    real_time::propagate(
        ions, electrons,
        [&](auto const& data) {
            rt.step(data);
            wp_mom.maybe_accumulate(data);
            wp_pos.maybe_accumulate(data);
            // Interior checkpoints: rolling `checkpoint` (what LJ_RESUME loads) plus
            // a RETAINED numbered snapshot, so >= 4 distinct resume points survive.
            if (data.iter() > 0 && data.iter() % CKPT_EVERY == 0) {
                electrons.save(CKPT);
                write_rt_state(data.iter());
                const std::string snap = OUT + "/ckpt_step" + tag6(data.iter());
                electrons.save(snap);
                std::cout << "  [ckpt] step " << data.iter() << " -> " << snap << "\n" << std::flush;
            }
        },
        options::theory{}.lda(),
        options::real_time{}.num_steps(N_STEPS).dt(DT * 1.0_atomictime),
        pert_with_cap, START);

    // ---- final checkpoint ---------------------------------------------------
    electrons.save(CKPT);
    write_rt_state(N_STEPS);

    const double wall = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - t_wall0).count();

    if (electrons.root()) {
        ix.close();
        std::ofstream s(OUT + "/run_summary.txt");
        s << std::setprecision(14)
          << "run = localised_jellium/wp_highdensity_sv/wp/" << env_s("LJ_OUT","wp") << "\n"
          << "run_type = wavepacket projectile, localised jellium slab TDDFT (ALDA)\n"
          << "campaign = classical-highdensity-sv (WP twin)\n"
          << "plan = docs/plans/wavepacket-highdensity-sv-twin.md\n"
          << "engine = inq-study\nxc = LDA (ALDA in TDDFT)\n"
          << "date_finished = " << iso_now() << "\nwall_time_s = " << wall << "\n"
          << "cell_bohr = " << Cfg::LX_BOHR << "x" << Cfg::LY_BOHR << "x" << Cfg::LZ_BOHR << "\n"
          << "periodicity = 2\nspacing_bohr = " << SPACING << "\n"
          << "slab_half_width = " << Cfg::SLAB_HALF_WIDTH
          << "  edge_width = " << Cfg::EDGE_WIDTH_BOHR << "\n"
          << "n_electrons = " << Cfg::N_ELECTRONS << "  n0_a0m3 = " << Cfg::N0
          << "  r_s = " << rs << "\n"
          << "n_states = " << n_states << "  extra_states = " << Cfg::EXTRA_STATES << "\n"
          << "wp_enabled = yes\nwp_state_index = " << wp_idx << "\n"
          << "wp_sigma_bohr = " << SIGMA_WP << "\n"
          << "wp_sigma_note = wavepacket sigma (psi width); density std = this/sqrt2 "
             "= classical sigma_pot\n"
          << "wp_sigma_density = " << SIGMA_WP/std::sqrt(2.0) << "\n"
          << "wp_k0_bohr_inv = " << K0 << "  wp_velocity = " << K0 << "\n"
          << "wp_drift_energy_ev = " << 0.5*K0*K0*HA << "\n"
          << "wp_localisation_energy_ev = " << 3.0/(4.0*SIGMA_WP*SIGMA_WP)*HA << "\n"
          << "launch_z = " << LAUNCH_Z << "\n"
          << "norm_after = " << report.norm_after << "  max_overlap = " << report.max_overlap << "\n"
          << "cap = " << (CAP_ON ? "on" : "off") << "  cap_eta_ha = " << CAP_ETA
          << "  cap_width_bohr = " << CAP_L << " per face\n"
          << "cap_mid_frac = " << CAP_MID_FRAC << "  cap_width_frac = " << CAP_WIDTH_FRAC << "\n"
          << "cap_band_hi_bohr = [" << z_cap_in << "," << Cfg::LZ_BOHR/2.0 << "]\n"
          << "cap_band_lo_bohr = [" << -Cfg::LZ_BOHR/2.0 << "," << -z_cap_in << "]\n"
          << "cap_note = energy_total NOT conserved (non-Hermitian); gate on norm instead\n"
          << "start_step = " << START << "  rt_num_steps = " << N_STEPS
          << "  dt_au = " << DT << "  total_time_au = " << (DT*N_STEPS) << "\n"
          << "save_every = " << SAVE_EVERY << "  wf_every = " << WF_EVERY
          << "  stats_every = 1  ckpt_every = " << CKPT_EVERY << "\n"
          << "spread_rate_bohr_per_au = " << 1.0/(std::sqrt(2.0)*SIGMA_WP) << "\n"
          << "t_transverse_overlap_au = 4.12  (6 sigma_d = L_xy = 35; step "
          << int(4.12/DT) << ")\n"
          << "analysis_note = fit S over the pre-overlap window; later steps are "
             "recorded but the packet is delocalised\n"
          << "gs_dir = " << GS_DIR << "\nrun_completed = true\n";
    }
    std::cout << "\nDone. Wall time " << wall << " s.\n";
    return 0;
}
